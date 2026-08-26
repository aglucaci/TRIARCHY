#!/usr/bin/env python3
"""TRIARCHY workflow adapters and control-aware integration utilities."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import os
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

sys.path.insert(0, str(Path(__file__).resolve().parent))
from triarchy_evidence import (  # noqa: E402
    evaluate as evaluate_evidence,
    fit_calibration,
    run_inference,
)


PROFILE_FIELDS = [
    "library_id",
    "assay",
    "feature_id",
    "original_name",
    "original_rank",
    "abundance",
]

RANK_LEVEL = {
    "unclassified": 0,
    "root": 0,
    "kingdom": 1,
    "domain": 1,
    "phylum": 2,
    "class": 3,
    "order": 4,
    "family": 5,
    "genus": 6,
    "species": 7,
    "strain": 8,
    "asv": 8,
}
RANK_SUPPORT = {rank: level / 8 for rank, level in RANK_LEVEL.items()}
TAXON_SYNONYMS = {
    "propionibacterium acnes": "Cutibacterium acnes",
    "lactobacillus fermentum": "Limosilactobacillus fermentum",
    "filobasidiella neoformans": "Cryptococcus neoformans",
}


def clamp(value: float) -> float:
    return max(0.0, min(1.0, value))


def fmt(value: float) -> str:
    return f"{value:.6g}"


def infer_kingdom(row: dict[str, str]) -> str:
    lineage = row.get("canonical_lineage", "").casefold()
    name = row.get("canonical_name", row.get("original_name", "")).casefold()
    if "fung" in lineage or any(
        token in name for token in ("aspergillus", "candida", "cladospor")
    ):
        return "Fungi"
    if "archaea" in lineage or "methan" in name:
        return "Archaea"
    return "Bacteria"


def normalized_taxon_name(value: str) -> str:
    normalized = " ".join(value.replace("_", " ").split())
    return TAXON_SYNONYMS.get(normalized.casefold(), normalized)


def expected_assays(kingdom: str) -> tuple[str, ...]:
    if kingdom == "Fungi":
        return ("WMS", "ITS")
    return ("WMS", "16S")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cmd_split(args: argparse.Namespace) -> None:
    rows = read_tsv(args.samples)
    if not rows:
        raise SystemExit("sample manifest is empty")
    fields = list(rows[0])
    base_dir = getattr(args, "base_dir", None) or args.samples.parent
    base_dir = Path(base_dir).resolve()
    by_assay = {"WMS": [], "16S": [], "ITS": []}
    for row in rows:
        for field in ("fastq_r1", "fastq_r2"):
            value = row.get(field, "")
            if value and not Path(value).is_absolute():
                row[field] = str((base_dir / value).resolve())
        by_assay[row["assay"]].append(row)
    args.outdir.mkdir(parents=True, exist_ok=True)
    for assay, filename in (("WMS", "wms.tsv"), ("16S", "16s.tsv"), ("ITS", "its.tsv")):
        write_tsv(args.outdir / filename, by_assay[assay], fields)


def cmd_bracken(args: argparse.Namespace) -> None:
    rows = []
    for row in read_tsv(args.input):
        rows.append(
            dict(
                library_id=args.library_id,
                assay="WMS",
                feature_id=row.get("taxonomy_id", ""),
                original_name=row.get("name", "Unclassified"),
                original_rank=row.get("taxonomy_lvl", "S").lower(),
                abundance=row.get("new_est_reads", "0"),
            )
        )
    write_tsv(args.output, rows, PROFILE_FIELDS)


def cmd_humann(args: argparse.Namespace) -> None:
    candidates = sorted(args.input_dir.glob("*_pathabundance.tsv"))
    rows = []
    if candidates:
        with candidates[0].open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
            for values in reader:
                if not values or values[0].startswith("#"):
                    continue
                rows.append(
                    {
                        "library_id": args.library_id,
                        "assay": "WMS_function",
                        "feature_id": values[0],
                        "original_name": values[0],
                        "abundance": values[1] if len(values) > 1 else "0",
                    }
                )
    write_tsv(
        args.output,
        rows,
        ["library_id", "assay", "feature_id", "original_name", "abundance"],
    )


def profile_files(paths: list[Path], filename: str) -> list[Path]:
    files = []
    for path in paths:
        candidate = path / filename if path.is_dir() else path
        if candidate.exists():
            files.append(candidate)
    return files


def cmd_summarize_qc(args: argparse.Namespace) -> None:
    """Standardize fastp, cutadapt, and host-removal QC into one assay table."""
    samples = {row["library_id"]: row for row in read_tsv(args.samples)}
    metrics: dict[str, dict[str, object]] = {}

    for root in args.qc_dirs:
        candidates = [root] if root.is_file() else list(root.rglob("*"))
        for path in [item for item in candidates if item.name.endswith(".fastp.json") or item.name == "fastp.json"]:
            library_id = (
                path.name.removesuffix(".fastp.json")
                if path.name != "fastp.json" else path.parent.name
            )
            payload = json.loads(path.read_text(encoding="utf-8"))
            before = payload.get("summary", {}).get("before_filtering", {})
            after = payload.get("summary", {}).get("after_filtering", {})
            input_reads = float(before.get("total_reads", 0))
            retained_reads = float(after.get("total_reads", input_reads))
            q30 = float(after.get("q30_rate", 1))
            metrics[library_id] = {
                "library_id": library_id, "assay": "WMS",
                "input_reads": fmt(input_reads), "retained_reads": fmt(retained_reads),
                "q30_rate": fmt(q30), "host_fraction": "0",
                "quality_factor": fmt(
                    clamp((retained_reads / input_reads if input_reads else 1) * q30)
                ),
                "qc_status": "failed" if input_reads and retained_reads == 0 else "pass",
                "qc_sources": str(path),
            }
        for path in [item for item in candidates if item.name.endswith(".cutadapt.json") or item.name == "cutadapt.json"]:
            library_id = (
                path.name.removesuffix(".cutadapt.json")
                if path.name != "cutadapt.json" else path.parent.name
            )
            payload = json.loads(path.read_text(encoding="utf-8"))
            counts = payload.get("read_counts", {})
            input_reads = float(counts.get("input", 0))
            retained_reads = float(counts.get("output", input_reads))
            factor = retained_reads / input_reads if input_reads else 1
            assay = samples.get(library_id, {}).get("assay", "")
            metrics[library_id] = {
                "library_id": library_id, "assay": assay,
                "input_reads": fmt(input_reads), "retained_reads": fmt(retained_reads),
                "q30_rate": "", "host_fraction": "0",
                "quality_factor": fmt(clamp(factor)),
                "qc_status": (
                    "failed" if input_reads and retained_reads == 0
                    else "low_depth" if retained_reads and retained_reads < 1000
                    else "pass"
                ),
                "qc_sources": str(path),
            }
        for path in [item for item in candidates if item.name.endswith(".bowtie2.log") or item.name == "bowtie2.log"]:
            library_id = (
                path.name.removesuffix(".bowtie2.log")
                if path.name != "bowtie2.log" else path.parent.name
            )
            text = path.read_text(encoding="utf-8", errors="replace")
            match = re.search(r"([0-9.]+)% overall alignment rate", text)
            host_fraction = float(match.group(1)) / 100 if match else 0
            row = metrics.setdefault(
                library_id,
                {
                    "library_id": library_id,
                    "assay": samples.get(library_id, {}).get("assay", "WMS"),
                    "input_reads": "", "retained_reads": "", "q30_rate": "",
                    "quality_factor": "1", "qc_status": "pass", "qc_sources": "",
                },
            )
            row["host_fraction"] = fmt(host_fraction)
            row["quality_factor"] = fmt(
                clamp(float(row.get("quality_factor") or 1) * (1 - host_fraction))
            )
            if host_fraction >= .8:
                row["qc_status"] = "host_dominated"
            row["qc_sources"] = ";".join(
                item for item in (str(row.get("qc_sources", "")), str(path)) if item
            )

    for library_id, sample in samples.items():
        metrics.setdefault(
            library_id,
            {
                "library_id": library_id, "assay": sample["assay"],
                "input_reads": "", "retained_reads": "", "q30_rate": "",
                "host_fraction": "", "quality_factor": "1",
                "qc_status": "not_observed", "qc_sources": "",
            },
        )
    write_tsv(
        args.output,
        [metrics[key] for key in sorted(metrics)],
        [
            "library_id", "assay", "input_reads", "retained_reads", "q30_rate",
            "host_fraction", "quality_factor", "qc_status", "qc_sources",
        ],
    )


def cmd_control(args: argparse.Namespace) -> None:
    samples = {row["library_id"]: row for row in read_tsv(args.samples)}
    edges = read_tsv(args.controls)
    linked: dict[str, list[str]] = defaultdict(list)
    linked_stages: dict[str, list[str]] = defaultdict(list)
    for edge in edges:
        linked[edge["library_id"]].append(edge["control_library_id"])
        linked_stages[edge["library_id"]].append(edge["control_stage"])

    tax_rows = [
        row
        for path in profile_files(args.taxa, "profile.tsv")
        for row in read_tsv(path)
    ]
    write_tsv(args.unfiltered_taxa, tax_rows, PROFILE_FIELDS)
    lookup: dict[tuple[str, str], float] = {}
    for row in tax_rows:
        lookup[(row["library_id"], row["original_name"].casefold())] = float(
            row.get("abundance") or 0
        )
    library_totals: dict[str, float] = defaultdict(float)
    library_values: dict[str, list[float]] = defaultdict(list)
    for row in tax_rows:
        value = max(0.0, float(row.get("abundance") or 0))
        library_totals[row["library_id"]] += value
        if value > 0:
            library_values[row["library_id"]].append(value)

    annotated = []
    high = []
    assay_biological = defaultdict(set)
    assay_controls = defaultdict(set)
    for library_id, sample in samples.items():
        target = assay_biological if sample["sample_role"] == "biological" else assay_controls
        target[sample["assay"]].add(library_id)

    for row in tax_rows:
        library_id = row["library_id"]
        assay = row["assay"]
        feature_key = row["original_name"].casefold()
        abundance = float(row.get("abundance") or 0)
        total = library_totals[library_id]
        relative_abundance = abundance / total if total else 0.0
        ordered = sorted(
            {
                float(candidate.get("abundance") or 0)
                for candidate in tax_rows
                if candidate["library_id"] == library_id
            },
            reverse=True,
        )
        abundance_rank = ordered.index(abundance) + 1
        positive = library_values[library_id]
        pseudocount = min(positive) / 2 if positive else 0.5
        log_values = [math.log(value + pseudocount) for value in positive]
        clr_center = sum(log_values) / len(log_values) if log_values else 0.0
        clr_abundance = math.log(abundance + pseudocount) - clr_center
        control_values = [
            lookup.get((control_id, feature_key), 0.0)
            for control_id in linked.get(library_id, [])
        ]
        stage_values: dict[str, list[float]] = defaultdict(list)
        for control_id, stage in zip(linked.get(library_id, []), linked_stages.get(library_id, [])):
            normalized_stage = stage.casefold().replace("-", "_")
            if normalized_stage in {"ext", "extraction_blank"}:
                normalized_stage = "extraction"
            elif normalized_stage in {"pcr_blank", "amplification"}:
                normalized_stage = "pcr"
            elif normalized_stage in {"lib", "library_blank"}:
                normalized_stage = "library"
            elif normalized_stage in {"sequencing", "seq", "run_blank"}:
                normalized_stage = "run"
            elif normalized_stage not in {"field", "extraction", "pcr", "library", "run"}:
                normalized_stage = "other"
            stage_values[normalized_stage].append(lookup.get((control_id, feature_key), 0.0))
        stage_scores = {
            stage: clamp(max(values, default=0.0) / max(abundance, 1e-12))
            for stage, values in stage_values.items()
        }
        max_control = max(control_values, default=0.0)
        ratio = max_control / max(abundance, 1e-12)
        control_ids = set(linked.get(library_id, []))
        biological_ids = assay_biological[assay]
        control_prevalence = (
            sum(lookup.get((item, feature_key), 0) > 0 for item in control_ids)
            / len(control_ids)
            if control_ids
            else 0.0
        )
        biological_prevalence = (
            sum(lookup.get((item, feature_key), 0) > 0 for item in biological_ids)
            / len(biological_ids)
            if biological_ids
            else 0.0
        )
        prevalence_signal = control_prevalence / max(biological_prevalence, 1e-12)
        batch_totals: dict[str, float] = defaultdict(float)
        for candidate_id in biological_ids:
            candidate = samples[candidate_id]
            batch_totals[candidate["extraction_batch"]] += lookup.get(
                (candidate_id, feature_key), 0.0
            )
        batch_total = sum(batch_totals.values())
        batch_score = max(batch_totals.values(), default=0.0) / batch_total if batch_total else 0.0
        abundance_signal = clamp(ratio)
        prevalence_evidence = clamp(prevalence_signal / 2)
        batch_evidence = batch_score if max_control > 0 else 0.0
        probability = clamp(
            0.45 * abundance_signal
            + 0.35 * prevalence_evidence
            + 0.20 * batch_evidence
        )
        decision = "exclude" if probability >= 0.5 else "retain"
        annotated_row = {
            **row,
            "native_count": fmt(abundance),
            "library_total_count": fmt(total),
            "native_normalized_value": fmt(relative_abundance),
            "normalization_method": "within_library_relative_abundance",
            "abundance_rank": abundance_rank,
            "clr_abundance": fmt(clr_abundance),
            "control_prevalence": fmt(control_prevalence),
            "biological_prevalence": fmt(biological_prevalence),
            "prevalence_contamination_score": fmt(prevalence_evidence),
            "control_abundance_ratio": fmt(ratio),
            "abundance_contamination_score": fmt(abundance_signal),
            "batch_localization_score": fmt(batch_score),
            "batch_contamination_score": fmt(batch_evidence),
            "field_contamination_score": fmt(stage_scores.get("field", 0.0)),
            "extraction_contamination_score": fmt(stage_scores.get("extraction", 0.0)),
            "pcr_contamination_score": fmt(stage_scores.get("pcr", 0.0)),
            "library_contamination_score": fmt(stage_scores.get("library", 0.0)),
            "run_contamination_score": fmt(stage_scores.get("run", 0.0)),
            "other_control_contamination_score": fmt(stage_scores.get("other", 0.0)),
            "probable_contamination_source": max(
                stage_scores,
                key=stage_scores.get,
                default="none",
            ) if max(stage_scores.values(), default=0.0) > 0 else "none",
            "contamination_probability": fmt(probability),
            "associated_controls": ";".join(linked.get(library_id, [])),
            "associated_control_types": ";".join(linked_stages.get(library_id, [])),
            "associated_batches": samples.get(library_id, {}).get("extraction_batch", ""),
            "decision": decision,
            "decision_reason": (
                "weighted prevalence, paired-control abundance, and batch evidence >= 0.5"
                if decision == "exclude"
                else "weighted contamination evidence < 0.5"
            ),
            "contamination_model": "paired-controls-v2",
        }
        annotated.append(annotated_row)
        if decision == "retain":
            high.append(annotated_row)

    annotated_fields = PROFILE_FIELDS + [
        "native_count",
        "library_total_count",
        "native_normalized_value",
        "normalization_method",
        "abundance_rank",
        "clr_abundance",
        "control_prevalence",
        "biological_prevalence",
        "prevalence_contamination_score",
        "control_abundance_ratio",
        "abundance_contamination_score",
        "batch_localization_score",
        "batch_contamination_score",
        "field_contamination_score",
        "extraction_contamination_score",
        "pcr_contamination_score",
        "library_contamination_score",
        "run_contamination_score",
        "other_control_contamination_score",
        "probable_contamination_source",
        "contamination_probability",
        "associated_controls",
        "associated_control_types",
        "associated_batches",
        "decision",
        "decision_reason",
        "contamination_model",
    ]
    write_tsv(args.annotated_taxa, annotated, annotated_fields)
    write_tsv(args.high_confidence, high, annotated_fields)

    function_rows = [
        row
        for path in profile_files(args.functions, "functions.tsv")
        for row in read_tsv(path)
    ]
    write_tsv(
        args.unfiltered_functions,
        function_rows,
        ["library_id", "assay", "feature_id", "original_name", "abundance"],
    )
    function_fields = [
        "library_id",
        "assay",
        "feature_id",
        "original_name",
        "abundance",
        "decision",
        "decision_reason",
    ]
    functions = [
        {
            **row,
            "decision": "retain",
            "decision_reason": "functional potential is reported separately from taxon confidence",
        }
        for row in function_rows
    ]
    write_tsv(args.annotated_functions, functions, function_fields)

    sample_rows = []
    biological_groups: dict[str, list[tuple[str, dict[str, str]]]] = defaultdict(list)
    for library_id, sample in samples.items():
        if sample["sample_role"] == "biological":
            biological_groups[sample["biological_sample_id"]].append((library_id, sample))
    for sample_id, libraries in sorted(biological_groups.items()):
        assays = {sample["assay"] for _, sample in libraries}
        missing = [assay for assay in ("WMS", "16S", "ITS") if assay not in assays]
        unlinked = [library_id for library_id, _ in libraries if not linked.get(library_id)]
        sample_rows.append(
            {
                "biological_sample_id": sample_id,
                "wms_qc_status": "pass" if "WMS" in assays else "not_measured",
                "16s_qc_status": "pass" if "16S" in assays else "not_measured",
                "its_qc_status": "pass" if "ITS" in assays else "not_measured",
                "control_status": "linked" if not unlinked else "incomplete",
                "cross_assay_status": "complete" if not missing else "partial",
                "recommended_use": "include" if not unlinked else "review",
                "failure_reasons": (
                    "; ".join(
                        ([f"missing assays: {','.join(missing)}"] if missing else [])
                        + ([f"unlinked controls: {','.join(unlinked)}"] if unlinked else [])
                    )
                ),
            }
        )
    write_tsv(
        args.sample_confidence,
        sample_rows,
        [
            "biological_sample_id",
            "wms_qc_status",
            "16s_qc_status",
            "its_qc_status",
            "control_status",
            "cross_assay_status",
            "recommended_use",
            "failure_reasons",
        ],
    )


def cmd_harmonize(args: argparse.Namespace) -> None:
    rows = read_tsv(args.input)
    supplied = {}
    if args.crosswalk:
        for row in read_tsv(args.crosswalk):
            source_database = row.get("source_database", "")
            source_taxon_id = row.get(
                "source_taxon_id", row.get("original_taxon_id", "")
            )
            if source_database and source_taxon_id:
                supplied[(source_database, source_taxon_id)] = row
    generated = {}
    harmonized = []
    for row in rows:
        canonical_name = normalized_taxon_name(row["original_name"])
        key = canonical_name.casefold()
        canonical_id = "TRI-" + hashlib.sha1(key.encode()).hexdigest()[:12].upper()
        original_rank = row["original_rank"].casefold()
        mapping = {
            "triarchy_taxon_id": canonical_id,
            "canonical_name": canonical_name,
            "canonical_rank": original_rank,
            "canonical_lineage": canonical_name,
            "source_database": row["assay"],
            "source_taxon_id": row["feature_id"],
            "original_name": row["original_name"],
            "original_rank": original_rank,
            "mapping_method": "normalized_exact_name",
            "mapping_confidence": "0.5",
            "database_version": "unrecorded",
            "database_coverage": "",
            "reference_status": "unrecorded",
        }
        override = supplied.get((row["assay"], row["feature_id"]))
        if override:
            mapping.update(override)
            mapping["triarchy_taxon_id"] = override.get(
                "triarchy_taxon_id",
                override.get("canonical_taxon_id", mapping["triarchy_taxon_id"]),
            )
            mapping["canonical_rank"] = override.get(
                "canonical_rank",
                override.get("harmonized_rank", mapping["canonical_rank"]),
            )
            mapping["database_coverage"] = override.get(
                "database_coverage",
                override.get("mapping_confidence", mapping["database_coverage"]),
            )
            mapping["reference_status"] = override.get(
                "reference_status", mapping["reference_status"]
            )
        mapping["canonical_rank"] = mapping["canonical_rank"].casefold()
        if mapping["canonical_rank"] not in RANK_LEVEL:
            mapping["mapping_confidence"] = "0"
            mapping["mapping_method"] = "unsupported_rank"
        mapping["taxonomic_resolution_score"] = fmt(
            RANK_SUPPORT.get(mapping["canonical_rank"], 0.0)
        )
        generated[(mapping["source_database"], mapping["source_taxon_id"])] = mapping
        harmonized.append({**row, **mapping})
    crosswalk_fields = [
        "triarchy_taxon_id",
        "canonical_name",
        "canonical_rank",
        "canonical_lineage",
        "source_database",
        "source_taxon_id",
        "original_name",
        "original_rank",
        "mapping_method",
        "mapping_confidence",
        "database_version",
        "database_coverage",
        "reference_status",
        "taxonomic_resolution_score",
    ]
    write_tsv(args.output_crosswalk, generated.values(), crosswalk_fields)
    write_tsv(args.output, harmonized, list(harmonized[0]) if harmonized else PROFILE_FIELDS + crosswalk_fields)


def cmd_classify(args: argparse.Namespace) -> None:
    """Run the formal assay-aware evidence-v3 inference model."""
    run_inference(
        harmonized=args.harmonized,
        samples=args.samples,
        output=args.output,
        model_output=getattr(args, "model_output", None),
        performance_path=getattr(args, "assay_performance", None),
        policy_path=getattr(args, "evidence_policy", None),
        qc_metrics_path=getattr(args, "qc_metrics", None),
    )


def cmd_fit_calibration(args: argparse.Namespace) -> None:
    fit_calibration(args.evidence, args.truth, args.output)


def cmd_evaluate_evidence(args: argparse.Namespace) -> None:
    evaluate_evidence(args.evidence, args.truth, args.output)


def pearson(values_a: list[float], values_b: list[float]) -> float:
    if len(values_a) < 2 or len(values_a) != len(values_b):
        return 0.0
    mean_a = sum(values_a) / len(values_a)
    mean_b = sum(values_b) / len(values_b)
    numerator = sum(
        (value_a - mean_a) * (value_b - mean_b)
        for value_a, value_b in zip(values_a, values_b)
    )
    denominator = math.sqrt(
        sum((value - mean_a) ** 2 for value in values_a)
        * sum((value - mean_b) ** 2 for value in values_b)
    )
    return numerator / denominator if denominator else 0.0


def cmd_concordance(args: argparse.Namespace) -> None:
    samples = {row["library_id"]: row for row in read_tsv(args.samples)}
    profiles: dict[tuple[str, str, str], dict[str, float]] = defaultdict(dict)
    for row in read_tsv(args.harmonized):
        sample = samples.get(row["library_id"])
        if not sample or sample["sample_role"] != "biological":
            continue
        kingdom = infer_kingdom(row)
        key = (sample["biological_sample_id"], row["assay"], kingdom)
        taxon = row["triarchy_taxon_id"]
        profiles[key][taxon] = profiles[key].get(taxon, 0.0) + float(
            row.get("native_normalized_value") or 0
        )
    output = []
    sample_ids = sorted({key[0] for key in profiles})
    comparisons = (
        ("Bacteria", "WMS", "16S"),
        ("Archaea", "WMS", "16S"),
        ("Fungi", "WMS", "ITS"),
    )
    for sample_id in sample_ids:
        for kingdom, assay_a, assay_b in comparisons:
            profile_a = profiles.get((sample_id, assay_a, kingdom), {})
            profile_b = profiles.get((sample_id, assay_b, kingdom), {})
            taxa_a, taxa_b = set(profile_a), set(profile_b)
            union = taxa_a | taxa_b
            shared = sorted(taxa_a & taxa_b)
            ranks_a = {
                taxon: rank
                for rank, taxon in enumerate(
                    sorted(profile_a, key=profile_a.get, reverse=True), start=1
                )
            }
            ranks_b = {
                taxon: rank
                for rank, taxon in enumerate(
                    sorted(profile_b, key=profile_b.get, reverse=True), start=1
                )
            }
            spearman = pearson(
                [float(ranks_a[taxon]) for taxon in shared],
                [float(ranks_b[taxon]) for taxon in shared],
            )
            rank_agreement = (
                sum(1 / (1 + abs(ranks_a[taxon] - ranks_b[taxon])) for taxon in shared)
                / len(shared)
                if shared
                else 0.0
            )
            output.append(
                {
                    "biological_sample_id": sample_id,
                    "kingdom": kingdom,
                    "assay_pair": f"{assay_a}_vs_{assay_b}",
                    "assay_a_taxa": len(taxa_a),
                    "assay_b_taxa": len(taxa_b),
                    "shared_taxa": len(shared),
                    "union_taxa": len(union),
                    "jaccard_presence": fmt(len(shared) / len(union) if union else 0),
                    "spearman_rank_correlation": fmt(spearman),
                    "rank_agreement_score": fmt(rank_agreement),
                    "scope": "detection_and_rank_only",
                }
            )
    write_tsv(
        args.output,
        output,
        [
            "biological_sample_id", "kingdom", "assay_pair", "assay_a_taxa",
            "assay_b_taxa", "shared_taxa", "union_taxa", "jaccard_presence",
            "spearman_rank_correlation", "rank_agreement_score", "scope",
        ],
    )


def cmd_diversity(args: argparse.Namespace) -> None:
    taxa = read_tsv(args.taxa)
    functions = read_tsv(args.functions)
    grouped: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    for row in taxa:
        grouped[(row["library_id"], row["assay"])][row["original_name"]] = float(
            row["abundance"]
        )
    for row in functions:
        grouped[(row["library_id"], "WMS_function")][row["feature_id"]] = float(
            row["abundance"]
        )
    alpha = []
    for (library_id, assay), profile in sorted(grouped.items()):
        values = list(profile.values())
        total = sum(values)
        proportions = [value / total for value in values if value > 0 and total > 0]
        shannon = -sum(value * math.log(value) for value in proportions)
        alpha.append(
            {
                "library_id": library_id,
                "assay": assay,
                "observed_features": len(values),
                "shannon": f"{shannon:.6g}",
                "scope": "within_assay",
            }
        )
    write_tsv(
        args.alpha_output,
        alpha,
        ["library_id", "assay", "observed_features", "shannon", "scope"],
    )
    beta = []
    assays = sorted({assay for _, assay in grouped})
    for assay in assays:
        libraries = sorted(
            library_id for library_id, candidate_assay in grouped if candidate_assay == assay
        )
        for index, library_a in enumerate(libraries):
            for library_b in libraries[index + 1 :]:
                profile_a = grouped[(library_a, assay)]
                profile_b = grouped[(library_b, assay)]
                features = set(profile_a) | set(profile_b)
                numerator = sum(
                    abs(profile_a.get(feature, 0.0) - profile_b.get(feature, 0.0))
                    for feature in features
                )
                denominator = sum(
                    profile_a.get(feature, 0.0) + profile_b.get(feature, 0.0)
                    for feature in features
                )
                beta.append(
                    {
                        "library_a": library_a,
                        "library_b": library_b,
                        "assay": assay,
                        "bray_curtis": f"{numerator / denominator if denominator else 0:.6g}",
                        "scope": "within_assay",
                    }
                )
    write_tsv(
        args.beta_output,
        beta,
        ["library_a", "library_b", "assay", "bray_curtis", "scope"],
    )


def cmd_provenance(args: argparse.Namespace) -> None:
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(args.project_dir), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    document = {
        "workflow": "TRIARCHY",
        "workflow_version": os.environ.get("NXF_VER", "unknown"),
        "git_commit": commit,
        "samples_sha256": checksum(args.samples),
        "controls_sha256": checksum(args.controls),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "run_mode": "production",
    }
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def cmd_report(args: argparse.Namespace) -> None:
    provenance = json.loads(args.provenance.read_text(encoding="utf-8"))
    evidence_model = json.loads(args.evidence_model.read_text(encoding="utf-8"))
    payload = {
        "evidence": read_tsv(args.evidence),
        "concordance": read_tsv(args.concordance),
        "alpha": read_tsv(args.alpha_diversity),
        "beta": read_tsv(args.beta_diversity),
        "confidence": read_tsv(args.sample_confidence),
        "qc": read_tsv(args.qc_metrics),
        "samples": read_tsv(args.samples),
        "controls": read_tsv(args.controls),
        "functions": read_tsv(args.functions),
        "model": evidence_model,
        "provenance": provenance,
    }
    serialized = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    document = DASHBOARD_TEMPLATE.replace("__TRIARCHY_DATA__", serialized)
    args.output.write_text(document, encoding="utf-8")


DASHBOARD_TEMPLATE = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>TRIARCHY · Scientific evidence dashboard</title>
<style>
:root{--navy:#082b58;--blue:#125a91;--teal:#078c97;--green:#58a62a;--ink:#17212b;--muted:#607080;--line:#dce5ea;--wash:#f5f8fa;--warn:#a65b00;--bad:#b42318;--good:#16794b;--white:#fff;--shadow:0 8px 24px rgba(8,43,88,.08)}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;font:15px/1.5 Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:var(--ink);background:var(--wash)}button,input,select{font:inherit}button{cursor:pointer}
.skip{position:absolute;left:-999px}.skip:focus{left:1rem;top:1rem;background:#fff;padding:.7rem;z-index:20}.app{display:grid;grid-template-columns:270px minmax(0,1fr);min-height:100vh}.sidebar{position:sticky;top:0;height:100vh;background:linear-gradient(165deg,#06244a,#073667 60%,#07516a);color:#fff;padding:1.35rem 1rem;overflow:auto}.brand{display:flex;align-items:center;gap:.75rem;padding:.3rem .55rem 1.2rem}.mark{width:42px;height:42px;border:2px solid #72d2c7;border-radius:50%;display:grid;place-items:center;font-weight:800;color:#b9ef85}.brand strong{font-size:1.18rem;letter-spacing:.14em}.brand small{display:block;opacity:.7;letter-spacing:0}
.nav-label{font-size:.68rem;text-transform:uppercase;letter-spacing:.13em;opacity:.6;margin:1rem .6rem .45rem}.nav button{width:100%;border:0;background:transparent;color:#d7e7f3;border-radius:10px;padding:.72rem .75rem;margin:.12rem 0;text-align:left;display:grid;grid-template-columns:26px 1fr;gap:.35rem;align-items:start}.nav button:hover,.nav button.active{background:rgba(255,255,255,.13);color:#fff}.nav b{font-size:.9rem}.nav small{display:block;font-size:.7rem;opacity:.66;margin-top:.1rem}.side-foot{margin:1.4rem .55rem;font-size:.72rem;color:#b7cad8}.main{min-width:0}.topbar{height:66px;background:rgba(255,255,255,.94);border-bottom:1px solid var(--line);display:flex;align-items:center;justify-content:space-between;padding:0 clamp(1rem,3vw,2.5rem);position:sticky;top:0;z-index:10;backdrop-filter:blur(10px)}.breadcrumb{font-weight:650;color:var(--navy)}.top-actions{display:flex;gap:.6rem}.btn{border:1px solid var(--line);background:#fff;border-radius:9px;padding:.55rem .8rem;color:var(--navy);font-weight:650}.btn.primary{background:var(--navy);color:#fff;border-color:var(--navy)}.content{padding:clamp(1.2rem,3vw,2.7rem);max-width:1500px;margin:auto}.page{display:none}.page.active{display:block}.eyebrow{text-transform:uppercase;letter-spacing:.12em;color:var(--teal);font-size:.74rem;font-weight:800}.hero{display:flex;justify-content:space-between;gap:2rem;align-items:end;margin-bottom:1.5rem}.hero h1{font-size:clamp(1.9rem,4vw,3.25rem);line-height:1.08;letter-spacing:-.035em;margin:.3rem 0;color:var(--navy)}.hero p{max-width:760px;color:var(--muted);font-size:1.04rem;margin:.4rem 0}.status{padding:.45rem .7rem;border-radius:99px;background:#e8f6ee;color:var(--good);font-weight:750;white-space:nowrap}
.kpis{display:grid;grid-template-columns:repeat(4,minmax(140px,1fr));gap:1rem;margin:1.2rem 0}.kpi,.card{background:#fff;border:1px solid var(--line);border-radius:15px;box-shadow:var(--shadow)}.kpi{padding:1rem 1.1rem;border-top:3px solid var(--teal)}.kpi .value{font-size:1.9rem;line-height:1.1;font-weight:780;color:var(--navy)}.kpi .label{font-size:.78rem;color:var(--muted);margin-top:.3rem}.grid{display:grid;grid-template-columns:repeat(12,1fr);gap:1rem;margin:1rem 0}.span-4{grid-column:span 4}.span-5{grid-column:span 5}.span-6{grid-column:span 6}.span-7{grid-column:span 7}.span-8{grid-column:span 8}.span-12{grid-column:span 12}.card{padding:1.15rem;min-width:0}.card h2,.card h3{color:var(--navy);margin:0 0 .35rem}.card h2{font-size:1.12rem}.card h3{font-size:.98rem}.sub{color:var(--muted);font-size:.82rem;margin:0 0 1rem}.attention{border-left:4px solid #ef9f27}.attention ul{padding-left:1.2rem;margin:.7rem 0}.attention li{margin:.35rem 0}.callout{padding:.9rem 1rem;background:#eef6f7;border-left:4px solid var(--teal);border-radius:5px;margin:1rem 0}.callout.warn{background:#fff7e8;border-color:#e29219}.question{font-size:1.08rem;color:var(--navy);font-weight:680;margin-bottom:1rem}.chart{min-height:210px;display:flex;flex-direction:column;justify-content:center}.bars{display:grid;gap:.55rem}.bar-row{display:grid;grid-template-columns:minmax(100px,160px) 1fr 42px;gap:.6rem;align-items:center;font-size:.8rem}.track{height:12px;background:#e7eef2;border-radius:9px;overflow:hidden}.fill{height:100%;border-radius:9px;background:linear-gradient(90deg,var(--teal),var(--green))}.legend{display:flex;gap:1rem;flex-wrap:wrap;font-size:.76rem;color:var(--muted)}.dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:.3rem}.matrix{overflow:auto}.matrix table{min-width:680px}.matrix td.center,.matrix th.center{text-align:center}.pill{display:inline-flex;padding:.2rem .48rem;border-radius:99px;background:#eef2f5;color:#405363;font-size:.72rem;font-weight:700}.pill.good{background:#e5f4ec;color:var(--good)}.pill.warn{background:#fff0d8;color:#925400}.pill.bad{background:#fee9e7;color:var(--bad)}.pill.info{background:#e6f3f6;color:#076d78}
.tools{display:flex;gap:.6rem;flex-wrap:wrap;align-items:center;margin:.8rem 0}.search,.select{border:1px solid #cdd9df;border-radius:9px;padding:.55rem .7rem;background:#fff;min-width:190px}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:11px;max-height:520px}table{border-collapse:separate;border-spacing:0;width:100%;font-size:.78rem;background:#fff}th,td{padding:.62rem .7rem;text-align:left;border-bottom:1px solid #e8eef1;white-space:nowrap}th{position:sticky;top:0;background:#edf4f5;color:#21445d;z-index:1;font-size:.7rem;text-transform:uppercase;letter-spacing:.035em}tbody tr:hover{background:#f3f9f8}.empty{text-align:center;padding:3rem 1rem;color:var(--muted)}.empty strong{display:block;color:var(--navy);font-size:1rem;margin-bottom:.25rem}.detail{background:linear-gradient(135deg,#f0f8f8,#f6f8fb);border:1px solid #d7e8e8;border-radius:12px;padding:1rem}.detail-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:.7rem}.metric{background:#fff;border-radius:9px;padding:.75rem;border:1px solid var(--line)}.metric b{display:block;color:var(--navy);font-size:1.1rem}.steps{display:grid;grid-template-columns:repeat(5,1fr);gap:.5rem;margin:1rem 0}.step{position:relative;background:#fff;border:1px solid var(--line);padding:.8rem;border-radius:10px;text-align:center;color:var(--navy);font-weight:700;font-size:.8rem}.step:not(:last-child):after{content:'›';position:absolute;right:-.48rem;color:var(--teal);z-index:2}.provenance{white-space:pre-wrap;word-break:break-word;background:#102b43;color:#d7ecf3;border-radius:10px;padding:1rem;font:12px/1.6 ui-monospace,SFMono-Regular,Menlo,monospace}.mobile-menu{display:none}
@media(max-width:980px){.app{grid-template-columns:1fr}.sidebar{position:fixed;z-index:30;left:-285px;width:270px;transition:.2s}.sidebar.open{left:0}.mobile-menu{display:inline-block}.kpis{grid-template-columns:repeat(2,1fr)}.span-4,.span-5,.span-6,.span-7,.span-8{grid-column:span 12}.steps{grid-template-columns:1fr 1fr}.hero{align-items:start;flex-direction:column;gap:.5rem}}@media(max-width:540px){.kpis{grid-template-columns:1fr 1fr}.top-actions .btn:not(.primary){display:none}.detail-grid{grid-template-columns:1fr}.content{padding:1rem}}
@media print{.sidebar,.topbar,.tools{display:none}.app{display:block}.page{display:block;break-before:page}.content{max-width:none}.card,.kpi{box-shadow:none}}
</style></head><body><a class="skip" href="#content">Skip to dashboard</a>
<div class="app"><aside class="sidebar" id="sidebar"><div class="brand"><div class="mark">△</div><div><strong>TRIARCHY</strong><small>Scientific evidence dashboard</small></div></div><div class="nav-label">Scientific journey</div><nav class="nav" id="nav">
<button class="active" data-page="overview"><span>01</span><span><b>Study overview</b><small>Did the experiment work?</small></span></button>
<button data-page="qc"><span>02</span><span><b>Sequencing QC</b><small>Which libraries are trustworthy?</small></span></button>
<button data-page="controls"><span>03</span><span><b>Controls</b><small>Is this contamination?</small></span></button>
<button data-page="taxonomy"><span>04</span><span><b>Taxonomic integration</b><small>What did each assay see?</small></span></button>
<button data-page="evidence"><span>05</span><span><b>Evidence & discordance</b><small>Why do assays disagree?</small></span></button>
<button data-page="ecology"><span>06</span><span><b>Ecology</b><small>Which conclusions survive?</small></span></button>
<button data-page="functions"><span>07</span><span><b>Functional biology</b><small>What functions were measured?</small></span></button>
<button data-page="provenance"><span>08</span><span><b>Export & provenance</b><small>Can I reproduce this?</small></span></button></nav><div class="side-foot">Assay-native measurements remain separate. Integration summarizes evidence, not interchangeable abundance.</div></aside>
<main class="main"><header class="topbar"><div><button class="btn mobile-menu" id="menu">☰</button> <span class="breadcrumb" id="breadcrumb">Study overview</span></div><div class="top-actions"><button class="btn" onclick="window.print()">Print / PDF</button><button class="btn primary" onclick="downloadAll()">Export data</button></div></header><div class="content" id="content">
<section class="page active" id="overview"><div class="hero"><div><div class="eyebrow">Run at a glance</div><h1>Did the experiment work?</h1><p>A decision-first summary of assay completeness, sample trustworthiness, and evidence requiring attention.</p></div><span class="status" id="run-status">Run assessed</span></div><div class="steps"><div class="step">Experiment</div><div class="step">Samples</div><div class="step">Taxa</div><div class="step">Discordance</div><div class="step">Conclusions</div></div><div class="kpis" id="overview-kpis"></div><div class="grid"><div class="card span-8"><h2>Sample × assay completeness</h2><p class="sub">A fast check that each biological sample has enough information for integration.</p><div class="matrix" id="completeness"></div></div><div class="card span-4 attention"><h2>Attention required</h2><p class="sub">Items worth reviewing before interpretation.</p><div id="attention"></div></div><div class="card span-12"><h2>Evidence landscape</h2><p class="sub">How TRIARCHY classified taxon-level evidence across all biological samples.</p><div class="chart" id="overview-evidence"></div></div></div></section>
<section class="page" id="qc"><div class="hero"><div><div class="eyebrow">Sequencing quality</div><h1>Which libraries are trustworthy?</h1><p>Compare retention, host burden, and QC decisions without losing the assay-specific context.</p></div></div><div class="tools"><select class="select" id="qc-assay"><option value="">All assays</option><option>WMS</option><option>16S</option><option>ITS</option></select><input class="search" id="qc-search" placeholder="Search libraries…"><button class="btn" onclick="downloadRows('qc')">Download QC TSV</button></div><div class="grid"><div class="card span-5"><h2>QC status by assay</h2><div class="chart" id="qc-chart"></div></div><div class="card span-7"><h2>Selected library</h2><p class="sub">Click a row below to inspect its key measurements.</p><div class="detail" id="qc-detail"></div></div><div class="card span-12"><div class="table-wrap" id="qc-table"></div></div></div></section>
<section class="page" id="controls"><div class="hero"><div><div class="eyebrow">Control-aware evidence</div><h1>Signal or contamination?</h1><p>Trace every biological library to its relevant field, extraction, PCR, and library controls.</p></div></div><div class="grid"><div class="card span-5"><h2>Sample–control links</h2><p class="sub">Select a library to reveal only the controls that matter to it.</p><select class="select" id="control-library"></select><div id="control-graph" class="chart"></div></div><div class="card span-7"><h2>Probable contamination sources</h2><div class="chart" id="contam-chart"></div></div><div class="card span-12"><h2>Contamination evidence</h2><div class="tools"><input class="search" id="control-search" placeholder="Search sample or taxon…"><button class="btn" onclick="downloadRows('evidence','contamination-evidence.tsv')">Download evidence</button></div><div class="table-wrap" id="control-table"></div></div></div></section>
<section class="page" id="taxonomy"><div class="hero"><div><div class="eyebrow">Cross-assay matrix</div><h1>What did each assay see?</h1><p>Detection is shown in assay scope; not-applicable evidence is never treated as absence.</p></div></div><div class="tools"><input class="search" id="tax-search" placeholder="Search taxon or sample…"><select class="select" id="tax-kingdom"><option value="">All kingdoms</option></select><button class="btn" onclick="downloadRows('evidence','taxonomic-integration.tsv')">Download matrix</button></div><div class="card"><div class="table-wrap" id="taxonomy-table"></div></div></section>
<section class="page" id="evidence"><div class="hero"><div><div class="eyebrow">TRIARCHY core</div><h1>Why do the assays disagree?</h1><p>Explore authenticity, contamination, detection power, reference coverage, and the decision reason together.</p></div></div><div class="grid"><div class="card span-5"><h2>Evidence-class distribution</h2><div class="chart" id="evidence-chart"></div></div><div class="card span-7"><h2>Taxon detail</h2><select class="select" id="evidence-taxon"></select><div class="detail" id="evidence-detail" style="margin-top:1rem"></div></div><div class="card span-12"><div class="tools"><input class="search" id="evidence-search" placeholder="Search evidence…"><select class="select" id="evidence-class"><option value="">All classes</option></select><button class="btn" onclick="downloadRows('evidence')">Download evidence TSV</button></div><div class="table-wrap" id="evidence-table"></div></div></div></section>
<section class="page" id="ecology"><div class="hero"><div><div class="eyebrow">Assay-native ecology</div><h1>Which conclusions survive?</h1><p>Diversity and concordance remain explicitly within assay scope. Cross-assay summaries compare detection and rank, never raw abundance.</p></div></div><div class="callout">Interpret each assay on its own measurement scale. Cross-assay concordance below reports detection overlap and rank agreement.</div><div class="grid"><div class="card span-6"><h2>Cross-assay presence agreement</h2><div class="chart" id="concordance-chart"></div></div><div class="card span-6"><h2>Alpha diversity</h2><div class="chart" id="alpha-chart"></div></div><div class="card span-12"><h2>Ecology tables</h2><div class="tools"><button class="btn" onclick="downloadRows('concordance')">Concordance TSV</button><button class="btn" onclick="downloadRows('alpha')">Alpha TSV</button><button class="btn" onclick="downloadRows('beta')">Beta TSV</button></div><div class="table-wrap" id="ecology-table"></div></div></div></section>
<section class="page" id="functions"><div class="hero"><div><div class="eyebrow">WMS functional evidence</div><h1>What functions were measured?</h1><p>Functional calls come from WMS. Complementary assays support organisms associated with functions—not the pathways themselves.</p></div></div><div class="callout warn">Do not describe 16S or ITS evidence as direct pathway validation.</div><div class="card"><div class="tools"><input class="search" id="function-search" placeholder="Search pathway or sample…"><button class="btn" onclick="downloadRows('functions')">Download functions TSV</button></div><div class="table-wrap" id="function-table"></div></div></section>
<section class="page" id="provenance"><div class="hero"><div><div class="eyebrow">Publication & reproducibility</div><h1>Export defensible evidence</h1><p>Download the tables behind every view, print publication-ready pages, and retain the exact run identity.</p></div></div><div class="grid"><div class="card span-5"><h2>Export center</h2><p class="sub">All exports contain the unrounded values used by the dashboard.</p><div class="tools"><button class="btn primary" onclick="downloadAll()">Download all tables</button><button class="btn" onclick="window.print()">Print report / PDF</button></div><div class="callout">Suggested citation metadata and model calibration status are included in the provenance bundle.</div></div><div class="card span-7"><h2>Run identity</h2><div class="detail" id="run-identity"></div></div><div class="card span-6"><h2>Run provenance</h2><pre class="provenance" id="provenance-json"></pre></div><div class="card span-6"><h2>Evidence model</h2><pre class="provenance" id="model-json"></pre></div></div></section>
</div></main></div><script id="triarchy-data" type="application/json">__TRIARCHY_DATA__</script>
<script>
const D=JSON.parse(document.getElementById('triarchy-data').textContent);const $=s=>document.querySelector(s),$$=s=>[...document.querySelectorAll(s)];
const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const num=v=>{const n=Number(v);return Number.isFinite(n)?n:0};const yes=v=>String(v).toLowerCase()==='true';
const title=s=>String(s||'').replaceAll('_',' ').replace(/\b\w/g,c=>c.toUpperCase());
const cls=s=>{s=String(s||'').toLowerCase();return /pass|complete|include|concordant|authentic/.test(s)?'good':/fail|contamin|discord/.test(s)?'bad':/limit|partial|restrict|indeterminate/.test(s)?'warn':'info'};
function pill(v){return `<span class="pill ${cls(v)}">${esc(title(v||'not available'))}</span>`}
function kpi(value,label){return `<div class="kpi"><div class="value">${esc(value)}</div><div class="label">${esc(label)}</div></div>`}
function bars(items){if(!items.length)return empty('No values available','This run did not emit data for this view.');const max=Math.max(...items.map(x=>x[1]),1);return `<div class="bars">${items.map(([l,v,c])=>`<div class="bar-row"><span title="${esc(l)}">${esc(l)}</span><div class="track"><div class="fill" style="width:${100*v/max}%;${c?`background:${c}`:''}"></div></div><b>${esc(v)}</b></div>`).join('')}</div>`}
function empty(a,b){return `<div class="empty"><strong>${esc(a)}</strong>${esc(b)}</div>`}
function table(rows,cols,id){if(!rows.length)return empty('No matching rows','Adjust the filters or inspect the upstream output.');cols=cols||Object.keys(rows[0]);return `<table data-table="${id||''}"><thead><tr>${cols.map(c=>`<th>${esc(title(c))}</th>`).join('')}</tr></thead><tbody>${rows.map((r,i)=>`<tr data-index="${i}">${cols.map(c=>`<td>${/status|classification|state/.test(c)?pill(r[c]):esc(r[c])}</td>`).join('')}</tr>`).join('')}</tbody></table>`}
function download(name,rows){if(!rows?.length)return;const cols=Object.keys(rows[0]);const clean=v=>String(v??'').replaceAll('\t',' ').replaceAll('\n',' ');const t=[cols.join('\t'),...rows.map(r=>cols.map(c=>clean(r[c])).join('\t'))].join('\n');const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([t],{type:'text/tab-separated-values'}));a.download=name;a.click();URL.revokeObjectURL(a.href)}
function downloadRows(key,name){download(name||`${key}.tsv`,D[key])}function downloadAll(){Object.entries(D).filter(([,v])=>Array.isArray(v)).forEach(([k,v],i)=>setTimeout(()=>download(`triarchy-${k}.tsv`,v),i*120))}
const biological=D.samples.filter(r=>r.sample_role==='biological');const bioIds=[...new Set(biological.map(r=>r.biological_sample_id))];const controls=D.samples.filter(r=>r.sample_role!=='biological');const assays=['WMS','16S','ITS'];
function counts(rows,key){return rows.reduce((a,r)=>(a[r[key]||'not_available']=(a[r[key]||'not_available']||0)+1,a),{})}
function nav(){ $$('#nav button').forEach(b=>b.onclick=()=>{$$('#nav button').forEach(x=>x.classList.remove('active'));$$('.page').forEach(x=>x.classList.remove('active'));b.classList.add('active');$('#'+b.dataset.page).classList.add('active');$('#breadcrumb').textContent=b.querySelector('b').textContent;$('#sidebar').classList.remove('open');scrollTo(0,0)})}
function overview(){const evc=counts(D.evidence,'classification');const failed=D.confidence.filter(r=>r.recommended_use!=='include');const contam=D.evidence.filter(r=>num(r.contamination_probability)>=.5||/contamin/i.test(r.classification));const discord=D.evidence.filter(r=>/discord/i.test(r.classification));const high=D.evidence.filter(r=>num(r.biological_authenticity)>=.8&&num(r.contamination_probability)<.5);$('#overview-kpis').innerHTML=[kpi(bioIds.length,'Biological samples'),...assays.map(a=>kpi(biological.filter(r=>r.assay===a).length,`${a} libraries`)),kpi(controls.length,'Control libraries'),kpi(D.confidence.filter(r=>r.recommended_use==='include').length,'Samples passing'),kpi(high.length,'High-confidence taxa'),kpi(contam.length,'Probable contaminants')].join('');const byId=Object.fromEntries(D.confidence.map(r=>[r.biological_sample_id,r]));$('#completeness').innerHTML=table(bioIds.map(id=>{const c=byId[id]||{};return{sample:id,WMS:c.wms_qc_status||'missing','16S':c['16s_qc_status']||'missing',ITS:c.its_qc_status||'missing',controls:c.control_status||'unknown',integration:c.cross_assay_status||'unknown',use:c.recommended_use||'review'}}),['sample','WMS','16S','ITS','controls','integration','use']);const issues=[];if(failed.length)issues.push(`${failed.length} sample(s) excluded or restricted`);const missing=D.confidence.filter(r=>['wms_qc_status','16s_qc_status','its_qc_status'].some(k=>!r[k]||r[k]==='missing'));if(missing.length)issues.push(`${missing.length} sample(s) have missing assay evidence`);if(contam.length)issues.push(`${contam.length} taxon calls show contamination concern`);if(discord.length)issues.push(`${discord.length} unresolved discordance(s)`);issues.push(`${D.evidence.filter(r=>/indeterminate/.test(r.classification)).length} indeterminate evidence call(s)`);$('#attention').innerHTML=issues.some(x=>!x.startsWith('0'))?`<ul>${issues.filter(x=>!x.startsWith('0')).map(x=>`<li>${esc(x)}</li>`).join('')}</ul>`:empty('No critical flags','All emitted checks passed.');$('#overview-evidence').innerHTML=bars(Object.entries(evc).sort((a,b)=>b[1]-a[1]));}
function qc(){const render=()=>{const q=$('#qc-search').value.toLowerCase(),a=$('#qc-assay').value;const rows=D.qc.filter(r=>(!a||r.assay===a)&&Object.values(r).join(' ').toLowerCase().includes(q));$('#qc-table').innerHTML=table(rows,['library_id','assay','input_reads','retained_reads','q30_rate','host_fraction','quality_factor','qc_status'],'qc');$('#qc-chart').innerHTML=bars(Object.entries(counts(rows,'assay')).map(([k,v])=>[k,v]));const trs=$$('#qc-table tbody tr');trs.forEach((tr,i)=>tr.onclick=()=>qcDetail(rows[i]));qcDetail(rows[0])};$('#qc-search').oninput=render;$('#qc-assay').onchange=render;render()}
function qcDetail(r){$('#qc-detail').innerHTML=r?`<h3>${esc(r.library_id)}</h3><div class="detail-grid"><div class="metric"><span>Status</span><b>${pill(r.qc_status)}</b></div><div class="metric"><span>Retention</span><b>${r.input_reads&&num(r.input_reads)>0?Math.round(100*num(r.retained_reads)/num(r.input_reads))+'%':'Not reported'}</b></div><div class="metric"><span>Host fraction</span><b>${r.host_fraction||'—'}</b></div></div><p class="sub">Sources: ${esc(r.qc_sources||'not reported')}</p>`:empty('Select a library','Choose a row to see its metrics')}
function controlPage(){const libs=[...new Set(D.controls.map(r=>r.library_id))];$('#control-library').innerHTML=libs.map(x=>`<option>${esc(x)}</option>`).join('');const graph=()=>{const id=$('#control-library').value,links=D.controls.filter(r=>r.library_id===id);$('#control-graph').innerHTML=links.length?`<div class="detail"><h3>${esc(id)}</h3>${links.map(x=>`<p>↳ ${pill(x.control_stage)} <b>${esc(x.control_library_id)}</b><br><small>${esc(x.reason)}</small></p>`).join('')}</div>`:empty('No linked controls','No edges were emitted for this library.')};$('#control-library').onchange=graph;graph();const sources=counts(D.evidence.filter(r=>r.probable_contamination_source&&r.probable_contamination_source!=='none'),'probable_contamination_source');$('#contam-chart').innerHTML=Object.keys(sources).length?bars(Object.entries(sources)):empty('No localized contamination','No probable source was assigned in this run.');const render=()=>{const q=$('#control-search').value.toLowerCase();const rows=D.evidence.filter(r=>num(r.contamination_probability)>0||/contamin/i.test(r.classification)).filter(r=>Object.values(r).join(' ').toLowerCase().includes(q));$('#control-table').innerHTML=table(rows,['biological_sample_id','canonical_name','field_contamination_score','extraction_contamination_score','pcr_contamination_score','library_contamination_score','contamination_probability','probable_contamination_source','classification','decision_reason'])};$('#control-search').oninput=render;render()}
function taxonomy(){const kingdoms=[...new Set(D.evidence.map(r=>r.kingdom))].filter(Boolean);$('#tax-kingdom').innerHTML+=[...kingdoms].map(x=>`<option>${esc(x)}</option>`).join('');const render=()=>{const q=$('#tax-search').value.toLowerCase(),k=$('#tax-kingdom').value;const rows=D.evidence.filter(r=>(!k||r.kingdom===k)&&Object.values(r).join(' ').toLowerCase().includes(q)).map(r=>({sample:r.biological_sample_id,taxon:r.canonical_name,kingdom:r.kingdom,rank:r.comparison_rank,WMS:yes(r.wms_measured)?(yes(r.wms_detected)?'detected':'not detected'):'N/A','16S':yes(r['16s_measured'])?(yes(r['16s_detected'])?'detected':'not detected'):'N/A',ITS:yes(r.its_measured)?(yes(r.its_detected)?'detected':'not detected'):'N/A',integration:r.classification,reason:r.decision_reason}));$('#taxonomy-table').innerHTML=table(rows,['sample','taxon','kingdom','rank','WMS','16S','ITS','integration','reason'])};$('#tax-search').oninput=render;$('#tax-kingdom').onchange=render;render()}
function evidence(){const classes=[...new Set(D.evidence.map(r=>r.classification))].filter(Boolean);$('#evidence-class').innerHTML+=classes.map(x=>`<option value="${esc(x)}">${esc(title(x))}</option>`).join('');const taxa=[...new Set(D.evidence.map(r=>r.canonical_name))].filter(Boolean);$('#evidence-taxon').innerHTML=taxa.map(x=>`<option>${esc(x)}</option>`).join('');$('#evidence-chart').innerHTML=bars(Object.entries(counts(D.evidence,'classification')).sort((a,b)=>b[1]-a[1]));const detail=()=>{const tax=$('#evidence-taxon').value,rows=D.evidence.filter(r=>r.canonical_name===tax);if(!rows.length)return;const auth=rows.reduce((a,r)=>a+num(r.biological_authenticity),0)/rows.length,cont=rows.reduce((a,r)=>a+num(r.contamination_probability),0)/rows.length;$('#evidence-detail').innerHTML=`<h3><i>${esc(tax)}</i></h3><div class="detail-grid"><div class="metric"><span>Mean authenticity</span><b>${auth.toFixed(2)}</b></div><div class="metric"><span>Mean contamination</span><b>${cont.toFixed(2)}</b></div><div class="metric"><span>Samples</span><b>${rows.length}</b></div></div><p><b>Why TRIARCHY reached this result</b><br>${esc(rows[0].decision_reason||'No explanation emitted.')}</p>`};$('#evidence-taxon').onchange=detail;detail();const render=()=>{const q=$('#evidence-search').value.toLowerCase(),c=$('#evidence-class').value;const rows=D.evidence.filter(r=>(!c||r.classification===c)&&Object.values(r).join(' ').toLowerCase().includes(q));$('#evidence-table').innerHTML=table(rows,['biological_sample_id','canonical_name','kingdom','assay_support','biological_authenticity','authenticity_ci_low','authenticity_ci_high','contamination_probability','classification','discordance_reason','decision_reason'])};$('#evidence-search').oninput=render;$('#evidence-class').onchange=render;render()}
function ecology(){const pairs={};D.concordance.forEach(r=>{const k=r.assay_pair||'comparison';(pairs[k]??=[]).push(num(r.jaccard_presence))});$('#concordance-chart').innerHTML=bars(Object.entries(pairs).map(([k,v])=>[k,Math.round(100*v.reduce((a,b)=>a+b,0)/v.length)]));const alphaBy={};D.alpha.forEach(r=>{const k=r.assay||'unknown';(alphaBy[k]??=[]).push(num(r.shannon_diversity||r.shannon||0))});$('#alpha-chart').innerHTML=Object.keys(alphaBy).length?bars(Object.entries(alphaBy).map(([k,v])=>[k,(v.reduce((a,b)=>a+b,0)/v.length).toFixed(2)])):empty('No alpha-diversity rows','The selected run emitted an empty high-confidence diversity table.');$('#ecology-table').innerHTML=table(D.concordance,['biological_sample_id','kingdom','assay_pair','shared_taxa','union_taxa','jaccard_presence','spearman_rank_correlation','rank_agreement_score'])}
function functions(){const render=()=>{const q=$('#function-search').value.toLowerCase();$('#function-table').innerHTML=table(D.functions.filter(r=>Object.values(r).join(' ').toLowerCase().includes(q)))};$('#function-search').oninput=render;render()}
function provenance(){const p=D.provenance;$('#run-identity').innerHTML=`<div class="detail-grid"><div class="metric"><span>Workflow</span><b>${esc(p.workflow||'TRIARCHY')}</b></div><div class="metric"><span>Nextflow</span><b>${esc(p.workflow_version||'unknown')}</b></div><div class="metric"><span>Run mode</span><b>${esc(title(p.run_mode||'production'))}</b></div></div><p class="sub">Git commit: ${esc(p.git_commit||'unknown')}</p>`;$('#provenance-json').textContent=JSON.stringify(p,null,2);$('#model-json').textContent=JSON.stringify(D.model,null,2)}
nav();overview();qc();controlPage();taxonomy();evidence();ecology();functions();provenance();$('#menu').onclick=()=>$('#sidebar').classList.toggle('open');
</script></body></html>'''


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser()
    sub = root.add_subparsers(dest="command", required=True)

    split = sub.add_parser("split-manifest")
    split.add_argument("--samples", type=Path, required=True)
    split.add_argument(
        "--base-dir", type=Path,
        help="Directory used to resolve relative FASTQ paths (default: manifest directory)",
    )
    split.add_argument("--outdir", type=Path, required=True)
    split.set_defaults(func=cmd_split)

    bracken = sub.add_parser("standardize-bracken")
    bracken.add_argument("--library-id", required=True)
    bracken.add_argument("--input", type=Path, required=True)
    bracken.add_argument("--output", type=Path, required=True)
    bracken.set_defaults(func=cmd_bracken)

    humann = sub.add_parser("standardize-humann")
    humann.add_argument("--library-id", required=True)
    humann.add_argument("--input-dir", type=Path, required=True)
    humann.add_argument("--output", type=Path, required=True)
    humann.set_defaults(func=cmd_humann)

    control = sub.add_parser("control-annotate")
    control.add_argument("--samples", type=Path, required=True)
    control.add_argument("--controls", type=Path, required=True)
    control.add_argument("--taxa", type=Path, nargs="+", required=True)
    control.add_argument("--functions", type=Path, nargs="+", required=True)
    control.add_argument("--unfiltered-taxa", type=Path, required=True)
    control.add_argument("--unfiltered-functions", type=Path, required=True)
    control.add_argument("--annotated-taxa", type=Path, required=True)
    control.add_argument("--annotated-functions", type=Path, required=True)
    control.add_argument("--high-confidence", type=Path, required=True)
    control.add_argument("--sample-confidence", type=Path, required=True)
    control.set_defaults(func=cmd_control)

    qc = sub.add_parser("summarize-qc")
    qc.add_argument("--samples", type=Path, required=True)
    qc.add_argument("--qc-dirs", type=Path, nargs="+", required=True)
    qc.add_argument("--output", type=Path, required=True)
    qc.set_defaults(func=cmd_summarize_qc)

    harmonize = sub.add_parser("harmonize")
    harmonize.add_argument("--input", type=Path, required=True)
    harmonize.add_argument("--crosswalk", type=Path)
    harmonize.add_argument("--output", type=Path, required=True)
    harmonize.add_argument("--output-crosswalk", type=Path, required=True)
    harmonize.set_defaults(func=cmd_harmonize)

    classify = sub.add_parser("classify")
    classify.add_argument("--harmonized", type=Path, required=True)
    classify.add_argument("--samples", type=Path, required=True)
    classify.add_argument("--assay-performance", type=Path)
    classify.add_argument("--evidence-policy", type=Path)
    classify.add_argument("--qc-metrics", type=Path)
    classify.add_argument("--model-output", type=Path)
    classify.add_argument("--output", type=Path, required=True)
    classify.set_defaults(func=cmd_classify)

    calibration = sub.add_parser("fit-calibration")
    calibration.add_argument("--evidence", type=Path, required=True)
    calibration.add_argument("--truth", type=Path, required=True)
    calibration.add_argument("--output", type=Path, required=True)
    calibration.set_defaults(func=cmd_fit_calibration)

    evaluation = sub.add_parser("evaluate-evidence")
    evaluation.add_argument("--evidence", type=Path, required=True)
    evaluation.add_argument("--truth", type=Path, required=True)
    evaluation.add_argument("--output", type=Path, required=True)
    evaluation.set_defaults(func=cmd_evaluate_evidence)

    concordance = sub.add_parser("concordance")
    concordance.add_argument("--harmonized", type=Path, required=True)
    concordance.add_argument("--samples", type=Path, required=True)
    concordance.add_argument("--output", type=Path, required=True)
    concordance.set_defaults(func=cmd_concordance)

    diversity = sub.add_parser("diversity")
    diversity.add_argument("--taxa", type=Path, required=True)
    diversity.add_argument("--functions", type=Path, required=True)
    diversity.add_argument("--alpha-output", type=Path, required=True)
    diversity.add_argument("--beta-output", type=Path, required=True)
    diversity.set_defaults(func=cmd_diversity)

    provenance = sub.add_parser("provenance")
    provenance.add_argument("--samples", type=Path, required=True)
    provenance.add_argument("--controls", type=Path, required=True)
    provenance.add_argument("--project-dir", type=Path, required=True)
    provenance.add_argument("--output", type=Path, required=True)
    provenance.set_defaults(func=cmd_provenance)

    report = sub.add_parser("report")
    report.add_argument("--evidence", type=Path, required=True)
    report.add_argument("--concordance", type=Path, required=True)
    report.add_argument("--alpha-diversity", type=Path, required=True)
    report.add_argument("--beta-diversity", type=Path, required=True)
    report.add_argument("--sample-confidence", type=Path, required=True)
    report.add_argument("--qc-metrics", type=Path, required=True)
    report.add_argument("--samples", type=Path, required=True)
    report.add_argument("--controls", type=Path, required=True)
    report.add_argument("--functions", type=Path, required=True)
    report.add_argument("--evidence-model", type=Path, required=True)
    report.add_argument("--provenance", type=Path, required=True)
    report.add_argument("--output", type=Path, required=True)
    report.set_defaults(func=cmd_report)
    return root


def main() -> int:
    args = parser().parse_args()
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

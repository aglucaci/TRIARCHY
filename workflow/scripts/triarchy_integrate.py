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


def cmd_mock_profile(args: argparse.Namespace) -> None:
    library_token = args.library_id.upper()
    control = (
        any(token in library_token for token in ("NEG", "BLANK", "NTC"))
        or library_token.startswith("EB")
        or "_EB" in library_token
    )
    if control:
        taxa = [("CTRL001", "Ralstonia", "genus", 80.0)]
    elif args.assay == "WMS":
        taxa = [
            ("WMS001", "Pseudomonas", "genus", 60.0),
            ("WMS002", "Aspergillus", "genus", 20.0),
        ]
    elif args.assay == "16S":
        taxa = [
            ("16S001", "Pseudomonas", "genus", 65.0),
            ("16S002", "Bacillus", "genus", 35.0),
        ]
    else:
        taxa = [
            ("ITS001", "Aspergillus", "genus", 70.0),
            ("ITS002", "Cladosporium", "genus", 30.0),
        ]
    rows = [
        dict(
            library_id=args.library_id,
            assay=args.assay,
            feature_id=feature_id,
            original_name=name,
            original_rank=rank,
            abundance=abundance,
        )
        for feature_id, name, rank, abundance in taxa
    ]
    write_tsv(args.output, rows, PROFILE_FIELDS)


def cmd_mock_functions(args: argparse.Namespace) -> None:
    rows = [
        {
            "library_id": args.library_id,
            "assay": "WMS_function",
            "feature_id": "PWY-0001",
            "original_name": "mock_pathway",
            "abundance": 1.0,
        }
    ]
    write_tsv(
        args.output,
        rows,
        ["library_id", "assay", "feature_id", "original_name", "abundance"],
    )


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
        "stub_run": os.environ.get("NXF_STUB", "").lower() == "true",
    }
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def table_html(path: Path) -> str:
    rows = read_tsv(path)
    if not rows:
        return "<p>No rows.</p>"
    fields = list(rows[0])
    head = "".join(f"<th>{html.escape(field)}</th>" for field in fields)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row.get(field, '')))}</td>" for field in fields) + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def cmd_report(args: argparse.Namespace) -> None:
    provenance = json.loads(args.provenance.read_text(encoding="utf-8"))
    evidence_model = json.loads(args.evidence_model.read_text(encoding="utf-8"))
    calibration_status = evidence_model["policy"].get("calibration", {}).get(
        "status", "uncalibrated"
    )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>TRIARCHY analytical report</title>
<style>
body{{font-family:system-ui,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#17212b}}
table{{border-collapse:collapse;width:100%;font-size:.85rem;margin:1rem 0 2rem}}
th,td{{border:1px solid #d9dfe7;padding:.4rem;text-align:left}}th{{background:#eef5f3}}
code,pre{{background:#f4f6f8;padding:.2rem .35rem}}.warning{{border-left:4px solid #d88900;padding:.7rem}}
</style></head><body>
<h1>TRIARCHY analytical report</h1>
<p class="warning">Cross-assay classifications summarize evidence; they do not make
raw relative abundances interchangeable.</p>
<p class="warning">Authenticity probability calibration status:
<strong>{html.escape(calibration_status)}</strong>. Synthetic calibration must
not be represented as empirical clinical or environmental validation.</p>
<h2>Evidence model</h2><pre>{html.escape(json.dumps(evidence_model, indent=2))}</pre>
<h2>Run provenance</h2><pre>{html.escape(json.dumps(provenance, indent=2))}</pre>
<h2>Which samples are trustworthy?</h2>{table_html(args.sample_confidence)}
<h2>Which taxa are supported?</h2>{table_html(args.evidence)}
<h2>Where do assays agree?</h2>{table_html(args.concordance)}
<h2>Within-assay alpha diversity</h2>{table_html(args.alpha_diversity)}
<h2>Within-assay beta diversity</h2>{table_html(args.beta_diversity)}
</body></html>
"""
    args.output.write_text(document, encoding="utf-8")


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

    mock = sub.add_parser("mock-profile")
    mock.add_argument("--library-id", required=True)
    mock.add_argument("--assay", choices=("WMS", "16S", "ITS"), required=True)
    mock.add_argument("--output", type=Path, required=True)
    mock.set_defaults(func=cmd_mock_profile)

    functions = sub.add_parser("mock-functions")
    functions.add_argument("--library-id", required=True)
    functions.add_argument("--output", type=Path, required=True)
    functions.set_defaults(func=cmd_mock_functions)

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

#!/usr/bin/env python3
"""Score a completed TRIARCHY run against ZymoBIOMICS mock-community truth.

Inputs are the run's results directory and the cohort truth table written by
the real-world downloader (for TRI-REAL-05: data/real_world/TRI-REAL-05/truth/
mock_truth.tsv). Comparisons are made at genus level, because 16S V4 amplicons
rarely resolve these organisms to species; genus synonyms used by SILVA, UNITE
and NCBI (for example Limosilactobacillus, Escherichia-Shigella) are accepted.

Outputs, written to --outdir (default <results>/mock_evaluation):
  mock_taxon_calls.tsv      one row per mock sample x expected genus
  mock_unexpected_taxa.tsv  evidence rows in mock samples for genera not in the standard
  mock_evaluation.json      summary metrics
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ASSAYS = ("WMS", "16S", "ITS")
ABUNDANCE_BINS = (
    (">=1%", 1.0, math.inf),
    ("0.1-1%", 0.1, 1.0),
    ("0.01-0.1%", 0.01, 0.1),
    ("<0.01%", 0.0, 0.01),
)
TRUE_VALUES = {"1", "true", "yes", "present"}

CALL_FIELDS = [
    "biological_sample_id", "standard", "dilution", "technical_replicate", "genus", "kingdom",
    "genomic_dna_pct", "abundance_bin", "expected_assays",
    "wms_detected", "16s_detected", "its_detected",
    "wms_relative_abundance", "16s_relative_abundance", "its_relative_abundance",
    "classification", "cross_assay_status", "authenticity_state", "biological_authenticity",
    "evidence_rank",
]
UNEXPECTED_FIELDS = [
    "biological_sample_id", "canonical_name", "comparison_rank", "kingdom",
    "wms_detected", "16s_detected", "its_detected",
    "classification", "cross_assay_status", "authenticity_state", "probable_contamination_source",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def genus_of(name: str) -> str:
    return name.strip().split(" ")[0] if name.strip() else ""


def is_true(value: str) -> bool:
    return str(value).strip().casefold() in TRUE_VALUES


def abundance_bin(pct: float) -> str:
    for label, low, high in ABUNDANCE_BINS:
        if low <= pct < high:
            return label
    return ABUNDANCE_BINS[-1][0]


def fraction(numerator: int, denominator: int) -> float | None:
    return round(numerator / denominator, 4) if denominator else None


def library_genus_abundance(harmonized: list[dict[str, str]]) -> dict[str, dict[str, float]]:
    """Relative abundance of each genus within each library (assay-native counts)."""
    totals: dict[str, float] = defaultdict(float)
    by_genus: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for row in harmonized:
        try:
            value = max(0.0, float(row.get("abundance") or 0))
        except ValueError:
            continue
        library = row["library_id"]
        totals[library] += value
        genus = genus_of(row.get("canonical_name") or row.get("original_name", ""))
        if genus:
            by_genus[library][genus] += value
    return {
        library: {genus: value / totals[library] for genus, value in genera.items() if totals[library]}
        for library, genera in by_genus.items()
    }


def evaluate(results: Path, truth_path: Path, samples_path: Path) -> tuple[dict, list, list]:
    evidence_path = results / "integrated" / "taxon_evidence.tsv"
    harmonized_path = results / "integrated" / "harmonized_taxa.tsv"
    for path in (evidence_path, harmonized_path, truth_path, samples_path):
        if not path.is_file():
            raise FileNotFoundError(f"required input is missing: {path}")

    truth = read_tsv(truth_path)
    evidence = read_tsv(evidence_path)
    abundance = library_genus_abundance(read_tsv(harmonized_path))
    libraries: dict[tuple[str, str], str] = {}
    for row in read_tsv(samples_path):
        libraries[(row["biological_sample_id"], row["assay"])] = row["library_id"]

    evidence_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in evidence:
        evidence_by_sample[row["biological_sample_id"]].append(row)

    mock_samples = sorted({row["biological_sample_id"] for row in truth})
    expected_aliases: dict[str, set[str]] = defaultdict(set)
    calls: list[dict[str, object]] = []
    for item in truth:
        sample = item["biological_sample_id"]
        aliases = set(item["genus_aliases"].split(";"))
        expected_aliases[sample] |= aliases
        matching = [row for row in evidence_by_sample[sample] if genus_of(row["canonical_name"]) in aliases]
        genus_rows = [row for row in matching if row.get("comparison_rank") == "genus"] or matching
        best = max(genus_rows, key=lambda row: float(row.get("biological_authenticity") or 0), default=None)
        call: dict[str, object] = {
            "biological_sample_id": sample,
            "standard": item["standard"],
            "dilution": item["dilution"],
            "technical_replicate": item["technical_replicate"],
            "genus": item["genus"],
            "kingdom": item["kingdom"],
            "genomic_dna_pct": item["genomic_dna_pct"],
            "abundance_bin": abundance_bin(float(item["genomic_dna_pct"])),
            "expected_assays": item["expected_assays"],
            "classification": best["classification"] if best else "no_evidence_row",
            "cross_assay_status": best["cross_assay_status"] if best else "no_evidence_row",
            "authenticity_state": best["authenticity_state"] if best else "no_evidence_row",
            "biological_authenticity": best.get("biological_authenticity", "") if best else "",
            "evidence_rank": best.get("comparison_rank", "") if best else "",
        }
        for assay in ASSAYS:
            key = assay.lower()
            library = libraries.get((sample, assay), "")
            observed = sum(abundance.get(library, {}).get(alias, 0.0) for alias in aliases)
            detected = observed > 0 or any(is_true(row.get(f"{key}_detected", "")) for row in matching)
            call[f"{key}_detected"] = str(detected).lower()
            call[f"{key}_relative_abundance"] = f"{observed:.6g}"
        calls.append(call)

    unexpected: list[dict[str, object]] = []
    for sample in mock_samples:
        for row in evidence_by_sample[sample]:
            if genus_of(row["canonical_name"]) not in expected_aliases[sample]:
                unexpected.append(row)

    summary = summarize(calls, unexpected, mock_samples, evidence, set(mock_samples))
    return summary, calls, unexpected


def summarize(
    calls: list[dict[str, object]],
    unexpected: list[dict[str, str]],
    mock_samples: list[str],
    evidence: list[dict[str, str]],
    mock_set: set[str],
) -> dict[str, object]:
    recall: dict[str, dict[str, object]] = {}
    for assay in ASSAYS:
        key = f"{assay.lower()}_detected"
        relevant = [call for call in calls if assay in str(call["expected_assays"]).split(";")]
        by_bin = {}
        for label, _, _ in ABUNDANCE_BINS:
            subset = [call for call in relevant if call["abundance_bin"] == label]
            by_bin[label] = {
                "expected": len(subset),
                "detected": sum(call[key] == "true" for call in subset),
                "recall": fraction(sum(call[key] == "true" for call in subset), len(subset)),
            }
        by_dilution = {}
        for dilution in sorted({str(call["dilution"]) for call in relevant}, key=lambda d: (len(d), d)):
            subset = [call for call in relevant if call["dilution"] == dilution]
            by_dilution[dilution or "undiluted"] = fraction(sum(call[key] == "true" for call in subset), len(subset))
        detected = sum(call[key] == "true" for call in relevant)
        recall[assay] = {
            "expected": len(relevant),
            "detected": detected,
            "recall": fraction(detected, len(relevant)),
            "by_genomic_dna_abundance": by_bin,
            "by_dilution": by_dilution,
        }

    supported = sum(call["authenticity_state"] == "biologically_supported" for call in calls)
    unexpected_supported = sum(row.get("authenticity_state") == "biologically_supported" for row in unexpected)
    non_mock = [row for row in evidence if row["biological_sample_id"] not in mock_set]
    return {
        "mock_samples": len(mock_samples),
        "expected_sample_genus_pairs": len(calls),
        "assay_recall": recall,
        "expected_taxa": {
            "biologically_supported": supported,
            "fraction_supported": fraction(supported, len(calls)),
            "cross_assay_status": dict(Counter(str(call["cross_assay_status"]) for call in calls)),
            "authenticity_state": dict(Counter(str(call["authenticity_state"]) for call in calls)),
        },
        "unexpected_taxa_in_mocks": {
            "evidence_rows": len(unexpected),
            "biologically_supported": unexpected_supported,
            "false_support_rate": fraction(unexpected_supported, len(unexpected)),
            "cross_assay_status": dict(Counter(row.get("cross_assay_status", "") for row in unexpected)),
        },
        "non_mock_samples": {
            "evidence_rows": len(non_mock),
            "contamination_associated": sum(
                row.get("cross_assay_status") == "contamination_associated" for row in non_mock
            ),
            "cross_assay_status": dict(Counter(row.get("cross_assay_status", "") for row in non_mock)),
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, required=True, help="TRIARCHY results directory")
    parser.add_argument("--truth", type=Path, required=True, help="mock_truth.tsv from the cohort")
    parser.add_argument("--samples", type=Path, required=True, help="the cohort samples.tsv used for the run")
    parser.add_argument("--outdir", type=Path, help="output directory (default: <results>/mock_evaluation)")
    args = parser.parse_args(argv)
    outdir = args.outdir or args.results / "mock_evaluation"
    try:
        summary, calls, unexpected = evaluate(args.results, args.truth, args.samples)
    except (FileNotFoundError, KeyError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    summary = {
        "evaluation": "TRIARCHY ZymoBIOMICS mock-community truth",
        "evaluated_at": datetime.now(timezone.utc).isoformat(),
        "results": str(args.results),
        "truth": str(args.truth),
        **summary,
    }
    write_tsv(outdir / "mock_taxon_calls.tsv", CALL_FIELDS, calls)
    write_tsv(outdir / "mock_unexpected_taxa.tsv", UNEXPECTED_FIELDS, unexpected)
    (outdir / "mock_evaluation.json").parent.mkdir(parents=True, exist_ok=True)
    (outdir / "mock_evaluation.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    recall = summary["assay_recall"]
    print(
        f"Mock evaluation: {summary['mock_samples']} samples; recall WMS {recall['WMS']['recall']}, "
        f"16S {recall['16S']['recall']}, ITS {recall['ITS']['recall']}; "
        f"expected taxa biologically supported {summary['expected_taxa']['fraction_supported']}; "
        f"unexpected taxa falsely supported {summary['unexpected_taxa_in_mocks']['false_support_rate']}"
    )
    print(f"Wrote {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

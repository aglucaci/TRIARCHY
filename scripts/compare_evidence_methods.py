#!/usr/bin/env python3
"""Compare TRIARCHY and baseline taxon-evidence predictions against truth."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


STATUS_ALIASES = {
    "concordant": "concordant_detection",
    "single_assay_supported": "single_assay_credible",
    "single-assay_supported": "single_assay_credible",
    "assay_limited": "assay_limited_nondetection",
    "contamination_associated": "probable_contaminant",
    "probable_contamination": "probable_contaminant",
    "discordant": "genuine_discordance",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def prediction_key(row: dict[str, str]) -> tuple[str, str]:
    return (
        row.get("biological_sample_id", row.get("sample_id", "")),
        row.get("canonical_name", row.get("taxon", "")),
    )


def normalize_status(value: str) -> str:
    normalized = value.strip().casefold().replace("-", "_")
    return STATUS_ALIASES.get(normalized, normalized)


def abundance_value(row: dict[str, str]) -> float | None:
    for field in (
        "truth_abundance",
        "relative_abundance",
        "genomic_dna_abundance",
        "genome_copy_abundance",
        "expected_abundance",
    ):
        value = row.get(field, "")
        if value == "":
            continue
        try:
            return float(value)
        except ValueError:
            continue
    return None


def compare(method: str, prediction_path: Path, truth: dict[tuple[str, str], dict[str, str]]) -> dict[str, object]:
    predictions = {prediction_key(row): row for row in read_tsv(prediction_path)}
    tp = fp = fn = contaminant_tp = contaminant_fp = status_correct = matched = 0
    false_concordance = false_discordance = rank_correct = rank_matched = 0
    low_abundance_total = low_abundance_retained = 0
    brier = []
    for key, expected in truth.items():
        row = predictions.get(key, {})
        actual = expected.get("truth_status", "").casefold() in {
            "authentic", "present", "true", "1",
        }
        predicted = row.get("biological_state", "").casefold() == "present"
        tp += int(actual and predicted)
        fp += int(not actual and predicted)
        fn += int(actual and not predicted)
        expected_contaminant = expected.get("truth_status", "").casefold() == "contaminant"
        predicted_contaminant = row.get("classification", "") in {
            "probable_contaminant", "cross_contamination_candidate",
        }
        contaminant_tp += int(expected_contaminant and predicted_contaminant)
        contaminant_fp += int(not expected_contaminant and predicted_contaminant)
        expected_status = expected.get(
            "expected_cross_assay_status", expected.get("classification", "")
        )
        if expected_status and row:
            expected_status = normalize_status(expected_status)
            predicted_status = normalize_status(row.get("classification", ""))
            matched += 1
            status_correct += int(predicted_status == expected_status)
            false_concordance += int(
                predicted_status == "concordant_detection"
                and expected_status != "concordant_detection"
            )
            false_discordance += int(
                predicted_status == "genuine_discordance"
                and expected_status != "genuine_discordance"
            )
        expected_rank = normalize_status(expected.get("supported_comparison_rank", ""))
        predicted_rank = normalize_status(row.get("comparison_rank", ""))
        if expected_rank and row:
            rank_matched += 1
            rank_correct += int(expected_rank == predicted_rank)
        abundance = abundance_value(expected)
        if actual and abundance is not None and abundance <= 0.01:
            low_abundance_total += 1
            low_abundance_retained += int(predicted)
        calibration_score = row.get("biological_authenticity", row.get("probability_authentic", ""))
        if calibration_score != "":
            brier.append((float(calibration_score) - int(actual)) ** 2)
    precision = safe_ratio(tp, tp + fp)
    recall = safe_ratio(tp, tp + fn)
    contaminant_precision = safe_ratio(contaminant_tp, contaminant_tp + contaminant_fp)
    contaminant_total = sum(
        row.get("truth_status", "").casefold() == "contaminant"
        for row in truth.values()
    )
    return {
        "method": method,
        "prediction_file": prediction_path,
        "truth_observations": len(truth),
        "presence_precision": f"{precision:.6g}",
        "presence_recall": f"{recall:.6g}",
        "presence_f1": f"{safe_ratio(2 * precision * recall, precision + recall):.6g}",
        "contaminant_precision": f"{contaminant_precision:.6g}",
        "contaminant_recall": f"{safe_ratio(contaminant_tp, contaminant_total):.6g}",
        "explanatory_status_accuracy": f"{safe_ratio(status_correct, matched):.6g}",
        "false_concordance_rate": f"{safe_ratio(false_concordance, matched):.6g}",
        "false_discordance_rate": f"{safe_ratio(false_discordance, matched):.6g}",
        "comparison_rank_accuracy": f"{safe_ratio(rank_correct, rank_matched):.6g}",
        "low_abundance_authentic_retention": (
            f"{safe_ratio(low_abundance_retained, low_abundance_total):.6g}"
            if low_abundance_total else ""
        ),
        "brier_score": f"{sum(brier) / len(brier):.6g}" if brier else "",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument(
        "--prediction", action="append", required=True, metavar="METHOD=FILE",
        help="Repeat for independent, naive, harmonized, decontam, global-control, Com-2seq, and TRIARCHY outputs.",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    truth = {prediction_key(row): row for row in read_tsv(args.truth)}
    results = []
    for value in args.prediction:
        if "=" not in value:
            parser.error("--prediction must use METHOD=FILE")
        method, filename = value.split("=", 1)
        results.append(compare(method, Path(filename), truth))
    fields = list(results[0])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

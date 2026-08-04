from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "scripts" / "compare_evidence_methods.py"
SPEC = importlib.util.spec_from_file_location("compare_evidence_methods", MODULE)
assert SPEC and SPEC.loader
comparison = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(comparison)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class EvidenceComparisonTests(unittest.TestCase):
    def test_perfect_predictions_have_unit_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            truth_path = root / "truth.tsv"
            prediction_path = root / "prediction.tsv"
            truth_rows = [
                {
                    "biological_sample_id": "S1",
                    "taxon": "TaxonA",
                    "truth_status": "authentic",
                    "expected_cross_assay_status": "concordant_detection",
                },
                {
                    "biological_sample_id": "S1",
                    "taxon": "TaxonB",
                    "truth_status": "contaminant",
                    "expected_cross_assay_status": "probable_contaminant",
                },
            ]
            prediction_rows = [
                {
                    "biological_sample_id": "S1",
                    "canonical_name": "TaxonA",
                    "biological_state": "present",
                    "classification": "concordant_detection",
                    "probability_authentic": "1",
                },
                {
                    "biological_sample_id": "S1",
                    "canonical_name": "TaxonB",
                    "biological_state": "indeterminate",
                    "classification": "probable_contaminant",
                    "probability_authentic": "0",
                },
            ]
            write_tsv(truth_path, truth_rows)
            write_tsv(prediction_path, prediction_rows)
            truth = {
                comparison.prediction_key(row): row
                for row in comparison.read_tsv(truth_path)
            }
            result = comparison.compare("TRIARCHY", prediction_path, truth)
            self.assertEqual("1", result["presence_precision"])
            self.assertEqual("1", result["presence_recall"])
            self.assertEqual("1", result["contaminant_precision"])
            self.assertEqual("1", result["contaminant_recall"])
            self.assertEqual("1", result["explanatory_status_accuracy"])
            self.assertEqual("0", result["brier_score"])


if __name__ == "__main__":
    unittest.main()

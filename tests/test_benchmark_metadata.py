from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALID = ROOT / "benchmarks" / "TRI-00" / "valid"
INVALID = ROOT / "benchmarks" / "TRI-00" / "invalid"
MODULE_PATH = ROOT / "scripts" / "validate_benchmark_metadata.py"
SPEC = importlib.util.spec_from_file_location("triarchy_metadata", MODULE_PATH)
assert SPEC and SPEC.loader
metadata = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(metadata)


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


class MetadataContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.sample_fields, self.samples = read_rows(VALID / "samples.tsv")
        self.control_fields, self.controls = read_rows(
            VALID / "control_relationships.tsv"
        )

    def write_case(
        self,
        samples: list[dict[str, str]] | None = None,
        controls: list[dict[str, str]] | None = None,
    ) -> tuple[Path, Path, tempfile.TemporaryDirectory[str]]:
        temp = tempfile.TemporaryDirectory()
        root = Path(temp.name)
        sample_path = root / "samples.tsv"
        control_path = root / "control_relationships.tsv"
        with sample_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, self.sample_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(samples if samples is not None else self.samples)
        with control_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, self.control_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(controls if controls is not None else self.controls)
        return sample_path, control_path, temp

    def assert_error(
        self,
        code: str,
        samples: list[dict[str, str]] | None = None,
        controls: list[dict[str, str]] | None = None,
    ) -> None:
        sample_path, control_path, temp = self.write_case(samples, controls)
        self.addCleanup(temp.cleanup)
        with self.assertRaises(metadata.MetadataError) as caught:
            metadata.validate(sample_path, control_path)
        self.assertEqual(code, caught.exception.code)

    def test_valid_graph(self) -> None:
        self.assertEqual(
            (16, 9),
            metadata.validate(
                VALID / "samples.tsv", VALID / "control_relationships.tsv"
            ),
        )

    def test_duplicate_library_id(self) -> None:
        rows = [row.copy() for row in self.samples]
        rows[1]["library_id"] = rows[0]["library_id"]
        self.assert_error("TRI-META-001", samples=rows)

    def test_missing_biological_sample_id(self) -> None:
        rows = [row.copy() for row in self.samples]
        rows[0]["biological_sample_id"] = ""
        self.assert_error("TRI-META-002", samples=rows)

    def test_invalid_assay(self) -> None:
        rows = [row.copy() for row in self.samples]
        rows[0]["assay"] = "WGS"
        self.assert_error("TRI-META-003", samples=rows)

    def test_paired_missing_r2(self) -> None:
        rows = [row.copy() for row in self.samples]
        rows[0]["fastq_r2"] = ""
        self.assert_error("TRI-META-004", samples=rows)

    def test_duplicate_fastq_path(self) -> None:
        rows = [row.copy() for row in self.samples]
        rows[1]["fastq_r1"] = rows[0]["fastq_r1"]
        self.assert_error("TRI-META-005", samples=rows)

    def test_unknown_control_id(self) -> None:
        rows = [row.copy() for row in self.controls]
        rows[0]["control_library_id"] = "DOES_NOT_EXIST"
        self.assert_error("TRI-META-010", controls=rows)

    def test_pcr_control_on_wms(self) -> None:
        rows = [row.copy() for row in self.controls]
        rows[0].update(
            control_library_id="SEQ_NEG_WMS",
            control_stage="pcr",
            relationship="run_level",
        )
        self.assert_error("TRI-META-011", controls=rows)

    def test_cross_batch_control(self) -> None:
        rows = [row.copy() for row in self.controls]
        rows[0]["control_library_id"] = "EXT_NEG_02"
        self.assert_error("TRI-META-014", controls=rows)

    def test_serialized_invalid_fixture_catalog(self) -> None:
        _, cases = read_rows(INVALID / "expected_errors.tsv")
        for case in cases:
            with self.subTest(case=case["case"]):
                fixture = INVALID / case["fixture"]
                if fixture.name.endswith(".samples.tsv"):
                    _, sample_rows = read_rows(fixture)
                    sample_path, control_path, temp = self.write_case(
                        samples=sample_rows, controls=[]
                    )
                else:
                    _, control_rows = read_rows(fixture)
                    sample_path, control_path, temp = self.write_case(
                        controls=control_rows
                    )
                self.addCleanup(temp.cleanup)
                with self.assertRaises(metadata.MetadataError) as caught:
                    metadata.validate(sample_path, control_path)
                self.assertEqual(case["error_code"], caught.exception.code)
                self.assertIn(
                    case["expected_message_fragment"], str(caught.exception)
                )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from importlib.machinery import SourceFileLoader
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "bin" / "triarchy"
SPEC = importlib.util.spec_from_loader(
    "triarchy_cli", SourceFileLoader("triarchy_cli", str(MODULE))
)
assert SPEC and SPEC.loader
cli = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cli)


class TriarchyCliTests(unittest.TestCase):
    def test_validate_accepts_metadata_when_fastq_check_is_skipped(self) -> None:
        result = cli.main([
            "validate",
            "--input",
            str(ROOT / "benchmarks" / "TRI-01" / "metadata" / "samples.tsv"),
            "--paired-controls",
            str(ROOT / "benchmarks" / "TRI-01" / "metadata" / "control_relationships.tsv"),
            "--skip-fastq-check",
        ])
        self.assertEqual(0, result)

    def test_database_manifest_reports_missing_required_reference(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            report = cli.check_databases(
                ROOT / "conf" / "databases.example.yml",
                Path(temp),
            )
            self.assertFalse(report["ok"])
            self.assertEqual("missing", report["references"][0]["status"])


if __name__ == "__main__":
    unittest.main()

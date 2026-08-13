from __future__ import annotations

import importlib.util
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "download_real_world_datasets.py"
SPEC = importlib.util.spec_from_file_location("triarchy_real_world", MODULE_PATH)
assert SPEC and SPEC.loader
real_world = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = real_world
SPEC.loader.exec_module(real_world)


def row(run: str, group: str, assay: str, role: str = "biological") -> dict[str, str]:
    values = {field: "" for field in real_world.ENA_FIELDS}
    values.update(
        run_accession=run,
        sample_alias=group,
        triarchy_dataset_id="TRI-REAL-TEST",
        triarchy_sample_group=group,
        triarchy_assay=assay,
        triarchy_sample_role=role,
        fastq_ftp=(
            f"ftp.sra.ebi.ac.uk/vol1/fastq/{run}/{run}_1.fastq.gz;"
            f"ftp.sra.ebi.ac.uk/vol1/fastq/{run}/{run}_2.fastq.gz"
        ),
        fastq_md5="a;b",
        fastq_bytes="10;20",
    )
    return values


class RealWorldDownloaderTests(unittest.TestCase):
    def test_assay_inference(self) -> None:
        self.assertEqual("WMS", real_world.infer_assay({"library_strategy": "WGS", "experiment_title": "shotgun metagenome"}))
        self.assertEqual("16S", real_world.infer_assay({"library_strategy": "AMPLICON", "experiment_title": "16S V4 amplicon"}))
        self.assertEqual("ITS", real_world.infer_assay({"library_strategy": "AMPLICON", "experiment_title": "fungal ITS1 amplicon"}))
        self.assertEqual("AMPLICON", real_world.infer_assay({"library_strategy": "AMPLICON", "experiment_title": "MiSeq sequencing"}))
        self.assertEqual("UNKNOWN", real_world.infer_assay({"library_strategy": "RNA-Seq", "experiment_title": "transcriptome"}))

    def test_role_inference(self) -> None:
        self.assertEqual("negative_control", real_world.infer_sample_role({"sample_title": "extraction blank"}))
        self.assertEqual("positive_control", real_world.infer_sample_role({"sample_title": "Zymo mock dilution"}))
        self.assertEqual("biological", real_world.infer_sample_role({"sample_title": "leaf sample"}))

    def test_ena_ftp_paths_are_downloadable_https(self) -> None:
        self.assertEqual(
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR/ERR_1.fastq.gz",
            real_world.ena_file_url("ftp.sra.ebi.ac.uk/vol1/fastq/ERR/ERR_1.fastq.gz"),
        )
        self.assertEqual(
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR/ERR_1.fastq.gz",
            real_world.ena_file_url("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR/ERR_1.fastq.gz"),
        )

    def test_complete_case_selection_requires_all_three_assays(self) -> None:
        rows = [
            row("RUN1", "S1", "WMS"),
            row("RUN2", "S1", "16S"),
            row("RUN3", "S1", "ITS"),
            row("RUN4", "S2", "WMS"),
            row("RUN5", "S2", "16S"),
        ]
        selected, summary = real_world.select_groups(rows, "complete_cases", 10, None, include_all_runs=False)
        self.assertEqual({"S1"}, {item["triarchy_sample_group"] for item in selected})
        self.assertEqual(1, summary["complete_groups"])
        self.assertEqual(0, summary["wms_amplicon_groups"])

    def test_wms_amplicon_groups_are_selectable_when_marker_is_unresolved(self) -> None:
        rows = [
            row("RUN1", "S1", "WMS"),
            row("RUN2", "S1", "AMPLICON"),
            row("RUN3", "S2", "AMPLICON"),
        ]
        selected, summary = real_world.select_groups(rows, "complete_cases", 10, None, include_all_runs=False)
        self.assertEqual({"S1"}, {item["triarchy_sample_group"] for item in selected})
        self.assertEqual(1, summary["wms_amplicon_groups"])

    def test_qiita_grouping_strips_extraction_suffix(self) -> None:
        self.assertEqual(
            "mock.community.d4",
            real_world.normalize_qiita_core("qiita_sid_12201:12201.mock.community.d4.homebrew.L"),
        )
        self.assertEqual(
            "BLANK.C8",
            real_world.normalize_qiita_core("12201.BLANK.C8.microbiome.H"),
        )

    def test_mock_control_selection_prioritizes_controls(self) -> None:
        rows = [
            row("BIO1", "BIO", "WMS"),
            row("BIO2", "BIO", "16S"),
            row("BIO3", "BIO", "ITS"),
            row("NEG1", "NEG", "WMS", "negative_control"),
            row("NEG2", "NEG", "16S", "negative_control"),
        ]
        selected, summary = real_world.select_groups(rows, "mock_controls", 1, None, include_all_runs=False)
        self.assertEqual({"NEG"}, {item["triarchy_sample_group"] for item in selected})
        self.assertEqual(1, summary["control_or_mock_groups"])

    def test_mock_control_selection_keeps_room_for_biological_groups(self) -> None:
        rows = []
        for index in range(1, 8):
            rows.extend([
                row(f"NEG{index}A", f"NEG{index}", "WMS", "negative_control"),
                row(f"NEG{index}B", f"NEG{index}", "AMPLICON", "negative_control"),
            ])
        rows.extend([
            row("BIO1A", "BIO1", "WMS"),
            row("BIO1B", "BIO1", "AMPLICON"),
        ])
        selected, summary = real_world.select_groups(rows, "mock_controls", 6, 2, include_all_runs=False)
        self.assertEqual(2, summary["selected_control_or_mock_groups"])
        self.assertIn("BIO1", {item["triarchy_sample_group"] for item in selected})

    def test_download_manifest_expands_paired_fastqs(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            dataset_dir = Path(temp) / "TRI-REAL-TEST"
            manifest = real_world.build_download_manifest(dataset_dir, [row("RUN1", "S1", "WMS")])
        self.assertEqual(2, len(manifest))
        self.assertEqual("reads/RUN1/RUN1_1.fastq.gz", manifest[0]["local_path"])
        self.assertEqual("a", manifest[0]["md5"])
        self.assertEqual("20", manifest[1]["bytes"])

    def test_draft_triarchy_metadata_maps_amplicon_and_fastq_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            dataset_dir = Path(temp) / "TRI-REAL-TEST"
            samples, controls, summary = real_world.build_triarchy_metadata(
                dataset_dir,
                [row("RUN1", "S1", "AMPLICON")],
                "16S",
            )
        self.assertEqual([], controls)
        self.assertEqual(1, summary["samples"])
        self.assertEqual("16S", samples[0]["assay"])
        self.assertEqual("../reads/RUN1/RUN1_1.fastq.gz", samples[0]["fastq_r1"])
        self.assertEqual("../reads/RUN1/RUN1_2.fastq.gz", samples[0]["fastq_r2"])

    def test_download_one_supports_quiet_local_transfer(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "source.fastq.gz"
            destination = root / "reads" / "copy.fastq.gz"
            payload = b"ACGT" * 1024
            source.write_bytes(payload)
            result = real_world.download_one(
                source.as_uri(),
                destination,
                hashlib.md5(payload).hexdigest(),
                str(len(payload)),
                force=False,
                verify_md5=True,
                timeout=10,
                retries=1,
                progress_label="test [1/1]",
                progress_interval=0,
                quiet=True,
            )
        self.assertEqual("downloaded", result["status"])
        self.assertEqual(len(payload), result["bytes"])

    def test_human_progress_formatters(self) -> None:
        self.assertEqual("1.0 KiB", real_world.human_bytes(1024))
        self.assertEqual("1m 05s", real_world.format_duration(65))


if __name__ == "__main__":
    unittest.main()

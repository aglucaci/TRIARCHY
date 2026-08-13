from __future__ import annotations

import csv
import gzip
import hashlib
import json
import unittest
from pathlib import Path

from scripts.validate_benchmark_metadata import validate


ROOT = Path(__file__).resolve().parents[1]
DATASET_IDS = (
    "TRI-01", "TRI-02", "TRI-03", "TRI-05", "TRI-06", "TRI-07",
    "TRI-10", "TRI-11", "TRI-12", "TRI-13", "TRI-14", "TRI-15", "TRI-16",
    "TRI-17", "TRI-18", "TRI-19",
    "TRI-20", "TRI-21", "TRI-30", "TRI-31",
)
EXPECTED_LIBRARY_COUNTS = {
    "TRI-01": 12,
    "TRI-02": 18,
    "TRI-03": 18,
    "TRI-05": 24,
    "TRI-06": 18,
    "TRI-07": 25,
    "TRI-10": 21,
    "TRI-11": 24,
    "TRI-12": 112,
    "TRI-13": 75,
    "TRI-14": 48,
    "TRI-15": 97,
    "TRI-16": 30,
    "TRI-17": 36,
    "TRI-18": 54,
    "TRI-19": 84,
    "TRI-20": 45,
    "TRI-21": 114,
    "TRI-30": 192,
    "TRI-31": 51,
}
REQUIRED_TRUTH = {
    "community_members.tsv",
    "sample_composition_genomic_dna.tsv",
    "sample_composition_genome_copy.tsv",
    "sample_composition_16s.tsv",
    "sample_composition_its.tsv",
    "read_origin_wms.tsv.gz",
    "true_16s_asvs.fasta",
    "true_its_asvs.fasta",
    "taxonomic_crosswalk.tsv",
    "contaminants.tsv",
    "control_relationships.tsv",
    "expected_sample_qc.tsv",
    "expected_taxon_status.tsv",
    "expected_pathways.tsv",
    "generation_manifest.json",
    "precomputed_profiles.tsv",
    "expected_assertions.json",
}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


class BenchmarkDatasetTests(unittest.TestCase):
    def test_metadata_fastqs_and_truth_are_complete(self) -> None:
        for dataset_id in DATASET_IDS:
            with self.subTest(dataset=dataset_id):
                root = ROOT / "benchmarks" / dataset_id
                sample_path = root / "metadata" / "samples.tsv"
                control_path = root / "metadata" / "control_relationships.tsv"
                sample_count, _ = validate(sample_path, control_path)
                samples = rows(sample_path)
                truth = rows(root / "truth" / "expected_taxa.tsv")
                origins = rows(root / "truth" / "read_origins.tsv")
                qc = {
                    row["library_id"]: row
                    for row in rows(root / "truth" / "expected_sample_qc.tsv")
                }
                self.assertEqual(sample_count, len(samples))
                self.assertEqual(EXPECTED_LIBRARY_COUNTS[dataset_id], len(samples))
                self.assertEqual({"WMS", "16S", "ITS"}, {row["assay"] for row in samples})
                sample_ids = {row["library_id"] for row in samples}
                truth_ids = {row["library_id"] for row in truth}
                self.assertTrue(truth_ids <= sample_ids)
                self.assertEqual(
                    {
                        library_id for library_id in truth_ids
                        if qc[library_id]["expected_qc_status"] not in {"failed", "empty"}
                    },
                    {row["library_id"] for row in origins},
                )
                for sample in samples:
                    fields = ["fastq_r1"]
                    if sample["library_layout"] == "paired":
                        fields.append("fastq_r2")
                    else:
                        self.assertEqual("", sample["fastq_r2"])
                    for field in fields:
                        fastq = (sample_path.parent / sample[field]).resolve()
                        self.assertTrue(fastq.is_file(), fastq)
                        with gzip.open(fastq, "rt", encoding="ascii") as handle:
                            lines = sum(1 for _ in handle)
                        if qc[sample["library_id"]]["expected_qc_status"] in {"failed", "empty", "no_target_template"}:
                            self.assertEqual(0, lines)
                        else:
                            self.assertGreater(lines, 0)
                        self.assertEqual(0, lines % 4)

    def test_required_truth_package_is_complete(self) -> None:
        for dataset_id in DATASET_IDS:
            with self.subTest(dataset=dataset_id):
                truth_dir = ROOT / "benchmarks" / dataset_id / "truth"
                self.assertTrue(REQUIRED_TRUTH <= {path.name for path in truth_dir.iterdir()})
                for name in REQUIRED_TRUTH:
                    self.assertTrue((truth_dir / name).is_file(), name)

    def test_manifests_match_generated_files(self) -> None:
        for dataset_id in DATASET_IDS:
            with self.subTest(dataset=dataset_id):
                root = ROOT / "benchmarks" / dataset_id
                manifest = json.loads(
                    (root / "truth" / "generation_manifest.json").read_text()
                )
                self.assertEqual(dataset_id, manifest["dataset_id"])
                self.assertEqual(EXPECTED_LIBRARY_COUNTS[dataset_id], manifest["library_count"])
                for relative, expected in manifest["checksums"].items():
                    path = root / relative
                    self.assertTrue(path.is_file(), path)
                    self.assertEqual(expected, digest(path), relative)

    def test_artificial_references_are_redistributable(self) -> None:
        required = {
            "artificial_bacterial.fasta", "artificial_archaeal.fasta",
            "artificial_fungal.fasta", "artificial_host.fasta",
        }
        for dataset_id in DATASET_IDS:
            with self.subTest(dataset=dataset_id):
                reference_dir = ROOT / "benchmarks" / dataset_id / "references"
                reference_names = {path.name for path in reference_dir.iterdir()}
                self.assertTrue(required <= reference_names)
                if dataset_id != "TRI-31":
                    self.assertEqual(required, reference_names)
                else:
                    self.assertTrue({"mini_function_catalog.fasta", "mini_function_mapping.tsv"} <= reference_names)
                for name in required:
                    text = (reference_dir / name).read_text(encoding="ascii")
                    if name == "artificial_host.fasta":
                        self.assertTrue(text.startswith(">TRIARCHY_ARTIFICIAL_HOST\n"))

    def test_designed_complexity_and_scale(self) -> None:
        self.assertGreater(EXPECTED_LIBRARY_COUNTS["TRI-02"], EXPECTED_LIBRARY_COUNTS["TRI-01"])
        self.assertGreater(EXPECTED_LIBRARY_COUNTS["TRI-12"], EXPECTED_LIBRARY_COUNTS["TRI-10"])
        self.assertGreaterEqual(EXPECTED_LIBRARY_COUNTS["TRI-30"], 144)
        self.assertEqual(51, EXPECTED_LIBRARY_COUNTS["TRI-31"])


if __name__ == "__main__":
    unittest.main()

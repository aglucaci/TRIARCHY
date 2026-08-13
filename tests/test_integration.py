from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "workflow" / "scripts" / "triarchy_integrate.py"
SPEC = importlib.util.spec_from_file_location("triarchy_integrate", MODULE_PATH)
assert SPEC and SPEC.loader
integration = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(integration)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class IntegrationTests(unittest.TestCase):
    def test_split_manifest_resolves_relative_fastqs(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = ROOT / "benchmarks" / "TRI-00" / "valid" / "samples.tsv"
            output = root / "split"
            integration.cmd_split(argparse.Namespace(samples=source, outdir=output))
            wms = read_tsv(output / "wms.tsv")
            self.assertTrue(wms)
            self.assertTrue(Path(wms[0]["fastq_r1"]).is_absolute())
            self.assertEqual({"WMS"}, {row["assay"] for row in wms})

    def test_split_manifest_uses_original_base_after_staging(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            original = root / "dataset" / "metadata"
            staged = root / "work"
            reads = root / "dataset" / "reads"
            original.mkdir(parents=True)
            staged.mkdir()
            reads.mkdir()
            (reads / "R1.fastq.gz").touch()
            manifest = staged / "samples.tsv"
            fields = ["library_id", "assay", "fastq_r1", "fastq_r2"]
            integration.write_tsv(
                manifest,
                [{"library_id": "L1", "assay": "WMS", "fastq_r1": "../reads/R1.fastq.gz", "fastq_r2": ""}],
                fields,
            )
            integration.cmd_split(
                argparse.Namespace(samples=manifest, base_dir=original, outdir=staged / "split")
            )
            row = read_tsv(staged / "split" / "wms.tsv")[0]
            self.assertEqual((reads / "R1.fastq.gz").resolve(), Path(row["fastq_r1"]))

    def test_mock_profiles_harmonize_by_canonical_name(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "profiles.tsv"
            rows = []
            for library_id, assay in (("S1_WMS", "WMS"), ("S1_16S", "16S")):
                path = root / f"{library_id}.tsv"
                integration.cmd_mock_profile(
                    argparse.Namespace(
                        library_id=library_id, assay=assay, output=path
                    )
                )
                rows.extend(read_tsv(path))
            integration.write_tsv(source, rows, integration.PROFILE_FIELDS)
            annotated_fields = integration.PROFILE_FIELDS + [
                "control_abundance_ratio",
                "contamination_probability",
                "associated_controls",
                "decision",
                "decision_reason",
            ]
            for row in rows:
                row.update(
                    control_abundance_ratio="0",
                    contamination_probability="0",
                    associated_controls="",
                    decision="retain",
                    decision_reason="test",
                )
            integration.write_tsv(source, rows, annotated_fields)
            harmonized = root / "harmonized.tsv"
            crosswalk = root / "crosswalk.tsv"
            integration.cmd_harmonize(
                argparse.Namespace(
                    input=source,
                    crosswalk=None,
                    output=harmonized,
                    output_crosswalk=crosswalk,
                )
            )
            calls = read_tsv(harmonized)
            pseudomonas = {
                row["triarchy_taxon_id"]
                for row in calls
                if row["canonical_name"] == "Pseudomonas"
            }
            self.assertEqual(1, len(pseudomonas))

    def test_diversity_emits_alpha_and_beta(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            taxa = root / "taxa.tsv"
            functions = root / "functions.tsv"
            integration.write_tsv(
                taxa,
                [
                    {
                        "library_id": library,
                        "assay": "16S",
                        "feature_id": feature,
                        "original_name": feature,
                        "original_rank": "genus",
                        "abundance": abundance,
                    }
                    for library, feature, abundance in (
                        ("A", "Taxon1", 8),
                        ("A", "Taxon2", 2),
                        ("B", "Taxon1", 2),
                        ("B", "Taxon2", 8),
                    )
                ],
                integration.PROFILE_FIELDS,
            )
            integration.write_tsv(
                functions,
                [],
                ["library_id", "assay", "feature_id", "original_name", "abundance"],
            )
            alpha = root / "alpha.tsv"
            beta = root / "beta.tsv"
            integration.cmd_diversity(
                argparse.Namespace(
                    taxa=taxa,
                    functions=functions,
                    alpha_output=alpha,
                    beta_output=beta,
                )
            )
            self.assertEqual(2, len(read_tsv(alpha)))
            beta_rows = read_tsv(beta)
            self.assertEqual(1, len(beta_rows))
            self.assertEqual("0.6", beta_rows[0]["bray_curtis"])

    def test_qc_summary_captures_retention_and_host_dominance(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            samples = root / "samples.tsv"
            integration.write_tsv(
                samples,
                [{
                    "library_id": "S1_WMS",
                    "biological_sample_id": "S1",
                    "sample_role": "biological",
                    "assay": "WMS",
                }],
                ["library_id", "biological_sample_id", "sample_role", "assay"],
            )
            fastp = root / "S1_WMS.fastp.json"
            fastp.write_text(json.dumps({
                "summary": {
                    "before_filtering": {"total_reads": 10000},
                    "after_filtering": {"total_reads": 9000, "q30_rate": .9},
                }
            }))
            bowtie = root / "S1_WMS.bowtie2.log"
            bowtie.write_text("85.00% overall alignment rate\n")
            output = root / "library_qc.tsv"
            integration.cmd_summarize_qc(
                argparse.Namespace(
                    samples=samples, qc_dirs=[fastp, bowtie], output=output
                )
            )
            row = read_tsv(output)[0]
            self.assertEqual("host_dominated", row["qc_status"])
            self.assertEqual("0.85", row["host_fraction"])
            self.assertAlmostEqual(.1215, float(row["quality_factor"]))

    def test_ranked_evidence_exposes_score_components(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            samples = root / "samples.tsv"
            integration.write_tsv(
                samples,
                [
                    {
                        "library_id": f"S1_{assay}",
                        "biological_sample_id": "S1",
                        "sample_role": "biological",
                    }
                    for assay in ("WMS", "16S")
                ],
                ["library_id", "biological_sample_id", "sample_role"],
            )
            harmonized = root / "harmonized.tsv"
            calls = []
            for assay, relative in (("WMS", 0.6), ("16S", 0.7)):
                calls.append(
                    {
                        "library_id": f"S1_{assay}",
                        "assay": assay,
                        "triarchy_taxon_id": "TRI-PSEUDO",
                        "canonical_name": "Pseudomonas",
                        "canonical_rank": "genus",
                        "canonical_lineage": "Bacteria;Pseudomonas",
                        "native_normalized_value": relative,
                        "abundance_rank": 1,
                        "mapping_confidence": 1,
                        "contamination_probability": 0.05,
                    }
                )
            integration.write_tsv(harmonized, calls, list(calls[0]))
            output = root / "evidence.tsv"
            integration.cmd_classify(
                argparse.Namespace(
                    harmonized=harmonized,
                    samples=samples,
                    assay_performance=None,
                    evidence_policy=None,
                    model_output=root / "model.json",
                    output=output,
                )
            )
            evidence = read_tsv(output)
            self.assertEqual(1, len(evidence))
            row = evidence[0]
            self.assertEqual("concordant_detection", row["classification"])
            self.assertEqual("concordant", row["cross_assay_status"])
            self.assertIn(
                row["authenticity_state"],
                {"biologically_supported", "indeterminate"},
            )
            self.assertEqual("exact_or_canonical_agreement", row["taxonomic_compatibility_status"])
            self.assertEqual("genus", row["comparison_rank"])
            self.assertEqual("1", row["sample_evidence_rank"])
            self.assertGreater(float(row["presence_score"]), 0)
            self.assertGreater(float(row["biological_authenticity"]), 0)
            self.assertGreaterEqual(float(row["authenticity_ci_high"]), float(row["authenticity_ci_low"]))
            self.assertEqual(
                "synthetic_only_not_empirically_calibrated",
                row["calibration_status"],
            )
            self.assertGreater(float(row["wms_evidence_score"]), 0)
            self.assertGreater(float(row["16s_evidence_score"]), 0)
            self.assertTrue((root / "model.json").is_file())

    def test_mismatched_supported_ranks_are_ambiguous(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            samples = root / "samples.tsv"
            integration.write_tsv(
                samples,
                [
                    {
                        "library_id": f"S1_{assay}",
                        "biological_sample_id": "S1",
                        "sample_role": "biological",
                    }
                    for assay in ("WMS", "16S")
                ],
                ["library_id", "biological_sample_id", "sample_role"],
            )
            harmonized = root / "harmonized.tsv"
            calls = [
                {
                    "library_id": f"S1_{assay}",
                    "assay": assay,
                    "triarchy_taxon_id": "TRI-X",
                    "canonical_name": "Pseudomonas",
                    "canonical_rank": rank,
                    "canonical_lineage": "Bacteria;Pseudomonas",
                    "native_normalized_value": 0.5,
                    "abundance_rank": 1,
                    "mapping_confidence": 1,
                    "contamination_probability": 0,
                }
                for assay, rank in (("WMS", "species"), ("16S", "genus"))
            ]
            integration.write_tsv(harmonized, calls, list(calls[0]))
            output = root / "evidence.tsv"
            integration.cmd_classify(
                argparse.Namespace(
                    harmonized=harmonized,
                    samples=samples,
                    assay_performance=None,
                    evidence_policy=None,
                    model_output=None,
                    output=output,
                )
            )
            row = read_tsv(output)[0]
            self.assertEqual("rank_limited_agreement", row["classification"])
            self.assertEqual("taxonomy_limited", row["cross_assay_status"])
            self.assertEqual("hierarchical_agreement", row["taxonomic_compatibility_status"])
            self.assertEqual("genus", row["comparison_rank"])

    def test_evidence_v3_explains_distinct_nondetections_and_contamination(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            samples = root / "samples.tsv"
            sample_rows = []
            designs = {
                "S_REF": ("WMS", "ITS"),
                "S_LOW": ("WMS", "16S"),
                "S_SINGLE": ("16S",),
                "S_CONT": ("WMS", "16S"),
                "S_CROSS": ("WMS", "16S"),
                "S_DISC": ("WMS", "16S"),
            }
            for sample_id, assays in designs.items():
                for assay in assays:
                    sample_rows.append({
                        "library_id": f"{sample_id}_{assay}",
                        "biological_sample_id": sample_id,
                        "sample_role": "biological",
                        "assay": assay,
                    })
            integration.write_tsv(
                samples, sample_rows,
                ["library_id", "biological_sample_id", "sample_role", "assay"],
            )

            calls = []

            def add_call(
                sample_id: str, assay: str, taxon_id: str, name: str,
                contamination: float = 0, batch: float = 0, native: int = 10000,
            ) -> None:
                calls.append({
                    "library_id": f"{sample_id}_{assay}",
                    "assay": assay,
                    "triarchy_taxon_id": taxon_id,
                    "canonical_name": name,
                    "canonical_rank": "genus",
                    "canonical_lineage": (
                        f"Fungi;{name}" if name == "Aspergillus" else f"Bacteria;{name}"
                    ),
                    "native_count": native,
                    "native_normalized_value": .6,
                    "abundance_rank": 1,
                    "mapping_confidence": 1,
                    "contamination_probability": contamination,
                    "batch_localization_score": batch,
                })

            add_call("S_REF", "ITS", "T_REF", "Aspergillus")
            add_call("S_LOW", "16S", "T_LOW", "Pseudomonas")
            add_call("S_SINGLE", "16S", "T_SINGLE", "Pseudomonas")
            add_call("S_CONT", "WMS", "T_CONT", "Pseudomonas", .8)
            add_call("S_CONT", "16S", "T_CONT", "Pseudomonas", .8)
            add_call("S_CROSS", "WMS", "T_CROSS", "Pseudomonas", .5, .8)
            add_call("S_CROSS", "16S", "T_CROSS", "Pseudomonas", .5, .8)
            add_call("S_DISC", "16S", "T_DISC", "Pseudomonas")
            add_call("S_DISC", "WMS", "T_BACKGROUND", "Bacillus", native=1_000_000)
            harmonized = root / "harmonized.tsv"
            integration.write_tsv(harmonized, calls, list(calls[0]))
            output = root / "evidence.tsv"
            integration.cmd_classify(
                argparse.Namespace(
                    harmonized=harmonized,
                    samples=samples,
                    assay_performance=None,
                    evidence_policy=None,
                    model_output=None,
                    output=output,
                )
            )
            observed = {
                (row["biological_sample_id"], row["triarchy_taxon_id"]): row
                for row in read_tsv(output)
            }
            self.assertEqual(
                "reference_limited", observed[("S_REF", "T_REF")]["classification"]
            )
            self.assertEqual(
                "assay_limited_nondetection",
                observed[("S_LOW", "T_LOW")]["classification"],
            )
            self.assertEqual(
                "single_assay_credible",
                observed[("S_SINGLE", "T_SINGLE")]["classification"],
            )
            self.assertEqual(
                "probable_contaminant",
                observed[("S_CONT", "T_CONT")]["classification"],
            )
            self.assertEqual(
                "contamination_associated",
                observed[("S_CONT", "T_CONT")]["cross_assay_status"],
            )
            self.assertEqual(
                "cross_contamination_candidate",
                observed[("S_CROSS", "T_CROSS")]["classification"],
            )
            self.assertEqual("batch", observed[("S_CROSS", "T_CROSS")]["probable_contamination_source"])
            self.assertEqual(
                "genuine_discordance",
                observed[("S_DISC", "T_DISC")]["classification"],
            )
            self.assertEqual("true", observed[("S_REF", "T_REF")]["wms_measured"])
            self.assertEqual("false", observed[("S_REF", "T_REF")]["wms_detected"])

    def test_cross_assay_concordance_reports_presence_and_rank(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            samples = root / "samples.tsv"
            integration.write_tsv(
                samples,
                [
                    {
                        "library_id": f"S1_{assay}",
                        "biological_sample_id": "S1",
                        "sample_role": "biological",
                    }
                    for assay in ("WMS", "16S")
                ],
                ["library_id", "biological_sample_id", "sample_role"],
            )
            harmonized = root / "harmonized.tsv"
            calls = [
                {
                    "library_id": f"S1_{assay}",
                    "assay": assay,
                    "triarchy_taxon_id": taxon,
                    "canonical_name": taxon,
                    "canonical_lineage": "Bacteria",
                    "native_normalized_value": abundance,
                }
                for assay, taxon, abundance in (
                    ("WMS", "A", 0.7),
                    ("WMS", "B", 0.3),
                    ("16S", "A", 0.6),
                    ("16S", "C", 0.4),
                )
            ]
            integration.write_tsv(harmonized, calls, list(calls[0]))
            output = root / "concordance.tsv"
            integration.cmd_concordance(
                argparse.Namespace(
                    harmonized=harmonized, samples=samples, output=output
                )
            )
            bacterial = next(
                row for row in read_tsv(output) if row["kingdom"] == "Bacteria"
            )
            self.assertEqual("0.333333", bacterial["jaccard_presence"])
            self.assertEqual("1", bacterial["rank_agreement_score"])

    def test_calibration_and_truth_evaluation_are_serialized(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            evidence = root / "evidence.tsv"
            truth = root / "truth.tsv"
            rows = [
                {
                    "biological_sample_id": f"S{index}",
                    "canonical_name": f"T{index}",
                    "presence_score": score,
                    "probability_authentic": probability,
                    "biological_state": state,
                }
                for index, score, probability, state in (
                    (1, .1, .05, "absent"),
                    (2, .2, .10, "absent"),
                    (3, .8, .90, "present"),
                    (4, .9, .95, "present"),
                )
            ]
            integration.write_tsv(evidence, rows, list(rows[0]))
            truth_rows = [
                {
                    "biological_sample_id": f"S{index}",
                    "taxon": f"T{index}",
                    "truth_status": status,
                }
                for index, status in (
                    (1, "contaminant"), (2, "absent"),
                    (3, "authentic"), (4, "present"),
                )
            ]
            integration.write_tsv(
                truth, truth_rows,
                ["biological_sample_id", "taxon", "truth_status"],
            )
            calibration = root / "calibration.json"
            metrics = root / "metrics.json"
            integration.cmd_fit_calibration(
                argparse.Namespace(evidence=evidence, truth=truth, output=calibration)
            )
            integration.cmd_evaluate_evidence(
                argparse.Namespace(evidence=evidence, truth=truth, output=metrics)
            )
            model = json.loads(calibration.read_text())
            self.assertEqual("isotonic", model["calibration"]["method"])
            self.assertEqual(4, model["calibration"]["training_observations"])
            observed = json.loads(metrics.read_text())
            self.assertEqual(1, observed["precision"])
            self.assertEqual(1, observed["recall"])
            self.assertEqual(1, observed["f1"])


if __name__ == "__main__":
    unittest.main()

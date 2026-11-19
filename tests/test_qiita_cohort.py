"""Tests for the Qiita-curated extraction-plate cohort and mock evaluation."""

from __future__ import annotations

import json
import zlib
from pathlib import Path

import pytest

import qiita_cohort
from conftest import write_tsv
from evaluate_mock_truth import main as evaluate_main
from validate_metadata import CONTROL_FIELDS, SAMPLE_FIELDS, validate

ASSAYS = ("WMS", "16S", "ITS")


def qiita(name: str, sample_type: str = "soil", plate: str = "high", kit: str = "PowerSoil r1", **extra) -> dict:
    row = {
        "sample_name": f"12201.{name}",
        "round": "1",
        "extraction_kit_round": kit,
        "biomass_plate": plate,
        "sample_type": sample_type,
        "sample_type_2": extra.pop("sample_type_2", sample_type),
        "empo_3": extra.pop("empo_3", "Soil (non-saline)"),
        "host_subject_id": extra.pop("host_subject_id", "not applicable"),
        "sample_technical_replicate": extra.pop("rep", "1"),
        "gdna_well_id": extra.pop("well", "A1"),
    }
    row.update(extra)
    return row


def runs_for(name: str, assays=ASSAYS, start: int = 0) -> list[dict]:
    rows = []
    for index, assay in enumerate(assays):
        accession = f"ERR{zlib.crc32(f'{name}:{assay}'.encode()) % 10_000_000:07d}"
        paired = assay == "WMS"
        files = f"ftp.sra.ebi.ac.uk/x/{accession}_1.fastq.gz;ftp.sra.ebi.ac.uk/x/{accession}_2.fastq.gz" if paired \
            else f"ftp.sra.ebi.ac.uk/x/{accession}.fastq.gz"
        rows.append({
            "run_accession": accession,
            "sample_title": f"12201.{name}",
            "experiment_title": f"Illumina sequencing: qiita_ptid_{1000 + index}:12201.{name}",
            "triarchy_assay": assay,
            "fastq_ftp": files,
        })
    return rows


def fastq_paths(run: dict) -> tuple[str, str, str]:
    files = run["fastq_ftp"].split(";")
    names = [f"../reads/{run['run_accession']}/{Path(f).name}" for f in files]
    return ("paired", names[0], names[1]) if len(names) == 2 else ("single", names[0], "")


@pytest.fixture
def cohort():
    samples = [
        qiita("soil.A.1.standard.H"),
        qiita("soil.A.2.standard.H", rep="2"),
        qiita("feces.human.B.1.standard.H", sample_type="feces", sample_type_2="human feces",
              host_subject_id="human.B"),
        qiita("BLANK.C4.standard.H", sample_type="control blank", well="C4"),
        qiita("mock.community.d1.standard.H", sample_type="control positive",
              sample_type_2="ZymoBIOMICS Microbial Community Standard II"),
        qiita("soil.Z.1.standard.L", plate="low"),
        qiita("BLANK.B4.standard.L", sample_type="control blank", plate="low", well="B4"),
        qiita("soil.X.1.pro.H", kit="PowerSoil Pro"),
        qiita("soil.incomplete.1.standard.H"),
    ]
    runs = []
    for sample in samples:
        name = sample["sample_name"].removeprefix("12201.")
        runs += runs_for(name, ("WMS", "16S") if "incomplete" in name else ASSAYS)
    return samples, runs


def build(cohort, kits=("PowerSoil r1",), plates=("high", "low"), cap=None):
    samples, runs = cohort
    selected, selected_runs, summary = qiita_cohort.select_cohort(samples, runs, list(kits), list(plates), cap)
    return selected, summary, qiita_cohort.build_cohort_metadata(selected, selected_runs, fastq_paths)


def test_selects_complete_samples_on_requested_plates(cohort) -> None:
    selected, summary, _ = build(cohort)
    names = {s["sample_name"] for s in selected}
    assert "12201.soil.incomplete.1.standard.H" not in names
    assert "12201.soil.X.1.pro.H" not in names
    assert summary["plates"]["R1_PowerSoil_r1_high"] == {
        "extraction_blanks": 1, "zymo_mock_samples": 1, "biological_samples": 3,
    }
    assert summary["libraries"] == len(selected) * 3


def test_metadata_passes_triarchy_validation(cohort, tmp_path: Path) -> None:
    _, _, (samples, controls, _) = build(cohort)
    samples_tsv = write_tsv(tmp_path / "samples.tsv", samples, sorted(SAMPLE_FIELDS))
    controls_tsv = write_tsv(tmp_path / "controls.tsv", controls, sorted(CONTROL_FIELDS))
    libraries, edges = validate(samples_tsv, controls_tsv)
    assert libraries == len(samples)
    assert edges == len(controls)


def test_blanks_link_only_within_their_plate(cohort) -> None:
    _, _, (samples, controls, _) = build(cohort)
    batch = {row["library_id"]: row["extraction_batch"] for row in samples}
    assert controls
    for edge in controls:
        assert batch[edge["library_id"]] == batch[edge["control_library_id"]]
        assert edge["library_id"].rsplit("_", 1)[1] == edge["control_library_id"].rsplit("_", 1)[1]
        assert edge["control_stage"] == "extraction"
    high_targets = {e["library_id"] for e in controls if "C4" in e["control_library_id"]}
    assert "mock.community.d1.standard.H_WMS" in high_targets
    assert not any("standard.L" in target for target in high_targets)


def test_mock_is_evaluated_as_biological_with_truth(cohort) -> None:
    _, _, (samples, _, truth) = build(cohort)
    mock = [row for row in samples if row["biological_sample_id"] == "mock.community.d1.standard.H"]
    assert {row["sample_role"] for row in mock} == {"biological"}
    assert {row["subject_id"] for row in mock} == {"ZYMO_D6310"}
    assert len(truth) == 10
    listeria = next(row for row in truth if row["genus"] == "Listeria")
    assert listeria["genomic_dna_pct"] == "89.1"
    yeast = next(row for row in truth if row["genus"] == "Saccharomyces")
    assert yeast["expected_assays"] == "WMS;ITS"
    assert "Limosilactobacillus" in next(r for r in truth if r["genus"] == "Lactobacillus")["genus_aliases"]


@pytest.mark.parametrize("standard", ["D6300", "D6310"])
def test_zymo_genomic_dna_sums_to_100(standard: str) -> None:
    # Values are transcribed as printed; Zymo rounds D6310 (sum 99.97%).
    total = sum(float(row[2]) for row in qiita_cohort.ZYMO_COMPOSITION[standard])
    assert total == pytest.approx(100, abs=0.05)


def test_cap_keeps_every_sample_type_and_all_controls(cohort) -> None:
    selected, summary, _ = build(cohort, plates=("high",), cap=2)
    names = {s["sample_name"] for s in selected}
    assert "12201.BLANK.C4.standard.H" in names
    assert "12201.mock.community.d1.standard.H" in names
    assert summary["plates"]["R1_PowerSoil_r1_high"]["biological_samples"] == 2
    types = {s["sample_type"] for s in selected if s["sample_type"] in {"soil", "feces"}}
    assert types == {"soil", "feces"}


def test_rejects_unknown_kit_and_plate_without_blank(cohort) -> None:
    with pytest.raises(qiita_cohort.CohortError, match="unknown extraction kit"):
        build(cohort, kits=("NotAKit",))
    with pytest.raises(qiita_cohort.CohortError, match="no complete extraction blank"):
        build(cohort, kits=("PowerSoil Pro",))


def test_mock_evaluator_scores_detection_and_false_support(tmp_path: Path) -> None:
    sample = "mock.d1"
    truth = [
        {"biological_sample_id": sample, "standard": "D6310", "dilution": "d1", "technical_replicate": "1",
         "species": species, "genus": species.split()[0],
         "genus_aliases": ";".join(qiita_cohort.GENUS_ALIASES[species.split()[0]]),
         "kingdom": kingdom, "expected_assays": "WMS;ITS" if kingdom == "Fungi" else "WMS;16S",
         "genomic_dna_pct": gdna, "rrna_16s_pct": "", "rrna_16s_18s_pct": "",
         "genome_copy_pct": "", "cell_number_pct": "", "source": ""}
        for species, kingdom, gdna, *_ in qiita_cohort.ZYMO_COMPOSITION["D6310"][:4]
    ]
    truth_tsv = write_tsv(tmp_path / "truth.tsv", truth, qiita_cohort.TRUTH_FIELDS)
    samples_tsv = write_tsv(tmp_path / "samples.tsv", [
        {"library_id": f"{sample}_{a}", "biological_sample_id": sample, "assay": a} for a in ASSAYS
    ], ["library_id", "biological_sample_id", "assay"])
    integrated = tmp_path / "results" / "integrated"
    integrated.mkdir(parents=True)
    write_tsv(integrated / "harmonized_taxa.tsv", [
        {"library_id": f"{sample}_WMS", "canonical_name": "Listeria monocytogenes", "abundance": "90"},
        {"library_id": f"{sample}_WMS", "canonical_name": "Pseudomonas", "abundance": "10"},
        {"library_id": f"{sample}_16S", "canonical_name": "Listeria", "abundance": "100"},
        {"library_id": f"{sample}_ITS", "canonical_name": "Saccharomyces", "abundance": "5"},
    ], ["library_id", "canonical_name", "abundance"])
    evidence_fields = ["biological_sample_id", "canonical_name", "comparison_rank", "wms_detected",
                       "16s_detected", "its_detected", "classification", "cross_assay_status",
                       "authenticity_state", "biological_authenticity"]
    write_tsv(integrated / "taxon_evidence.tsv", [
        {"biological_sample_id": sample, "canonical_name": "Listeria", "comparison_rank": "genus",
         "wms_detected": "true", "16s_detected": "true", "its_detected": "false",
         "classification": "concordant_detection", "cross_assay_status": "concordant",
         "authenticity_state": "biologically_supported", "biological_authenticity": "0.95"},
        {"biological_sample_id": sample, "canonical_name": "Ralstonia", "comparison_rank": "genus",
         "wms_detected": "false", "16s_detected": "true", "its_detected": "false",
         "classification": "probable_contaminant", "cross_assay_status": "contamination_associated",
         "authenticity_state": "not_supported", "biological_authenticity": "0.05"},
    ], evidence_fields)

    outdir = tmp_path / "eval"
    assert evaluate_main(["--results", str(tmp_path / "results"), "--truth", str(truth_tsv),
                          "--samples", str(samples_tsv), "--outdir", str(outdir)]) == 0
    summary = json.loads((outdir / "mock_evaluation.json").read_text())
    # WMS expects Listeria, Pseudomonas, Bacillus, Saccharomyces; two were seen.
    assert (summary["assay_recall"]["WMS"]["expected"], summary["assay_recall"]["WMS"]["detected"]) == (4, 2)
    # 16S expects the three bacteria; ITS expects only the yeast.
    assert (summary["assay_recall"]["16S"]["expected"], summary["assay_recall"]["16S"]["detected"]) == (3, 1)
    assert (summary["assay_recall"]["ITS"]["expected"], summary["assay_recall"]["ITS"]["detected"]) == (1, 1)
    assert summary["assay_recall"]["WMS"]["by_genomic_dna_abundance"][">=1%"]["recall"] == 1.0
    assert summary["unexpected_taxa_in_mocks"]["evidence_rows"] == 1
    assert summary["unexpected_taxa_in_mocks"]["false_support_rate"] == 0
    assert summary["expected_taxa"]["biologically_supported"] == 1

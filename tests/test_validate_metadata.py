"""Tests for the sample manifest and control-graph contract."""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import write_tsv
from validate_metadata import (
    CONTROL_FIELDS,
    SAMPLE_FIELDS,
    MetadataError,
    validate,
    validate_controls,
    validate_samples,
)

SAMPLE_ORDER = sorted(SAMPLE_FIELDS)
CONTROL_ORDER = sorted(CONTROL_FIELDS)


def library(library_id: str, **overrides: str) -> dict[str, str]:
    row = {
        "library_id": library_id,
        "biological_sample_id": f"S_{library_id}",
        "specimen_id": f"SP_{library_id}",
        "subject_id": "SUBJ1",
        "assay": "16S",
        "sample_role": "biological",
        "sample_material": "soil",
        "environment_category": "terrestrial",
        "library_layout": "paired",
        "fastq_r1": f"/data/{library_id}_R1.fastq.gz",
        "fastq_r2": f"/data/{library_id}_R2.fastq.gz",
        "extraction_batch": "EX1",
        "pcr_plate": "P1",
        "library_batch": "LB1",
        "sequencing_run": "RUN1",
        "technical_replicate": "1",
    }
    row.update(overrides)
    return row


def edge(source: str, control: str, **overrides: str) -> dict[str, str]:
    row = {
        "library_id": source,
        "control_library_id": control,
        "control_stage": "extraction",
        "relationship": "same_batch",
        "priority": "1",
        "reason": "extraction blank",
    }
    row.update(overrides)
    return row


def blank(library_id: str, **overrides: str) -> dict[str, str]:
    return library(library_id, sample_role="negative_control", **overrides)


def expect_code(code: str, func, *args) -> None:
    with pytest.raises(MetadataError) as caught:
        func(*args)
    assert caught.value.code == code


def test_valid_study_round_trip(tmp_path: Path) -> None:
    samples = write_tsv(tmp_path / "samples.tsv", [library("L1"), blank("B1")], SAMPLE_ORDER)
    controls = write_tsv(tmp_path / "controls.tsv", [edge("L1", "B1")], CONTROL_ORDER)
    assert validate(samples, controls) == (2, 1)


def test_missing_column_is_reported(tmp_path: Path) -> None:
    fields = [field for field in SAMPLE_ORDER if field != "assay"]
    rows = [{key: value for key, value in library("L1").items() if key != "assay"}]
    samples = write_tsv(tmp_path / "samples.tsv", rows, fields)
    controls = write_tsv(tmp_path / "controls.tsv", [], CONTROL_ORDER)
    expect_code("TRI-META-000", validate, samples, controls)


@pytest.mark.parametrize(
    ("rows", "code"),
    [
        ([library("L1"), library("L1", fastq_r1="/x1", fastq_r2="/x2")], "TRI-META-001"),
        ([library("")], "TRI-META-002"),
        ([library("L1", assay="RNA")], "TRI-META-003"),
        ([library("L1", fastq_r2="")], "TRI-META-004"),
        ([library("L1"), library("L2", fastq_r1="/data/L1_R1.fastq.gz")], "TRI-META-005"),
        ([library("L1", sample_role="unknown")], "TRI-META-006"),
        ([library("L1", library_layout="interleaved")], "TRI-META-007"),
        ([library("L1", sequencing_run="")], "TRI-META-008"),
        ([library("L1", pcr_plate="")], "TRI-META-008"),
        (
            [library("L1"), library("L2", biological_sample_id="S_L1", subject_id="OTHER")],
            "TRI-META-009",
        ),
        (
            [library("L1"), library("L2", biological_sample_id="S_L1", specimen_id="SP_L1",
                                     sample_role="mock")],
            "TRI-META-016",
        ),
    ],
)
def test_sample_violations(rows: list[dict[str, str]], code: str) -> None:
    expect_code(code, validate_samples, rows)


def test_wms_library_does_not_need_pcr_plate() -> None:
    libraries = validate_samples([library("W1", assay="WMS", pcr_plate="")])
    assert "W1" in libraries


@pytest.mark.parametrize(
    ("libraries", "edges", "code"),
    [
        ([library("L1"), blank("B1")], [edge("L9", "B1")], "TRI-META-010"),
        ([library("L1"), blank("B1")], [edge("L1", "B9")], "TRI-META-010"),
        ([library("L1"), blank("B1")], [edge("L1", "B1", control_stage="storage")], "TRI-META-012"),
        ([library("L1"), blank("B1")], [edge("L1", "B1", relationship="nearby")], "TRI-META-012"),
        ([library("L1"), blank("B1")], [edge("L1", "B1", priority="0")], "TRI-META-012"),
        ([library("L1"), blank("B1")], [edge("L1", "B1", priority="high")], "TRI-META-012"),
        ([library("L1"), library("L2")], [edge("L1", "L2")], "TRI-META-013"),
        ([library("L1"), blank("B1", assay="ITS")], [edge("L1", "B1")], "TRI-META-013"),
        (
            [library("W1", assay="WMS"), blank("B1", assay="WMS")],
            [edge("W1", "B1", control_stage="pcr")],
            "TRI-META-011",
        ),
        ([library("L1"), blank("B1", extraction_batch="EX2")], [edge("L1", "B1")], "TRI-META-014"),
        (
            [blank("B1"), blank("B2")],
            [edge("B1", "B2"), edge("B2", "B1")],
            "TRI-META-015",
        ),
    ],
)
def test_control_violations(
    libraries: list[dict[str, str]], edges: list[dict[str, str]], code: str
) -> None:
    expect_code(code, validate_controls, edges, validate_samples(libraries))


def test_run_level_control_may_cross_batches() -> None:
    libraries = validate_samples([library("L1"), blank("B1", extraction_batch="EX2")])
    validate_controls([edge("L1", "B1", relationship="run_level")], libraries)

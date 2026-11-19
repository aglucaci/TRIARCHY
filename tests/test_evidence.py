"""Tests for the TRIARCHY evidence model building blocks."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from conftest import write_tsv
from triarchy_evidence import (
    authenticity_interval,
    authenticity_state,
    calibrated_probability,
    classify_reason,
    expected_assays,
    fit_calibration,
    fit_isotonic,
    infer_kingdom,
    library_depths,
    load_performance,
    load_policy,
    parent_cross_assay_status,
    sigmoid,
    strongest_contamination_source,
)


@pytest.fixture
def policy(root: Path) -> dict:
    return load_policy(root / "conf" / "evidence_policy.json")


def component(detected=False, measured=True, rank=7, reference=0.9, power=0.9) -> dict:
    return {
        "detected": detected,
        "measured": measured,
        "rank_level": rank,
        "reference": reference,
        "power": power,
    }


@pytest.mark.parametrize(
    ("row", "kingdom"),
    [
        ({"canonical_name": "Escherichia coli"}, "Bacteria"),
        ({"canonical_name": "Candida albicans"}, "Fungi"),
        ({"kingdom": "Fungi", "canonical_name": "Unknown"}, "Fungi"),
        ({"canonical_lineage": "Archaea;Euryarchaeota"}, "Archaea"),
        ({"canonical_name": "Methanobrevibacter smithii"}, "Archaea"),
    ],
)
def test_infer_kingdom(row: dict, kingdom: str) -> None:
    assert infer_kingdom(row) == kingdom


def test_expected_assays_pair_by_kingdom() -> None:
    assert expected_assays("Fungi") == ("WMS", "ITS")
    assert expected_assays("Bacteria") == ("WMS", "16S")
    assert expected_assays("Archaea") == ("WMS", "16S")


def test_sigmoid_is_stable_and_symmetric() -> None:
    assert sigmoid(0) == pytest.approx(0.5)
    assert sigmoid(1000) == pytest.approx(1.0)
    assert sigmoid(-1000) == pytest.approx(0.0)
    assert sigmoid(2) + sigmoid(-2) == pytest.approx(1.0)


def test_calibration_interpolates_and_clamps(policy: dict) -> None:
    calibration = policy["calibration"]
    assert calibrated_probability(0.0, calibration)[0] == pytest.approx(0.02)
    assert calibrated_probability(0.5, calibration)[0] == pytest.approx((0.25 + 0.62) / 2)
    assert calibrated_probability(1.5, calibration)[0] == pytest.approx(0.98)
    assert calibrated_probability(0.5, None) == (0.5, "uncalibrated")
    bad = {"breakpoints": [0, 1], "values": [0.1]}
    assert calibrated_probability(0.3, bad) == (0.3, "invalid_calibration")


def test_calibration_is_monotonic(policy: dict) -> None:
    grid = [i / 50 for i in range(51)]
    values = [calibrated_probability(x, policy["calibration"])[0] for x in grid]
    assert values == sorted(values)


def test_policy_file_overrides_defaults(tmp_path: Path) -> None:
    custom = tmp_path / "policy.json"
    custom.write_text(json.dumps({"contaminant_threshold": 0.9}), encoding="utf-8")
    policy = load_policy(custom)
    assert policy["contaminant_threshold"] == 0.9
    assert "present_score_threshold" in policy


def test_performance_table_loads_shipped_values(root: Path) -> None:
    performance = load_performance(root / "conf" / "assay_performance.tsv")
    assert performance[("ITS", "Fungi")]["reliability"] == pytest.approx(0.90)
    assert performance[("WMS", "Fungi")]["minimum_depth"] == pytest.approx(100000)


def test_interval_narrows_with_evidence() -> None:
    _, _, wide = authenticity_interval(0.7, 0)
    _, _, narrow = authenticity_interval(0.7, 50)
    assert narrow < wide
    low, high, _ = authenticity_interval(0.99, 0)
    assert 0.0 <= low <= high <= 1.0


def test_authenticity_state(policy: dict) -> None:
    assert authenticity_state(0.95, 0.85, 0.99, policy) == "biologically_supported"
    assert authenticity_state(0.05, 0.01, 0.15, policy) == "not_supported"
    assert authenticity_state(0.95, 0.5, 0.99, policy) == "indeterminate"


def test_strongest_contamination_source() -> None:
    assert strongest_contamination_source({"pcr": 0.4, "field": 0.1}, 0.2) == "pcr"
    assert strongest_contamination_source({"pcr": 0.1}, 0.6) == "batch"
    assert strongest_contamination_source({"pcr": 0.0}, 0.0) == "none"


def test_library_depths_flags_missing_native_counts() -> None:
    rows = [
        {"library_id": "L1", "native_count": "10"},
        {"library_id": "L1", "native_count": "5"},
        {"library_id": "L2", "native_count": ""},
    ]
    depths = library_depths(rows)
    assert depths["L1"] == 15
    assert depths["L2"] == -1


EXPECTED = ("WMS", "16S")


@pytest.mark.parametrize(
    ("components", "contamination", "batch", "calls", "classification"),
    [
        ({"WMS": component(True), "16S": component(True)}, 0.0, 0.0, [], "concordant_detection"),
        ({"WMS": component(True, rank=7), "16S": component(True, rank=6)}, 0.0, 0.0, [],
         "rank_limited_agreement"),
        ({"WMS": component(True), "16S": component(True)}, 0.7, 0.0, [], "probable_contaminant"),
        ({"WMS": component(True), "16S": component(True)}, 0.5, 0.7, [],
         "cross_contamination_candidate"),
        ({"WMS": component(True), "16S": component(True)}, 0.0, 0.0,
         [{"cross_assay_compatibility": "conflict"}], "genuine_discordance"),
        ({"WMS": component(True), "16S": component(measured=False)}, 0.0, 0.0, [],
         "single_assay_credible"),
        ({"WMS": component(True), "16S": component(reference=0.3)}, 0.0, 0.0, [],
         "reference_limited"),
        ({"WMS": component(True), "16S": component(power=0.3)}, 0.0, 0.0, [],
         "assay_limited_nondetection"),
        ({"WMS": component(True), "16S": component(power=0.95)}, 0.0, 0.0, [],
         "genuine_discordance"),
        ({"WMS": component(True), "16S": component(power=0.6)}, 0.0, 0.0, [],
         "single_assay_credible"),
        ({"WMS": component(), "16S": component()}, 0.0, 0.0, [], "indeterminate"),
    ],
)
def test_classify_reason(components, contamination, batch, calls, classification, policy) -> None:
    result, reason = classify_reason(calls, components, EXPECTED, contamination, batch, policy)
    assert result == classification
    assert reason


@pytest.mark.parametrize(
    ("classification", "status"),
    [
        ("concordant_detection", "concordant"),
        ("assay_limited_nondetection", "assay_limited"),
        ("reference_limited", "taxonomy_limited"),
        ("probable_contaminant", "contamination_associated"),
        ("genuine_discordance", "genuine_discordance"),
        ("anything_else", "indeterminate"),
    ],
)
def test_parent_cross_assay_status(classification: str, status: str) -> None:
    assert parent_cross_assay_status(classification) == status


def test_isotonic_fit_is_monotonic() -> None:
    model = fit_isotonic([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], [0, 1, 0, 1, 1, 1])
    assert model["values"] == sorted(model["values"])
    assert model["training_observations"] == 6


def test_fit_calibration_requires_both_labels(tmp_path: Path) -> None:
    fields = ["biological_sample_id", "canonical_name", "biological_authenticity"]
    evidence = write_tsv(tmp_path / "evidence.tsv", [
        {"biological_sample_id": "S1", "canonical_name": "A", "biological_authenticity": "0.9"},
    ], fields)
    truth = write_tsv(tmp_path / "truth.tsv", [
        {"biological_sample_id": "S1", "canonical_name": "A", "truth_status": "present"},
    ], ["biological_sample_id", "canonical_name", "truth_status"])
    with pytest.raises(ValueError):
        fit_calibration(evidence, truth, tmp_path / "out.json")

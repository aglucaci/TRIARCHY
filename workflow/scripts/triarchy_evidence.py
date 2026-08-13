#!/usr/bin/env python3
"""Formal taxon-evidence inference, calibration, and evaluation for TRIARCHY."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


ASSAYS = ("WMS", "16S", "ITS")
RANK_LEVEL = {
    "unclassified": 0, "root": 0, "kingdom": 1, "domain": 1, "phylum": 2,
    "class": 3, "order": 4, "family": 5, "genus": 6, "species": 7,
    "strain": 8, "asv": 8,
}
LEVEL_RANK = {
    0: "unclassified", 1: "kingdom", 2: "phylum", 3: "class", 4: "order",
    5: "family", 6: "genus", 7: "species", 8: "strain",
}
DEFAULT_PERFORMANCE = {
    ("WMS", "Bacteria"): dict(reliability=.85, minimum_depth=10_000, reference_coverage=.95, false_positive_rate=.01),
    ("WMS", "Archaea"): dict(reliability=.80, minimum_depth=20_000, reference_coverage=.80, false_positive_rate=.02),
    ("WMS", "Fungi"): dict(reliability=.60, minimum_depth=100_000, reference_coverage=.45, false_positive_rate=.03),
    ("16S", "Bacteria"): dict(reliability=.90, minimum_depth=5_000, reference_coverage=.95, false_positive_rate=.01),
    ("16S", "Archaea"): dict(reliability=.85, minimum_depth=10_000, reference_coverage=.80, false_positive_rate=.02),
    ("ITS", "Fungi"): dict(reliability=.90, minimum_depth=5_000, reference_coverage=.90, false_positive_rate=.01),
}
DEFAULT_POLICY: dict[str, Any] = {
    "model_version": "triarchy-authenticity-v1.0",
    "prior_presence": .25,
    "contaminant_threshold": .65,
    "cross_contamination_threshold": .45,
    "batch_localization_threshold": .65,
    "low_detection_power": .50,
    "adequate_detection_power": .80,
    "minimum_reference_coverage": .55,
    "minimum_mapping_confidence": .50,
    "present_score_threshold": .75,
    "absent_score_threshold": .20,
    "authenticity_supported_threshold": .80,
    "authenticity_absent_threshold": .20,
    "indeterminate_ci_overlap_low": .20,
    "indeterminate_ci_overlap_high": .80,
    "calibration": {
        "method": "piecewise_linear",
        "status": "synthetic_only_not_empirically_calibrated",
        "breakpoints": [0, .2, .4, .6, .8, 1],
        "values": [.02, .08, .25, .62, .88, .98],
    },
}
EVIDENCE_FIELDS = [
    "biological_sample_id", "triarchy_taxon_id", "canonical_name", "kingdom",
    "comparison_rank", "lowest_accurately_resolved_rank", "assay_support",
    "wms_measured", "wms_detected", "wms_library_depth", "wms_detection_power",
    "wms_reference_coverage", "wms_evidence_score",
    "16s_measured", "16s_detected", "16s_library_depth", "16s_detection_power",
    "16s_reference_coverage", "16s_evidence_score",
    "its_measured", "its_detected", "its_library_depth", "its_detection_power",
    "its_reference_coverage", "its_evidence_score",
    "abundance_support", "taxonomic_resolution_support", "mapping_support",
    "assay_reliability", "control_evidence", "field_contamination_score",
    "extraction_contamination_score", "pcr_contamination_score",
    "library_contamination_score", "run_contamination_score",
    "other_control_contamination_score", "probable_contamination_source",
    "contamination_probability", "batch_localization_score",
    "taxonomic_compatibility_score", "taxonomic_compatibility_status",
    "biological_authenticity", "authenticity_ci_low", "authenticity_ci_high",
    "authenticity_uncertainty", "authenticity_state", "cross_assay_status",
    "presence_score", "probability_authentic", "calibration_status",
    "biological_state", "sample_evidence_rank", "global_evidence_rank",
    "classification", "discordance_reason", "decision_reason", "model_version",
]


def clamp(value: float) -> float:
    return max(0.0, min(1.0, value))


def number(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def flag(value: Any) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes", "pass", "present"}


def fmt(value: float) -> str:
    return f"{value:.6g}"


def sigmoid(value: float) -> float:
    if value >= 0:
        z = math.exp(-value)
        return 1 / (1 + z)
    z = math.exp(value)
    return z / (1 + z)


def infer_kingdom(row: dict[str, str]) -> str:
    text = " ".join(
        row.get(key, "") for key in ("kingdom", "canonical_lineage", "canonical_name")
    ).casefold()
    if "fung" in text or any(word in text for word in ("aspergillus", "candida", "cladospor", "saccharomyces", "cryptococcus")):
        return "Fungi"
    if "archaea" in text or "methan" in text or "halobacter" in text:
        return "Archaea"
    return "Bacteria"


def expected_assays(kingdom: str) -> tuple[str, ...]:
    return ("WMS", "ITS") if kingdom == "Fungi" else ("WMS", "16S")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_performance(path: Path | None) -> dict[tuple[str, str], dict[str, float]]:
    result = {key: dict(values) for key, values in DEFAULT_PERFORMANCE.items()}
    if not path:
        return result
    for row in read_tsv(path):
        key = (row["assay"], row["kingdom"])
        values = result.setdefault(key, {})
        for field in ("reliability", "minimum_depth", "reference_coverage", "false_positive_rate"):
            if row.get(field, ""):
                values[field] = number(row[field])
    return result


def load_policy(path: Path | None) -> dict[str, Any]:
    policy = json.loads(json.dumps(DEFAULT_POLICY))
    if path:
        supplied = json.loads(path.read_text(encoding="utf-8"))
        policy.update({key: value for key, value in supplied.items() if key != "calibration"})
        if "calibration" in supplied:
            policy["calibration"] = supplied["calibration"]
    return policy


def calibrated_probability(score: float, calibration: dict[str, Any] | None) -> tuple[float, str]:
    if not calibration:
        return score, "uncalibrated"
    breaks = [number(value) for value in calibration.get("breakpoints", [])]
    values = [number(value) for value in calibration.get("values", [])]
    if len(breaks) != len(values) or not breaks:
        return score, "invalid_calibration"
    if score <= breaks[0]:
        return clamp(values[0]), str(calibration.get("status", "unspecified"))
    for left, right, low, high in zip(breaks, breaks[1:], values, values[1:]):
        if score <= right:
            fraction = (score - left) / max(right - left, 1e-12)
            return clamp(low + fraction * (high - low)), str(calibration.get("status", "unspecified"))
    return clamp(values[-1]), str(calibration.get("status", "unspecified"))


def stage_contamination_scores(calls: list[dict[str, str]]) -> dict[str, float]:
    fields = {
        "field": "field_contamination_score",
        "extraction": "extraction_contamination_score",
        "pcr": "pcr_contamination_score",
        "library": "library_contamination_score",
        "run": "run_contamination_score",
        "other": "other_control_contamination_score",
    }
    return {
        stage: max((number(row.get(field)) for row in calls), default=0.0)
        for stage, field in fields.items()
    }


def strongest_contamination_source(stage_scores: dict[str, float], batch_score: float) -> str:
    candidates = dict(stage_scores)
    candidates["batch"] = batch_score
    stage, score = max(candidates.items(), key=lambda item: item[1], default=("none", 0.0))
    return stage if score > 0 else "none"


def taxonomic_compatibility(
    calls: list[dict[str, str]],
    components: dict[str, dict[str, float | bool]],
    expected: tuple[str, ...],
) -> tuple[float, str]:
    expected_calls = [row for row in calls if row.get("assay") in expected]
    if not expected_calls:
        return 0.0, "no_comparable_assignments"
    if any(
        row.get("cross_assay_compatibility", "").casefold() == "conflict"
        or row.get("mapping_method", "").casefold() in {"incompatible", "conflict"}
        for row in expected_calls
    ):
        return 0.0, "incompatible"
    detected_expected = [assay for assay in expected if components[assay]["detected"]]
    if len(detected_expected) < 2:
        return 0.5, "single_assay_only"
    taxon_ids = {row.get("triarchy_taxon_id", "") for row in expected_calls}
    names = {row.get("canonical_name", "").casefold() for row in expected_calls}
    levels = {
        int(components[assay]["rank_level"])
        for assay in detected_expected
    }
    mapping = min(
        (number(row.get("mapping_confidence"), 1) for row in expected_calls),
        default=1,
    )
    if len(taxon_ids) == 1 or len(names) == 1:
        if len(levels) == 1:
            return clamp(mapping), "exact_or_canonical_agreement"
        return clamp(0.85 * mapping), "hierarchical_agreement"
    return clamp(0.25 * mapping), "ambiguous_compatibility"


def authenticity_interval(probability: float, evidence_units: float) -> tuple[float, float, float]:
    n_eff = max(2.0, 2.0 + evidence_units)
    stderr = math.sqrt(max(probability * (1 - probability), 1e-9) / n_eff)
    half_width = 1.96 * stderr
    low = clamp(probability - half_width)
    high = clamp(probability + half_width)
    return low, high, half_width


def parent_cross_assay_status(classification: str) -> str:
    if classification == "concordant_detection":
        return "concordant"
    if classification in {"single_assay_credible", "assay_limited_nondetection"}:
        return "assay_limited"
    if classification in {"rank_limited_agreement", "reference_limited"}:
        return "taxonomy_limited"
    if classification in {"probable_contaminant", "cross_contamination_candidate"}:
        return "contamination_associated"
    if classification == "genuine_discordance":
        return "genuine_discordance"
    return "indeterminate"


def authenticity_state(probability: float, low: float, high: float, policy: dict[str, Any]) -> str:
    supported = number(policy["authenticity_supported_threshold"])
    absent = number(policy["authenticity_absent_threshold"])
    overlap_low = number(policy["indeterminate_ci_overlap_low"])
    overlap_high = number(policy["indeterminate_ci_overlap_high"])
    if probability >= supported and low >= overlap_high:
        return "biologically_supported"
    if probability <= absent and high <= overlap_low:
        return "not_supported"
    return "indeterminate"


def library_depths(rows: list[dict[str, str]]) -> dict[str, float]:
    result: dict[str, float] = defaultdict(float)
    native_available: dict[str, bool] = defaultdict(bool)
    for row in rows:
        library = row["library_id"]
        if row.get("native_count", ""):
            native_available[library] = True
            result[library] += max(0, number(row["native_count"]))
    for row in rows:
        library = row["library_id"]
        if not native_available[library]:
            result[library] = -1
    return result


def classify_reason(
    calls: list[dict[str, str]],
    components: dict[str, dict[str, float | bool]],
    expected: tuple[str, ...],
    contamination: float,
    batch_score: float,
    policy: dict[str, Any],
) -> tuple[str, str]:
    detected = [assay for assay in expected if components[assay]["detected"]]
    conflict = any(
        row.get("cross_assay_compatibility", "").casefold() == "conflict"
        or row.get("mapping_method", "").casefold() in {"incompatible", "conflict"}
        for row in calls
    )
    if contamination >= number(policy["contaminant_threshold"]):
        return "probable_contaminant", "linked-control evidence exceeds the contaminant threshold"
    if (
        contamination >= number(policy["cross_contamination_threshold"])
        and batch_score >= number(policy["batch_localization_threshold"])
    ):
        return "cross_contamination_candidate", "control-associated signal is strongly localized to one batch"
    if conflict:
        return "genuine_discordance", "technically supported assignments are taxonomically incompatible"
    if len(detected) == len(expected):
        levels = {
            int(components[assay]["rank_level"])
            for assay in expected if components[assay]["detected"]
        }
        if len(levels) > 1:
            return "rank_limited_agreement", "assays agree only after roll-up to their shared defensible rank"
        return "concordant_detection", "all kingdom-appropriate assays support the canonical taxon"
    if len(detected) == 1:
        other = next(assay for assay in expected if assay not in detected)
        component = components[other]
        if not component["measured"]:
            return "single_assay_credible", f"{other} was not measured and supplies no contradictory evidence"
        if number(component["reference"]) < number(policy["minimum_reference_coverage"]):
            return "reference_limited", f"{other} reference coverage is inadequate for a defensible negative call"
        if number(component["power"]) < number(policy["low_detection_power"]):
            return "assay_limited_nondetection", f"{other} estimated detection power is insufficient"
        if number(component["power"]) >= number(policy["adequate_detection_power"]):
            return "genuine_discordance", f"{other} was technically adequate but did not support the detection"
        return "single_assay_credible", f"one expected assay supports the taxon; {other} is not decisively contradictory"
    if any(component["detected"] for component in components.values()):
        return "indeterminate", "support occurs outside the kingdom-appropriate assay pair"
    return "indeterminate", "no usable assay evidence is available"


def infer_evidence(
    harmonized_rows: list[dict[str, str]],
    sample_rows: list[dict[str, str]],
    performance: dict[tuple[str, str], dict[str, float]],
    policy: dict[str, Any],
    qc_rows: list[dict[str, str]] | None = None,
) -> list[dict[str, Any]]:
    samples = {row["library_id"]: row for row in sample_rows}
    qc_metrics = {
        row["library_id"]: row for row in (qc_rows or [])
    }
    biological_libraries: dict[str, dict[str, str]] = defaultdict(dict)
    for library_id, sample in samples.items():
        if sample.get("sample_role") == "biological":
            assay = sample.get("assay", library_id.rsplit("_", 1)[-1])
            biological_libraries[sample["biological_sample_id"]][assay] = library_id
    depths = library_depths(harmonized_rows)
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    catalog: dict[str, dict[str, str]] = {}
    for row in harmonized_rows:
        catalog.setdefault(row["triarchy_taxon_id"], row)
        sample = samples.get(row["library_id"])
        if sample and sample.get("sample_role") == "biological":
            grouped[(sample["biological_sample_id"], row["triarchy_taxon_id"])].append(row)

    evidence: list[dict[str, Any]] = []
    prior = clamp(number(policy["prior_presence"], .25))
    prior_log_odds = math.log(max(prior, 1e-9) / max(1 - prior, 1e-9))
    candidate_keys = [
        (sample_id, taxon_id)
        for sample_id in sorted(biological_libraries)
        for taxon_id in sorted(catalog)
    ]
    for sample_id, taxon_id in candidate_keys:
        calls = grouped.get((sample_id, taxon_id), [])
        representative = catalog[taxon_id]
        name = representative["canonical_name"]
        kingdom = infer_kingdom(representative)
        expected = expected_assays(kingdom)
        components: dict[str, dict[str, float | bool]] = {}
        max_contamination = max(
            (number(row.get("contamination_probability")) for row in calls),
            default=0,
        )
        max_batch = max(
            (number(row.get("batch_localization_score")) for row in calls),
            default=0,
        )
        max_batch_contamination = max(
            (number(row.get("batch_contamination_score")) for row in calls),
            default=0,
        )
        stage_scores = stage_contamination_scores(calls)
        source_batch_score = max_batch_contamination
        if not source_batch_score and (
            max_contamination >= number(policy["cross_contamination_threshold"])
            and max_batch >= number(policy["batch_localization_threshold"])
        ):
            source_batch_score = max_batch
        strongest_source = strongest_contamination_source(stage_scores, source_batch_score)
        log_odds = prior_log_odds
        evidence_units = 0.0
        for assay in ASSAYS:
            assay_calls = [row for row in calls if row["assay"] == assay]
            detected = bool(assay_calls)
            library = biological_libraries.get(sample_id, {}).get(assay, "")
            measured = bool(library)
            settings = performance.get((assay, kingdom), {
                "reliability": 0, "minimum_depth": 1, "reference_coverage": 0,
                "false_positive_rate": .05,
            })
            raw_depth = depths.get(library, 0)
            depth = number(settings["minimum_depth"]) if raw_depth < 0 and measured else max(0, raw_depth)
            depth_factor = 1 - math.exp(-depth / max(number(settings["minimum_depth"]), 1))
            quality = 1.0
            sample = samples.get(library, {})
            qc_row = qc_metrics.get(library, {})
            qc = qc_row.get(
                "qc_status", sample.get("qc_status", "pass")
            ).casefold()
            if qc in {"failed", "empty", "unusable"}:
                quality = 0
            elif qc in {"low_depth", "host_dominated", "review"}:
                quality = .4
            quality *= clamp(number(qc_row.get("quality_factor"), 1))
            call = max(
                assay_calls,
                key=lambda item: number(item.get("native_normalized_value")),
                default=None,
            )
            mapping = clamp(number(call.get("mapping_confidence"), 0) if call else 1)
            reference_value = (
                call.get("database_coverage") or call.get("reference_coverage")
                if call else ""
            )
            reference = clamp(number(
                reference_value if reference_value != "" else settings["reference_coverage"]
            ))
            reliability = clamp(number(settings["reliability"]))
            power = clamp(reliability * depth_factor * reference * quality) if measured else 0
            relative = clamp(number(call.get("native_normalized_value")) if call else 0)
            rank = (call.get("canonical_rank", "unclassified") if call else "unclassified").casefold()
            rank_level = RANK_LEVEL.get(rank, 0)
            abundance_rank = max(1, int(number(call.get("abundance_rank"), 1))) if call else 0
            detection_support = clamp(1 - math.exp(-5 * relative)) if call else 0
            abundance_support = 1 / abundance_rank if abundance_rank else 0
            resolution = rank_level / 8
            contamination = clamp(number(call.get("contamination_probability")) if call else 0)
            primer = clamp(number(
                call.get("primer_compatibility_score", call.get("primer_compatibility", 1))
            ) if call else 1)
            assay_evidence = (
                detection_support * math.sqrt(abundance_support)
                * resolution * mapping * reliability * primer * (1 - contamination)
            )
            false_positive = clamp(number(settings["false_positive_rate"], .01))
            false_positive = clamp(
                false_positive
                + 0.25 * contamination
                + 0.20 * max(stage_scores.values(), default=0.0)
            )
            if assay in expected:
                evidence_units += int(measured) + int(detected) + power + mapping + reference
                if detected:
                    log_odds += 4 * assay_evidence + math.log((power + .05) / (false_positive + .05))
                elif measured:
                    log_odds += math.log(max(1 - power, .01) / max(1 - false_positive, .01))
            components[assay] = {
                "measured": measured, "detected": detected, "depth": depth,
                "power": power, "reference": reference, "reliability": reliability,
                "rank_level": rank_level, "mapping": mapping,
                "abundance": abundance_support, "resolution": resolution,
                "detection": detection_support, "evidence": assay_evidence,
                "primer": primer,
            }
        log_odds -= 3 * max_contamination
        log_odds -= 2 * max(stage_scores.values(), default=0.0)
        biological_authenticity = sigmoid(log_odds)
        ci_low, ci_high, ci_half_width = authenticity_interval(
            biological_authenticity, evidence_units
        )
        auth_state = authenticity_state(
            biological_authenticity, ci_low, ci_high, policy
        )
        probability, calibration_status = calibrated_probability(
            biological_authenticity, policy.get("calibration")
        )
        classification, reason = classify_reason(
            calls, components, expected, max_contamination, max_batch, policy
        )
        compatibility_score, compatibility_status = taxonomic_compatibility(
            calls, components, expected
        )
        cross_assay_status = parent_cross_assay_status(classification)
        detected_levels = [
            int(components[assay]["rank_level"])
            for assay in expected if components[assay]["detected"]
        ]
        comparison_level = min(detected_levels, default=0)
        if auth_state == "biologically_supported":
            state = "present"
        elif auth_state == "not_supported":
            state = "absent"
        else:
            state = "indeterminate"
        evidence.append({
            "biological_sample_id": sample_id,
            "triarchy_taxon_id": taxon_id,
            "canonical_name": name,
            "kingdom": kingdom,
            "comparison_rank": LEVEL_RANK[comparison_level],
            "lowest_accurately_resolved_rank": LEVEL_RANK[comparison_level],
            "assay_support": ";".join(assay for assay in ASSAYS if components[assay]["detected"]),
            **{
                f"{assay.lower()}_{field}": (
                    str(value).lower() if isinstance(value, bool) else fmt(number(value))
                )
                for assay in ASSAYS
                for field, value in (
                    ("measured", components[assay]["measured"]),
                    ("detected", components[assay]["detected"]),
                    ("library_depth", components[assay]["depth"]),
                    ("detection_power", components[assay]["power"]),
                    ("reference_coverage", components[assay]["reference"]),
                    ("evidence_score", components[assay]["evidence"]),
                )
            },
            "abundance_support": fmt(max(number(item["abundance"]) for item in components.values())),
            "taxonomic_resolution_support": fmt(comparison_level / 8),
            "mapping_support": fmt(min(
                (number(components[assay]["mapping"]) for assay in ASSAYS if components[assay]["detected"]),
                default=0,
            )),
            "assay_reliability": ";".join(
                f"{assay}:{fmt(number(components[assay]['reliability']))}" for assay in expected
            ),
            "control_evidence": "linked_stage_specific",
            "field_contamination_score": fmt(stage_scores["field"]),
            "extraction_contamination_score": fmt(stage_scores["extraction"]),
            "pcr_contamination_score": fmt(stage_scores["pcr"]),
            "library_contamination_score": fmt(stage_scores["library"]),
            "run_contamination_score": fmt(stage_scores["run"]),
            "other_control_contamination_score": fmt(stage_scores["other"]),
            "probable_contamination_source": strongest_source,
            "contamination_probability": fmt(max_contamination),
            "batch_localization_score": fmt(max_batch),
            "taxonomic_compatibility_score": fmt(compatibility_score),
            "taxonomic_compatibility_status": compatibility_status,
            "biological_authenticity": fmt(biological_authenticity),
            "authenticity_ci_low": fmt(ci_low),
            "authenticity_ci_high": fmt(ci_high),
            "authenticity_uncertainty": fmt(ci_half_width),
            "authenticity_state": auth_state,
            "cross_assay_status": cross_assay_status,
            "presence_score": fmt(biological_authenticity),
            "probability_authentic": fmt(probability),
            "calibration_status": calibration_status,
            "biological_state": state,
            "classification": classification,
            "discordance_reason": "" if classification == "concordant_detection" else reason,
            "decision_reason": reason,
            "model_version": policy["model_version"],
        })
    evidence.sort(key=lambda row: (
        row["biological_sample_id"], -number(row["probability_authentic"]), row["canonical_name"]
    ))
    sample_rank: dict[str, int] = defaultdict(int)
    for global_rank, row in enumerate(
        sorted(evidence, key=lambda item: -number(item["probability_authentic"])), start=1
    ):
        row["global_evidence_rank"] = global_rank
    for row in evidence:
        sample_rank[row["biological_sample_id"]] += 1
        row["sample_evidence_rank"] = sample_rank[row["biological_sample_id"]]
    return evidence


def run_inference(
    harmonized: Path,
    samples: Path,
    output: Path,
    model_output: Path | None = None,
    performance_path: Path | None = None,
    policy_path: Path | None = None,
    qc_metrics_path: Path | None = None,
) -> None:
    performance = load_performance(performance_path)
    policy = load_policy(policy_path)
    rows = infer_evidence(
        read_tsv(harmonized),
        read_tsv(samples),
        performance,
        policy,
        read_tsv(qc_metrics_path) if qc_metrics_path else None,
    )
    write_tsv(output, rows, EVIDENCE_FIELDS)
    if model_output:
        payload = {
            "model_version": policy["model_version"],
            "policy": policy,
            "assay_performance": {
                f"{assay}:{kingdom}": values
                for (assay, kingdom), values in sorted(performance.items())
            },
            "output_rows": len(rows),
            "interpretation": (
                "biological_authenticity is the stage-I posterior-style authenticity "
                "estimate. cross_assay_status is the stage-II explanation. "
                "probability_authentic is retained for compatibility and is calibrated "
                "only to the status named in calibration.status; synthetic calibration "
                "is not empirical validation."
            ),
        }
        model_output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def fit_isotonic(scores: list[float], labels: list[int]) -> dict[str, Any]:
    blocks = [
        {"low": score, "high": score, "sum": float(label), "count": 1}
        for score, label in sorted(zip(scores, labels))
    ]
    index = 0
    while index < len(blocks) - 1:
        left = blocks[index]["sum"] / blocks[index]["count"]
        right = blocks[index + 1]["sum"] / blocks[index + 1]["count"]
        if left <= right:
            index += 1
            continue
        blocks[index:index + 2] = [{
            "low": blocks[index]["low"],
            "high": blocks[index + 1]["high"],
            "sum": blocks[index]["sum"] + blocks[index + 1]["sum"],
            "count": blocks[index]["count"] + blocks[index + 1]["count"],
        }]
        index = max(0, index - 1)
    return {
        "method": "isotonic",
        "status": "fitted_requires_external_validation",
        "breakpoints": [block["high"] for block in blocks],
        "values": [block["sum"] / block["count"] for block in blocks],
        "training_observations": len(scores),
    }


def fit_calibration(evidence_path: Path, truth_path: Path, output: Path) -> None:
    truth: dict[tuple[str, str], int] = {}
    for row in read_tsv(truth_path):
        sample = row.get("biological_sample_id", row.get("sample_id", ""))
        taxon = row.get("canonical_name", row.get("taxon", ""))
        label = row.get("truth_status", row.get("label", "")).casefold()
        truth[(sample, taxon)] = int(label in {"authentic", "present", "true", "1"})
    scores, labels = [], []
    for row in read_tsv(evidence_path):
        key = (row["biological_sample_id"], row["canonical_name"])
        if key in truth:
            scores.append(number(row.get("biological_authenticity") or row.get("presence_score")))
            labels.append(truth[key])
    if not scores or len(set(labels)) < 2:
        raise ValueError("calibration requires matched positive and negative truth observations")
    model = fit_isotonic(scores, labels)
    output.write_text(
        json.dumps({"calibration": model}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def evaluate(evidence_path: Path, truth_path: Path, output: Path) -> None:
    truth: dict[tuple[str, str], str] = {}
    for row in read_tsv(truth_path):
        key = (
            row.get("biological_sample_id", row.get("sample_id", "")),
            row.get("canonical_name", row.get("taxon", "")),
        )
        truth[key] = row.get("truth_status", row.get("label", "")).casefold()
    tp = fp = fn = tn = 0
    brier_values = []
    for row in read_tsv(evidence_path):
        key = (row["biological_sample_id"], row["canonical_name"])
        if key not in truth:
            continue
        actual = truth[key] in {"authentic", "present", "true", "1"}
        predicted = row["biological_state"] == "present"
        tp += int(actual and predicted)
        fp += int(not actual and predicted)
        fn += int(actual and not predicted)
        tn += int(not actual and not predicted)
        calibration_score = row.get("biological_authenticity", row.get("probability_authentic", ""))
        brier_values.append((number(calibration_score) - int(actual)) ** 2)
    precision = tp / (tp + fp) if tp + fp else 0
    recall = tp / (tp + fn) if tp + fn else 0
    result = {
        "matched_observations": tp + fp + fn + tn,
        "true_positive": tp, "false_positive": fp, "false_negative": fn, "true_negative": tn,
        "precision": precision, "recall": recall,
        "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0,
        "brier_score": sum(brier_values) / len(brier_values) if brier_values else None,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

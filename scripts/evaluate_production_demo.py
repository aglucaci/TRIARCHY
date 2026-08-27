#!/usr/bin/env python3
"""Evaluate a completed production demo against its designed outcomes."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalized(value: str) -> str:
    return " ".join(value.replace("_", " ").casefold().split())


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    expectations = read_tsv(args.expectations)
    evidence_path = args.results / "integrated" / "taxon_evidence.tsv"
    qc_path = args.results / "qc" / "library_qc.tsv"
    profiles_path = args.results / "profiles" / "control_annotated" / "taxa.tsv"
    report_path = args.results / "triarchy_report.html"
    required_paths = (evidence_path, qc_path, profiles_path, report_path)
    missing = [str(path) for path in required_paths if not path.is_file()]
    if missing:
        raise RuntimeError("missing production outputs: " + ", ".join(missing))

    evidence = read_tsv(evidence_path)
    qc = read_tsv(qc_path)
    profiles = read_tsv(profiles_path)
    checks = []
    for expectation in expectations:
        check_type = expectation["check_type"]
        target = expectation["target_id"]
        name = normalized(expectation["canonical_name"])
        observed = ""
        passed = False
        detail = ""
        if check_type == "classification":
            candidates = [
                row for row in evidence
                if row.get("biological_sample_id") == target
                and normalized(row.get("canonical_name", "")) == name
            ]
            observed = candidates[0].get("classification", "") if candidates else "missing"
            accepted = set(expectation["expected_values"].split("|"))
            passed = observed in accepted
            detail = candidates[0].get("decision_reason", "") if candidates else "taxon evidence row not found"
        elif check_type == "host_fraction":
            candidates = [row for row in qc if row.get("library_id") == target]
            observed = candidates[0].get("host_fraction", "") if candidates else "missing"
            try:
                value = float(observed)
                passed = float(expectation["minimum"]) <= value <= float(expectation["maximum"])
            except (TypeError, ValueError):
                passed = False
            detail = f"expected {expectation['minimum']}–{expectation['maximum']}"
        elif check_type == "profile_detection":
            candidates = [
                row for row in profiles
                if row.get("library_id") == target
                and name in normalized(row.get("original_name", ""))
                and float(row.get("abundance") or 0) > 0
            ]
            observed = "detected" if candidates else "missing"
            passed = bool(candidates)
            detail = candidates[0].get("original_name", "") if candidates else "profile row not found"
        else:
            observed = "unsupported_check_type"
            detail = check_type

        required = expectation.get("required", "true").casefold() == "true"
        checks.append({
            "check_id": expectation["check_id"],
            "check_type": check_type,
            "target_id": target,
            "canonical_name": expectation["canonical_name"],
            "expected": expectation["expected_values"] or (
                f"{expectation['minimum']}–{expectation['maximum']}"
            ),
            "observed": observed,
            "passed": passed,
            "required": required,
            "detail": detail,
        })

    required_failures = [row for row in checks if row["required"] and not row["passed"]]
    return {
        "evaluation": "TRIARCHY Production Demo v1.0",
        "evaluated_at": datetime.now(timezone.utc).isoformat(),
        "results": str(args.results.resolve()),
        "ok": not required_failures,
        "passed": sum(bool(row["passed"]) for row in checks),
        "total": len(checks),
        "required_failures": len(required_failures),
        "checks": checks,
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--results", type=Path, required=True)
    result.add_argument("--expectations", type=Path, required=True)
    result.add_argument("--output", type=Path)
    return result


def main() -> int:
    args = parser().parse_args()
    try:
        report = evaluate(args)
    except (OSError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    output = args.output or args.results / "production_demo_evaluation.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for check in report["checks"]:
        status = "PASS" if check["passed"] else "FAIL"
        print(
            f"{status} {check['check_id']}: expected {check['expected']}; "
            f"observed {check['observed']}"
        )
    print(f"Production demo evaluation: {report['passed']}/{report['total']} checks passed")
    print(f"Evaluation report: {output}")
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

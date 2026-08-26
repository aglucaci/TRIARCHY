#!/usr/bin/env python3
"""Validate TRIARCHY sample and sample-control metadata contracts."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


SAMPLE_FIELDS = {
    "library_id",
    "biological_sample_id",
    "specimen_id",
    "subject_id",
    "assay",
    "sample_role",
    "sample_material",
    "environment_category",
    "library_layout",
    "fastq_r1",
    "fastq_r2",
    "extraction_batch",
    "pcr_plate",
    "library_batch",
    "sequencing_run",
    "technical_replicate",
}
CONTROL_FIELDS = {
    "library_id",
    "control_library_id",
    "control_stage",
    "relationship",
    "priority",
    "reason",
}
ASSAYS = {"WMS", "16S", "ITS"}
ROLES = {"biological", "negative_control", "positive_control", "mock"}
CONTROL_STAGES = {"field", "extraction", "pcr", "library", "sequencing"}
RELATIONSHIPS = {"paired", "same_batch", "run_level"}


class MetadataError(ValueError):
    """A deterministic metadata contract violation."""

    def __init__(self, code: str, message: str):
        super().__init__(f"ERROR {code}: {message}")
        self.code = code


def read_tsv(path: Path, required_fields: set[str], kind: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(required_fields - fields)
        if missing:
            raise MetadataError(
                "TRI-META-000",
                f"{kind} is missing required columns: {', '.join(missing)}.",
            )
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def validate_samples(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    libraries: dict[str, dict[str, str]] = {}
    fastqs: dict[str, str] = {}
    identities: dict[str, tuple[str, str]] = {}
    roles: dict[str, str] = {}

    for line, row in enumerate(rows, start=2):
        library_id = row["library_id"]
        if not library_id:
            raise MetadataError("TRI-META-002", f"line {line} has an empty library_id.")
        if library_id in libraries:
            raise MetadataError("TRI-META-001", f"duplicate library_id {library_id}.")
        if not row["biological_sample_id"]:
            raise MetadataError(
                "TRI-META-002",
                f"library {library_id} has an empty biological_sample_id.",
            )
        if row["assay"] not in ASSAYS:
            raise MetadataError(
                "TRI-META-003",
                f"library {library_id} has invalid assay {row['assay']!r}.",
            )
        if row["sample_role"] not in ROLES:
            raise MetadataError(
                "TRI-META-006",
                f"library {library_id} has invalid sample_role {row['sample_role']!r}.",
            )
        if row["library_layout"] not in {"paired", "single"}:
            raise MetadataError(
                "TRI-META-007",
                f"library {library_id} has invalid library_layout {row['library_layout']!r}.",
            )
        if row["library_layout"] == "paired" and not row["fastq_r2"]:
            raise MetadataError(
                "TRI-META-004",
                f"paired library {library_id} is missing fastq_r2.",
            )
        for field in ("fastq_r1", "fastq_r2"):
            path = row[field]
            if not path:
                continue
            if path in fastqs:
                raise MetadataError(
                    "TRI-META-005",
                    f"FASTQ path {path} is used by both {fastqs[path]} and {library_id}.",
                )
            fastqs[path] = library_id

        required_values = (
            "specimen_id",
            "subject_id",
            "sample_material",
            "environment_category",
            "fastq_r1",
            "extraction_batch",
            "library_batch",
            "sequencing_run",
            "technical_replicate",
        )
        empty = [field for field in required_values if not row[field]]
        if empty:
            raise MetadataError(
                "TRI-META-008",
                f"library {library_id} has empty required fields: {', '.join(empty)}.",
            )
        if row["assay"] in {"16S", "ITS"} and not row["pcr_plate"]:
            raise MetadataError(
                "TRI-META-008",
                f"amplicon library {library_id} has an empty pcr_plate.",
            )

        sample_id = row["biological_sample_id"]
        identity = (row["specimen_id"], row["subject_id"])
        if sample_id in identities and identities[sample_id] != identity:
            raise MetadataError(
                "TRI-META-009",
                f"biological_sample_id {sample_id} maps to inconsistent specimen or subject identifiers.",
            )
        role = row["sample_role"]
        if sample_id in roles and roles[sample_id] != role:
            raise MetadataError(
                "TRI-META-016",
                f"biological_sample_id {sample_id} is assigned conflicting sample roles "
                f"{roles[sample_id]} and {role}.",
            )
        identities[sample_id] = identity
        roles[sample_id] = role
        libraries[library_id] = row

    return libraries


def _batch_field(stage: str) -> str | None:
    return {
        "extraction": "extraction_batch",
        "pcr": "pcr_plate",
        "library": "library_batch",
        "sequencing": "sequencing_run",
    }.get(stage)


def validate_controls(
    rows: list[dict[str, str]], libraries: dict[str, dict[str, str]]
) -> None:
    graph: dict[str, set[str]] = defaultdict(set)

    for line, edge in enumerate(rows, start=2):
        source_id = edge["library_id"]
        control_id = edge["control_library_id"]
        if source_id not in libraries:
            raise MetadataError(
                "TRI-META-010",
                f"line {line} references unknown library_id {source_id}.",
            )
        if control_id not in libraries:
            raise MetadataError(
                "TRI-META-010",
                f"line {line} references unknown control_id {control_id}.",
            )
        stage = edge["control_stage"]
        if stage not in CONTROL_STAGES:
            raise MetadataError(
                "TRI-META-012",
                f"control edge {source_id}->{control_id} has invalid stage {stage!r}.",
            )
        if edge["relationship"] not in RELATIONSHIPS:
            raise MetadataError(
                "TRI-META-012",
                f"control edge {source_id}->{control_id} has invalid relationship.",
            )

        source = libraries[source_id]
        control = libraries[control_id]
        if control["sample_role"] not in {"negative_control", "positive_control", "mock"}:
            raise MetadataError(
                "TRI-META-013",
                f"control_id {control_id} is not a control library.",
            )
        if source["assay"] != control["assay"]:
            raise MetadataError(
                "TRI-META-013",
                f"{source_id} ({source['assay']}) cannot link to {control_id} ({control['assay']}).",
            )
        if stage == "pcr" and source["assay"] == "WMS":
            raise MetadataError(
                "TRI-META-011",
                f"PCR control {control_id} cannot be attached to WMS library {source_id}.",
            )
        field = _batch_field(stage)
        if (
            edge["relationship"] in {"paired", "same_batch"}
            and field
            and source[field] != control[field]
        ):
            raise MetadataError(
                "TRI-META-014",
                f"control_id {control_id} belongs to {field} {control[field]}, "
                f"but sample {source_id} belongs to {field} {source[field]}.",
            )
        try:
            if int(edge["priority"]) < 1:
                raise ValueError
        except ValueError:
            raise MetadataError(
                "TRI-META-012",
                f"control edge {source_id}->{control_id} has invalid priority {edge['priority']!r}.",
            ) from None
        graph[source_id].add(control_id)

    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(node: str) -> None:
        if node in visiting:
            raise MetadataError("TRI-META-015", f"cyclic control relationship includes {node}.")
        if node in visited:
            return
        visiting.add(node)
        for target in graph[node]:
            visit(target)
        visiting.remove(node)
        visited.add(node)

    for node in list(graph):
        visit(node)


def validate(sample_path: Path, control_path: Path) -> tuple[int, int]:
    samples = read_tsv(sample_path, SAMPLE_FIELDS, "sample manifest")
    controls = read_tsv(control_path, CONTROL_FIELDS, "control relationships")
    libraries = validate_samples(samples)
    validate_controls(controls, libraries)
    return len(samples), len(controls)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samples", type=Path)
    parser.add_argument("controls", type=Path)
    args = parser.parse_args(argv)
    try:
        sample_count, edge_count = validate(args.samples, args.controls)
    except (OSError, MetadataError) as error:
        print(error, file=sys.stderr)
        return 2
    print(f"VALID: {sample_count} libraries; {edge_count} control relationships.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

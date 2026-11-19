"""Repository-level consistency checks: versions, schema, and CLI."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path


def version(root: Path) -> str:
    return (root / "VERSION").read_text(encoding="utf-8").strip()


def test_version_is_semver(root: Path) -> None:
    assert re.fullmatch(r"\d+\.\d+\.\d+", version(root))


def test_versions_agree(root: Path) -> None:
    expected = version(root)
    config = (root / "nextflow.config").read_text(encoding="utf-8")
    assert f"version = '{expected}'" in config
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    assert f'version: "{expected}"' in citation
    for name in ("docker.config", "apptainer.config", "singularity.config"):
        text = (root / "conf" / name).read_text(encoding="utf-8")
        assert f"ghcr.io/aglucaci/triarchy:{expected}" in text, name


def test_cli_reports_version(root: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(root / "triarchy"), "--version"],
        capture_output=True, text=True, check=True,
    )
    assert result.stdout.strip() == f"triarchy {version(root)}"


def config_params(root: Path) -> set[str]:
    names: set[str] = set()
    for path in (root / "nextflow.config", root / "conf" / "primers.config"):
        text = path.read_text(encoding="utf-8")
        block = re.search(r"^params\s*\{(.*?)^\}", text, re.S | re.M)
        assert block, path
        names.update(re.findall(r"^\s+([a-z0-9_]+)\s*=", block.group(1), re.M))
    return names


def schema_params(root: Path) -> dict[str, dict]:
    schema = json.loads((root / "nextflow_schema.json").read_text(encoding="utf-8"))
    params: dict[str, dict] = {}
    for group in schema["$defs"].values():
        params.update(group["properties"])
    return params


def test_schema_documents_every_parameter(root: Path) -> None:
    documented = set(schema_params(root))
    declared = config_params(root)
    assert declared - documented == set(), "parameters missing from nextflow_schema.json"
    assert documented - declared == set(), "schema documents undeclared parameters"


def test_schema_defaults_match_primers(root: Path) -> None:
    params = schema_params(root)
    text = (root / "conf" / "primers.config").read_text(encoding="utf-8")
    for name, value in re.findall(r"(\w+)\s*=\s*'([^']*)'", text):
        assert params[name]["default"] == value, name


def test_provenance_records_triarchy_version(root: Path, tmp_path: Path) -> None:
    samples = tmp_path / "samples.tsv"
    controls = tmp_path / "controls.tsv"
    samples.write_text("library_id\n", encoding="utf-8")
    controls.write_text("library_id\n", encoding="utf-8")
    output = tmp_path / "provenance.json"
    subprocess.run(
        [
            sys.executable, str(root / "workflow" / "scripts" / "triarchy_integrate.py"),
            "provenance", "--samples", str(samples), "--controls", str(controls),
            "--project-dir", str(root), "--output", str(output),
        ],
        check=True,
    )
    document = json.loads(output.read_text(encoding="utf-8"))
    assert document["triarchy_version"] == version(root)


def test_parameters_doc_is_current(root: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(root / "scripts" / "render_parameters_doc.py"), "--check"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

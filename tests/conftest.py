"""Shared fixtures for TRIARCHY unit tests."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


@pytest.fixture
def root() -> Path:
    return ROOT

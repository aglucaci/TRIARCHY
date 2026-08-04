#!/usr/bin/env python3
"""Run Kraken2's k2 CLI with compatibility for trailing-slash NCBI paths."""

from __future__ import annotations

import os
import re
import runpy
import sys
import urllib.parse
from pathlib import Path


def basename_without_trailing_separator(path: os.PathLike[str] | str) -> str:
    """Match File::Basename behavior for NCBI directory URLs ending in '/'."""
    value = os.fspath(path)
    if value not in {"", "/", "\\"}:
        value = value.rstrip("/\\")
    return ORIGINAL_BASENAME(value)


def urlsplit_without_duplicate_path_separators(
    url: str, *args: object, **kwargs: object
) -> urllib.parse.SplitResult:
    """Normalize doubled separators introduced by trailing NCBI URL slashes."""
    result = ORIGINAL_URLSPLIT(url, *args, **kwargs)
    return result._replace(path=re.sub(r"/{2,}", "/", result.path))


if len(sys.argv) < 2:
    raise SystemExit("usage: run_k2_compat.py /path/to/k2 [k2 arguments]")

K2_SCRIPT = Path(sys.argv[1]).resolve()
if not K2_SCRIPT.is_file():
    raise SystemExit(f"k2 script not found: {K2_SCRIPT}")

ORIGINAL_BASENAME = os.path.basename
ORIGINAL_URLSPLIT = urllib.parse.urlsplit
os.path.basename = basename_without_trailing_separator
urllib.parse.urlsplit = urlsplit_without_duplicate_path_separators
sys.argv = [str(K2_SCRIPT), *sys.argv[2:]]
runpy.run_path(str(K2_SCRIPT), run_name="__main__")

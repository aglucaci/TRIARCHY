"""Downloader fallback when ENA serves a directory instead of a processed FASTQ."""

from __future__ import annotations

import argparse
import functools
import gzip
import hashlib
import http.server
import threading
from pathlib import Path

import pytest

import download_real_world_datasets as downloader


@pytest.fixture
def server(tmp_path: Path):
    root = tmp_path / "server"
    root.mkdir()
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(root))
    handler.log_message = lambda *args, **kwargs: None  # type: ignore[attr-defined]
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield root, f"http://127.0.0.1:{httpd.server_address[1]}"
    httpd.shutdown()


def fastq(path: Path, mate: int) -> tuple[str, str]:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(f"@read1/{mate}\nACGT\n+\nIIII\n")
    return hashlib.md5(path.read_bytes()).hexdigest(), str(path.stat().st_size)


def item(local: str, url: str, md5: str, size: str, fb_url: str, fb_md5: str, fb_size: str) -> dict:
    return {
        "run_accession": "ERR1", "local_path": local, "source_url": url, "md5": md5, "bytes": size,
        "fallback_url": fb_url, "fallback_md5": fb_md5, "fallback_bytes": fb_size,
    }


def run(tmp_path: Path, manifest: list[dict]) -> tuple[list, dict]:
    args = argparse.Namespace(
        force=False, verify_md5=True, timeout=10, retries=3, progress_interval=60, quiet=True
    )
    summary = {"downloaded": 0, "skipped_existing": 0, "bytes": 0, "seconds": 0.0}
    config = downloader.DATASETS["TRI-REAL-05"]
    failures = downloader.download_manifest_files(args, config, tmp_path / "data", manifest, summary)
    return failures, summary


def test_directory_on_server_falls_back_to_submitted_pair(server, tmp_path: Path) -> None:
    root, base = server
    md5_1, size_1 = fastq(root / "fastq/ERR1_1.fastq.gz", 1)
    (root / "fastq/ERR1_2.fastq.gz").mkdir()  # ENA defect: the path is a directory
    sub_md5_1, sub_size_1 = fastq(root / "run/S.R1.fastq.gz", 1)
    sub_md5_2, sub_size_2 = fastq(root / "run/S.R2.fastq.gz", 2)
    manifest = [
        item("reads/ERR1/ERR1_1.fastq.gz", f"{base}/fastq/ERR1_1.fastq.gz", md5_1, size_1,
             f"{base}/run/S.R1.fastq.gz", sub_md5_1, sub_size_1),
        item("reads/ERR1/ERR1_2.fastq.gz", f"{base}/fastq/ERR1_2.fastq.gz", "0" * 32, "100",
             f"{base}/run/S.R2.fastq.gz", sub_md5_2, sub_size_2),
    ]
    failures, summary = run(tmp_path, manifest)
    assert failures == []
    assert summary["downloaded"] == 2
    reads = tmp_path / "data/reads/ERR1"
    # Both mates come from the submitted upload so the pair is consistent.
    assert hashlib.md5((reads / "ERR1_1.fastq.gz").read_bytes()).hexdigest() == sub_md5_1
    assert hashlib.md5((reads / "ERR1_2.fastq.gz").read_bytes()).hexdigest() == sub_md5_2

    failures, summary = run(tmp_path, manifest)
    assert failures == [] and summary["skipped_existing"] == 2


def test_failure_is_reported_and_other_runs_continue(server, tmp_path: Path) -> None:
    root, base = server
    (root / "fastq/BAD_1.fastq.gz").mkdir(parents=True)
    good_md5, good_size = fastq(root / "fastq/GOOD_1.fastq.gz", 1)
    manifest = [
        {**item("reads/BAD/BAD_1.fastq.gz", f"{base}/fastq/BAD_1.fastq.gz", "0" * 32, "10", "", "", ""),
         "run_accession": "BAD"},
        {**item("reads/GOOD/GOOD_1.fastq.gz", f"{base}/fastq/GOOD_1.fastq.gz", good_md5, good_size, "", "", ""),
         "run_accession": "GOOD"},
    ]
    failures, summary = run(tmp_path, manifest)
    assert [failure["run_accession"] for failure in failures] == ["BAD"]
    assert "not a downloadable file" in failures[0]["error"]
    assert summary["downloaded"] == 1

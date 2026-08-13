#!/usr/bin/env python3
"""Run TRIARCHY benchmark datasets and measure wall time, memory, and ETA."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import re
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUNDLE = ROOT / "benchmarks" / "bundles" / "TRIARCHY_benchmark_bundle_v0.1.json"
DEFAULT_REPORT_DIR = ROOT / "output" / "benchmarks"
DEFAULT_RESULTS_DIR = ROOT / "results" / "resource-benchmarks"


def dataset_directory(dataset: str) -> Path:
    """Resolve a benchmark ID or an explicit benchmark directory."""
    supplied = Path(dataset).expanduser()
    if supplied.is_absolute() or len(supplied.parts) > 1:
        return supplied.resolve()
    return ROOT / "benchmarks" / dataset


def validate_dataset(dataset: str) -> Path:
    directory = dataset_directory(dataset)
    required = (
        directory / "metadata" / "samples.tsv",
        directory / "metadata" / "control_relationships.tsv",
        directory / "reads",
    )
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        locations = "\n  - ".join(missing)
        raise ValueError(
            f"benchmark dataset {dataset!r} is not installed or is incomplete; "
            f"missing:\n  - {locations}\n"
            "Extract the runtime benchmark so it contains metadata/ and reads/, "
            "then pass its directory to --dataset (or install it under benchmarks/<ID>)."
        )
    return directory


def manifest_library_count(dataset: str) -> int:
    path = dataset_directory(dataset) / "metadata" / "samples.tsv"
    with path.open(encoding="utf-8", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


def input_fastq_bytes(dataset: str) -> int:
    reads = dataset_directory(dataset) / "reads"
    return sum(path.stat().st_size for path in reads.glob("*.fastq.gz"))


def parse_duration(value: str) -> float | None:
    if not value or value == "-":
        return None
    units = {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}
    matches = re.findall(r"([0-9.]+)\s*(ms|[smhd])", value)
    return sum(float(number) * units[unit] for number, unit in matches) if matches else None


def parse_bytes(value: str) -> int | None:
    if not value or value == "-":
        return None
    match = re.fullmatch(r"([0-9.]+)\s*([KMGT]?B)", value.strip(), re.IGNORECASE)
    if not match:
        return None
    scales = {"B": 1, "KB": 1024, "MB": 1024**2, "GB": 1024**3, "TB": 1024**4}
    return round(float(match.group(1)) * scales[match.group(2).upper()])


def trace_metrics(trace_path: Path) -> dict[str, object]:
    if not trace_path.is_file():
        return {"task_count": 0, "failed_task_count": None, "trace_peak_rss_bytes": None,
                "trace_cpu_seconds": None, "summed_task_realtime_seconds": None}
    task_count = failed = 0
    peak_values: list[int] = []
    task_realtime = 0.0
    cpu_seconds = 0.0
    cpu_observed = False
    with trace_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            task_count += 1
            failed += row.get("status") != "COMPLETED"
            peak = parse_bytes(row.get("peak_rss", ""))
            if peak is not None:
                peak_values.append(peak)
            realtime = parse_duration(row.get("realtime", ""))
            if realtime is not None:
                task_realtime += realtime
                cpu = row.get("%cpu", "").rstrip("%")
                if cpu and cpu != "-":
                    cpu_seconds += realtime * float(cpu) / 100
                    cpu_observed = True
    return {
        "task_count": task_count,
        "failed_task_count": failed,
        "trace_peak_rss_bytes": max(peak_values) if peak_values else None,
        "trace_cpu_seconds": round(cpu_seconds, 3) if cpu_observed else None,
        "summed_task_realtime_seconds": round(task_realtime, 3),
    }


def process_tree_rss(root_pid: int) -> int:
    """Return process-tree RSS in bytes using the POSIX ps command."""
    try:
        completed = subprocess.run(
            ["ps", "-axo", "pid=,ppid=,rss="], check=False,
            capture_output=True, text=True,
        )
    except (OSError, PermissionError):
        return 0
    children: dict[int, list[int]] = {}
    rss_kib: dict[int, int] = {}
    for line in completed.stdout.splitlines():
        try:
            pid, parent, rss = (int(value) for value in line.split())
        except ValueError:
            continue
        children.setdefault(parent, []).append(pid)
        rss_kib[pid] = rss
    descendants = [root_pid]
    for pid in descendants:
        descendants.extend(children.get(pid, ()))
    return sum(rss_kib.get(pid, 0) for pid in descendants) * 1024


def child_rusage_peak_bytes() -> int:
    """Return the best available child-process max RSS estimate."""
    value = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return int(value if platform.system() == "Darwin" else value * 1024)


def run_dataset(dataset_id: str, mode: str, results_root: Path, extra_args: list[str]) -> dict[str, object]:
    dataset_dir = dataset_directory(dataset_id)
    dataset_label = dataset_dir.name
    outdir = results_root / mode / dataset_label
    log_dir = results_root / mode / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    command = [str(ROOT / "run_triarchy.sh")]
    if mode == "stub":
        command.append("--stub")
    command += [
        "--samples", str(dataset_dir / "metadata" / "samples.tsv"),
        "--controls", str(dataset_dir / "metadata" / "control_relationships.tsv"),
        "--outdir", str(outdir),
        *extra_args,
    ]
    started = time.monotonic()
    peak_rss = 0
    memory_method = "process_tree_ps"
    with (log_dir / f"{dataset_label}.log").open("w", encoding="utf-8") as log:
        process = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        while process.poll() is None:
            peak_rss = max(peak_rss, process_tree_rss(process.pid))
            time.sleep(0.2)
        peak_rss = max(peak_rss, process_tree_rss(process.pid))
    wall = time.monotonic() - started
    if peak_rss == 0:
        peak_rss = child_rusage_peak_bytes()
        memory_method = "child_rusage_fallback"
    trace = trace_metrics(outdir / "pipeline_info" / "trace.tsv")
    trace_cpu_seconds = trace["trace_cpu_seconds"]
    return {
        "dataset_id": dataset_label,
        "dataset_dir": str(dataset_dir),
        "mode": mode,
        "exit_code": process.returncode,
        "library_count": manifest_library_count(dataset_id),
        "input_fastq_bytes": input_fastq_bytes(dataset_id),
        "wall_seconds": round(wall, 3),
        "wall_hours": round(wall / 3600, 6),
        "libraries_per_wall_hour": round(manifest_library_count(dataset_id) / wall * 3600, 3),
        "observed_process_tree_peak_rss_bytes": peak_rss,
        "observed_peak_rss_gib": round(peak_rss / 1024**3, 4),
        "peak_memory_gib_hours": round(peak_rss / 1024**3 * wall / 3600, 6),
        "memory_measurement_method": memory_method,
        **trace,
        "trace_cpu_hours": round(float(trace_cpu_seconds) / 3600, 6) if trace_cpu_seconds is not None else None,
        "results_dir": str(outdir.relative_to(ROOT)),
        "log_file": str((log_dir / f"{dataset_label}.log").relative_to(ROOT)),
    }


def write_reports(records: list[dict[str, object]], report_dir: Path, mode: str) -> tuple[Path, Path]:
    report_dir.mkdir(parents=True, exist_ok=True)
    csv_path = report_dir / f"resource-profile-{mode}.csv"
    json_path = report_dir / f"resource-profile-{mode}.json"
    fields = list(records[0]) if records else []
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": mode,
        "warning": (
            "Stub measurements characterize orchestration only and must not be used to size "
            "production Kraken2, HUMAnN, DADA2, or reference-database runs."
            if mode == "stub" else
            "Production measurements are specific to the supplied reads, references, tools, hardware, and executor."
        ),
        "dataset_count": len(records),
        "library_count": sum(int(row["library_count"]) for row in records),
        "total_wall_seconds": round(sum(float(row["wall_seconds"]) for row in records), 3),
        "total_wall_hours": round(sum(float(row["wall_hours"]) for row in records), 6),
        "peak_process_tree_rss_bytes": max((int(row["observed_process_tree_peak_rss_bytes"]) for row in records), default=0),
        "peak_process_tree_rss_gib": round(max((float(row["observed_peak_rss_gib"]) for row in records), default=0), 4),
        "total_peak_memory_gib_hours": round(sum(float(row["peak_memory_gib_hours"]) for row in records), 6),
        "runs": records,
    }
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return csv_path, json_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument(
        "--dataset", action="append",
        help="Benchmark ID or dataset directory; overrides bundle and may be repeated",
    )
    parser.add_argument("--mode", choices=("stub", "production"), default="stub")
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--nextflow-arg", action="append", default=[], help="Extra argument passed to run_triarchy.sh")
    args = parser.parse_args()
    bundle_path = args.bundle if args.bundle.is_absolute() else ROOT / args.bundle
    datasets = args.dataset or json.loads(bundle_path.read_text())["datasets"]
    results_dir = args.results_dir if args.results_dir.is_absolute() else ROOT / args.results_dir
    report_dir = args.report_dir if args.report_dir.is_absolute() else ROOT / args.report_dir
    try:
        for item in datasets:
            validate_dataset(item)
    except ValueError as error:
        parser.error(str(error))
    total_libraries = sum(manifest_library_count(item) for item in datasets)
    records: list[dict[str, object]] = []
    completed_libraries = 0
    started = time.monotonic()
    for index, dataset_id in enumerate(datasets, 1):
        elapsed_before = time.monotonic() - started
        eta_before = (
            elapsed_before / completed_libraries * (total_libraries - completed_libraries)
            if completed_libraries else None
        )
        print(f"[{index}/{len(datasets)}] running {dataset_id} ({args.mode})", flush=True)
        record = run_dataset(dataset_id, args.mode, results_dir, args.nextflow_arg)
        records.append(record)
        completed_libraries += int(record["library_count"])
        elapsed = time.monotonic() - started
        remaining = total_libraries - completed_libraries
        eta = elapsed / completed_libraries * remaining if completed_libraries else 0
        record["eta_before_run_seconds"] = round(eta_before, 3) if eta_before is not None else None
        record["eta_after_run_seconds"] = round(eta, 3)
        print(
            f"  wall={record['wall_seconds']:.3f}s "
            f"peak_rss={int(record['observed_process_tree_peak_rss_bytes']) / 1024**2:.1f} MiB "
            f"remaining_eta={eta:.1f}s exit={record['exit_code']}",
            flush=True,
        )
        if record["exit_code"] != 0:
            csv_path, json_path = write_reports(records, report_dir, args.mode)
            print(f"partial reports: {csv_path} {json_path}", file=sys.stderr)
            return int(record["exit_code"])
    csv_path, json_path = write_reports(records, report_dir, args.mode)
    print(f"reports: {csv_path} {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

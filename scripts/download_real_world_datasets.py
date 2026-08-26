#!/usr/bin/env python3
"""Plan and optionally download public real-world datasets for TRIARCHY.

The default mode fetches ENA run metadata and writes manifests only. Add
--download when you are ready to retrieve the FASTQ files listed in the
manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT / "data" / "real_world"
ENA_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
ASSAYS = ("WMS", "16S", "ITS", "AMPLICON")
COMPLETE_TRIARCHY_ASSAYS = frozenset(("WMS", "16S", "ITS"))
CHUNK_SIZE = 1024 * 1024
DEFAULT_PROGRESS_INTERVAL = 2.0

ENA_FIELDS = [
    "run_accession",
    "study_accession",
    "secondary_study_accession",
    "sample_accession",
    "secondary_sample_accession",
    "experiment_accession",
    "experiment_alias",
    "sample_alias",
    "experiment_title",
    "study_title",
    "sample_title",
    "scientific_name",
    "tax_id",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_name",
    "library_construction_protocol",
    "library_layout",
    "target_gene",
    "instrument_platform",
    "instrument_model",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "submitted_ftp",
    "sra_ftp",
]

SELECTED_FIELDS = [
    *ENA_FIELDS,
    "triarchy_dataset_id",
    "triarchy_sample_group",
    "triarchy_assay",
    "triarchy_sample_role",
]

DOWNLOAD_FIELDS = [
    "dataset_id",
    "run_accession",
    "sample_group",
    "assay",
    "sample_role",
    "file_index",
    "source_url",
    "local_path",
    "md5",
    "bytes",
]

SAMPLE_FIELDS = [
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
]

CONTROL_FIELDS = [
    "library_id",
    "control_library_id",
    "control_stage",
    "relationship",
    "priority",
    "reason",
]


@dataclass(frozen=True)
class DatasetConfig:
    dataset_id: str
    study_accession: str
    title: str
    purpose: str
    default_max_groups: int
    selection: str
    qiita_id: str = ""
    note: str = ""
    group_mode: str = "sample_alias"


DATASETS = {
    config.dataset_id: config
    for config in (
        DatasetConfig(
            dataset_id="TRI-REAL-01",
            study_accession="ERP124610",
            qiita_id="12201",
            title="Shaffer mixed bacterial/fungal extraction and mock cohort",
            purpose=(
                "Primary empirical TRIARCHY validation target with physical "
                "mixed-community material, dilution structure, blanks, and "
                "matched WMS/16S/ITS runs."
            ),
            default_max_groups=30,
            selection="mock_controls",
            group_mode="qiita_core",
        ),
        DatasetConfig(
            dataset_id="TRI-REAL-02",
            study_accession="PRJEB31530",
            title="Arabidopsis phyllosphere matched WMS/16S/ITS cohort",
            purpose=(
                "Secondary empirical environmental validation target for "
                "complete matched WMS/16S/ITS cases across microbial-load "
                "variation."
            ),
            default_max_groups=36,
            selection="complete_cases",
        ),
        DatasetConfig(
            dataset_id="TRI-REAL-03",
            study_accession="PRJEB42019",
            qiita_id="13114",
            title="EMP500 matched multi-environment subset",
            purpose=(
                "Optional breadth cohort: select complete matched cases "
                "across several environment categories."
            ),
            default_max_groups=24,
            selection="complete_cases",
        ),
        DatasetConfig(
            dataset_id="TRI-REAL-04",
            study_accession="ERP117287",
            qiita_id="11666",
            title="HCHS/SOL matched microbiome subset",
            purpose=(
                "Optional stress-test cohort for split development and "
                "validation checks."
            ),
            default_max_groups=50,
            selection="complete_cases",
            note="Human-associated dataset; use as a supplemental stress test, not environmental validation.",
        ),
    )
}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize_path(path: Path) -> str:
    return str(path).replace("\\", "/")


def human_bytes(value: int | float | None) -> str:
    if value is None:
        return "unknown"
    size = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024 or unit == "TiB":
            return f"{size:.1f} {unit}" if unit != "B" else f"{int(size)} B"
        size /= 1024
    return f"{size:.1f} TiB"


def format_duration(seconds: float | None) -> str:
    if seconds is None:
        return "unknown"
    seconds = max(0, int(round(seconds)))
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours:d}h {minutes:02d}m {secs:02d}s"
    if minutes:
        return f"{minutes:d}m {secs:02d}s"
    return f"{secs:d}s"


def parse_int(value: str) -> int | None:
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def split_ena_list(value: str) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(";") if item.strip()]


def text_blob(row: dict[str, str]) -> str:
    fields = (
        "target_gene",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_name",
        "library_construction_protocol",
        "experiment_alias",
        "experiment_title",
        "sample_title",
        "sample_alias",
        "study_title",
    )
    return " ".join(row.get(field, "") for field in fields).lower()


def infer_assay(row: dict[str, str]) -> str:
    text = text_blob(row)
    strategy = row.get("library_strategy", "").upper()
    if any(token in text for token in ("internal transcribed", " its ", "its1", "its2", "fungal")):
        return "ITS"
    if any(token in text for token in (" 16s", "16s ", "v3-v4", "v4", "rrna amplicon")):
        return "16S"
    if strategy in {"WGS", "AMPLICON"} and "16s" in text:
        return "16S"
    if strategy == "AMPLICON":
        return "AMPLICON"
    if strategy == "WGS" or any(
        token in text
        for token in (
            "shotgun",
            "metagenomic",
            "whole genome",
            "whole-genome",
            "wms",
            "wgs",
        )
    ):
        return "WMS"
    return "UNKNOWN"


def infer_sample_role(row: dict[str, str]) -> str:
    text = text_blob(row)
    if any(token in text for token in ("blank", "negative control", "ntc", "no template", "water control")):
        return "negative_control"
    if any(token in text for token in ("mock", "zymo", "positive control", "standard", "dilution")):
        return "positive_control"
    if "control" in text:
        return "control"
    return "biological"


def normalize_qiita_core(value: str) -> str:
    if not value:
        return value
    if ":" in value:
        value = value.rsplit(":", 1)[1]
    if value.startswith("12201."):
        value = value[len("12201.") :]
    parts = value.split(".")
    extraction_tokens = {
        "homebrew",
        "magmax",
        "microbiome",
        "nucleomag",
        "powersoil",
        "zymo",
        "zymoMA",
    }
    if len(parts) >= 3 and parts[-1] in {"H", "L"} and parts[-2] in extraction_tokens:
        parts = parts[:-2]
    return ".".join(parts)


def sample_group(row: dict[str, str], group_mode: str = "sample_alias") -> str:
    for field in ("sample_alias", "secondary_sample_accession", "sample_accession", "experiment_accession"):
        value = row.get(field, "").strip()
        if value:
            return normalize_qiita_core(value) if group_mode == "qiita_core" else value
    return row.get("run_accession", "").strip()


def annotate_row(config: DatasetConfig, row: dict[str, str]) -> dict[str, str]:
    annotated = {field: row.get(field, "") for field in ENA_FIELDS}
    annotated.update(
        triarchy_dataset_id=config.dataset_id,
        triarchy_sample_group=sample_group(row, config.group_mode),
        triarchy_assay=infer_assay(row),
        triarchy_sample_role=infer_sample_role(row),
    )
    return annotated


def ena_query_url(accession: str) -> str:
    query = f'study_accession="{accession}" OR secondary_study_accession="{accession}"'
    params = {
        "result": "read_run",
        "query": query,
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
        "limit": "0",
    }
    return f"{ENA_SEARCH_URL}?{urllib.parse.urlencode(params)}"


def fetch_text(url: str, timeout: int, retries: int) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "TRIARCHY-real-world-downloader/1.0"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            last_error = error
            if attempt < retries:
                time.sleep(min(2**attempt, 30))
    raise RuntimeError(f"failed to fetch {url}: {last_error}") from last_error


def fetch_ena_runs(config: DatasetConfig, timeout: int, retries: int) -> tuple[str, list[dict[str, str]]]:
    url = ena_query_url(config.study_accession)
    text = fetch_text(url, timeout=timeout, retries=retries)
    if not text.strip():
        return url, []
    lines = text.splitlines()
    reader = csv.DictReader(lines, delimiter="\t")
    rows = [{field: row.get(field, "") for field in ENA_FIELDS} for row in reader]
    return url, rows


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def group_rows(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    groups: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(row["triarchy_sample_group"], []).append(row)
    return groups


def group_assays(rows: list[dict[str, str]]) -> set[str]:
    return {row["triarchy_assay"] for row in rows if row["triarchy_assay"] in {*ASSAYS, "UNKNOWN"}}


def is_complete_triarchy_group(rows: list[dict[str, str]]) -> bool:
    return COMPLETE_TRIARCHY_ASSAYS <= group_assays(rows)


def is_wms_amplicon_group(rows: list[dict[str, str]]) -> bool:
    assays = group_assays(rows)
    return "WMS" in assays and "AMPLICON" in assays


def group_role(rows: list[dict[str, str]]) -> str:
    roles = {row["triarchy_sample_role"] for row in rows}
    if "negative_control" in roles:
        return "negative_control"
    if "positive_control" in roles:
        return "positive_control"
    if "control" in roles:
        return "control"
    return "biological"


def select_groups(
    rows: list[dict[str, str]],
    selection: str,
    max_groups: int,
    max_control_groups: int | None,
    include_all_runs: bool,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    if include_all_runs:
        selected = sorted(rows, key=lambda row: (row["triarchy_sample_group"], row["run_accession"]))
        return selected, {
            "selection": "all_runs",
            "selected_groups": len(group_rows(selected)),
            "selected_runs": len(selected),
        }

    groups = group_rows(rows)
    complete = [
        group
        for group, group_rows_ in groups.items()
        if is_complete_triarchy_group(group_rows_)
    ]
    wms_amplicon = [
        group
        for group, group_rows_ in groups.items()
        if group not in complete and is_wms_amplicon_group(group_rows_)
    ]
    controls = [
        group
        for group, group_rows_ in groups.items()
        if group_role(group_rows_) != "biological"
    ]
    if selection == "mock_controls":
        control_set = set(controls)
        matched = sorted(complete) + sorted(wms_amplicon)
        control_limit = min(max(max_groups // 3, 1), 10, max_groups) if max_control_groups is None else max_control_groups
        control_selection = sorted(control_set)[:control_limit]
        biological_selection = [
            group for group in matched
            if group not in control_set
        ][: max(max_groups - len(control_selection), 0)]
        ordered = control_selection + biological_selection
    else:
        ordered = sorted(complete) + sorted(wms_amplicon)
    wanted = ordered[:max_groups]
    wanted_set = set(wanted)
    selected = [
        row for row in sorted(rows, key=lambda item: (item["triarchy_sample_group"], item["run_accession"]))
        if row["triarchy_sample_group"] in wanted_set
    ]
    return selected, {
        "selection": selection,
        "available_groups": len(groups),
        "complete_groups": len(complete),
        "wms_amplicon_groups": len(wms_amplicon),
        "control_or_mock_groups": len(controls),
        "selected_control_or_mock_groups": len([group for group in wanted if group in set(controls)]),
        "max_groups": max_groups,
        "max_control_groups": max_control_groups,
        "selected_groups": len(wanted),
        "selected_runs": len(selected),
    }


def ena_file_url(value: str) -> str:
    if value.startswith("ftp://"):
        value = value[len("ftp://") :]
    if value.startswith("https://"):
        return value
    if value.startswith("http://"):
        return "https://" + value[len("http://") :]
    if value.startswith("ftp.sra.ebi.ac.uk/"):
        return "https://" + value
    return "https://" + value


def build_download_manifest(dataset_dir: Path, selected_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    manifest: list[dict[str, object]] = []
    for row in selected_rows:
        urls = split_ena_list(row.get("fastq_ftp", ""))
        md5s = split_ena_list(row.get("fastq_md5", ""))
        sizes = split_ena_list(row.get("fastq_bytes", ""))
        for index, url in enumerate(urls, 1):
            source_url = ena_file_url(url)
            filename = Path(urllib.parse.urlparse(source_url).path).name
            run = row["run_accession"] or f"run_{len(manifest) + 1}"
            local_path = dataset_dir / "reads" / run / filename
            manifest.append(
                {
                    "dataset_id": row["triarchy_dataset_id"],
                    "run_accession": run,
                    "sample_group": row["triarchy_sample_group"],
                    "assay": row["triarchy_assay"],
                    "sample_role": row["triarchy_sample_role"],
                    "file_index": index,
                    "source_url": source_url,
                    "local_path": normalize_path(local_path.relative_to(dataset_dir)),
                    "md5": md5s[index - 1] if index <= len(md5s) else "",
                    "bytes": sizes[index - 1] if index <= len(sizes) else "",
                }
            )
    return manifest


def safe_id(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    value = re.sub(r"_+", "_", value).strip("_.-")
    return value or "unknown"


def qiita_ptid(row: dict[str, str]) -> str:
    match = re.search(r"qiita_ptid_([0-9]+)", row.get("experiment_title", ""))
    return match.group(1) if match else "unknown"


def infer_extraction_batch(row: dict[str, str]) -> str:
    text = " ".join((row.get("sample_alias", ""), row.get("sample_title", "")))
    for token in (
        "homebrew",
        "magmax",
        "microbiome2",
        "microbiome20",
        "microbiome",
        "nucleomag",
        "powersoil",
        "qiagen",
        "zymoMA",
        "zymo",
    ):
        if re.search(rf"(^|[. _-]){re.escape(token)}([. _-]|$)", text, flags=re.IGNORECASE):
            return f"EXT_{safe_id(token).upper()}"
    return f"EXT_{safe_id(qiita_ptid(row)).upper()}"


def infer_material(row: dict[str, str]) -> str:
    text = text_blob(row)
    for token, label in (
        ("blank", "blank"),
        ("mock", "mixed_mock"),
        ("dust", "dust"),
        ("soil", "soil"),
        ("stool", "stool"),
        ("oral", "oral"),
        ("food", "food"),
        ("built", "built_environment"),
        ("skin", "skin"),
    ):
        if token in text:
            return label
    return safe_id(row.get("scientific_name", "environmental_sample")).lower()


def infer_environment(row: dict[str, str]) -> str:
    scientific = row.get("scientific_name", "").strip()
    return safe_id(scientific).lower() if scientific else "environmental_microbiome"


def assay_for_metadata(row: dict[str, str], amplicon_assay: str) -> str:
    assay = row["triarchy_assay"]
    return amplicon_assay if assay == "AMPLICON" else assay


def fastq_paths_for_sample(dataset_dir: Path, metadata_dir: Path, row: dict[str, str]) -> tuple[str, str, str]:
    paths = []
    for url in split_ena_list(row.get("fastq_ftp", "")):
        source_url = ena_file_url(url)
        filename = Path(urllib.parse.urlparse(source_url).path).name
        paths.append(dataset_dir / "reads" / row["run_accession"] / filename)
    layout = "paired" if len(paths) >= 2 else "single"
    r1 = normalize_path(Path(os.path.relpath(paths[0], metadata_dir))) if paths else ""
    r2 = normalize_path(Path(os.path.relpath(paths[1], metadata_dir))) if len(paths) >= 2 else ""
    return layout, r1, r2


def build_triarchy_metadata(
    dataset_dir: Path,
    selected_rows: list[dict[str, str]],
    amplicon_assay: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    metadata_dir = dataset_dir / "metadata"
    samples: list[dict[str, object]] = []
    replicate_counts: dict[tuple[str, str], int] = {}
    seen_libraries: set[str] = set()

    for row in selected_rows:
        assay = assay_for_metadata(row, amplicon_assay)
        if assay not in {"WMS", "16S", "ITS"}:
            continue
        group = safe_id(row["triarchy_sample_group"])
        replicate_key = (group, assay)
        replicate_counts[replicate_key] = replicate_counts.get(replicate_key, 0) + 1
        replicate = replicate_counts[replicate_key]
        base_library = safe_id(f"{group}_{assay}_{row['run_accession']}")
        library_id = base_library
        suffix = 2
        while library_id in seen_libraries:
            library_id = f"{base_library}_{suffix}"
            suffix += 1
        seen_libraries.add(library_id)

        layout, fastq_r1, fastq_r2 = fastq_paths_for_sample(dataset_dir, metadata_dir, row)
        role = row["triarchy_sample_role"]
        if role == "control":
            role = "negative_control"
        samples.append(
            {
                "library_id": library_id,
                "biological_sample_id": group,
                "specimen_id": group,
                "subject_id": "CONTROL" if role != "biological" else group,
                "assay": assay,
                "sample_role": role,
                "sample_material": infer_material(row),
                "environment_category": infer_environment(row),
                "library_layout": layout,
                "fastq_r1": fastq_r1,
                "fastq_r2": fastq_r2,
                "extraction_batch": infer_extraction_batch(row),
                "pcr_plate": "" if assay == "WMS" else f"PCR_{assay}_{safe_id(qiita_ptid(row))}",
                "library_batch": f"LIB_{safe_id(qiita_ptid(row))}",
                "sequencing_run": f"SEQ_{safe_id(qiita_ptid(row))}",
                "technical_replicate": replicate,
            }
        )

    controls: list[dict[str, object]] = []
    summary = {
        "samples": len(samples),
        "controls": len(controls),
        "amplicon_assay_mapping": amplicon_assay,
        "warning": (
            "Draft metadata maps unresolved ENA AMPLICON runs to the requested "
            "TRIARCHY assay and leaves control relationships empty. Curate "
            "samples.tsv and control_relationships.tsv before biological claims."
        ),
    }
    return samples, controls, summary


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def existing_file_ok(path: Path, expected_md5: str, expected_bytes: str, verify_md5: bool) -> bool:
    if not path.is_file():
        return False
    if expected_bytes:
        try:
            if path.stat().st_size != int(expected_bytes):
                return False
        except ValueError:
            return False
    return not (verify_md5 and expected_md5) or md5sum(path) == expected_md5


def progress_message(label: str, done: int, total: int | None, started: float) -> str:
    elapsed = max(time.monotonic() - started, 1e-9)
    speed = done / elapsed
    if total:
        percent = done / total * 100
        eta = (total - done) / speed if speed > 0 else None
        return (
            f"{label} {percent:6.2f}% "
            f"{human_bytes(done)}/{human_bytes(total)} "
            f"at {human_bytes(speed)}/s ETA {format_duration(eta)}"
        )
    return f"{label} {human_bytes(done)} at {human_bytes(speed)}/s elapsed {format_duration(elapsed)}"


def emit_progress(message: str, tty: bool, final: bool = False) -> None:
    if tty:
        ending = "\n" if final else ""
        print("\r" + message.ljust(120), end=ending, file=sys.stderr, flush=True)
    else:
        print(message, file=sys.stderr, flush=True)


def download_one(
    source_url: str,
    destination: Path,
    expected_md5: str,
    expected_bytes: str,
    force: bool,
    verify_md5: bool,
    timeout: int,
    retries: int,
    progress_label: str,
    progress_interval: float,
    quiet: bool,
) -> dict[str, object]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not force and existing_file_ok(destination, expected_md5, expected_bytes, verify_md5):
        return {"status": "skipped_existing", "bytes": destination.stat().st_size, "seconds": 0.0}
    part = destination.with_suffix(destination.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            expected_size = parse_int(expected_bytes)
            started = time.monotonic()
            last_progress = 0.0
            downloaded = 0
            tty = sys.stderr.isatty()
            if not quiet:
                size_label = human_bytes(expected_size) if expected_size else "unknown size"
                print(
                    f"{progress_label} downloading {destination.name} ({size_label})",
                    file=sys.stderr,
                    flush=True,
                )
            request = urllib.request.Request(source_url, headers={"User-Agent": "TRIARCHY-real-world-downloader/1.0"})
            with urllib.request.urlopen(request, timeout=timeout) as response, part.open("wb") as handle:
                while True:
                    block = response.read(CHUNK_SIZE)
                    if not block:
                        break
                    handle.write(block)
                    downloaded += len(block)
                    now = time.monotonic()
                    if not quiet and now - last_progress >= progress_interval:
                        emit_progress(
                            progress_message(progress_label, downloaded, expected_size, started),
                            tty=tty,
                        )
                        last_progress = now
            seconds = time.monotonic() - started
            if expected_size is not None and part.stat().st_size != expected_size:
                raise RuntimeError(
                    f"size mismatch for {source_url}: got {part.stat().st_size}, expected {expected_size}"
                )
            if verify_md5 and expected_md5:
                if not quiet:
                    emit_progress(
                        f"{progress_label} verifying md5 for {destination.name}",
                        tty=tty,
                        final=True,
                    )
                actual = md5sum(part)
                if actual != expected_md5:
                    raise RuntimeError(f"md5 mismatch for {source_url}: got {actual}, expected {expected_md5}")
            part.replace(destination)
            if not quiet:
                emit_progress(
                    progress_message(progress_label, destination.stat().st_size, expected_size, started),
                    tty=tty,
                    final=True,
                )
            return {"status": "downloaded", "bytes": destination.stat().st_size, "seconds": seconds}
        except (urllib.error.URLError, TimeoutError, OSError, RuntimeError, ValueError) as error:
            last_error = error
            if part.exists():
                part.unlink()
            if attempt < retries:
                if not quiet:
                    print(
                        f"{progress_label} attempt {attempt}/{retries} failed: {error}; retrying",
                        file=sys.stderr,
                        flush=True,
                    )
                time.sleep(min(2**attempt, 30))
    raise RuntimeError(f"failed to download {source_url}: {last_error}") from last_error


def summarize(rows: list[dict[str, str]]) -> dict[str, dict[str, int]]:
    summary: dict[str, dict[str, int]] = {
        "assay_counts": {assay: 0 for assay in (*ASSAYS, "UNKNOWN")},
        "role_counts": {},
    }
    for row in rows:
        assay = row["triarchy_assay"] if row["triarchy_assay"] in summary["assay_counts"] else "UNKNOWN"
        summary["assay_counts"][assay] += 1
        role = row["triarchy_sample_role"]
        summary["role_counts"][role] = summary["role_counts"].get(role, 0) + 1
    return summary


def display_data_dir(dataset_dir: Path) -> str:
    try:
        return normalize_path(dataset_dir.relative_to(ROOT))
    except ValueError:
        return str(dataset_dir)


def process_dataset(args: argparse.Namespace, config: DatasetConfig) -> dict[str, object]:
    max_groups = args.max_groups if args.max_groups is not None else config.default_max_groups
    dataset_dir = args.data_dir / config.dataset_id
    metadata_dir = dataset_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    query_url, ena_rows = fetch_ena_runs(config, timeout=args.timeout, retries=args.retries)
    annotated = [annotate_row(config, row) for row in ena_rows]
    selected, selection_summary = select_groups(
        annotated,
        selection=config.selection,
        max_groups=max_groups,
        max_control_groups=args.max_control_groups,
        include_all_runs=args.include_all_runs,
    )
    download_manifest = build_download_manifest(dataset_dir, selected)
    samples, controls, metadata_summary = build_triarchy_metadata(
        dataset_dir,
        selected,
        amplicon_assay=args.amplicon_assay,
    )

    write_tsv(metadata_dir / "ena_runs.tsv", ENA_FIELDS, ena_rows)
    write_tsv(metadata_dir / "selected_runs.tsv", SELECTED_FIELDS, selected)
    write_tsv(metadata_dir / "download_manifest.tsv", DOWNLOAD_FIELDS, download_manifest)
    write_tsv(metadata_dir / "samples.tsv", SAMPLE_FIELDS, samples)
    write_tsv(metadata_dir / "control_relationships.tsv", CONTROL_FIELDS, controls)

    manifest_bytes = sum(parse_int(str(item["bytes"])) or 0 for item in download_manifest)
    download_summary = {
        "enabled": args.download,
        "downloaded": 0,
        "skipped_existing": 0,
        "bytes": 0,
        "seconds": 0.0,
        "manifest_bytes": manifest_bytes,
    }
    if args.download:
        started_downloads = time.monotonic()
        if not args.quiet:
            print(
                f"{config.dataset_id}: downloading {len(download_manifest)} FASTQ files "
                f"({human_bytes(manifest_bytes)} expected) into {display_data_dir(dataset_dir)}/reads",
                file=sys.stderr,
                flush=True,
            )
        for index, item in enumerate(download_manifest, 1):
            local_path = dataset_dir / str(item["local_path"])
            progress_label = f"{config.dataset_id} [{index}/{len(download_manifest)}]"
            result = download_one(
                str(item["source_url"]),
                local_path,
                str(item["md5"]),
                str(item["bytes"]),
                force=args.force,
                verify_md5=args.verify_md5,
                timeout=args.timeout,
                retries=args.retries,
                progress_label=progress_label,
                progress_interval=args.progress_interval,
                quiet=args.quiet,
            )
            status = str(result["status"])
            if status in download_summary:
                download_summary[status] = int(download_summary[status]) + 1
            download_summary["bytes"] = int(download_summary["bytes"]) + int(result["bytes"])
            download_summary["seconds"] = float(download_summary["seconds"]) + float(result["seconds"])
            if not args.quiet:
                print(
                    f"{progress_label} {status} {item['local_path']} "
                    f"({human_bytes(int(result['bytes']))})",
                    file=sys.stderr,
                    flush=True,
                )
        elapsed = time.monotonic() - started_downloads
        if not args.quiet:
            speed = int(download_summary["bytes"]) / elapsed if elapsed > 0 else 0
            print(
                f"{config.dataset_id}: download pass complete in {format_duration(elapsed)}; "
                f"downloaded {download_summary['downloaded']}, "
                f"skipped {download_summary['skipped_existing']}, "
                f"transferred {human_bytes(int(download_summary['bytes']))} "
                f"at {human_bytes(speed)}/s",
                file=sys.stderr,
                flush=True,
            )

    payload = {
        "dataset": asdict(config),
        "generated_at": now_utc(),
        "script": "scripts/download_real_world_datasets.py",
        "ena_query_url": query_url,
        "data_dir": display_data_dir(dataset_dir),
        "fetched_runs": len(ena_rows),
        "selected_runs": len(selected),
        "manifest_files": len(download_manifest),
        "selection_summary": selection_summary,
        "selected_summary": summarize(selected),
        "draft_triarchy_metadata": metadata_summary,
        "download_summary": download_summary,
        "triarchy_metadata_note": (
            "This script creates ENA run and file manifests. Curate "
            "TRIARCHY samples.tsv and control_relationships.tsv from the "
            "study metadata before running biological validation."
        ),
    }
    (metadata_dir / "dataset.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def list_datasets() -> None:
    for config in DATASETS.values():
        qiita = f", Qiita {config.qiita_id}" if config.qiita_id else ""
        print(f"{config.dataset_id}: {config.study_accession}{qiita} - {config.title}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        choices=[*DATASETS.keys(), "all"],
        help="Dataset to fetch. May be repeated. Default: TRI-REAL-01.",
    )
    parser.add_argument("--list-datasets", action="store_true", help="Print configured datasets and exit")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR, help="Output root; default data/real_world")
    parser.add_argument("--max-groups", type=int, help="Override selected sample/control groups for each dataset")
    parser.add_argument(
        "--max-control-groups",
        type=int,
        help="For mock-control selection, cap control/mock groups before adding biological groups",
    )
    parser.add_argument(
        "--amplicon-assay",
        choices=("16S", "ITS"),
        default="16S",
        help="TRIARCHY assay label for unresolved ENA AMPLICON runs in draft samples.tsv",
    )
    parser.add_argument("--include-all-runs", action="store_true", help="Select every ENA run returned by the study query")
    parser.add_argument("--download", action="store_true", help="Download FASTQ files listed in the selected manifest")
    parser.add_argument("--verify-md5", action="store_true", help="Verify downloaded files against ENA MD5 values when present")
    parser.add_argument("--force", action="store_true", help="Replace existing downloaded files")
    parser.add_argument("--timeout", type=int, default=120, help="HTTP timeout in seconds")
    parser.add_argument("--retries", type=int, default=3, help="Network retry attempts")
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=DEFAULT_PROGRESS_INTERVAL,
        help="Seconds between download progress updates; default 2",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress download progress output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.list_datasets:
        list_datasets()
        return 0
    requested = args.dataset or ["TRI-REAL-01"]
    dataset_ids = list(DATASETS) if "all" in requested else requested
    args.data_dir = args.data_dir if args.data_dir.is_absolute() else ROOT / args.data_dir
    try:
        for dataset_id in dataset_ids:
            payload = process_dataset(args, DATASETS[dataset_id])
            print(
                f"{dataset_id}: fetched {payload['fetched_runs']} runs, "
                f"selected {payload['selected_runs']} runs, "
                f"manifested {payload['manifest_files']} FASTQ files",
                flush=True,
            )
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    if not args.download:
        dataset_flags = " ".join(f"--dataset {dataset_id}" for dataset_id in dataset_ids)
        print(
            "Metadata and download manifests are ready. To fetch FASTQs, rerun:\n"
            f"python scripts/download_real_world_datasets.py {dataset_flags} --download --verify-md5"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

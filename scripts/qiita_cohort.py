#!/usr/bin/env python3
"""Build curated TRIARCHY cohorts from Qiita sample information plus ENA runs.

The Shaffer et al. (2022) DNA-extraction comparison (Qiita 12201, ENA
ERP124610) records, for every extraction, the kit, extraction round,
extraction plate (high- or low-biomass) and well. That makes it possible to
link each extraction blank to exactly the samples extracted on the same plate,
and to identify ZymoBIOMICS mock-community dilution series with a published
composition. This module turns that metadata into TRIARCHY sample manifests,
control relationships and a mock-community truth table.
"""

from __future__ import annotations

import csv
import io
import re
import time
import urllib.error
import urllib.request
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

QIITA_SAMPLE_INFO_URL = "https://qiita.ucsd.edu/public_download/?data=sample_information&study_id={study_id}"

COMPLETE_ASSAYS = ("WMS", "16S", "ITS")

# Theoretical compositions printed in Table 1 of the Zymo Research instruction
# manuals for D6300 (Microbial Community Standard) and D6310 (Microbial
# Community Standard II, Log Distribution). Columns: genomic DNA %, 16S-only %,
# 16S & 18S %, genome copy %, cell number %. "NA" = not applicable (yeasts in
# the 16S-only column).
ZYMO_SOURCES = {
    "D6300": "https://files.zymoresearch.com/protocols/_d6300_zymobiomics_microbial_community_standard.pdf",
    "D6310": (
        "https://files.zymoresearch.com/protocols/"
        "_d6310_zymobiomics_microbial_community_standard_ii_(log_distribution).pdf"
    ),
}
ZYMO_COMPOSITION = {
    "D6300": [
        ("Pseudomonas aeruginosa", "Bacteria", "12", "4.2", "3.6", "6.1", "6.1"),
        ("Escherichia coli", "Bacteria", "12", "10.1", "8.9", "8.5", "8.5"),
        ("Salmonella enterica", "Bacteria", "12", "10.4", "9.1", "8.7", "8.8"),
        ("Lactobacillus fermentum", "Bacteria", "12", "18.4", "16.1", "21.6", "21.9"),
        ("Enterococcus faecalis", "Bacteria", "12", "9.9", "8.7", "14.6", "14.6"),
        ("Staphylococcus aureus", "Bacteria", "12", "15.5", "13.6", "15.2", "15.3"),
        ("Listeria monocytogenes", "Bacteria", "12", "14.1", "12.4", "13.9", "13.9"),
        ("Bacillus subtilis", "Bacteria", "12", "17.4", "15.3", "10.3", "10.3"),
        ("Saccharomyces cerevisiae", "Fungi", "2", "NA", "9.3", "0.57", "0.29"),
        ("Cryptococcus neoformans", "Fungi", "2", "NA", "3.3", "0.37", "0.18"),
    ],
    "D6310": [
        ("Listeria monocytogenes", "Bacteria", "89.1", "95.9", "91.9", "94.8", "94.9"),
        ("Pseudomonas aeruginosa", "Bacteria", "8.9", "2.8", "2.7", "4.2", "4.2"),
        ("Bacillus subtilis", "Bacteria", "0.89", "1.2", "1.1", "0.7", "0.7"),
        ("Saccharomyces cerevisiae", "Fungi", "0.89", "NA", "4.1", "0.23", "0.12"),
        ("Escherichia coli", "Bacteria", "0.089", "0.069", "0.066", "0.058", "0.058"),
        ("Salmonella enterica", "Bacteria", "0.089", "0.07", "0.067", "0.059", "0.059"),
        ("Lactobacillus fermentum", "Bacteria", "0.0089", "0.012", "0.012", "0.015", "0.015"),
        ("Enterococcus faecalis", "Bacteria", "0.00089", "0.00067", "0.00064", "0.001", "0.001"),
        ("Cryptococcus neoformans", "Fungi", "0.00089", "NA", "0.0014", "0.00015", "0.00007"),
        ("Staphylococcus aureus", "Bacteria", "0.000089", "0.0001", "0.0001", "0.0001", "0.0001"),
    ],
}
ZYMO_PRODUCTS = {
    "ZymoBIOMICS Microbial Community Standard I": "D6300",
    "ZymoBIOMICS Microbial Community Standard II": "D6310",
}

# Genus names each reference database may use for the mock organisms.
GENUS_ALIASES = {
    "Pseudomonas": ("Pseudomonas",),
    "Escherichia": ("Escherichia", "Escherichia-Shigella", "Shigella"),
    "Salmonella": ("Salmonella",),
    "Lactobacillus": ("Lactobacillus", "Limosilactobacillus"),
    "Enterococcus": ("Enterococcus",),
    "Staphylococcus": ("Staphylococcus",),
    "Listeria": ("Listeria",),
    "Bacillus": ("Bacillus",),
    "Saccharomyces": ("Saccharomyces",),
    "Cryptococcus": ("Cryptococcus", "Filobasidiella"),
}

TRUTH_FIELDS = [
    "biological_sample_id",
    "standard",
    "dilution",
    "technical_replicate",
    "species",
    "genus",
    "genus_aliases",
    "kingdom",
    "expected_assays",
    "genomic_dna_pct",
    "rrna_16s_pct",
    "rrna_16s_18s_pct",
    "genome_copy_pct",
    "cell_number_pct",
    "source",
]

CONTROL_SAMPLE_TYPES = {"control blank": "negative_control"}
POSITIVE_SAMPLE_TYPE = "control positive"


class CohortError(RuntimeError):
    """The requested cohort cannot be built from the available metadata."""


def fetch_qiita_sample_info(
    study_id: str, cache_path: Path, timeout: int, retries: int, refresh: bool = False
) -> list[dict[str, str]]:
    """Download (or reuse) the public Qiita sample-information template."""
    if refresh or not cache_path.is_file():
        url = QIITA_SAMPLE_INFO_URL.format(study_id=study_id)
        last_error: Exception | None = None
        payload = b""
        for attempt in range(1, retries + 1):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "TRIARCHY-real-world-downloader/1.0"})
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    payload = response.read()
                break
            except (urllib.error.URLError, TimeoutError, OSError) as error:
                last_error = error
                if attempt < retries:
                    time.sleep(min(2**attempt, 30))
        else:
            raise CohortError(f"failed to fetch Qiita sample information: {last_error}") from last_error
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            names = [name for name in archive.namelist() if name.endswith(".txt")]
            if not names:
                raise CohortError("Qiita sample-information archive contains no template")
            text = archive.read(sorted(names)[-1]).decode("utf-8", errors="replace")
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(text, encoding="utf-8")
    with cache_path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_id(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    value = re.sub(r"_+", "_", value).strip("_.-")
    return value or "unknown"


def plate_key(sample: dict[str, str]) -> str:
    return safe_id(
        f"R{sample['round']}_{sample['extraction_kit_round']}_{sample['biomass_plate']}"
    ).replace(".", "_")


def prep_template(run: dict[str, str]) -> str:
    match = re.search(r"qiita_ptid_([0-9]+)", run.get("experiment_title", ""))
    return match.group(1) if match else "unknown"


def technical_replicate(sample: dict[str, str]) -> str:
    replicate = sample.get("sample_technical_replicate", "")
    if replicate.isdigit():
        return replicate
    match = re.search(r"\.([0-9]+)\.[^.]+\.[HL]$", sample["sample_name"])
    if match and not re.search(r"\.d[0-9]+\.[^.]+\.[HL]$", sample["sample_name"]):
        return match.group(1)
    return "1"


def mock_dilution(sample_name: str) -> str:
    match = re.search(r"\.d([0-9]+)\.", sample_name)
    return f"d{match.group(1)}" if match else ""


def sample_role(sample: dict[str, str]) -> str:
    if sample.get("sample_type") in CONTROL_SAMPLE_TYPES:
        return CONTROL_SAMPLE_TYPES[sample["sample_type"]]
    return "biological"


def zymo_standard(sample: dict[str, str]) -> str:
    return ZYMO_PRODUCTS.get(sample.get("sample_type_2", ""), "")


def available_plates(
    qiita_rows: list[dict[str, str]], runs_by_sample: dict[str, dict[str, dict[str, str]]]
) -> Counter[tuple[str, str, str]]:
    plates: Counter[tuple[str, str, str]] = Counter()
    for sample in qiita_rows:
        if set(COMPLETE_ASSAYS) <= set(runs_by_sample.get(sample["sample_name"], {})):
            plates[(sample["round"], sample["extraction_kit_round"], sample["biomass_plate"])] += 1
    return plates


def select_cohort(
    qiita_rows: list[dict[str, str]],
    ena_rows: list[dict[str, str]],
    kits: list[str],
    plates: list[str],
    max_biological_per_plate: int | None = None,
) -> tuple[list[dict[str, str]], dict[str, dict[str, str]], dict[str, object]]:
    """Select complete WMS/16S/ITS samples on the requested extraction plates.

    `kits` are Qiita `extraction_kit_round` values (for example "PowerSoil r1").
    Every blank and every ZymoBIOMICS standard on a selected plate is kept;
    other biological samples are optionally capped per plate, stratified by
    sample type so every material is represented.
    """
    runs_by_sample: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for run in ena_rows:
        assay = run["triarchy_assay"]
        if assay not in COMPLETE_ASSAYS:
            continue
        existing = runs_by_sample[run["sample_title"]].get(assay)
        if existing is None or run["run_accession"] < existing["run_accession"]:
            runs_by_sample[run["sample_title"]][assay] = run

    known = {sample["extraction_kit_round"] for sample in qiita_rows}
    unknown = [kit for kit in kits if kit not in known]
    if unknown:
        raise CohortError(
            f"unknown extraction kit(s) {', '.join(unknown)}; choose from {', '.join(sorted(known))}"
        )

    by_plate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for sample in qiita_rows:
        if sample["extraction_kit_round"] not in kits or sample["biomass_plate"] not in plates:
            continue
        if not set(COMPLETE_ASSAYS) <= set(runs_by_sample.get(sample["sample_name"], {})):
            continue
        by_plate[plate_key(sample)].append(sample)
    if not by_plate:
        raise CohortError("no complete WMS/16S/ITS samples on the requested extraction plates")

    selected: list[dict[str, str]] = []
    plate_summary: dict[str, dict[str, int]] = {}
    for plate, samples in sorted(by_plate.items()):
        blanks = [s for s in samples if sample_role(s) == "negative_control"]
        mocks = [s for s in samples if zymo_standard(s)]
        others = [s for s in samples if s not in blanks and s not in mocks]
        if not blanks:
            raise CohortError(f"extraction plate {plate} has no complete extraction blank")
        if max_biological_per_plate is not None:
            others = stratified_cap(others, max_biological_per_plate)
        chosen = sorted(blanks + mocks + others, key=lambda s: s["sample_name"])
        selected.extend(chosen)
        plate_summary[plate] = {
            "extraction_blanks": len(blanks),
            "zymo_mock_samples": len(mocks),
            "biological_samples": len(others),
        }

    selected_runs = {
        sample["sample_name"]: runs_by_sample[sample["sample_name"]] for sample in selected
    }
    summary = {
        "extraction_kits": kits,
        "biomass_plates": plates,
        "max_biological_per_plate": max_biological_per_plate,
        "plates": plate_summary,
        "samples": len(selected),
        "libraries": len(selected) * len(COMPLETE_ASSAYS),
    }
    return selected, selected_runs, summary


def stratified_cap(samples: list[dict[str, str]], limit: int) -> list[dict[str, str]]:
    """Round-robin across sample types so a cap keeps every material."""
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for sample in sorted(samples, key=lambda s: s["sample_name"]):
        groups[sample.get("sample_type_2") or sample.get("sample_type", "")].append(sample)
    chosen: list[dict[str, str]] = []
    while len(chosen) < limit and any(groups.values()):
        for key in sorted(groups):
            if groups[key] and len(chosen) < limit:
                chosen.append(groups[key].pop(0))
    return chosen


def build_cohort_metadata(
    selected: list[dict[str, str]],
    selected_runs: dict[str, dict[str, dict[str, str]]],
    fastq_paths,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Return TRIARCHY samples, control relationships and mock truth rows.

    `fastq_paths(run)` returns (layout, fastq_r1, fastq_r2) for an ENA run.
    """
    samples: list[dict[str, object]] = []
    library_ids: dict[tuple[str, str], str] = {}
    for sample in selected:
        name = sample["sample_name"]
        sample_id = safe_id(name.removeprefix("12201."))
        role = sample_role(sample)
        standard = zymo_standard(sample)
        subject = sample.get("host_subject_id", "")
        if role != "biological":
            subject = "CONTROL"
        elif standard:
            subject = f"ZYMO_{standard}"
        elif not subject or subject.casefold() in {"not applicable", "not provided", "not collected"}:
            subject = sample_id
        material = (
            f"zymobiomics_{standard.lower()}" if standard
            else safe_id(sample.get("sample_type_2") or sample.get("sample_type") or "unknown").lower()
        )
        environment = safe_id(sample.get("empo_3") or sample.get("env_material") or "unknown").lower()
        replicate = technical_replicate(sample)
        plate = plate_key(sample)
        for assay in COMPLETE_ASSAYS:
            run = selected_runs[name][assay]
            layout, r1, r2 = fastq_paths(run)
            library_id = f"{sample_id}_{assay}"
            library_ids[(name, assay)] = library_id
            ptid = prep_template(run)
            samples.append({
                "library_id": library_id,
                "biological_sample_id": sample_id,
                "specimen_id": sample_id,
                "subject_id": safe_id(subject),
                "assay": assay,
                "sample_role": role,
                "sample_material": material,
                "environment_category": environment,
                "library_layout": layout,
                "fastq_r1": r1,
                "fastq_r2": r2,
                "extraction_batch": plate,
                "pcr_plate": "" if assay == "WMS" else f"{plate}_{assay}",
                "library_batch": f"QIITA_PREP_{ptid}",
                "sequencing_run": f"QIITA_PREP_{ptid}",
                "technical_replicate": replicate,
            })

    controls: list[dict[str, object]] = []
    blanks_by_plate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for sample in selected:
        if sample_role(sample) == "negative_control":
            blanks_by_plate[plate_key(sample)].append(sample)
    for sample in selected:
        if sample_role(sample) != "biological":
            continue
        for blank in blanks_by_plate[plate_key(sample)]:
            for assay in COMPLETE_ASSAYS:
                controls.append({
                    "library_id": library_ids[(sample["sample_name"], assay)],
                    "control_library_id": library_ids[(blank["sample_name"], assay)],
                    "control_stage": "extraction",
                    "relationship": "same_batch",
                    "priority": 1,
                    "reason": (
                        f"extraction blank in well {blank.get('gdna_well_id', '?')} "
                        f"of extraction plate {plate_key(sample)}"
                    ),
                })

    truth: list[dict[str, object]] = []
    for sample in selected:
        standard = zymo_standard(sample)
        if not standard:
            continue
        sample_id = safe_id(sample["sample_name"].removeprefix("12201."))
        for species, kingdom, gdna, s16, s16_18, genome, cells in ZYMO_COMPOSITION[standard]:
            genus = species.split()[0]
            truth.append({
                "biological_sample_id": sample_id,
                "standard": standard,
                "dilution": mock_dilution(sample["sample_name"]),
                "technical_replicate": technical_replicate(sample),
                "species": species,
                "genus": genus,
                "genus_aliases": ";".join(GENUS_ALIASES[genus]),
                "kingdom": kingdom,
                "expected_assays": "WMS;ITS" if kingdom == "Fungi" else "WMS;16S",
                "genomic_dna_pct": gdna,
                "rrna_16s_pct": s16,
                "rrna_16s_18s_pct": s16_18,
                "genome_copy_pct": genome,
                "cell_number_pct": cells,
                "source": ZYMO_SOURCES[standard],
            })
    return samples, controls, truth


def subset_qiita_rows(qiita_rows: list[dict[str, str]], selected: list[dict[str, str]]) -> list[dict[str, str]]:
    wanted = {sample["sample_name"] for sample in selected}
    return [row for row in qiita_rows if row["sample_name"] in wanted]

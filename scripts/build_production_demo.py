#!/usr/bin/env python3
"""Build the deterministic TRIARCHY Production Demo v1.0.

All reads are generated from version-pinned public NCBI nucleotide
records. The resulting manifests are valid production inputs; the normal
TRIARCHY workflow and installed production databases process them.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import random
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data" / "production_demo"
DATASET_NAME = "TRIARCHY Production Demo"
GENERATOR_VERSION = "1.0"
DEFAULT_SEED = 20260827
USER_AGENT = "TRIARCHY-production-demo/1.0"

SAMPLE_FIELDS = [
    "library_id", "biological_sample_id", "specimen_id", "subject_id",
    "assay", "sample_role", "sample_material", "environment_category",
    "library_layout", "fastq_r1", "fastq_r2", "extraction_batch",
    "pcr_plate", "library_batch", "sequencing_run", "technical_replicate",
]
CONTROL_FIELDS = [
    "library_id", "control_library_id", "control_stage", "relationship",
    "priority", "reason",
]
COMPOSITION_FIELDS = [
    "library_id", "assay", "reference_key", "organism", "accession",
    "target_fraction", "read_pairs",
]
REFERENCE_FIELDS = [
    "reference_key", "organism", "ncbi_taxid", "accession", "molecule",
    "sequence_length", "sha256", "source_url",
]
TRUTH_FIELDS = [
    "biological_sample_id", "taxon", "ncbi_taxid", "wms_expected",
    "16s_expected", "its_expected", "expected_interpretation", "rationale",
]
CROSSWALK_FIELDS = [
    "source_database", "source_name", "triarchy_taxon_id", "canonical_name",
    "canonical_rank", "canonical_lineage", "mapping_method",
    "mapping_confidence", "database_version", "database_coverage",
    "reference_status",
]
EXPECTATION_FIELDS = [
    "check_id", "check_type", "target_id", "canonical_name",
    "expected_values", "minimum", "maximum", "required", "rationale",
]
ASSAY_PERFORMANCE_FIELDS = [
    "assay", "kingdom", "reliability", "minimum_depth",
    "reference_coverage", "false_positive_rate",
]


@dataclass(frozen=True)
class ReferenceSpec:
    key: str
    organism: str
    taxid: int
    accession: str
    molecule: str
    minimum_length: int


@dataclass(frozen=True)
class LibrarySpec:
    library_id: str
    biological_sample_id: str
    subject_id: str
    assay: str
    sample_role: str
    sample_material: str
    environment_category: str
    extraction_batch: str
    pcr_plate: str
    library_batch: str
    read_pairs: int
    composition: tuple[tuple[str, float], ...]


REFERENCES = {
    item.key: item
    for item in (
        ReferenceSpec("ecoli_genome", "Escherichia coli K-12 MG1655", 562, "NC_000913.3", "chromosome", 4_000_000),
        ReferenceSpec("pseudomonas_genome", "Pseudomonas aeruginosa PAO1", 287, "NC_002516.2", "chromosome", 6_000_000),
        ReferenceSpec("staphylococcus_genome", "Staphylococcus epidermidis ATCC 12228", 1282, "NC_004461.1", "chromosome", 2_000_000),
        ReferenceSpec("saccharomyces_genome", "Saccharomyces cerevisiae S288C chromosome I", 4932, "NC_001133.9", "chromosome", 200_000),
        ReferenceSpec("human_mitochondrion", "Homo sapiens mitochondrion", 9606, "NC_012920.1", "mitochondrion", 16_000),
        ReferenceSpec("ecoli_16s", "Escherichia coli", 562, "NR_114042.1", "16S rRNA", 1_400),
        ReferenceSpec("pseudomonas_16s", "Pseudomonas aeruginosa", 287, "NR_026078.1", "16S rRNA", 1_500),
        ReferenceSpec("staphylococcus_16s", "Staphylococcus epidermidis", 1282, "NR_113957.1", "16S rRNA", 1_400),
        ReferenceSpec("saccharomyces_its", "Saccharomyces cerevisiae S288C chromosome XII ITS locus", 4932, "NC_001144.5", "chromosome XII ITS locus", 1_000_000),
        ReferenceSpec("candida_its", "Candida albicans CBS 8758 ITS2", 5476, "FJ662402.1", "ITS2", 500),
    )
}


def compositions(*items: tuple[str, float]) -> tuple[tuple[str, float], ...]:
    total = sum(weight for _, weight in items)
    if not math.isclose(total, 1.0, abs_tol=1e-9):
        raise ValueError(f"composition weights sum to {total}, not 1")
    return items


LIBRARIES = (
    LibrarySpec("BIO_A_WMS", "BIO_A", "SUBJ_A", "WMS", "biological", "indoor_air", "built_environment", "EB01", "", "LIB_WMS01", 80_000, compositions(("ecoli_genome", .45), ("pseudomonas_genome", .35), ("saccharomyces_genome", .18), ("staphylococcus_genome", .02))),
    LibrarySpec("BIO_A_16S", "BIO_A", "SUBJ_A", "16S", "biological", "indoor_air", "built_environment", "EB01", "PCR16_01", "LIB_16S01", 8_000, compositions(("ecoli_16s", .55), ("pseudomonas_16s", .43), ("staphylococcus_16s", .02))),
    LibrarySpec("BIO_A_ITS", "BIO_A", "SUBJ_A", "ITS", "biological", "indoor_air", "built_environment", "EB01", "PCRITS_01", "LIB_ITS01", 8_000, compositions(("saccharomyces_its", .75), ("candida_its", .25))),
    LibrarySpec("BIO_B_WMS", "BIO_B", "SUBJ_B", "WMS", "biological", "settled_dust", "built_environment", "EB02", "", "LIB_WMS01", 80_000, compositions(("ecoli_genome", .20), ("pseudomonas_genome", .20), ("saccharomyces_genome", .43), ("staphylococcus_genome", .02), ("human_mitochondrion", .15))),
    LibrarySpec("BIO_B_16S", "BIO_B", "SUBJ_B", "16S", "biological", "settled_dust", "built_environment", "EB02", "PCR16_01", "LIB_16S01", 8_000, compositions(("ecoli_16s", .45), ("pseudomonas_16s", .50), ("staphylococcus_16s", .05))),
    LibrarySpec("BIO_B_ITS", "BIO_B", "SUBJ_B", "ITS", "biological", "settled_dust", "built_environment", "EB02", "PCRITS_01", "LIB_ITS01", 8_000, compositions(("saccharomyces_its", .80), ("candida_its", .20))),
    LibrarySpec("BIO_LOW_WMS", "BIO_LOW", "SUBJ_LOW", "WMS", "biological", "low_biomass_air", "built_environment", "EB01", "", "LIB_WMS01", 30_000, compositions(("ecoli_genome", .45), ("pseudomonas_genome", .18), ("saccharomyces_genome", .10), ("staphylococcus_genome", .27))),
    LibrarySpec("BIO_LOW_16S", "BIO_LOW", "SUBJ_LOW", "16S", "biological", "low_biomass_air", "built_environment", "EB01", "PCR16_01", "LIB_16S01", 5_000, compositions(("ecoli_16s", .20), ("pseudomonas_16s", .10), ("staphylococcus_16s", .70))),
    LibrarySpec("BIO_LOW_ITS", "BIO_LOW", "SUBJ_LOW", "ITS", "biological", "low_biomass_air", "built_environment", "EB01", "PCRITS_01", "LIB_ITS01", 5_000, compositions(("saccharomyces_its", .65), ("candida_its", .35))),
    LibrarySpec("EXT_NEG_EB01_WMS", "EXT_NEG_EB01", "CONTROL", "WMS", "negative_control", "extraction_blank", "laboratory_control", "EB01", "", "LIB_WMS01", 10_000, compositions(("staphylococcus_genome", .95), ("ecoli_genome", .05))),
    LibrarySpec("EXT_NEG_EB01_16S", "EXT_NEG_EB01", "CONTROL", "16S", "negative_control", "extraction_blank", "laboratory_control", "EB01", "PCR16_01", "LIB_16S01", 4_000, compositions(("staphylococcus_16s", .95), ("ecoli_16s", .05))),
    LibrarySpec("EXT_NEG_EB01_ITS", "EXT_NEG_EB01", "CONTROL", "ITS", "negative_control", "extraction_blank", "laboratory_control", "EB01", "PCRITS_01", "LIB_ITS01", 1_000, compositions(("saccharomyces_its", .80), ("candida_its", .20))),
    LibrarySpec("EXT_NEG_EB02_WMS", "EXT_NEG_EB02", "CONTROL", "WMS", "negative_control", "extraction_blank", "laboratory_control", "EB02", "", "LIB_WMS01", 10_000, compositions(("staphylococcus_genome", .95), ("saccharomyces_genome", .05))),
    LibrarySpec("EXT_NEG_EB02_16S", "EXT_NEG_EB02", "CONTROL", "16S", "negative_control", "extraction_blank", "laboratory_control", "EB02", "PCR16_01", "LIB_16S01", 4_000, compositions(("staphylococcus_16s", .95), ("pseudomonas_16s", .05))),
    LibrarySpec("EXT_NEG_EB02_ITS", "EXT_NEG_EB02", "CONTROL", "ITS", "negative_control", "extraction_blank", "laboratory_control", "EB02", "PCRITS_01", "LIB_ITS01", 1_000, compositions(("saccharomyces_its", .85), ("candida_its", .15))),
    LibrarySpec("PCR16_NEG_16S", "PCR16_NEG", "CONTROL", "16S", "negative_control", "pcr_blank", "laboratory_control", "EB_CONTROL", "PCR16_01", "LIB_16S01", 4_000, compositions(("staphylococcus_16s", .98), ("ecoli_16s", .02))),
    LibrarySpec("PCRITS_NEG_ITS", "PCRITS_NEG", "CONTROL", "ITS", "negative_control", "pcr_blank", "laboratory_control", "EB_CONTROL", "PCRITS_01", "LIB_ITS01", 1_000, compositions(("saccharomyces_its", .90), ("candida_its", .10))),
    LibrarySpec("LIB_NEG_WMS", "LIB_NEG", "CONTROL", "WMS", "negative_control", "library_blank", "laboratory_control", "EB_CONTROL", "", "LIB_WMS01", 10_000, compositions(("staphylococcus_genome", .98), ("ecoli_genome", .02))),
    LibrarySpec("POS_MIXED_WMS", "POS_MIXED", "CONTROL", "WMS", "positive_control", "defined_mixed_community", "laboratory_control", "EB01", "", "LIB_WMS01", 30_000, compositions(("ecoli_genome", .40), ("pseudomonas_genome", .35), ("saccharomyces_genome", .25))),
    LibrarySpec("POS_MIXED_16S", "POS_MIXED", "CONTROL", "16S", "positive_control", "defined_mixed_community", "laboratory_control", "EB01", "PCR16_01", "LIB_16S01", 5_000, compositions(("ecoli_16s", .50), ("pseudomonas_16s", .50))),
    LibrarySpec("POS_MIXED_ITS", "POS_MIXED", "CONTROL", "ITS", "positive_control", "defined_mixed_community", "laboratory_control", "EB01", "PCRITS_01", "LIB_ITS01", 5_000, compositions(("saccharomyces_its", .60), ("candida_its", .40))),
)


IUPAC = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"),
    "T": frozenset("T"), "R": frozenset("AG"), "Y": frozenset("CT"),
    "S": frozenset("CG"), "W": frozenset("AT"), "K": frozenset("GT"),
    "M": frozenset("AC"), "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"), "N": frozenset("ACGT"),
}
COMPLEMENT = str.maketrans("ACGTRYMKWSBDHVN", "TGCAYRKMWSVHDBN")


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(COMPLEMENT)[::-1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_fasta(text: str) -> tuple[str, str]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError("response is not FASTA")
    sequence = "".join(line for line in lines[1:] if not line.startswith(">"))
    sequence = re.sub(r"\s+", "", sequence).upper().replace("U", "T")
    if not sequence or re.search(r"[^ACGTRYSWKMBDHVN]", sequence):
        raise ValueError("FASTA contains an invalid or empty nucleotide sequence")
    return lines[0][1:], sequence


def reference_urls(accession: str) -> tuple[str, str]:
    encoded = urllib.parse.quote(accession)
    return (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={encoded}&rettype=fasta&retmode=text",
        "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi"
        f"?id={encoded}&db=nuccore&report=fasta&retmode=text",
    )


def fetch_reference(
    spec: ReferenceSpec,
    destination: Path,
    timeout: int,
    retries: int,
    offline: bool,
) -> tuple[str, str, str]:
    if destination.is_file():
        header, sequence = parse_fasta(destination.read_text(encoding="utf-8"))
        if spec.accession in header and len(sequence) >= spec.minimum_length:
            return header, sequence, reference_urls(spec.accession)[0]
    if offline:
        raise RuntimeError(f"offline mode requires cached reference {destination}")

    errors: list[str] = []
    destination.parent.mkdir(parents=True, exist_ok=True)
    curl = shutil.which("curl")
    for mirror_index, url in enumerate(reference_urls(spec.accession), 1):
        if curl:
            partial = destination.with_suffix(destination.suffix + f".part.{mirror_index}")
            command = [
                curl,
                "--fail", "--location", "--http1.1", "--silent", "--show-error",
                "--retry", str(max(retries - 1, 0)), "--retry-all-errors",
                "--retry-delay", "3", "--connect-timeout", str(timeout),
                "--speed-limit", "1024", "--speed-time", "120",
                "--continue-at", "-", "--output", str(partial), url,
            ]
            completed = subprocess.run(command, check=False)
            try:
                header, sequence = parse_fasta(partial.read_text(encoding="utf-8"))
                if spec.accession not in header:
                    raise ValueError(f"response header does not contain {spec.accession}")
                if len(sequence) < spec.minimum_length:
                    raise ValueError(
                        f"{spec.accession} has {len(sequence)} bases; expected at least {spec.minimum_length}"
                    )
                normalized = f">{header}\n" + "\n".join(
                    sequence[index:index + 80] for index in range(0, len(sequence), 80)
                ) + "\n"
                destination.parent.mkdir(parents=True, exist_ok=True)
                normalized_partial = destination.with_suffix(destination.suffix + ".normalized")
                normalized_partial.write_text(normalized, encoding="utf-8")
                normalized_partial.replace(destination)
                partial.unlink(missing_ok=True)
                time.sleep(.34)
                return header, sequence, url
            except (OSError, ValueError) as error:
                errors.append(f"{url} curl exit {completed.returncode}: {error}")
                partial.unlink(missing_ok=True)
            continue

        for attempt in range(1, retries + 1):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    text = response.read().decode("utf-8")
                header, sequence = parse_fasta(text)
                if spec.accession not in header:
                    raise ValueError(f"response header does not contain {spec.accession}")
                if len(sequence) < spec.minimum_length:
                    raise ValueError(
                        f"{spec.accession} has {len(sequence)} bases; expected at least {spec.minimum_length}"
                    )
                normalized = f">{header}\n" + "\n".join(
                    sequence[index:index + 80] for index in range(0, len(sequence), 80)
                ) + "\n"
                destination.parent.mkdir(parents=True, exist_ok=True)
                partial = destination.with_suffix(destination.suffix + ".part")
                partial.write_text(normalized, encoding="utf-8")
                partial.replace(destination)
                time.sleep(.34)
                return header, sequence, url
            except (urllib.error.URLError, TimeoutError, OSError, ValueError) as error:
                errors.append(f"{url} attempt {attempt}: {error}")
                if attempt < retries:
                    time.sleep(min(2**attempt, 15))
    raise RuntimeError(
        f"unable to download NCBI reference {spec.accession}: " + "; ".join(errors)
    )


def load_primers(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8")
    primers = {}
    for key in (
        "cutadapt_16s_forward", "cutadapt_16s_reverse",
        "cutadapt_its_forward", "cutadapt_its_reverse",
    ):
        match = re.search(rf"\b{key}\s*=\s*['\"]([^'\"]+)['\"]", text, flags=re.IGNORECASE)
        if not match:
            raise RuntimeError(f"missing {key} in {path}")
        primers[key] = match.group(1).upper()
    return primers


def primer_hits(sequence: str, primer: str, max_mismatches: int = 2) -> list[tuple[int, int]]:
    hits = []
    width = len(primer)
    for start in range(0, len(sequence) - width + 1):
        observed = sequence[start:start + width]
        mismatches = sum(base not in IUPAC[code] for base, code in zip(observed, primer))
        if mismatches <= max_mismatches:
            hits.append((start, mismatches))
    return hits


def extract_amplicon(
    sequence: str,
    forward: str,
    reverse: str,
    target_length: int,
    minimum_length: int,
    maximum_length: int,
) -> tuple[str, dict[str, object]]:
    reverse_binding = reverse_complement(reverse)
    candidates = []
    for orientation, oriented in (("forward", sequence), ("reverse_complement", reverse_complement(sequence))):
        for forward_start, forward_mismatches in primer_hits(oriented, forward):
            for reverse_start, reverse_mismatches in primer_hits(oriented, reverse_binding):
                end = reverse_start + len(reverse_binding)
                length = end - forward_start
                if reverse_start > forward_start and minimum_length <= length <= maximum_length:
                    score = (
                        forward_mismatches + reverse_mismatches,
                        abs(length - target_length),
                        forward_start,
                    )
                    candidates.append((score, oriented[forward_start:end], orientation, forward_start, end))
    if not candidates:
        raise RuntimeError(
            f"no {minimum_length}-{maximum_length} bp amplicon matches configured primers"
        )
    score, amplicon, orientation, start, end = min(candidates, key=lambda item: item[0])
    return amplicon, {
        "orientation": orientation,
        "start": start + 1,
        "end": end,
        "length": len(amplicon),
        "primer_mismatches": score[0],
    }


def allocated_counts(total: int, composition: tuple[tuple[str, float], ...]) -> dict[str, int]:
    raw = [(key, total * weight) for key, weight in composition]
    result = {key: int(math.floor(value)) for key, value in raw}
    remainder = total - sum(result.values())
    order = sorted(raw, key=lambda item: (-(item[1] - math.floor(item[1])), item[0]))
    for key, _ in order[:remainder]:
        result[key] += 1
    return result


@contextmanager
def deterministic_gzip_writer(path: Path) -> Iterator[TextIO]:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".part")
    with partial.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="ascii", newline="\n") as text:
                yield text
    partial.replace(path)


def mutate(sequence: str, rng: random.Random, error_rate: float) -> str:
    if error_rate <= 0:
        return sequence
    bases = "ACGT"
    output = []
    for base in sequence:
        if base in bases and rng.random() < error_rate:
            output.append(rng.choice(bases.replace(base, "")))
        else:
            output.append(base)
    return "".join(output)


def valid_fragment(sequence: str, length: int, rng: random.Random) -> str:
    for _ in range(1000):
        start = rng.randrange(0, len(sequence) - length + 1)
        fragment = sequence[start:start + length]
        if set(fragment) <= set("ACGT"):
            return fragment
    raise RuntimeError("unable to sample an unambiguous fragment from reference")


def write_library(
    library: LibrarySpec,
    references: dict[str, str],
    amplicons: dict[str, str],
    reads_dir: Path,
    seed: int,
    scale: float,
    error_rate: float,
) -> tuple[Path, Path, list[dict[str, object]]]:
    total_pairs = max(1, int(round(library.read_pairs * scale)))
    counts = allocated_counts(total_pairs, library.composition)
    schedule = [key for key, count in counts.items() for _ in range(count)]
    library_seed = int.from_bytes(
        hashlib.sha256(f"{seed}:{library.library_id}".encode()).digest()[:8], "big"
    )
    rng = random.Random(library_seed)
    rng.shuffle(schedule)
    r1_path = reads_dir / f"{library.library_id}_R1.fastq.gz"
    r2_path = reads_dir / f"{library.library_id}_R2.fastq.gz"
    read_length = 150 if library.assay == "WMS" else 250

    with deterministic_gzip_writer(r1_path) as r1_handle, deterministic_gzip_writer(r2_path) as r2_handle:
        for pair_number, reference_key in enumerate(schedule, 1):
            if library.assay == "WMS":
                source = references[reference_key]
                fragment_length = rng.randint(320, 420)
                fragment = valid_fragment(source, fragment_length, rng)
                if rng.random() < .5:
                    fragment = reverse_complement(fragment)
                read1 = fragment[:read_length]
                read2 = reverse_complement(fragment[-read_length:])
            else:
                source = amplicons[reference_key]
                read1 = source[:read_length]
                read2 = reverse_complement(source[-read_length:])
            read1 = mutate(read1, rng, error_rate)
            read2 = mutate(read2, rng, error_rate)
            name = f"{library.library_id}:{pair_number:08d}:{reference_key}"
            r1_handle.write(f"@{name}/1\n{read1}\n+\n{'I' * len(read1)}\n")
            r2_handle.write(f"@{name}/2\n{read2}\n+\n{'I' * len(read2)}\n")

    composition_rows = [
        {
            "library_id": library.library_id,
            "assay": library.assay,
            "reference_key": key,
            "organism": REFERENCES[key].organism,
            "accession": REFERENCES[key].accession,
            "target_fraction": f"{weight:.6g}",
            "read_pairs": counts[key],
        }
        for key, weight in library.composition
    ]
    return r1_path, r2_path, composition_rows


def sample_rows(reads_dir: Path) -> list[dict[str, object]]:
    rows = []
    for library in LIBRARIES:
        r1 = reads_dir / f"{library.library_id}_R1.fastq.gz"
        r2 = reads_dir / f"{library.library_id}_R2.fastq.gz"
        rows.append({
            "library_id": library.library_id,
            "biological_sample_id": library.biological_sample_id,
            "specimen_id": library.biological_sample_id,
            "subject_id": library.subject_id,
            "assay": library.assay,
            "sample_role": library.sample_role,
            "sample_material": library.sample_material,
            "environment_category": library.environment_category,
            "library_layout": "paired",
            "fastq_r1": Path("..") / "reads" / r1.name,
            "fastq_r2": Path("..") / "reads" / r2.name,
            "extraction_batch": library.extraction_batch,
            "pcr_plate": library.pcr_plate,
            "library_batch": library.library_batch,
            "sequencing_run": "SEQ_DEMO_01",
            "technical_replicate": "1",
        })
    return rows


def control_rows() -> list[dict[str, object]]:
    rows = []
    sources = [library for library in LIBRARIES if library.sample_role in {"biological", "positive_control"}]
    for library in sources:
        extraction_control = f"EXT_NEG_{library.extraction_batch}_{library.assay}"
        rows.append({
            "library_id": library.library_id,
            "control_library_id": extraction_control,
            "control_stage": "extraction",
            "relationship": "same_batch",
            "priority": 1,
            "reason": f"Matched extraction blank for {library.extraction_batch}",
        })
        if library.assay == "16S":
            rows.append({
                "library_id": library.library_id,
                "control_library_id": "PCR16_NEG_16S",
                "control_stage": "pcr",
                "relationship": "same_batch",
                "priority": 2,
                "reason": "16S-specific PCR blank from PCR16_01",
            })
        elif library.assay == "ITS":
            rows.append({
                "library_id": library.library_id,
                "control_library_id": "PCRITS_NEG_ITS",
                "control_stage": "pcr",
                "relationship": "same_batch",
                "priority": 2,
                "reason": "ITS-specific PCR blank from PCRITS_01",
            })
        else:
            rows.append({
                "library_id": library.library_id,
                "control_library_id": "LIB_NEG_WMS",
                "control_stage": "library",
                "relationship": "same_batch",
                "priority": 2,
                "reason": "WMS-specific library blank from LIB_WMS01",
            })
    return rows


def truth_rows() -> list[dict[str, object]]:
    return [
        {"biological_sample_id": "BIO_A", "taxon": "Escherichia coli", "ncbi_taxid": 562, "wms_expected": "detected", "16s_expected": "detected", "its_expected": "not_applicable", "expected_interpretation": "concordant_detection", "rationale": "Authentic bacterial signal in WMS and 16S."},
        {"biological_sample_id": "BIO_A", "taxon": "Pseudomonas aeruginosa", "ncbi_taxid": 287, "wms_expected": "detected", "16s_expected": "detected", "its_expected": "not_applicable", "expected_interpretation": "concordant_detection", "rationale": "Authentic bacterial signal in WMS and 16S."},
        {"biological_sample_id": "BIO_A", "taxon": "Saccharomyces cerevisiae", "ncbi_taxid": 4932, "wms_expected": "detected", "16s_expected": "not_applicable", "its_expected": "detected", "expected_interpretation": "concordant_detection", "rationale": "Authentic fungal signal in WMS and ITS."},
        {"biological_sample_id": "BIO_A", "taxon": "Candida albicans", "ncbi_taxid": 5476, "wms_expected": "not_detected", "16s_expected": "not_applicable", "its_expected": "detected", "expected_interpretation": "assay_limited_nondetection", "rationale": "ITS-positive fungal target with inadequate WMS fungal detection power."},
        {"biological_sample_id": "BIO_B", "taxon": "Homo sapiens", "ncbi_taxid": 9606, "wms_expected": "host_removed", "16s_expected": "not_applicable", "its_expected": "not_applicable", "expected_interpretation": "host_removed", "rationale": "Fifteen percent of generated WMS pairs originate from human mitochondrion."},
        {"biological_sample_id": "BIO_LOW", "taxon": "Staphylococcus epidermidis", "ncbi_taxid": 1282, "wms_expected": "detected", "16s_expected": "detected", "its_expected": "not_applicable", "expected_interpretation": "probable_contaminant", "rationale": "High relative abundance in a low-biomass sample and stronger linked blank signal."},
        {"biological_sample_id": "POS_MIXED", "taxon": "defined positive-control community", "ncbi_taxid": "", "wms_expected": "E. coli;P. aeruginosa;S. cerevisiae", "16s_expected": "E. coli;P. aeruginosa", "its_expected": "S. cerevisiae;C. albicans", "expected_interpretation": "expected_recovery", "rationale": "Positive-control recovery across all applicable assays."},
    ]


def crosswalk_rows() -> list[dict[str, object]]:
    taxa = (
        ("Escherichia", "TRI-NCBI-561", "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia", ("Escherichia coli", "Escherichia", "Escherichia-Shigella")),
        ("Pseudomonas", "TRI-NCBI-286", "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas", ("Pseudomonas aeruginosa", "Pseudomonas")),
        ("Staphylococcus", "TRI-NCBI-1279", "Bacteria;Bacillati;Bacillota;Bacilli;Caryophanales;Staphylococcaceae;Staphylococcus", ("Staphylococcus epidermidis", "Staphylococcus")),
        ("Saccharomyces", "TRI-NCBI-4930", "Eukaryota;Fungi;Ascomycota;Saccharomycetes;Saccharomycetales;Saccharomycetaceae;Saccharomyces", ("Saccharomyces cerevisiae", "Saccharomyces")),
        ("Candida", "TRI-NCBI-5475", "Eukaryota;Fungi;Ascomycota;Saccharomycetes;Saccharomycetales;Debaryomycetaceae;Candida", ("Candida albicans", "Candida")),
    )
    rows = []
    for canonical_name, taxon_id, lineage, aliases in taxa:
        for assay in ("WMS", "16S", "ITS"):
            for alias in aliases:
                rows.append({
                    "source_database": assay,
                    "source_name": alias,
                    "triarchy_taxon_id": taxon_id,
                    "canonical_name": canonical_name,
                    "canonical_rank": "genus",
                    "canonical_lineage": lineage,
                    "mapping_method": "demo_truth_shared_genus",
                    "mapping_confidence": "1",
                    "database_version": "production_installed",
                    "database_coverage": "0.95" if assay != "WMS" or canonical_name != "Saccharomyces" else "0.45",
                    "reference_status": "version_pinned_public_reference",
                })
    return rows


def expectation_rows() -> list[dict[str, object]]:
    return [
        {"check_id": "BIO_A_ECOLI", "check_type": "classification", "target_id": "BIO_A", "canonical_name": "Escherichia", "expected_values": "concordant_detection|rank_limited_agreement", "minimum": "", "maximum": "", "required": "true", "rationale": "WMS and 16S concordance."},
        {"check_id": "BIO_A_PSEUDOMONAS", "check_type": "classification", "target_id": "BIO_A", "canonical_name": "Pseudomonas", "expected_values": "concordant_detection|rank_limited_agreement", "minimum": "", "maximum": "", "required": "true", "rationale": "WMS and 16S concordance."},
        {"check_id": "BIO_A_SACCHAROMYCES", "check_type": "classification", "target_id": "BIO_A", "canonical_name": "Saccharomyces", "expected_values": "concordant_detection|rank_limited_agreement", "minimum": "", "maximum": "", "required": "true", "rationale": "WMS and ITS concordance."},
        {"check_id": "BIO_A_CANDIDA", "check_type": "classification", "target_id": "BIO_A", "canonical_name": "Candida", "expected_values": "assay_limited_nondetection", "minimum": "", "maximum": "", "required": "true", "rationale": "ITS detection with underpowered WMS fungal profiling."},
        {"check_id": "BIO_LOW_STAPH", "check_type": "classification", "target_id": "BIO_LOW", "canonical_name": "Staphylococcus", "expected_values": "probable_contaminant", "minimum": "", "maximum": "", "required": "true", "rationale": "Linked negative-control abundance dominates low-biomass evidence."},
        {"check_id": "BIO_B_HOST", "check_type": "host_fraction", "target_id": "BIO_B_WMS", "canonical_name": "Homo sapiens", "expected_values": "", "minimum": ".08", "maximum": ".25", "required": "true", "rationale": "Approximately 15% human mitochondrial pairs should be removed."},
        {"check_id": "POS_WMS_ECOLI", "check_type": "profile_detection", "target_id": "POS_MIXED_WMS", "canonical_name": "Escherichia", "expected_values": "detected", "minimum": "", "maximum": "", "required": "true", "rationale": "Positive-control WMS recovery."},
        {"check_id": "POS_WMS_PSEUDOMONAS", "check_type": "profile_detection", "target_id": "POS_MIXED_WMS", "canonical_name": "Pseudomonas", "expected_values": "detected", "minimum": "", "maximum": "", "required": "true", "rationale": "Positive-control WMS recovery."},
        {"check_id": "POS_WMS_SACCHAROMYCES", "check_type": "profile_detection", "target_id": "POS_MIXED_WMS", "canonical_name": "Saccharomyces", "expected_values": "detected", "minimum": "", "maximum": "", "required": "true", "rationale": "Positive-control WMS recovery."},
        {"check_id": "POS_16S_ECOLI", "check_type": "profile_detection", "target_id": "POS_MIXED_16S", "canonical_name": "Escherichia", "expected_values": "detected", "minimum": "", "maximum": "", "required": "true", "rationale": "Positive-control 16S recovery."},
        {"check_id": "POS_16S_PSEUDOMONAS", "check_type": "profile_detection", "target_id": "POS_MIXED_16S", "canonical_name": "Pseudomonas", "expected_values": "detected", "minimum": "", "maximum": "", "required": "true", "rationale": "Positive-control 16S recovery."},
        {"check_id": "POS_ITS_SACCHAROMYCES", "check_type": "profile_detection", "target_id": "POS_MIXED_ITS", "canonical_name": "Saccharomyces", "expected_values": "detected", "minimum": "", "maximum": "", "required": "true", "rationale": "Positive-control ITS recovery."},
        {"check_id": "POS_ITS_CANDIDA", "check_type": "profile_detection", "target_id": "POS_MIXED_ITS", "canonical_name": "Candida", "expected_values": "detected", "minimum": "", "maximum": "", "required": "true", "rationale": "Positive-control ITS recovery."},
    ]


def assay_performance_rows() -> list[dict[str, object]]:
    """Use production defaults with an explicit demo-study fungal prior."""
    rows = []
    with (ROOT / "conf" / "assay_performance.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        for source in csv.DictReader(handle, delimiter="\t"):
            row = dict(source)
            if row["assay"] == "WMS" and row["kingdom"] == "Fungi":
                row["reference_coverage"] = "0.65"
            rows.append(row)
    return rows


def validate_fastq(path: Path) -> int:
    lines = 0
    with gzip.open(path, "rt", encoding="ascii") as handle:
        while True:
            record = [handle.readline() for _ in range(4)]
            if not record[0]:
                break
            if any(not line for line in record) or not record[0].startswith("@") or not record[2].startswith("+"):
                raise RuntimeError(f"invalid FASTQ record in {path}")
            sequence = record[1].strip()
            quality = record[3].strip()
            if len(sequence) != len(quality):
                raise RuntimeError(f"sequence/quality length mismatch in {path}")
            lines += 1
    if not lines:
        raise RuntimeError(f"empty FASTQ: {path}")
    return lines


def ready_dataset(output: Path, seed: int, scale: float, error_rate: float) -> bool:
    manifest = output / "dataset.json"
    if not manifest.is_file():
        return False
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    expected = payload.get("parameters", {})
    if payload.get("dataset") != DATASET_NAME or payload.get("version") != GENERATOR_VERSION:
        return False
    if expected != {"seed": seed, "scale": scale, "error_rate": error_rate}:
        return False
    reads_ready = all(
        (output / "reads" / f"{library.library_id}_{mate}.fastq.gz").is_file()
        for library in LIBRARIES for mate in ("R1", "R2")
    )
    metadata_ready = all(
        (output / "metadata" / filename).is_file()
        for filename in (
            "samples.tsv", "control_relationships.tsv", "library_composition.tsv",
            "references.tsv", "truth.tsv", "taxonomy_crosswalk.tsv",
            "expected_outcomes.tsv", "assay_performance.tsv",
        )
    )
    references_ready = all(
        (output / "references" / f"{spec.accession}.fasta").is_file()
        for spec in REFERENCES.values()
    )
    return reads_ready and metadata_ready and references_ready and (
        output / "checksums.sha256"
    ).is_file()


def build(args: argparse.Namespace) -> dict[str, object]:
    output = args.output.resolve()
    if not args.force and ready_dataset(output, args.seed, args.scale, args.error_rate):
        print(f"TRIARCHY production demo dataset already exists: {output}")
        return json.loads((output / "dataset.json").read_text(encoding="utf-8"))

    references_dir = output / "references"
    reads_dir = output / "reads"
    metadata_dir = output / "metadata"
    references: dict[str, str] = {}
    reference_rows: list[dict[str, object]] = []
    for index, spec in enumerate(REFERENCES.values(), 1):
        destination = references_dir / f"{spec.accession}.fasta"
        print(f"[{index:02d}/{len(REFERENCES):02d}] NCBI {spec.accession} — {spec.organism}")
        _, sequence, source_url = fetch_reference(
            spec, destination, timeout=args.timeout, retries=args.retries, offline=args.offline
        )
        references[spec.key] = sequence
        reference_rows.append({
            "reference_key": spec.key,
            "organism": spec.organism,
            "ncbi_taxid": spec.taxid,
            "accession": spec.accession,
            "molecule": spec.molecule,
            "sequence_length": len(sequence),
            "sha256": sha256(destination),
            "source_url": source_url,
        })

    primers = load_primers(ROOT / "conf" / "primers.config")
    amplicons: dict[str, str] = {}
    amplicon_metadata: dict[str, dict[str, object]] = {}
    for key in ("ecoli_16s", "pseudomonas_16s", "staphylococcus_16s"):
        amplicons[key], amplicon_metadata[key] = extract_amplicon(
            references[key], primers["cutadapt_16s_forward"],
            primers["cutadapt_16s_reverse"], 460, 350, 600,
        )
    for key in ("saccharomyces_its", "candida_its"):
        amplicons[key], amplicon_metadata[key] = extract_amplicon(
            references[key], primers["cutadapt_its_forward"],
            primers["cutadapt_its_reverse"], 350, 180, 650,
        )

    composition_output: list[dict[str, object]] = []
    fastq_checks = []
    for index, library in enumerate(LIBRARIES, 1):
        pairs = max(1, int(round(library.read_pairs * args.scale)))
        print(f"[{index:02d}/{len(LIBRARIES):02d}] {library.library_id}: {pairs:,} paired reads")
        r1, r2, rows = write_library(
            library, references, amplicons, reads_dir, args.seed, args.scale,
            args.error_rate,
        )
        count_r1 = validate_fastq(r1)
        count_r2 = validate_fastq(r2)
        if count_r1 != pairs or count_r2 != pairs:
            raise RuntimeError(f"FASTQ count mismatch for {library.library_id}")
        composition_output.extend(rows)
        fastq_checks.extend((r1, r2))

    samples = sample_rows(reads_dir)
    controls = control_rows()
    write_tsv(metadata_dir / "samples.tsv", SAMPLE_FIELDS, samples)
    write_tsv(metadata_dir / "control_relationships.tsv", CONTROL_FIELDS, controls)
    write_tsv(metadata_dir / "library_composition.tsv", COMPOSITION_FIELDS, composition_output)
    write_tsv(metadata_dir / "references.tsv", REFERENCE_FIELDS, reference_rows)
    write_tsv(metadata_dir / "truth.tsv", TRUTH_FIELDS, truth_rows())
    write_tsv(metadata_dir / "taxonomy_crosswalk.tsv", CROSSWALK_FIELDS, crosswalk_rows())
    write_tsv(metadata_dir / "expected_outcomes.tsv", EXPECTATION_FIELDS, expectation_rows())
    write_tsv(
        metadata_dir / "assay_performance.tsv",
        ASSAY_PERFORMANCE_FIELDS,
        assay_performance_rows(),
    )

    sys.path.insert(0, str(ROOT / "scripts"))
    from validate_metadata import validate  # pylint: disable=import-outside-toplevel

    library_count, control_count = validate(
        metadata_dir / "samples.tsv", metadata_dir / "control_relationships.tsv"
    )
    checksum_paths = sorted(
        fastq_checks
        + list(metadata_dir.glob("*.tsv"))
        + list(references_dir.glob("*.fasta"))
    )
    checksum_text = "".join(
        f"{sha256(path)}  {path.relative_to(output)}\n" for path in checksum_paths
    )
    (output / "checksums.sha256").write_text(checksum_text, encoding="utf-8")

    payload = {
        "dataset": DATASET_NAME,
        "version": GENERATOR_VERSION,
        "parameters": {
            "seed": args.seed,
            "scale": args.scale,
            "error_rate": args.error_rate,
        },
        "libraries": library_count,
        "control_relationships": control_count,
        "biological_samples": ["BIO_A", "BIO_B", "BIO_LOW"],
        "reference_accessions": [spec.accession for spec in REFERENCES.values()],
        "amplicons": amplicon_metadata,
        "read_pairs": sum(max(1, int(round(item.read_pairs * args.scale))) for item in LIBRARIES),
        "production_mode": True,
        "limitations": [
            "Reads are synthetic observations generated from authentic public reference sequences.",
            "The dataset tests pipeline behavior and evidence logic, not biological accuracy or limit of detection.",
        ],
    }
    (output / "dataset.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    (output / "README.md").write_text(
        "# TRIARCHY Production Demo v1.0\n\n"
        "Deterministic paired reads generated from version-pinned NCBI nucleotide records. "
        "These are normal production inputs and require the complete TRIARCHY production databases.\n\n"
        "Run from the TRIARCHY repository with:\n\n"
        "```bash\ntriarchy demo --profile laptop --resume\n```\n\n"
        "See `metadata/truth.tsv`, `metadata/library_composition.tsv`, and "
        "`metadata/expected_outcomes.tsv` for the designed evidence patterns.\n",
        encoding="utf-8",
    )
    print(
        f"VALID: {library_count} production libraries; {control_count} control relationships; "
        f"{payload['read_pairs']:,} paired reads."
    )
    print(f"Dataset: {output}")
    return payload


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    result.add_argument("--seed", type=int, default=DEFAULT_SEED)
    result.add_argument("--scale", type=float, default=1.0, help="Scale all library read counts; default 1.0")
    result.add_argument("--error-rate", type=float, default=.001, help="Per-base substitution rate; default 0.001")
    result.add_argument("--timeout", type=int, default=60)
    result.add_argument("--retries", type=int, default=3)
    result.add_argument("--offline", action="store_true", help="Require cached NCBI references")
    result.add_argument("--force", action="store_true", help="Regenerate FASTQs and metadata")
    return result


def main() -> int:
    args = parser().parse_args()
    if args.scale <= 0:
        raise SystemExit("--scale must be positive")
    if not 0 <= args.error_rate <= .05:
        raise SystemExit("--error-rate must be between 0 and 0.05")
    if args.retries < 1 or args.timeout < 1:
        raise SystemExit("--retries and --timeout must be positive")
    try:
        build(args)
    except (OSError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

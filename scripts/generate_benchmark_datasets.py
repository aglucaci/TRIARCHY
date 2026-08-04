#!/usr/bin/env python3
"""Generate deterministic, repository-sized TRIARCHY benchmark datasets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import random
from dataclasses import dataclass, field
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "benchmarks"
ASSAYS = ("WMS", "16S", "ITS")
SAMPLE_FIELDS = [
    "library_id", "biological_sample_id", "specimen_id", "subject_id", "assay",
    "sample_role", "sample_material", "environment_category", "library_layout",
    "fastq_r1", "fastq_r2", "extraction_batch", "pcr_plate", "library_batch",
    "sequencing_run", "technical_replicate", "biomass_measure", "biomass_unit",
]
CONTROL_FIELDS = [
    "library_id", "control_library_id", "control_stage", "relationship",
    "priority", "reason",
]
TRUTH_FIELDS = [
    "library_id", "biological_sample_id", "assay", "taxon", "kingdom",
    "genus", "expected_rank", "expected_fraction", "truth_status",
]


@dataclass(frozen=True)
class Taxon:
    name: str
    kingdom: str
    genus: str
    accession: str
    genome_size_mb: float
    marker_16s_copies: float = 0
    marker_its_copies: float = 0
    recovery_16s: float = 1
    recovery_its: float = 1
    reference_status: str = "reference_seen"
    synonym: str = ""


TAXA_LIST = (
    Taxon("Listeria_monocytogenes", "Bacteria", "Listeria", "GCF_000196035.1", 3.0, 6),
    Taxon("Pseudomonas_aeruginosa", "Bacteria", "Pseudomonas", "GCF_000006765.1", 6.3, 4),
    Taxon("Bacillus_subtilis", "Bacteria", "Bacillus", "GCF_000009045.1", 4.2, 10),
    Taxon("Escherichia_coli", "Bacteria", "Escherichia", "GCF_000005845.2", 4.6, 7),
    Taxon("Salmonella_enterica", "Bacteria", "Salmonella", "GCF_000006945.2", 4.8, 7),
    Taxon("Lactobacillus_fermentum", "Bacteria", "Limosilactobacillus", "GCF_000010145.1", 2.0, 5, synonym="Limosilactobacillus_fermentum"),
    Taxon("Enterococcus_faecalis", "Bacteria", "Enterococcus", "GCF_000007785.1", 3.2, 4),
    Taxon("Staphylococcus_aureus", "Bacteria", "Staphylococcus", "GCF_000013425.1", 2.8, 5),
    Taxon("Ralstonia_pickettii", "Bacteria", "Ralstonia", "GCF_000020205.1", 5.1, 4),
    Taxon("Cutibacterium_acnes", "Bacteria", "Cutibacterium", "GCF_000008345.1", 2.6, 3, synonym="Propionibacterium_acnes"),
    Taxon("Sphingomonas_paucimobilis", "Bacteria", "Sphingomonas", "GCF_000021825.1", 4.0, 2),
    Taxon("Burkholderia_cepacia", "Bacteria", "Burkholderia", "GCF_000009485.1", 8.7, 6),
    Taxon("Shigella_flexneri", "Bacteria", "Shigella", "GCF_000006925.2", 4.6, 7, recovery_16s=0.45),
    Taxon("Bacillus_cereus", "Bacteria", "Bacillus", "GCF_000007825.1", 5.4, 12),
    Taxon("Deinococcus_radiodurans", "Bacteria", "Deinococcus", "GCF_000008565.1", 3.3, 3, recovery_16s=0.2, reference_status="close_relative_only"),
    Taxon("Novelus_aeris", "Bacteria", "Novelus", "HELDOUT_BACTERIUM_01", 3.8, 2, reference_status="genus_held_out"),
    Taxon("Methanobrevibacter_smithii", "Archaea", "Methanobrevibacter", "GCF_000016525.1", 1.9, 2),
    Taxon("Halobacterium_salinarum", "Archaea", "Halobacterium", "GCF_000006805.1", 2.6, 1),
    Taxon("Saccharomyces_cerevisiae", "Fungi", "Saccharomyces", "GCF_000146045.2", 12.1, marker_its_copies=150),
    Taxon("Cryptococcus_neoformans", "Fungi", "Cryptococcus", "GCF_000149245.1", 19.1, marker_its_copies=55, synonym="Filobasidiella_neoformans"),
    Taxon("Aspergillus_fumigatus", "Fungi", "Aspergillus", "GCF_000002655.1", 29.4, marker_its_copies=45),
    Taxon("Candida_albicans", "Fungi", "Candida", "GCF_000182965.3", 14.3, marker_its_copies=80, synonym="Candida_africana_complex"),
    Taxon("Cladosporium_cladosporioides", "Fungi", "Cladosporium", "GCF_003614015.1", 33.5, marker_its_copies=35),
    Taxon("Penicillium_chrysogenum", "Fungi", "Penicillium", "GCF_000710275.1", 32.2, marker_its_copies=50),
    Taxon("Alternaria_alternata", "Fungi", "Alternaria", "GCF_001642055.1", 33.3, marker_its_copies=30, recovery_its=0.5),
    Taxon("Fusarium_oxysporum", "Fungi", "Fusarium", "GCF_000149955.1", 61.4, marker_its_copies=60),
    Taxon("Malassezia_restricta", "Fungi", "Malassezia", "GCF_003290485.1", 7.4, marker_its_copies=25, recovery_its=0.3),
    Taxon("Novelomyces_aerophilus", "Fungi", "Novelomyces", "HELDOUT_FUNGUS_01", 28.0, marker_its_copies=20, reference_status="genus_held_out"),
)
TAXA = {taxon.name: taxon for taxon in TAXA_LIST}
CONTAMINANTS = {
    "Ralstonia_pickettii": ("universal_reagent", "extraction"),
    "Cutibacterium_acnes": ("skin_associated", "field"),
    "Sphingomonas_paucimobilis": ("water_associated", "extraction"),
    "Burkholderia_cepacia": ("batch_B02", "extraction"),
    "Penicillium_chrysogenum": ("fungal_reagent", "pcr"),
    "Malassezia_restricta": ("index_bleed", "sequencing"),
}
ZYMO = [
    "Listeria_monocytogenes", "Pseudomonas_aeruginosa", "Bacillus_subtilis",
    "Escherichia_coli", "Salmonella_enterica", "Lactobacillus_fermentum",
    "Enterococcus_faecalis", "Staphylococcus_aureus",
    "Saccharomyces_cerevisiae", "Cryptococcus_neoformans",
]
ENVIRONMENTAL = [
    "Pseudomonas_aeruginosa", "Bacillus_subtilis", "Cutibacterium_acnes",
    "Sphingomonas_paucimobilis", "Deinococcus_radiodurans",
    "Methanobrevibacter_smithii", "Aspergillus_fumigatus",
    "Cladosporium_cladosporioides", "Penicillium_chrysogenum",
    "Alternaria_alternata", "Malassezia_restricta",
]


@dataclass(frozen=True)
class Specimen:
    name: str
    role: str
    batch: str
    biomass: float
    environment: str = "benchmark"
    material: str = "synthetic_community"
    assays: tuple[str, ...] = ASSAYS
    cohort: str = ""
    qc: tuple[tuple[str, str], ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class Dataset:
    dataset_id: str
    title: str
    purpose: str
    specimens: tuple[Specimen, ...]
    mode: str
    provenance_class: str = "synthetic"


def specimens_for(dataset_id: str) -> tuple[Specimen, ...]:
    if dataset_id == "TRI-01":
        return (
            Specimen("BIO_A", "biological", "B01", 100, "bacteria_dominated"),
            Specimen("BIO_B", "biological", "B01", 70, "fungal_rich"),
            Specimen("NEG_EXT", "negative_control", "B01", 0.1),
            Specimen("POS_MOCK", "positive_control", "B01", 100),
        )
    if dataset_id == "TRI-02":
        return (
            Specimen("WMS_STRESS", "biological", "B01", 100, qc=(("WMS", "host_20pct_adapter_low_quality_duplicates"),)),
            Specimen("AMP_STRESS", "biological", "B01", 80, qc=(("16S", "primer_mismatch_and_chimeras"),)),
            Specimen("FUNGAL_STRESS", "biological", "B01", 60, qc=(("ITS", "variable_length_intragenomic_variants"),)),
            Specimen("NEG_EXT", "negative_control", "B01", 0.1),
            Specimen("NEG_PCR16", "negative_control", "B01", 0.1, assays=("16S",)),
            Specimen("NEG_PCRITS", "negative_control", "B01", 0.1, assays=("ITS",)),
            Specimen("NEG_LIBWMS", "negative_control", "B01", 0.1, assays=("WMS",)),
            Specimen("POS_MOCK", "positive_control", "B01", 100),
        )
    if dataset_id == "TRI-03":
        return tuple(
            [Specimen(f"BALANCED_{i:02d}", "biological", "B01", 100) for i in range(1, 4)]
            + [Specimen(f"NEG_EXT_{i}", "negative_control", "B01", 0.1) for i in (1, 2)]
            + [Specimen("POS_MIXED", "positive_control", "B01", 100)]
        )
    if dataset_id == "TRI-05":
        items = [
            Specimen(f"CONTAM_{batch}_{level}", "biological", batch, biomass)
            for batch, level, biomass in (
                ("B01", "LOW", 1e3), ("B01", "MED", 1e5),
                ("B02", "LOW", 1e3), ("B02", "HIGH", 1e7), ("B02", "MAX", 1e8),
            )
        ]
        items += [Specimen(f"NEG_EXT_{batch}", "negative_control", batch, 0) for batch in ("B01", "B02")]
        items += [Specimen("POS_MIXED", "positive_control", "B01", 100)]
        return tuple(items)
    if dataset_id == "TRI-06":
        items = [
            Specimen(f"TAX_{challenge}", "biological", "B01", 100, challenge)
            for challenge in ("synonym", "unequal_rank", "historical_fungal_name", "genus_held_out", "one_to_many")
        ]
        items += [Specimen("NEG_EXT", "negative_control", "B01", 0)]
        return tuple(items)
    if dataset_id == "TRI-07":
        return (
            Specimen("F01_COMPLETE", "biological", "B01", 100),
            Specimen("F02_WMS_16S", "biological", "B01", 80, assays=("WMS", "16S")),
            Specimen("F03_WMS_ITS", "biological", "B01", 80, assays=("WMS", "ITS")),
            Specimen("F04_16S_ITS", "biological", "B01", 80, assays=("16S", "ITS")),
            Specimen("F05_FAILED_ITS", "biological", "B02", 60, qc=(("ITS", "failed"),)),
            Specimen("F06_LOW_16S", "biological", "B02", 50, qc=(("16S", "low_depth"),)),
            Specimen("F07_HOST_WMS", "biological", "B02", 40, qc=(("WMS", "host_dominated"),)),
            Specimen("NEG_EMPTY", "negative_control", "B01", 0, qc=tuple((a, "empty") for a in ASSAYS)),
            Specimen("POS_PARTIAL", "positive_control", "B01", 100, qc=(("ITS", "failed"),)),
            Specimen("NEG_PCR16_B02", "negative_control", "B02", 0, assays=("16S",)),
        )
    if dataset_id == "TRI-10":
        return tuple(
            [Specimen(f"EVEN_{i:02d}", "biological", "B01", 100) for i in range(1, 5)]
            + [Specimen(f"NEG_EXT_{i}", "negative_control", "B01", 0.1) for i in range(1, 3)]
            + [Specimen("POS_MOCK", "positive_control", "B01", 100)]
        )
    if dataset_id == "TRI-11":
        return tuple(
            [
                Specimen(f"LOG_{depth}_{rep}", "biological", "B01", biomass)
                for depth, biomass in (("LOW", 1), ("MED", 100), ("HIGH", 10000))
                for rep in (1, 2)
            ]
            + [Specimen("NEG_EXT", "negative_control", "B01", 0.1),
               Specimen("POS_MOCK", "positive_control", "B01", 100)]
        )
    if dataset_id == "TRI-12":
        items = [
            Specimen(f"LOW_{batch}_{level}_{rep}", "biological", batch, float(level))
            for batch in ("B01", "B02")
            for level in (10**3, 10**4, 10**5, 10**6, 10**7, 10**8)
            for rep in range(1, 3)
        ]
        for batch in ("B01", "B02"):
            items += [Specimen(f"FIELD_{batch}_{i}", "negative_control", batch, 0, assays=ASSAYS) for i in (1, 2)]
            items += [Specimen(f"NEG_EXT_{batch}_{i}", "negative_control", batch, 0, assays=ASSAYS) for i in (1, 2)]
            items += [Specimen(f"NEG_PCR16_{batch}_{i}", "negative_control", batch, 0, assays=("16S",)) for i in (1, 2)]
            items += [Specimen(f"NEG_PCRITS_{batch}_{i}", "negative_control", batch, 0, assays=("ITS",)) for i in (1, 2)]
            items += [Specimen(f"NEG_LIBWMS_{batch}", "negative_control", batch, 0, assays=("WMS",))]
            items += [Specimen(f"POS_{batch}", "positive_control", batch, 100)]
        return tuple(items)
    if dataset_id == "TRI-13":
        items = [
            Specimen(f"{habitat}_{rep:02d}", "biological", f"B{1 + (rep % 2):02d}", 10 ** (rep % 4 + 2), habitat)
            for habitat in ("outdoor_air", "supply_air", "return_air", "settled_dust")
            for rep in range(1, 6)
        ]
        items += [
            Specimen(f"NEG_EXT_B{batch:02d}_{rep}", "negative_control", f"B{batch:02d}", 0)
            for batch in (1, 2) for rep in (1, 2)
        ]
        items += [Specimen("POS_MOCK", "positive_control", "B01", 100)]
        return tuple(items)
    if dataset_id == "TRI-14":
        items = [
            Specimen(f"HARM_{challenge}_{rep}", "biological", f"B{rep:02d}", 100, challenge)
            for challenge in ("synonym", "unequal_rank", "close_relative", "held_out")
            for rep in (1, 2, 3)
        ]
        items += [Specimen(f"NEG_EXT_B{i:02d}", "negative_control", f"B{i:02d}", 0) for i in (1, 2, 3)]
        items += [Specimen("POS_MOCK", "positive_control", "B01", 100)]
        return tuple(items)
    if dataset_id == "TRI-15":
        items = [
            Specimen(f"BAL_{batch}_{group}", "biological", batch, 100, group, cohort="balanced")
            for batch in ("B01", "B02", "B03", "B04") for group in ("A", "B", "C", "D")
        ]
        items += [
            Specimen(f"CONF_{batch}_{rep}", "biological", batch, 100, group, cohort="confounded")
            for batch, group in zip(("B01", "B02", "B03", "B04"), ("A", "B", "C", "D"))
            for rep in (1, 2, 3)
        ]
        items += [Specimen(f"NEG_EXT_{batch}", "negative_control", batch, 0) for batch in ("B01", "B02", "B03", "B04")]
        items += [Specimen("NEG_PCR16_B02", "negative_control", "B02", 0, assays=("16S",), qc=(("16S", "carryover_anomalous"),))]
        return tuple(items)
    if dataset_id == "TRI-16":
        return (
            Specimen("S01_COMPLETE", "biological", "B01", 100),
            Specimen("S02_WMS_16S", "biological", "B01", 80, assays=("WMS", "16S")),
            Specimen("S03_WMS_ITS", "biological", "B01", 80, assays=("WMS", "ITS")),
            Specimen("S04_16S_ITS", "biological", "B01", 80, assays=("16S", "ITS")),
            Specimen("S05_WMS_ONLY", "biological", "B01", 80, assays=("WMS",)),
            Specimen("S06_FAILED_ITS", "biological", "B02", 60, qc=(("ITS", "failed"),)),
            Specimen("S07_LOW_16S", "biological", "B02", 50, qc=(("16S", "low_depth"),)),
            Specimen("S08_HOST_WMS", "biological", "B02", 40, qc=(("WMS", "host_dominated"),)),
            Specimen("S09_FWD_ONLY", "biological", "B02", 30, assays=("WMS",), qc=(("WMS", "forward_only"),)),
            Specimen("S10_MISMATCHED", "biological", "B02", 30, assays=("16S",), qc=(("16S", "mismatched_pairs"),)),
            Specimen("NEG_EMPTY", "negative_control", "B01", 0, qc=tuple((a, "empty") for a in ASSAYS)),
            Specimen("POS_PARTIAL", "positive_control", "B01", 100, qc=(("ITS", "failed"),)),
            Specimen("NEG_EXT_B02", "negative_control", "B02", 0),
        )
    if dataset_id == "TRI-17":
        items = []
        for depth, biomass in (("D01", 1), ("D02", 10), ("D03", 100), ("D04", 1000), ("D05", 10000)):
            for replicate in (1, 2):
                qc = tuple((assay, "low_depth") for assay in ASSAYS) if depth == "D01" else ()
                items.append(Specimen(
                    f"SENS_{depth}_R{replicate}", "biological", "B01", biomass,
                    depth, "sensitivity_matrix", cohort="depth_primer_reference", qc=qc,
                ))
        items += [
            Specimen("NEG_EXT", "negative_control", "B01", 0),
            Specimen("POS_MOCK", "positive_control", "B01", 100),
        ]
        return tuple(items)
    if dataset_id == "TRI-18":
        items = [
            Specimen(
                f"CARRY_{batch}_R{replicate}", "biological", batch,
                10 if replicate == 1 else 100, "carryover_challenge",
                "synthetic_community", cohort="sequencing_carryover",
            )
            for batch in ("B01", "B02") for replicate in range(1, 7)
        ]
        for batch in ("B01", "B02"):
            items += [
                Specimen(f"NEG_EXT_{batch}", "negative_control", batch, 0),
                Specimen(f"NEG_PCR16_{batch}", "negative_control", batch, 0, assays=("16S",)),
                Specimen(f"NEG_PCRITS_{batch}", "negative_control", batch, 0, assays=("ITS",)),
                Specimen(f"NEG_LIBWMS_{batch}", "negative_control", batch, 0, assays=("WMS",)),
                Specimen(f"POS_{batch}", "positive_control", batch, 100),
            ]
        return tuple(items)
    if dataset_id == "TRI-19":
        items = [
            Specimen(
                f"REPEAT_{community}_{batch}_R{replicate}", "biological", batch,
                100, community, "synthetic_community", cohort="repeatability",
            )
            for community in ("BACTERIAL", "FUNGAL", "MIXED", "LOW_ABUNDANCE")
            for batch in ("B01", "B02", "B03")
            for replicate in (1, 2)
        ]
        items += [
            Specimen(f"NEG_EXT_{batch}", "negative_control", batch, 0)
            for batch in ("B01", "B02", "B03")
        ]
        items += [Specimen("POS_REPEAT", "positive_control", "B01", 100)]
        return tuple(items)
    if dataset_id == "TRI-20":
        items = [
            Specimen(f"{standard}_REP{rep}", "biological", "B01" if standard != "WHOLE_CELL" else "B02", 100, standard, "commercial_standard")
            for standard in ("EVEN_DNA", "WHOLE_CELL", "LOG_DNA") for rep in (1, 2, 3)
        ]
        items += [Specimen(f"NEG_EXT_{i}", "negative_control", "B02", 0) for i in (1, 2, 3)]
        items += [Specimen(f"NEG_PCR16_{i}", "negative_control", "B01", 0, assays=("16S",)) for i in (1, 2)]
        items += [Specimen(f"NEG_PCRITS_{i}", "negative_control", "B01", 0, assays=("ITS",)) for i in (1, 2)]
        items += [Specimen(f"NEG_LIB_{i}", "negative_control", "B01", 0, assays=("WMS",)) for i in (1, 2)]
        items += [Specimen("POS_REPEAT", "positive_control", "B01", 100)]
        return tuple(items)
    if dataset_id == "TRI-21":
        items = [
            Specimen(f"{profile}_{biomass}_R{rep}", "biological", f"B{1 + rep % 2:02d}", value, profile, "custom_environmental_mock")
            for profile in ("EVEN", "LOG", "FUNGAL_RICH")
            for biomass, value in (("LOW", 1), ("MED", 100), ("HIGH", 10000))
            for rep in (1, 2, 3, 4)
        ]
        items += [Specimen(f"NEG_EXT_B{i:02d}", "negative_control", f"B{i:02d}", 0) for i in (1, 2)]
        return tuple(items)
    if dataset_id == "TRI-30":
        items = [
            Specimen(
                f"{material}_{environment}_R{rep}", "biological",
                f"B{1 + ((rep + material_i + env_i) % 4):02d}",
                1000 if material == "settled_dust" else 10,
                environment, material, cohort="matched_environmental",
            )
            for material_i, material in enumerate(("air_filter", "active_air", "settled_dust", "surface_swab"))
            for env_i, environment in enumerate(("outdoor_intake", "indoor_supply", "indoor_return"))
            for rep in (1, 2, 3, 4)
        ]
        for batch in ("B01", "B02", "B03", "B04"):
            items += [
                Specimen(f"FIELD_{batch}", "negative_control", batch, 0),
                Specimen(f"NEG_EXT_{batch}", "negative_control", batch, 0),
                Specimen(f"NEG_PCR16_{batch}", "negative_control", batch, 0, assays=("16S",)),
                Specimen(f"NEG_PCRITS_{batch}", "negative_control", batch, 0, assays=("ITS",)),
                Specimen(f"NEG_LIBWMS_{batch}", "negative_control", batch, 0, assays=("WMS",)),
                Specimen(f"POS_{batch}", "positive_control", batch, 100),
            ]
        return tuple(items)
    raise ValueError(dataset_id)


DATASET_META = (
    ("TRI-01", "Tiny end-to-end smoke test", "Complete three-assay execution with biological samples and controls.", "smoke", "synthetic"),
    ("TRI-02", "Assay-specific branch tests", "Host, primer, chimera, ITS-length, and unknown-feature stressors.", "branch", "synthetic"),
    ("TRI-03", "Balanced mixed-kingdom truth community", "Replicated bacterial, archaeal, and fungal truth with negative and positive controls.", "even", "synthetic"),
    ("TRI-05", "Compact contamination gradient", "Five-bundle low-biomass and batch-specific contamination challenge with assay-specific truth.", "contamination", "synthetic"),
    ("TRI-06", "Compact taxonomy challenge", "Synonyms, unequal ranks, historical names, ambiguous mappings, and held-out references.", "harmonization", "synthetic"),
    ("TRI-07", "Compact failure challenge", "Missing modalities, low depth, host domination, empty data, and failed controls.", "missing", "synthetic"),
    ("TRI-10", "Even mixed-kingdom community", "Assay-appropriate recovery of the ten-member even Zymo-like community.", "even", "synthetic"),
    ("TRI-11", "Log-abundance detection challenge", "Detection behavior over six orders of taxon abundance and three depths.", "log", "synthetic"),
    ("TRI-12", "Low-biomass contamination ladder", "Paired-control filtering over six biomass levels and two batches.", "contamination", "synthetic"),
    ("TRI-13", "Structured environmental gradient", "Recovery of known habitat, diversity, taxon, and pathway gradients.", "habitat", "synthetic"),
    ("TRI-14", "Taxonomic harmonization stress test", "Synonyms, unequal ranks, close relatives, and held-out references.", "harmonization", "synthetic"),
    ("TRI-15", "Batch and paired-control challenge", "Balanced and deliberately confounded batch layouts.", "batch", "synthetic"),
    ("TRI-16", "Missing-modality and failure challenge", "Graceful integration of missing, failed, low-depth, and unreliable assays.", "missing", "synthetic"),
    ("TRI-17", "Depth, primer, and reference factorial", "Jointly test depth-limited detection, marker recovery bias, copy number, and held-out references.", "sensitivity_matrix", "synthetic"),
    ("TRI-18", "Sequencing carryover and index leakage", "Separate run-local carryover from extraction and assay-specific contamination using linked controls.", "carryover", "synthetic"),
    ("TRI-19", "Technical and extraction repeatability", "Measure within-run, between-run, and extraction-batch stability of evidence classifications.", "repeatability", "synthetic"),
    ("TRI-20", "Commercial mixed-community standards", "Executable design surrogate for even DNA, whole-cell, and log DNA standards.", "physical", "physical_design_surrogate"),
    ("TRI-21", "Custom environmental mixed-kingdom mock", "Executable design surrogate spanning profile, biomass, and kingdom diversity.", "envmock", "physical_design_surrogate"),
    ("TRI-30", "Matched built-environment microbiome study", "Executable balanced design surrogate for a 48-sample matched study.", "environmental", "environmental_design_surrogate"),
)
DATASETS = tuple(
    Dataset(dataset_id, title, purpose, specimens_for(dataset_id), mode, provenance)
    for dataset_id, title, purpose, mode, provenance in DATASET_META
)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            text = "\t".join(fields) + "\n"
            text += "".join(
                "\t".join(str(row.get(field, "")) for field in fields) + "\n"
                for row in rows
            )
            zipped.write(text.encode())


def normalize(weights: dict[str, float]) -> dict[str, float]:
    positive = {name: value for name, value in weights.items() if value > 0}
    total = sum(positive.values())
    return {name: value / total for name, value in positive.items()} if total else {}


def base_composition(dataset: Dataset, specimen: Specimen) -> dict[str, float]:
    if specimen.role == "negative_control":
        weights = {"Ralstonia_pickettii": 5, "Sphingomonas_paucimobilis": 2}
        if specimen.batch == "B02":
            weights["Burkholderia_cepacia"] = 4
        if "PCRITS" in specimen.name:
            weights = {"Penicillium_chrysogenum": 1}
        elif "PCR16" in specimen.name:
            weights = {"Ralstonia_pickettii": 1, "Sphingomonas_paucimobilis": 1}
        elif "FIELD" in specimen.name:
            weights["Cutibacterium_acnes"] = 3
        if any(status == "empty" for _, status in specimen.qc):
            return {}
        return normalize(weights)
    if specimen.role == "positive_control":
        return normalize({name: 1 for name in ZYMO})
    if dataset.mode in {"even", "physical"} and "LOG" not in specimen.environment:
        return normalize({name: (12 if TAXA[name].kingdom != "Fungi" else 2) for name in ZYMO})
    if dataset.mode == "log" or "LOG" in specimen.environment:
        values = (89.1, 8.9, .89, .89, .089, .089, .0089, .00089, .00089, .000089)
        return normalize(dict(zip(ZYMO, values)))
    if dataset.mode == "branch":
        names = ZYMO[:5] + ["Methanobrevibacter_smithii", "Aspergillus_fumigatus", "Candida_albicans", "Novelus_aeris"]
        return normalize({name: 1 for name in names})
    if dataset.mode == "contamination":
        authentic = ZYMO[:8] + ["Methanobrevibacter_smithii", "Aspergillus_fumigatus", "Candida_albicans", "Cladosporium_cladosporioides"]
        fraction = min(0.8, 800 / max(specimen.biomass, 1000))
        weights = {name: (1 - fraction) / len(authentic) for name in authentic}
        weights["Ralstonia_pickettii"] = fraction * 0.55
        weights["Sphingomonas_paucimobilis"] = fraction * 0.25
        weights["Penicillium_chrysogenum"] = fraction * 0.2
        if specimen.batch == "B02":
            weights["Burkholderia_cepacia"] = fraction * 0.3
        return normalize(weights)
    if dataset.mode == "sensitivity_matrix":
        values = (80, 10, 5, 2, 1, .5, .25, .1, .05, .025)
        names = [
            "Listeria_monocytogenes", "Pseudomonas_aeruginosa", "Bacillus_subtilis",
            "Saccharomyces_cerevisiae", "Escherichia_coli", "Salmonella_enterica",
            "Deinococcus_radiodurans", "Alternaria_alternata", "Novelus_aeris",
            "Novelomyces_aerophilus",
        ]
        return normalize(dict(zip(names, values)))
    if dataset.mode == "carryover":
        weights = {
            "Pseudomonas_aeruginosa": 5, "Bacillus_subtilis": 5,
            "Methanobrevibacter_smithii": 1, "Aspergillus_fumigatus": 2,
            "Candida_albicans": 2, "Ralstonia_pickettii": 0.08,
            "Sphingomonas_paucimobilis": 0.04, "Malassezia_restricta": 0.01,
        }
        if specimen.batch == "B02":
            weights["Burkholderia_cepacia"] = 0.12
        return normalize(weights)
    if dataset.mode == "repeatability":
        profiles = {
            "BACTERIAL": {"Pseudomonas_aeruginosa": 5, "Bacillus_subtilis": 4, "Escherichia_coli": 2},
            "FUNGAL": {"Aspergillus_fumigatus": 5, "Candida_albicans": 4, "Cladosporium_cladosporioides": 2},
            "MIXED": {"Pseudomonas_aeruginosa": 3, "Bacillus_subtilis": 3, "Aspergillus_fumigatus": 2, "Candida_albicans": 2},
            "LOW_ABUNDANCE": {"Listeria_monocytogenes": 100, "Deinococcus_radiodurans": 1, "Alternaria_alternata": .5},
        }
        return normalize(profiles[specimen.environment])
    if dataset.mode == "harmonization":
        return normalize({
            "Escherichia_coli": 2, "Shigella_flexneri": 2,
            "Bacillus_subtilis": 2, "Bacillus_cereus": 2,
            "Lactobacillus_fermentum": 2, "Cryptococcus_neoformans": 2,
            "Candida_albicans": 2, "Novelus_aeris": 1, "Novelomyces_aerophilus": 1,
        })
    if dataset.mode == "envmock":
        names = ENVIRONMENTAL + ["Halobacterium_salinarum", "Fusarium_oxysporum"]
        weights = {name: 1 for name in names}
        if specimen.environment == "FUNGAL_RICH":
            weights = {name: (8 if TAXA[name].kingdom == "Fungi" else 1) for name in names}
        elif specimen.environment == "LOG":
            weights = {name: 10 ** (-i / 2) for i, name in enumerate(names)}
        return normalize(weights)
    if dataset.mode in {"habitat", "environmental"}:
        weights = {name: 1 for name in ENVIRONMENTAL}
        env = specimen.environment.lower()
        material = specimen.material.lower()
        if "outdoor" in env:
            weights["Alternaria_alternata"] = 8
            weights["Cladosporium_cladosporioides"] = 7
        if "supply" in env:
            weights["Pseudomonas_aeruginosa"] = 4
            weights["Aspergillus_fumigatus"] = 0.4
        if "return" in env:
            weights["Cutibacterium_acnes"] = 8
            weights["Malassezia_restricta"] = 5
        if "dust" in env or "settled_dust" in material:
            weights.update({"Bacillus_subtilis": 5, "Penicillium_chrysogenum": 7, "Fusarium_oxysporum": 3})
        return normalize(weights)
    if dataset.mode == "batch":
        groups = {
            "a": ["Pseudomonas_aeruginosa", "Aspergillus_fumigatus"],
            "b": ["Bacillus_subtilis", "Candida_albicans"],
            "c": ["Escherichia_coli", "Cladosporium_cladosporioides"],
            "d": ["Staphylococcus_aureus", "Alternaria_alternata"],
        }
        names = groups.get(specimen.environment.lower(), ENVIRONMENTAL[:4])
        weights = {name: 5 for name in names}
        weights["Listeria_monocytogenes"] = 1
        weights["Ralstonia_pickettii"] = 0.1
        if specimen.batch == "B02":
            weights["Burkholderia_cepacia"] = 0.5
        return normalize(weights)
    if specimen.environment == "fungal_rich":
        return normalize({"Pseudomonas_aeruginosa": 1, "Bacillus_subtilis": 1, "Aspergillus_fumigatus": 5, "Candida_albicans": 5})
    return normalize({"Pseudomonas_aeruginosa": 3, "Bacillus_subtilis": 3, "Methanobrevibacter_smithii": 1, "Aspergillus_fumigatus": 1, "Candida_albicans": 1})


def scale_composition(base: dict[str, float], scale: str) -> dict[str, float]:
    weights: dict[str, float] = {}
    for name, fraction in base.items():
        taxon = TAXA[name]
        if scale == "genomic_dna":
            factor = 1
        elif scale == "genome_copy":
            factor = 1 / taxon.genome_size_mb
        elif scale == "16s":
            if taxon.kingdom == "Fungi":
                continue
            factor = taxon.marker_16s_copies * taxon.recovery_16s
        elif scale == "its":
            if taxon.kingdom != "Fungi":
                continue
            factor = taxon.marker_its_copies * taxon.recovery_its
        else:
            raise ValueError(scale)
        weights[name] = fraction * factor
    return normalize(weights)


def assay_profile(base: dict[str, float], assay: str) -> dict[str, float]:
    return scale_composition(base, {"WMS": "genome_copy", "16S": "16s", "ITS": "its"}[assay])


def sequence_for(taxon: str, assay: str, mate: int, length: int = 100) -> str:
    seed = int(hashlib.sha256(f"{taxon}:{assay}:{mate}".encode()).hexdigest()[:16], 16)
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def allocations(profile: dict[str, float], total: int) -> list[tuple[str, int]]:
    if not profile or total == 0:
        return []
    result = [(name, max(1, round(fraction * total))) for name, fraction in profile.items()]
    return result


def write_fastq(
    path: Path, library_id: str, assay: str, profile: dict[str, float],
    mate: int, total: int,
) -> list[dict[str, object]]:
    path.parent.mkdir(parents=True, exist_ok=True)
    origins: list[dict[str, object]] = []
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            read_index = 0
            for taxon, count in allocations(profile, total):
                sequence = sequence_for(taxon, assay, mate)
                status = "contaminant" if taxon in CONTAMINANTS else "authentic"
                for _ in range(count):
                    read_index += 1
                    read_id = f"{library_id}:{read_index:05d}"
                    zipped.write(
                        f"@{read_id}/{mate} taxon={taxon}\n{sequence}\n+\n{'I' * len(sequence)}\n".encode()
                    )
                    if mate == 1:
                        origins.append({
                            "library_id": library_id, "read_id": read_id,
                            "assay": assay, "taxon": taxon, "truth_status": status,
                        })
    return origins


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def control_stage(name: str) -> str:
    if "FIELD" in name:
        return "field"
    if "PCR" in name:
        return "pcr"
    if "LIB" in name:
        return "library"
    return "extraction"


def expected_qc(specimen: Specimen, assay: str) -> str:
    return dict(specimen.qc).get(assay, "pass")


def generate(dataset: Dataset) -> None:
    root = BENCHMARKS / dataset.dataset_id
    metadata_dir, reads_dir, truth_dir = root / "metadata", root / "reads", root / "truth"
    reference_dir = root / "references"
    if reads_dir.is_dir():
        for generated_fastq in reads_dir.glob("*.fastq.gz"):
            generated_fastq.unlink()
    samples: list[dict[str, object]] = []
    controls: list[dict[str, object]] = []
    expected_taxa: list[dict[str, object]] = []
    origins: list[dict[str, object]] = []
    sample_qc: list[dict[str, object]] = []
    modality_status: list[dict[str, object]] = []
    taxon_status: list[dict[str, object]] = []
    pathways: list[dict[str, object]] = []
    scale_rows: dict[str, list[dict[str, object]]] = {
        scale: [] for scale in ("genomic_dna", "genome_copy", "16s", "its")
    }

    controls_by_batch_assay: dict[tuple[str, str], list[Specimen]] = {}
    for specimen in dataset.specimens:
        if specimen.role == "negative_control":
            for assay in specimen.assays:
                controls_by_batch_assay.setdefault((specimen.batch, assay), []).append(specimen)

    for specimen in dataset.specimens:
        base = base_composition(dataset, specimen)
        for scale in scale_rows:
            for taxon, fraction in scale_composition(base, scale).items():
                scale_rows[scale].append({
                    "biological_sample_id": specimen.name,
                    "taxon": taxon,
                    "expected_fraction": f"{fraction:.10f}",
                    "truth_scale": scale,
                })
        for planned_assay in ASSAYS:
            state = "present" if planned_assay in specimen.assays else "missing_by_design"
            modality_status.append({
                "biological_sample_id": specimen.name,
                "assay": planned_assay,
                "modality_status": state if state != "present" else expected_qc(specimen, planned_assay),
                "integration_interpretation": (
                    "not_evidence_of_absence" if state == "missing_by_design"
                    else "uninterpretable" if expected_qc(specimen, planned_assay) in {"failed", "empty"}
                    else "sensitivity_limited" if expected_qc(specimen, planned_assay) in {"low_depth", "host_dominated"}
                    else "usable"
                ),
            })
        for assay in specimen.assays:
            library_id = f"{specimen.name}_{assay}"
            qc = expected_qc(specimen, assay)
            single = qc == "forward_only"
            r1 = reads_dir / f"{library_id}_R1.fastq.gz"
            r2 = reads_dir / f"{library_id}_R2.fastq.gz"
            profile = assay_profile(base, assay)
            if not profile and qc == "pass":
                qc = "no_target_template"
            total = 0 if qc in {"failed", "empty"} else 4 if qc == "low_depth" else 24
            if qc == "no_target_template":
                total = 0
            r2_total = max(0, total - 3) if qc == "mismatched_pairs" else total
            origins.extend(write_fastq(r1, library_id, assay.lower(), profile, 1, total))
            if not single:
                write_fastq(r2, library_id, assay.lower(), profile, 2, r2_total)
            samples.append({
                "library_id": library_id,
                "biological_sample_id": specimen.name,
                "specimen_id": specimen.name,
                "subject_id": specimen.name if specimen.role == "biological" else "CONTROL",
                "assay": assay,
                "sample_role": specimen.role,
                "sample_material": specimen.material,
                "environment_category": specimen.environment,
                "library_layout": "single" if single else "paired",
                "fastq_r1": f"../reads/{r1.name}",
                "fastq_r2": "" if single else f"../reads/{r2.name}",
                "extraction_batch": specimen.batch,
                "pcr_plate": "" if assay == "WMS" else f"P_{assay}_{specimen.batch}",
                "library_batch": f"LB_{specimen.batch}",
                "sequencing_run": f"SR_{1 + int(specimen.batch[-1]) % 2}",
                "technical_replicate": 1,
                "biomass_measure": specimen.biomass,
                "biomass_unit": "genome_equivalents_relative",
            })
            sample_qc.append({
                "library_id": library_id,
                "expected_qc_status": qc,
                "expected_read_pairs": total if not single else 0,
                "expected_host_fraction": "0.90" if qc == "host_dominated" else "0.20" if "host_20pct" in qc else "0.00",
                "expected_action": "exclude" if qc in {"failed", "empty", "mismatched_pairs"} else "qualify" if qc in {"low_depth", "host_dominated", "forward_only", "no_target_template"} else "retain",
            })
            for taxon, fraction in profile.items():
                info = TAXA[taxon]
                status = "probable_contaminant" if taxon in CONTAMINANTS else (
                    "reference_held_out" if info.reference_status != "reference_seen" else "concordant"
                )
                expected_taxa.append({
                    "library_id": library_id, "biological_sample_id": specimen.name,
                    "assay": assay, "taxon": taxon, "kingdom": info.kingdom,
                    "genus": info.genus, "expected_rank": "genus" if taxon in {"Escherichia_coli", "Shigella_flexneri"} and assay == "16S" else "species",
                    "expected_fraction": f"{fraction:.10f}",
                    "truth_status": "contaminant" if taxon in CONTAMINANTS else "authentic",
                })
                taxon_status.append({
                    "library_id": library_id, "taxon": taxon,
                    "expected_cross_assay_status": status,
                    "supported_comparison_rank": "genus" if assay in {"16S", "ITS"} and dataset.mode == "harmonization" else "species",
                    "reason": info.reference_status if status == "reference_held_out" else status,
                })
            if assay == "WMS" and specimen.role == "biological":
                pathway = "PWY_RETURN_OCCUPANCY" if "return" in specimen.environment else "PWY_PLANT_ASSOCIATED" if "outdoor" in specimen.environment else "PWY_CORE_METABOLISM"
                pathways.append({
                    "library_id": library_id, "pathway_id": pathway,
                    "expected_direction": "enriched" if pathway != "PWY_CORE_METABOLISM" else "present",
                    "expected_fraction": "0.10",
                })

            if specimen.role == "biological":
                candidates = controls_by_batch_assay.get((specimen.batch, assay), [])
                used_stages: set[str] = set()
                for control in candidates:
                    stage = control_stage(control.name)
                    if stage in used_stages or (stage == "pcr" and assay == "WMS"):
                        continue
                    used_stages.add(stage)
                    controls.append({
                        "library_id": library_id,
                        "control_library_id": f"{control.name}_{assay}",
                        "control_stage": stage,
                        "relationship": "paired",
                        "priority": len(used_stages),
                        "reason": f"batch-local {stage} control",
                    })

    write_tsv(metadata_dir / "samples.tsv", SAMPLE_FIELDS, samples)
    write_tsv(metadata_dir / "control_relationships.tsv", CONTROL_FIELDS, controls)
    write_tsv(truth_dir / "control_relationships.tsv", CONTROL_FIELDS, controls)
    write_tsv(truth_dir / "expected_taxa.tsv", TRUTH_FIELDS, expected_taxa)
    write_tsv(
        truth_dir / "precomputed_profiles.tsv",
        ["library_id", "biological_sample_id", "assay", "feature_id",
         "canonical_taxon", "count", "relative_abundance", "truth_status"],
        [{
            "library_id": row["library_id"],
            "biological_sample_id": row["biological_sample_id"],
            "assay": row["assay"],
            "feature_id": f"SYN:{row['taxon']}",
            "canonical_taxon": row["taxon"],
            "count": max(1, round(float(row["expected_fraction"]) * 10000)),
            "relative_abundance": row["expected_fraction"],
            "truth_status": row["truth_status"],
        } for row in expected_taxa],
    )
    write_tsv(
        truth_dir / "read_origins.tsv",
        ["library_id", "read_id", "assay", "taxon", "truth_status"],
        origins,
    )
    write_tsv_gz(
        truth_dir / "read_origin_wms.tsv.gz",
        ["library_id", "read_id", "assay", "taxon", "truth_status"],
        [row for row in origins if row["assay"] == "WMS"],
    )
    for scale, rows in scale_rows.items():
        filename = {
            "genomic_dna": "sample_composition_genomic_dna.tsv",
            "genome_copy": "sample_composition_genome_copy.tsv",
            "16s": "sample_composition_16s.tsv",
            "its": "sample_composition_its.tsv",
        }[scale]
        write_tsv(
            truth_dir / filename,
            ["biological_sample_id", "taxon", "expected_fraction", "truth_scale"],
            rows,
        )
    members = []
    seen_taxa = sorted({row["taxon"] for row in expected_taxa})
    for name in seen_taxa:
        taxon = TAXA[name]
        members.append({
            "taxon": name, "kingdom": taxon.kingdom, "genus": taxon.genus,
            "source_genome_accession": taxon.accession,
            "source_genome_sha256": hashlib.sha256(taxon.accession.encode()).hexdigest(),
            "genome_size_mb": taxon.genome_size_mb,
            "marker_16s_copies": taxon.marker_16s_copies,
            "marker_its_copies": taxon.marker_its_copies,
            "reference_status": taxon.reference_status,
        })
    write_tsv(
        truth_dir / "community_members.tsv",
        ["taxon", "kingdom", "genus", "source_genome_accession",
         "source_genome_sha256", "genome_size_mb", "marker_16s_copies",
         "marker_its_copies", "reference_status"],
        members,
    )
    write_tsv(
        truth_dir / "taxonomic_crosswalk.tsv",
        ["original_taxon_name", "original_taxon_id", "source_database", "source_rank",
         "canonical_taxon_id", "canonical_name", "harmonized_rank", "mapping_method",
         "mapping_confidence", "cross_assay_status"],
        [{
            "original_taxon_name": taxon.synonym or taxon.name,
            "original_taxon_id": f"SRC:{index + 1}",
            "source_database": "frozen_synthetic_2026_01",
            "source_rank": "species",
            "canonical_taxon_id": f"TRI:{index + 1}",
            "canonical_name": taxon.name,
            "harmonized_rank": "genus" if taxon.name in {"Escherichia_coli", "Shigella_flexneri"} else "species",
            "mapping_method": "synonym" if taxon.synonym else "exact" if taxon.reference_status == "reference_seen" else "held_out",
            "mapping_confidence": "1.0" if taxon.reference_status == "reference_seen" else "0.5",
            "cross_assay_status": "ambiguous" if taxon.reference_status != "reference_seen" else "supported",
        } for index, taxon in enumerate(TAXA_LIST)],
    )
    write_tsv(
        truth_dir / "contaminants.tsv",
        ["taxon", "contamination_pattern", "stage", "expected_batches", "expected_modalities"],
        [{
            "taxon": name, "contamination_pattern": pattern, "stage": stage,
            "expected_batches": "B02" if pattern == "batch_B02" else "all",
            "expected_modalities": "ITS" if TAXA[name].kingdom == "Fungi" else "WMS,16S",
        } for name, (pattern, stage) in CONTAMINANTS.items()],
    )
    write_tsv(
        truth_dir / "expected_sample_qc.tsv",
        ["library_id", "expected_qc_status", "expected_read_pairs",
         "expected_host_fraction", "expected_action"],
        sample_qc,
    )
    write_tsv(
        truth_dir / "expected_modality_status.tsv",
        ["biological_sample_id", "assay", "modality_status", "integration_interpretation"],
        modality_status,
    )
    write_tsv(
        truth_dir / "expected_taxon_status.tsv",
        ["library_id", "taxon", "expected_cross_assay_status",
         "supported_comparison_rank", "reason"],
        taxon_status,
    )
    write_tsv(
        truth_dir / "expected_pathways.tsv",
        ["library_id", "pathway_id", "expected_direction", "expected_fraction"],
        pathways,
    )
    (truth_dir / "true_16s_asvs.fasta").write_text(
        "".join(f">{name}|copy=1\n{sequence_for(name, '16s_asv', 1, 250)}\n" for name in seen_taxa if TAXA[name].kingdom != "Fungi"),
        encoding="ascii",
    )
    (truth_dir / "true_its_asvs.fasta").write_text(
        "".join(f">{name}|variant=1\n{sequence_for(name, 'its_asv', 1, 300)}\n" for name in seen_taxa if TAXA[name].kingdom == "Fungi"),
        encoding="ascii",
    )
    reference_dir.mkdir(parents=True, exist_ok=True)
    for kingdom, filename in (
        ("Bacteria", "artificial_bacterial.fasta"),
        ("Archaea", "artificial_archaeal.fasta"),
        ("Fungi", "artificial_fungal.fasta"),
    ):
        (reference_dir / filename).write_text(
            "".join(
                f">{name}|{TAXA[name].accession}\n{sequence_for(name, 'reference', 1, 1000)}\n"
                for name in seen_taxa if TAXA[name].kingdom == kingdom
            ),
            encoding="ascii",
        )
    (reference_dir / "artificial_host.fasta").write_text(
        f">TRIARCHY_ARTIFICIAL_HOST\n{sequence_for('artificial_host', 'host', 1, 1000)}\n",
        encoding="ascii",
    )
    assertions = {
        "dataset": dataset.dataset_id,
        "schema_version": "1.0",
        "assertions": [
            {"id": "library_count", "operator": "equals", "expected": len(samples)},
            {"id": "metadata_valid", "operator": "equals", "expected": True},
            {"id": "required_assays_present", "operator": "set_equals", "expected": list(ASSAYS)},
            {"id": "unexpected_workflow_failure", "operator": "equals", "expected": 0},
            {"id": "taxon_status_truth", "operator": "table_match_with_tolerance",
             "truth_file": "expected_taxon_status.tsv", "relative_abundance_tolerance": 0.001},
        ],
    }
    (truth_dir / "expected_assertions.json").write_text(
        json.dumps(assertions, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    checksums = {
        str(path.relative_to(root)): sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.name not in {"generation_manifest.json", "README.md"}
    }
    manifest = {
        "dataset_id": dataset.dataset_id,
        "schema_version": "2.0",
        "generator": "scripts/generate_benchmark_datasets.py",
        "generator_version": "2.0.0",
        "provenance_class": dataset.provenance_class,
        "random_seed": "SHA-256-derived per taxon/assay/mate",
        "simulator": {"name": "TRIARCHY deterministic CI simulator", "version": "2.0.0"},
        "primer_sequences": {"16S_forward": "CCTACGGGNGGCWGCAG", "16S_reverse": "GACTACHVGGGTATCTAATCC", "ITS_forward": "GTGAATCATCGAATCTTTGAA", "ITS_reverse": "TCCTCCGCTTATTGATATGC"},
        "read_length": 100,
        "insert_size_model": "fixed deterministic paired surrogate",
        "sequencing_error_model": "perfect CI reads with scenario flags in expected_sample_qc.tsv",
        "copy_number_assumptions": "community_members.tsv",
        "reference_database_versions": {"Kraken2": "standard_frozen_2026_01", "SILVA": "138.2", "EUKARYOME": "2.0-99pct", "HUMAnN": "3.9"},
        "library_count": len(samples),
        "purpose": dataset.purpose,
        "limitations": [
            "Repository-sized reads exercise software contracts; publication-scale reads must be generated from frozen genomes.",
            "Physical and environmental design surrogates do not claim laboratory or field observations.",
        ],
        "checksums": checksums,
    }
    (truth_dir / "generation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    qualifier = (
        "This is an executable design surrogate, not empirical laboratory or field data. "
        if dataset.provenance_class != "synthetic" else ""
    )
    (root / "README.md").write_text(
        f"# {dataset.dataset_id} — {dataset.title}\n\n{dataset.purpose}\n\n"
        f"{qualifier}The deterministic fixture contains {len(samples)} libraries, "
        f"{len(controls)} batch-local control links, paired or explicitly single-end FASTQs, "
        "assay-specific abundance scales, read origins, ASV truth, artificial references, "
        "precomputed profiles, machine-readable assertions, taxonomy mappings, contamination "
        "truth, QC expectations, pathways, checksums, and provenance.\n\n"
        "Regenerate from the repository root:\n\n```bash\n"
        f"python3 scripts/generate_benchmark_datasets.py --dataset {dataset.dataset_id}\n"
        "```\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset", action="append",
        choices=[item.dataset_id for item in DATASETS],
        help="Dataset to generate; repeat as needed. Default: all.",
    )
    args = parser.parse_args()
    wanted = set(args.dataset or [item.dataset_id for item in DATASETS])
    for dataset in DATASETS:
        if dataset.dataset_id in wanted:
            generate(dataset)
            print(f"generated {dataset.dataset_id}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

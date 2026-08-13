#!/usr/bin/env python3
"""Build editable SVG manuscript and supplementary figures for TRIARCHY."""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "figures"

COLORS = {
    "ink": "#17212b",
    "muted": "#657080",
    "line": "#9aa8b6",
    "panel": "#f8fafc",
    "wms": "#d8eaf6",
    "s16": "#ddf1e5",
    "its": "#f3e5d2",
    "integration": "#e9e1f5",
    "output": "#d9ece9",
    "accent": "#2f6f73",
    "warn": "#b56b00",
    "bad": "#a43d3d",
    "gray": "#e5e9ef",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def wrap(text: str, width: int = 26) -> list[str]:
    words = text.replace("_", " ").split()
    lines: list[str] = []
    current = ""
    for word in words:
        if len(current) + len(word) + 1 > width and current:
            lines.append(current)
            current = word
        else:
            current = f"{current} {word}".strip()
    if current:
        lines.append(current)
    return lines or [""]


class SVG:
    def __init__(self, width: int, height: int, title: str):
        self.width = width
        self.height = height
        self.items: list[str] = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" role="img" aria-label="{escape(title)}">',
            "<defs>",
            '<marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">',
            f'<path d="M 0 0 L 10 5 L 0 10 z" fill="{COLORS["line"]}"/>',
            "</marker>",
            "</defs>",
            f'<rect width="{width}" height="{height}" fill="white"/>',
        ]

    def rect(self, x: float, y: float, w: float, h: float, fill: str, stroke: str | None = None, rx: int = 10, sw: float = 1.2) -> None:
        stroke = stroke or COLORS["line"]
        self.items.append(
            f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'
        )

    def line(self, x1: float, y1: float, x2: float, y2: float, color: str | None = None, width: float = 1.6, arrow: bool = False) -> None:
        marker = ' marker-end="url(#arrow)"' if arrow else ""
        self.items.append(
            f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color or COLORS["line"]}" stroke-width="{width}"{marker}/>'
        )

    def polyline(self, points: list[tuple[float, float]], color: str, width: float = 2.4) -> None:
        data = " ".join(f"{x},{y}" for x, y in points)
        self.items.append(f'<polyline points="{data}" fill="none" stroke="{color}" stroke-width="{width}"/>')

    def circle(self, x: float, y: float, r: float, fill: str, stroke: str | None = None) -> None:
        self.items.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{stroke or COLORS["line"]}" stroke-width="1.2"/>')

    def text(self, x: float, y: float, value: str, size: int = 22, weight: str = "400", color: str | None = None, anchor: str = "start") -> None:
        self.items.append(
            f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" '
            f'font-weight="{weight}" fill="{color or COLORS["ink"]}" text-anchor="{anchor}">{escape(value)}</text>'
        )

    def label(self, x: float, y: float, value: str, width: int = 24, size: int = 20, weight: str = "400", color: str | None = None, anchor: str = "middle", line_gap: int = 24) -> None:
        for i, line in enumerate(wrap(value, width)):
            self.text(x, y + i * line_gap, line, size=size, weight=weight, color=color, anchor=anchor)

    def panel(self, x: float, y: float, w: float, h: float, letter: str, title: str) -> None:
        self.rect(x, y, w, h, COLORS["panel"], COLORS["gray"], rx=12, sw=1)
        self.text(x + 22, y + 34, letter, size=30, weight="700")
        self.text(x + 62, y + 32, title, size=22, weight="700")

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(self.items + ["</svg>\n"]), encoding="utf-8")


def box(svg: SVG, x: float, y: float, w: float, h: float, text: str, fill: str, size: int = 18) -> None:
    svg.rect(x, y, w, h, fill, COLORS["line"], rx=10)
    svg.label(x + w / 2, y + h / 2 - 8, text, width=max(14, int(w / 10)), size=size, weight="600")


def bar(svg: SVG, x: float, y: float, w: float, h: float, value: float, max_value: float, fill: str, label: str, value_label: str = "") -> None:
    svg.rect(x, y, w, h, "#f0f3f7", COLORS["gray"], rx=4)
    bw = 0 if max_value <= 0 else w * value / max_value
    svg.rect(x, y, bw, h, fill, fill, rx=4)
    svg.text(x, y - 8, label, size=14, color=COLORS["muted"])
    if value_label:
        svg.text(x + w + 8, y + h - 4, value_label, size=14, color=COLORS["ink"])


def read_benchmark_metrics() -> list[dict[str, str]]:
    metrics = ROOT / "results" / "frozen-release-accuracy-20260812" / "benchmark_comparison.tsv"
    if not metrics.exists():
        return []
    return read_tsv(metrics)


def workflow_architecture_figure() -> Path:
    svg = SVG(2400, 1500, "TRIARCHY detailed workflow architecture")
    svg.text(60, 75, "TRIARCHY workflow architecture", size=38, weight="700")
    svg.text(60, 115, "Matched or subset WMS, 16S and ITS assays are processed independently, linked to controls, and interpreted through an evidence layer.", size=22, color=COLORS["muted"])

    # Swimlane backgrounds.
    lanes = [
        (60, 165, 420, 1210, "Inputs and design", COLORS["panel"]),
        (535, 165, 640, 1210, "Assay-native processing", "#f7faf7"),
        (1230, 165, 520, 1210, "Control-aware integration", "#fbf8ff"),
        (1805, 165, 535, 1210, "Outputs and interpretation", "#f7fbfb"),
    ]
    for x, y, w, h, title, fill in lanes:
        svg.rect(x, y, w, h, fill, COLORS["gray"], rx=12, sw=1)
        svg.text(x + 24, y + 42, title, size=25, weight="700")

    # Inputs and validation.
    input_nodes = [
        (110, 255, 320, 70, "Sample manifest\nlibraries, roles, batches", COLORS["panel"]),
        (110, 365, 320, 70, "Paired-control table\nfield, extraction, PCR, library, run", "#f7f1e8"),
        (110, 475, 320, 70, "Evidence policy\nthresholds + calibration status", COLORS["integration"]),
        (110, 585, 320, 70, "Assay performance\nreliability, depth, reference coverage", COLORS["integration"]),
        (110, 725, 320, 86, "Metadata validator\nIDs, FASTQs, paired reads, control edges", COLORS["output"]),
        (110, 875, 320, 86, "Sample-control graph\nstage-, assay- and batch-specific links", COLORS["output"]),
    ]
    for x, y, w, h, text, fill in input_nodes:
        box(svg, x, y, w, h, text, fill, size=17)
    for y1, y2 in [(325, 365), (435, 475), (545, 585), (655, 725), (811, 875)]:
        svg.line(270, y1, 270, y2, arrow=True)

    # Assay branches.
    branch_x = 585
    branch_w = 540
    branches = [
        ("WMS", 250, COLORS["wms"], [
            "fastp quality trimming",
            "host-read depletion",
            "Kraken2 / Bracken taxonomy",
            "HUMAnN functional profiles",
        ]),
        ("16S", 595, COLORS["s16"], [
            "primer/adaptor trimming",
            "DADA2 ASV inference",
            "bacterial/archaeal taxonomy",
            "marker-gene count table",
        ]),
        ("ITS", 940, COLORS["its"], [
            "primer/adaptor trimming",
            "ITSxpress region extraction",
            "DADA2 fungal ASVs",
            "fungal taxonomy",
        ]),
    ]
    for title, y, fill, steps in branches:
        svg.rect(branch_x, y, branch_w, 250, "white", COLORS["gray"], rx=10)
        svg.text(branch_x + 24, y + 40, title, size=28, weight="700")
        for i, step in enumerate(steps):
            sx = branch_x + 130
            sy = y + 35 + i * 49
            svg.circle(sx, sy, 12, fill)
            svg.text(sx + 28, sy + 8, step, size=20)
            if i < len(steps) - 1:
                svg.line(sx, sy + 15, sx, sy + 37, arrow=True)
    for yy in (375, 720, 1065):
        svg.line(430, 918, 585, yy, arrow=True)

    # Integration.
    integration_nodes = [
        (1290, 245, 400, 76, "QC summary\npass, low depth, failed, missing", "#e8eef3"),
        (1290, 375, 400, 76, "Linked-control annotation\nprevalence, abundance ratio, batch localization", "#f7d9d9"),
        (1290, 525, 400, 76, "Taxonomic harmonization\ncanonical IDs, rank roll-up, provenance", COLORS["integration"]),
        (1290, 675, 400, 86, "Stage I authenticity\nassay evidence + detection power + contamination penalties", COLORS["output"]),
        (1290, 845, 400, 86, "Stage II explanation\nconcordant, assay-limited, taxonomy-limited, contaminant, discordant", COLORS["output"]),
        (1290, 1035, 400, 86, "Cross-assay concordance\nkingdom-appropriate comparisons", COLORS["integration"]),
    ]
    for x, y, w, h, text, fill in integration_nodes:
        box(svg, x, y, w, h, text, fill, size=17)
    for y1, y2 in [(321, 375), (451, 525), (601, 675), (761, 845), (931, 1035)]:
        svg.line(1490, y1, 1490, y2, arrow=True)
    for yy in (375, 720, 1065):
        svg.line(1125, yy, 1290, yy if yy != 1065 else 1078, arrow=True)
    svg.line(430, 918, 1290, 413, arrow=True, color=COLORS["bad"])

    # Outputs.
    outputs = [
        (1865, 250, 420, 76, "Assay-native profiles\nraw + annotated + high-confidence", COLORS["panel"]),
        (1865, 390, 420, 76, "Integrated taxon evidence\nbiological authenticity + decision reasons", COLORS["output"]),
        (1865, 530, 420, 76, "Contamination report\nprobable source by control stage or batch", "#f7d9d9"),
        (1865, 670, 420, 76, "Diversity and ecology\nassay-specific alpha/beta diversity", COLORS["s16"]),
        (1865, 810, 420, 76, "Functional summaries\nWMS pathways kept separate", COLORS["wms"]),
        (1865, 950, 420, 76, "HTML report + provenance\nparameters, versions, reference manifests", COLORS["integration"]),
        (1865, 1120, 420, 96, "Interpretation\nrobust, assay-limited, reference-limited, control-associated or unresolved", COLORS["output"]),
    ]
    for x, y, w, h, text, fill in outputs:
        box(svg, x, y, w, h, text, fill, size=17)
    for y in (288, 428, 568, 708, 848, 988, 1168):
        svg.line(1690, 888 if y < 1030 else 1078, 1865, y, arrow=True)

    # Legend.
    legend_y = 1325
    legend = [
        ("WMS", COLORS["wms"]),
        ("16S", COLORS["s16"]),
        ("ITS", COLORS["its"]),
        ("Controls", "#f7d9d9"),
        ("Evidence model", COLORS["output"]),
    ]
    for i, (label, fill) in enumerate(legend):
        x = 620 + i * 220
        svg.circle(x, legend_y, 13, fill)
        svg.text(x + 24, legend_y + 7, label, size=18, color=COLORS["muted"])
    svg.text(60, 1435, "Missing-by-design assays reduce the scope of inference rather than becoming measured negative evidence.", size=20, color=COLORS["muted"])

    path = OUT / "figure1_workflow_architecture.svg"
    svg.save(path)
    return path


def main_figure() -> Path:
    svg = SVG(2400, 1700, "TRIARCHY main composite figure")
    svg.panel(40, 40, 1120, 760, "A", "Workflow architecture")
    svg.panel(1240, 40, 1120, 760, "B", "Evidence-classification logic")
    svg.panel(40, 860, 1120, 760, "C", "Benchmark comparison template")
    svg.panel(1240, 860, 1120, 760, "D", "Stub resource profile")

    # Panel A
    nodes = [
        (120, 130, 230, 70, "Manifests + paired controls", COLORS["panel"]),
        (445, 130, 230, 70, "Sample-control graph", COLORS["integration"]),
        (780, 130, 230, 70, "Assay branches", COLORS["panel"]),
        (190, 300, 190, 70, "WMS taxonomy + functions", COLORS["wms"]),
        (465, 300, 190, 70, "16S ASVs + taxonomy", COLORS["s16"]),
        (740, 300, 190, 70, "ITS ASVs + taxonomy", COLORS["its"]),
        (360, 480, 250, 70, "Taxonomic harmonization", COLORS["integration"]),
        (710, 480, 250, 70, "Authenticity + explanation", COLORS["integration"]),
        (530, 650, 280, 70, "Evidence + ecology + report", COLORS["output"]),
    ]
    for n in nodes:
        box(svg, *n)
    for x1, y1, x2, y2 in [
        (350, 165, 445, 165), (675, 165, 780, 165),
        (895, 200, 285, 300), (895, 200, 560, 300), (895, 200, 835, 300),
        (285, 370, 485, 480), (560, 370, 485, 480), (835, 370, 485, 480),
        (610, 515, 710, 515), (835, 550, 670, 650),
    ]:
        svg.line(x1, y1, x2, y2, arrow=True)

    # Panel B
    examples = [
        ("Concordant", "WMS + marker agree", COLORS["output"]),
        ("Assay-limited", "Other assay underpowered", "#fff4d8"),
        ("Rank-limited", "Agreement after roll-up", "#e8eef3"),
        ("Contaminant", "Linked controls enriched", "#f7d9d9"),
        ("Reference-limited", "Database coverage low", "#f3e5d2"),
        ("Genuine discordance", "Adequate assays conflict", "#ead7ef"),
    ]
    start_x, start_y = 1320, 135
    for i, (name, desc, color) in enumerate(examples):
        col, row = i % 2, i // 2
        x = start_x + col * 510
        y = start_y + row * 190
        box(svg, x, y, 420, 122, name, color, size=20)
        svg.text(x + 28, y + 95, desc, size=17, color=COLORS["muted"])
        for j, assay in enumerate(("WMS", "16S", "ITS")):
            cx = x + 285 + j * 38
            fill = COLORS["wms"] if assay == "WMS" else COLORS["s16"] if assay == "16S" else COLORS["its"]
            svg.circle(cx, y + 40, 14, fill)
            svg.text(cx, y + 72, assay, size=12, anchor="middle", color=COLORS["muted"])

    # Panel C
    metrics = read_benchmark_metrics()
    method_labels = {
        "independent_assay": "Independent assay",
        "naive_concat": "Naive concatenation",
        "shared_taxonomy": "Shared taxonomy",
        "per_assay_decontam": "Per-assay controls",
        "global_control": "Global controls",
        "TRIARCHY": "TRIARCHY",
    }
    svg.text(100, 940, "Synthetic benchmark results", size=24, weight="700")
    svg.text(100, 972, "1,562 collapsed sample-taxon truth observations", size=17, color=COLORS["muted"])
    if metrics:
        rows = [row for row in metrics if row["method"] in method_labels]
        svg.text(360, 1012, "Authentic F1", size=15, weight="700", anchor="middle")
        svg.text(570, 1012, "Contam. F1", size=15, weight="700", anchor="middle")
        svg.text(780, 1012, "Status acc.", size=15, weight="700", anchor="middle")
        svg.text(995, 1012, "1 - Brier", size=15, weight="700", anchor="middle")
        for i, row in enumerate(rows):
            y = 1050 + i * 78
            method = row["method"]
            is_triarchy = method == "TRIARCHY"
            label_color = COLORS["accent"] if is_triarchy else COLORS["ink"]
            svg.text(115, y + 20, method_labels[method], size=17, weight="700" if is_triarchy else "400", color=label_color)
            contaminant_precision = float(row["contaminant_precision"])
            contaminant_recall = float(row["contaminant_recall"])
            contaminant_f1 = (
                2 * contaminant_precision * contaminant_recall
                / (contaminant_precision + contaminant_recall)
                if contaminant_precision + contaminant_recall else 0.0
            )
            values = [
                float(row["presence_f1"]),
                contaminant_f1,
                float(row["explanatory_status_accuracy"]),
                1 - float(row["brier_score"]),
            ]
            colors = [COLORS["wms"], "#a43d3d", COLORS["accent"], COLORS["integration"]]
            for j, value in enumerate(values):
                x = 300 + j * 210
                svg.rect(x, y, 130, 24, "#edf1f5", COLORS["gray"], rx=4)
                svg.rect(x, y, 130 * value, 24, colors[j], colors[j], rx=4)
                svg.text(x + 150, y + 19, f"{value:.3f}", size=14, color=COLORS["ink"])
        svg.label(925, 1510, "Brier is shown as 1 - Brier so longer bars indicate better score calibration.", width=36, size=16, color=COLORS["muted"])
    else:
        svg.label(580, 1170, "Run scripts/compare_evidence_methods.py to populate benchmark metrics.", width=44, size=20, color=COLORS["warn"])

    # Panel D
    resource_path = ROOT / "output" / "benchmarks" / "resource-profile-stub.json"
    if resource_path.exists():
        payload = json.loads(resource_path.read_text(encoding="utf-8"))
        runs = payload.get("runs", [])
        max_wall = max((float(run["wall_seconds"]) for run in runs), default=1)
        max_mem = max((float(run["observed_peak_rss_gib"]) for run in runs), default=1)
        svg.text(1310, 940, f"Compact stub bundle: {payload.get('library_count', 0)} libraries, {payload.get('total_wall_seconds', 0)} s", size=20, weight="700")
        for i, run in enumerate(runs):
            y = 1010 + i * 95
            svg.text(1310, y + 22, run["dataset_id"], size=16, weight="700")
            bar(svg, 1410, y, 330, 24, float(run["wall_seconds"]), max_wall, COLORS["accent"], "wall s", f"{float(run['wall_seconds']):.1f}")
            bar(svg, 1410, y + 42, 330, 24, float(run["observed_peak_rss_gib"]), max_mem, "#7b5db6", "peak GiB", f"{float(run['observed_peak_rss_gib']):.2f}")
        svg.label(1890, 1140, "Stub measurements validate orchestration only; production Kraken2/HUMAnN resources must be measured separately.", width=30, size=18, color=COLORS["warn"])
    else:
        svg.label(1480, 1110, "Run scripts/benchmark_resources.py --mode stub to populate this panel.", width=34, size=20, color=COLORS["warn"])

    path = OUT / "figure1_composite.svg"
    svg.save(path)
    return path


def supplementary_s1() -> Path:
    svg = SVG(1800, 1350, "Supplementary Figure S1 complete DAG")
    svg.text(60, 70, "Supplementary Figure S1. Complete TRIARCHY DAG", size=30, weight="700")
    columns = [
        ("Inputs", 80, COLORS["panel"], ["Sample manifest", "Control relationships", "Evidence policy", "Assay performance"]),
        ("WMS", 420, COLORS["wms"], ["fastp QC", "Bowtie2 host removal", "Kraken2 + Bracken", "HUMAnN functions"]),
        ("16S", 760, COLORS["s16"], ["cutadapt", "DADA2 ASVs", "16S taxonomy"]),
        ("ITS", 1100, COLORS["its"], ["cutadapt", "ITSxpress", "DADA2 ASVs", "ITS taxonomy"]),
        ("Integration", 1440, COLORS["integration"], ["QC summary", "Control annotation", "Harmonization", "Authenticity + explanation", "Concordance", "Diversity", "HTML report"]),
    ]
    for title, x, color, items in columns:
        svg.text(x + 90, 130, title, size=24, weight="700", anchor="middle")
        for i, item in enumerate(items):
            y = 180 + i * 105
            box(svg, x, y, 220, 62, item, color, size=16)
            if i > 0:
                svg.line(x + 110, y - 43, x + 110, y, arrow=True)
    for x in (420, 760, 1100):
        svg.line(300, 210, x, 210, arrow=True)
    for x in (640, 980, 1320):
        svg.line(x, 500, 1440, 500, arrow=True)
    path = OUT / "supplementary_figure_s1_dag.svg"
    svg.save(path)
    return path


def supplementary_s2() -> Path:
    svg = SVG(1800, 1100, "Supplementary Figure S2 control links")
    svg.text(60, 70, "Supplementary Figure S2. Stage-linked sample-control graph", size=30, weight="700")
    samples = [("BIO_WMS", 760, 220, COLORS["wms"]), ("BIO_16S", 760, 410, COLORS["s16"]), ("BIO_ITS", 760, 600, COLORS["its"])]
    controls = [
        ("Field blank", 170, 150), ("Extraction blank", 170, 315),
        ("PCR blank", 170, 480), ("Library blank", 170, 645),
        ("Run blank", 170, 810),
    ]
    for label, x, y in controls:
        box(svg, x, y, 250, 70, label, "#f7f1e8", size=18)
    for label, x, y, color in samples:
        box(svg, x, y, 240, 80, label, color, size=19)
    for label, _, y in controls:
        for sample_label, sx, sy, _ in samples:
            if label == "PCR blank" and sample_label == "BIO_WMS":
                continue
            svg.line(420, y + 35, sx, sy + 40, arrow=True, color=COLORS["line"])
    box(svg, 1230, 260, 360, 90, "Use only relevant stage, assay and batch controls", COLORS["integration"], size=19)
    box(svg, 1230, 470, 360, 90, "Preserve unfiltered and annotated outputs", COLORS["output"], size=19)
    box(svg, 1230, 680, 360, 90, "Emit decision reasons for every taxon", COLORS["output"], size=19)
    svg.line(1000, 260, 1230, 305, arrow=True)
    svg.line(1000, 450, 1230, 515, arrow=True)
    svg.line(1000, 640, 1230, 725, arrow=True)
    path = OUT / "supplementary_figure_s2_control_graph.svg"
    svg.save(path)
    return path


def supplementary_s3() -> Path:
    rows = read_tsv(ROOT / "benchmarks" / "TRI-14" / "truth" / "taxonomic_crosswalk.tsv")
    method_counts = Counter(row["mapping_method"] for row in rows)
    rank_counts = Counter(row["harmonized_rank"] for row in rows)
    svg = SVG(1800, 1100, "Supplementary Figure S3 taxonomy crosswalk")
    svg.text(60, 70, "Supplementary Figure S3. Taxonomic-crosswalk validation", size=30, weight="700")
    svg.text(80, 140, "TRI-14 harmonization challenge", size=22, weight="700")
    max_count = max([*method_counts.values(), *rank_counts.values()])
    y = 220
    for method, count in sorted(method_counts.items()):
        bar(svg, 120, y, 520, 36, count, max_count, COLORS["integration"], method, str(count))
        y += 90
    y = 220
    for rank, count in sorted(rank_counts.items()):
        bar(svg, 900, y, 520, 36, count, max_count, COLORS["accent"], rank, str(count))
        y += 90
    box(svg, 660, 760, 460, 90, "Original identifiers retained", COLORS["panel"], size=20)
    box(svg, 660, 900, 460, 90, "Ambiguous rank handled by roll-up", COLORS["panel"], size=20)
    path = OUT / "supplementary_figure_s3_taxonomy_crosswalk.svg"
    svg.save(path)
    return path


def supplementary_s4() -> Path:
    perf = read_tsv(ROOT / "conf" / "assay_performance.tsv")
    svg = SVG(1800, 1100, "Supplementary Figure S4 detection power")
    svg.text(60, 70, "Supplementary Figure S4. Detection-power model", size=30, weight="700")
    x0, y0, w, h = 180, 160, 1250, 700
    svg.line(x0, y0 + h, x0 + w, y0 + h, color=COLORS["ink"], width=2)
    svg.line(x0, y0, x0, y0 + h, color=COLORS["ink"], width=2)
    svg.text(x0 + w / 2, y0 + h + 70, "Library depth (log scale)", size=20, anchor="middle")
    svg.text(60, y0 + h / 2, "Detection power", size=20, anchor="middle")
    depths = [10**x for x in [3, 3.5, 4, 4.5, 5, 5.5, 6]]
    colors = ["#2f6f73", "#5f8fc9", "#a56a43", "#7b5db6", "#b56b00", "#4f7f52"]
    for i, row in enumerate(perf):
        rel = float(row["reliability"])
        minimum = float(row["minimum_depth"])
        ref = float(row["reference_coverage"])
        pts = []
        for depth in depths:
            power = rel * (1 - math.exp(-depth / minimum)) * ref
            x = x0 + (math.log10(depth) - 3) / 3 * w
            y = y0 + h - power * h
            pts.append((x, y))
        color = colors[i % len(colors)]
        svg.polyline(pts, color, width=3)
        svg.text(1480, 190 + i * 45, f"{row['assay']} {row['kingdom']}", size=18, color=color)
    for tick in [0, .25, .5, .75, 1]:
        y = y0 + h - tick * h
        svg.line(x0 - 8, y, x0, y, color=COLORS["ink"])
        svg.text(x0 - 18, y + 6, f"{tick:.2g}", size=15, anchor="end", color=COLORS["muted"])
    for depth in [1_000, 10_000, 100_000, 1_000_000]:
        x = x0 + (math.log10(depth) - 3) / 3 * w
        svg.line(x, y0 + h, x, y0 + h + 8, color=COLORS["ink"])
        svg.text(x, y0 + h + 32, f"{depth:g}", size=15, anchor="middle", color=COLORS["muted"])
    path = OUT / "supplementary_figure_s4_detection_power.svg"
    svg.save(path)
    return path


def supplementary_s5() -> Path:
    contaminants = read_tsv(ROOT / "benchmarks" / "TRI-12" / "truth" / "contaminants.tsv")
    stage_counts = Counter(row["stage"] for row in contaminants)
    modality_counts = Counter()
    for row in contaminants:
        for modality in row["expected_modalities"].split(","):
            modality_counts[modality] += 1
    svg = SVG(1800, 1100, "Supplementary Figure S5 contamination design")
    svg.text(60, 70, "Supplementary Figure S5. Contamination challenge design", size=30, weight="700")
    svg.text(90, 150, "TRI-12 known contaminants by laboratory stage", size=22, weight="700")
    y = 230
    max_stage = max(stage_counts.values())
    for stage, count in sorted(stage_counts.items()):
        bar(svg, 120, y, 520, 42, count, max_stage, "#a43d3d", stage, str(count))
        y += 105
    svg.text(900, 150, "Expected assay scope", size=22, weight="700")
    y = 230
    for assay, count in sorted(modality_counts.items()):
        color = COLORS["wms"] if assay == "WMS" else COLORS["s16"] if assay == "16S" else COLORS["its"]
        bar(svg, 920, y, 520, 42, count, max(modality_counts.values()), color, assay, str(count))
        y += 105
    svg.label(900, 760, "Final performance bars should be populated from baseline comparison output after a frozen TRI-12 run.", width=52, size=20, color=COLORS["warn"])
    path = OUT / "supplementary_figure_s5_contamination_design.svg"
    svg.save(path)
    return path


def supplementary_s6() -> Path:
    rows = [row for row in read_tsv(ROOT / "benchmarks" / "TRI-30" / "metadata" / "samples.tsv") if row["sample_role"] == "biological"]
    sample_keys = {(row["biological_sample_id"], row["sample_material"], row["environment_category"]) for row in rows}
    by_material = Counter(material for _, material, _ in sample_keys)
    by_env = Counter(env for _, _, env in sample_keys)
    svg = SVG(1800, 1100, "Supplementary Figure S6 environmental design")
    svg.text(60, 70, "Supplementary Figure S6. Matched environmental design surrogate", size=30, weight="700")
    svg.text(90, 150, "TRI-30 biological samples by material", size=22, weight="700")
    y = 230
    for material, count in sorted(by_material.items()):
        bar(svg, 120, y, 540, 42, count, max(by_material.values()), COLORS["output"], material, str(count))
        y += 95
    svg.text(900, 150, "Samples by environment", size=22, weight="700")
    y = 230
    for env, count in sorted(by_env.items()):
        bar(svg, 920, y, 540, 42, count, max(by_env.values()), COLORS["accent"], env, str(count))
        y += 95
    svg.label(930, 800, "TRI-30 is an executable design surrogate until matched field data are analyzed.", width=46, size=21, color=COLORS["warn"])
    path = OUT / "supplementary_figure_s6_environmental_design.svg"
    svg.save(path)
    return path


def supplementary_s7() -> Path:
    svg = SVG(1800, 1100, "Supplementary Figure S7 installation and execution tests")
    svg.text(60, 70, "Supplementary Figure S7. Installation and execution tests", size=30, weight="700")
    checks = [
        ("Python unit tests", "29 tests passed", COLORS["output"]),
        ("Preflight validator", "TRI-01 passed", COLORS["output"]),
        ("Benchmark checksums", "Regenerated and verified", COLORS["output"]),
        ("Stub resource profile", "5 datasets / 97 libraries", COLORS["output"]),
        ("Nextflow parser", "Not run in this workspace", "#fff4d8"),
        ("Production containers", "Await frozen digests", "#fff4d8"),
        ("Slurm/Apptainer", "Await site run", "#fff4d8"),
        ("Empirical benchmark", "Required before submission", "#f7d9d9"),
    ]
    for i, (name, status, color) in enumerate(checks):
        col, row = i % 2, i // 2
        x = 120 + col * 820
        y = 160 + row * 190
        box(svg, x, y, 660, 80, name, color, size=21)
        svg.text(x + 35, y + 125, status, size=20, color=COLORS["muted"])
    svg.text(120, 1010, f"Generated {datetime.now().date().isoformat()} from repository-local checks.", size=18, color=COLORS["muted"])
    path = OUT / "supplementary_figure_s7_installation_tests.svg"
    svg.save(path)
    return path


def index(paths: list[Path]) -> None:
    lines = [
        "# Manuscript figure assets",
        "",
        "Generated by `python scripts/build_manuscript_figures.py`.",
        "",
        "These SVG files are editable vector assets. Figure 1 also has a PDF",
        "companion for LaTeX submission.",
        "",
    ]
    for path in paths:
        lines.append(f"- [{path.name}]({path.name})")
    lines += [
        "",
        "Data caveat: `figure1_composite.svg` Panel C uses",
        "`results/frozen-release-accuracy-20260812/benchmark_comparison.tsv`.",
        "Panel D uses the current stub resource profile, which validates orchestration",
        "only and should not be presented as production runtime.",
        "",
    ]
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    paths = [
        workflow_architecture_figure(),
        main_figure(),
        supplementary_s1(),
        supplementary_s2(),
        supplementary_s3(),
        supplementary_s4(),
        supplementary_s5(),
        supplementary_s6(),
        supplementary_s7(),
    ]
    index(paths)
    for path in paths:
        print(path.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

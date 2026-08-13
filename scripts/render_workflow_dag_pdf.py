#!/usr/bin/env python3
"""Render the overall TRIARCHY workflow DAG as a vector PDF.

The matching editable Mermaid source is docs/overall-workflow-dag.mmd.
"""

import argparse
import os
from pathlib import Path

try:
    from reportlab.lib.colors import HexColor
    from reportlab.lib.pagesizes import TABLOID, landscape
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.pdfbase.pdfmetrics import stringWidth
    from reportlab.pdfgen import canvas
except ModuleNotFoundError as error:
    if error.name == "reportlab":
        raise SystemExit(
            "reportlab is required; run ./install.sh or install environment.yml, "
            "then use .tools/env/bin/python"
        ) from error
    raise


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "output" / "pdf" / "triarchy-overall-workflow-dag.pdf"


def configure_fonts() -> tuple[str, str]:
    """Use embedded TrueType fonts when available, otherwise PDF core fonts."""
    candidates = []
    if os.environ.get("TRIARCHY_DAG_FONT_DIR"):
        candidates.append(Path(os.environ["TRIARCHY_DAG_FONT_DIR"]))
    candidates.extend((
        Path("/usr/share/fonts/truetype/dejavu"),
        Path("/Library/Fonts"),
    ))
    for directory in candidates:
        regular = directory / "DejaVuSans.ttf"
        bold = directory / "DejaVuSans-Bold.ttf"
        if regular.is_file() and bold.is_file():
            pdfmetrics.registerFont(TTFont("TRIARCHY-Regular", str(regular)))
            pdfmetrics.registerFont(TTFont("TRIARCHY-Bold", str(bold)))
            return "TRIARCHY-Regular", "TRIARCHY-Bold"
    return "Helvetica", "Helvetica-Bold"


REGULAR, BOLD = configure_fonts()

PAGE_W, PAGE_H = landscape(TABLOID)
NAVY = HexColor("#18324A")
INK = HexColor("#17222D")
MUTED = HexColor("#556574")
LINE = HexColor("#8293A1")
INPUT = HexColor("#E8EEF3")
WMS = HexColor("#D8EAF6")
S16 = HexColor("#DDF1E5")
ITS = HexColor("#F3E5D2")
INTEGRATION = HexColor("#E9E1F5")
OUTPUT_FILL = HexColor("#D9ECE9")
WHITE = HexColor("#FFFFFF")


def wrap(text: str, font: str, size: float, width: float) -> list[str]:
    lines = []
    for paragraph in text.split("\n"):
        words = paragraph.split()
        current = ""
        for word in words:
            trial = f"{current} {word}".strip()
            if current and stringWidth(trial, font, size) > width:
                lines.append(current)
                current = word
            else:
                current = trial
        lines.append(current)
    return lines


def node(c, x, y, w, h, text, fill, border=NAVY, font_size=9.2):
    c.setFillColor(fill)
    c.setStrokeColor(border)
    c.setLineWidth(1.1)
    c.roundRect(x, y, w, h, 8, fill=1, stroke=1)
    lines = wrap(text, BOLD, font_size, w - 16)
    total = len(lines) * (font_size + 2)
    baseline = y + (h + total) / 2 - font_size
    c.setFillColor(INK)
    c.setFont(BOLD, font_size)
    for line in lines:
        c.drawCentredString(x + w / 2, baseline, line)
        baseline -= font_size + 2


def arrow(c, x1, y1, x2, y2, color=LINE, width=1.35):
    c.setStrokeColor(color)
    c.setFillColor(color)
    c.setLineWidth(width)
    c.line(x1, y1, x2, y2)
    import math
    angle = math.atan2(y2 - y1, x2 - x1)
    length = 7
    spread = 0.48
    p1 = (x2 - length * math.cos(angle - spread), y2 - length * math.sin(angle - spread))
    p2 = (x2 - length * math.cos(angle + spread), y2 - length * math.sin(angle + spread))
    path = c.beginPath()
    path.moveTo(x2, y2)
    path.lineTo(*p1)
    path.lineTo(*p2)
    path.close()
    c.drawPath(path, fill=1, stroke=0)


def poly_arrow(c, points, color=LINE, width=1.35):
    c.setStrokeColor(color)
    c.setLineWidth(width)
    path = c.beginPath()
    path.moveTo(*points[0])
    for point in points[1:-1]:
        path.lineTo(*point)
    c.drawPath(path, stroke=1, fill=0)
    arrow(c, *points[-2], *points[-1], color=color, width=width)


def label(c, x, y, text, color=MUTED, size=8.5):
    c.setFillColor(color)
    c.setFont(REGULAR, size)
    c.drawString(x, y, text)


def render(output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    c = canvas.Canvas(str(output), pagesize=(PAGE_W, PAGE_H))
    c.setTitle("TRIARCHY overall workflow DAG")
    c.setAuthor("TRIARCHY")

    c.setFillColor(WHITE)
    c.rect(0, 0, PAGE_W, PAGE_H, fill=1, stroke=0)
    c.setFillColor(NAVY)
    c.rect(0, PAGE_H - 58, PAGE_W, 58, fill=1, stroke=0)
    c.setFillColor(WHITE)
    c.setFont(BOLD, 21)
    c.drawString(42, PAGE_H - 37, "TRIARCHY overall workflow DAG")
    c.setFont(REGULAR, 9.5)
    c.drawRightString(PAGE_W - 42, PAGE_H - 34, "Control-aware WMS + 16S + ITS integration")

    # Phase labels
    for x, title in ((42, "INPUT AND VALIDATION"), (250, "ASSAY-SPECIFIC PROCESSING"),
                     (718, "CONTROL-AWARE INTEGRATION"), (1034, "REPORTING")):
        c.setFillColor(MUTED)
        c.setFont(BOLD, 9)
        c.drawString(x, PAGE_H - 84, title)

    # Nodes
    node(c, 42, 615, 160, 58, "Sample manifest\nControl relationships", INPUT)
    node(c, 42, 520, 160, 54, "Validate metadata", INPUT)
    node(c, 42, 425, 160, 54, "Split manifests by assay", INPUT)
    node(c, 42, 190, 160, 54, "Run provenance", INPUT)

    node(c, 250, 610, 124, 52, "WMS quality control\nfastp", WMS)
    node(c, 402, 610, 124, 52, "Host removal\nBowtie2", WMS)
    node(c, 554, 655, 130, 52, "Taxonomy\nKraken 2 + Bracken", WMS)
    node(c, 554, 565, 130, 52, "Functions\nHUMAnN", WMS)

    node(c, 250, 430, 124, 52, "16S primer/QC\ncutadapt", S16)
    node(c, 472, 430, 150, 52, "16S ASVs + taxonomy\nDADA2", S16)

    node(c, 250, 250, 124, 52, "ITS primer/QC\ncutadapt", ITS)
    node(c, 402, 250, 124, 52, "ITS extraction\nITSxpress", ITS)
    node(c, 554, 250, 130, 52, "ITS ASVs + taxonomy\nDADA2", ITS)

    node(c, 472, 345, 150, 48, "Library QC summary", INPUT)

    node(c, 718, 510, 168, 62, "Linked-control annotation\nTaxa + functions", INTEGRATION)
    node(c, 718, 400, 168, 58, "Taxonomic harmonization\nCanonical shared ranks", INTEGRATION)
    node(c, 922, 485, 170, 58, "Authenticity + explanation\nPolicy + assay reliability", INTEGRATION)
    node(c, 922, 390, 170, 58, "Cross-assay concordance", INTEGRATION)
    node(c, 718, 285, 168, 58, "Assay-specific diversity\nAlpha + beta", INTEGRATION)
    node(c, 718, 190, 168, 54, "Sample confidence", INTEGRATION)

    node(c, 1034, 260, 160, 78, "Interactive report\nEvidence + concordance\nDiversity + provenance", OUTPUT_FILL)
    node(c, 1034, 135, 160, 74, "Published outputs\nQC / profiles / integrated\ndiversity / pipeline info", OUTPUT_FILL)

    # Input flow
    arrow(c, 122, 615, 122, 574)
    arrow(c, 122, 520, 122, 479)
    poly_arrow(c, [(122, 425), (122, 370), (226, 370), (226, 636), (250, 636)])
    poly_arrow(c, [(122, 425), (226, 425), (226, 456), (250, 456)])
    poly_arrow(c, [(122, 425), (226, 425), (226, 276), (250, 276)])
    poly_arrow(c, [(42, 642), (20, 642), (20, 217), (42, 217)])

    # Assay branches
    arrow(c, 374, 636, 402, 636)
    poly_arrow(c, [(526, 636), (540, 636), (540, 681), (554, 681)])
    poly_arrow(c, [(526, 636), (540, 636), (540, 591), (554, 591)])
    arrow(c, 374, 456, 472, 456)
    arrow(c, 374, 276, 402, 276)
    arrow(c, 526, 276, 554, 276)

    # QC aggregation
    poly_arrow(c, [(312, 610), (312, 369), (472, 369)])
    poly_arrow(c, [(464, 610), (464, 369), (472, 369)])
    poly_arrow(c, [(312, 430), (312, 369), (472, 369)])
    poly_arrow(c, [(312, 250), (312, 369), (472, 369)])

    # Profiles to control annotation
    poly_arrow(c, [(684, 681), (700, 681), (700, 541), (718, 541)])
    poly_arrow(c, [(684, 591), (700, 591), (700, 541), (718, 541)])
    poly_arrow(c, [(622, 456), (700, 456), (700, 541), (718, 541)])
    poly_arrow(c, [(684, 276), (700, 276), (700, 541), (718, 541)])

    # Integration
    arrow(c, 802, 510, 802, 458)
    poly_arrow(c, [(886, 429), (904, 429), (904, 514), (922, 514)])
    arrow(c, 886, 429, 922, 419)
    poly_arrow(c, [(802, 510), (696, 510), (696, 314), (718, 314)])
    poly_arrow(c, [(802, 510), (696, 510), (696, 260), (802, 260), (802, 244)])
    poly_arrow(c, [(622, 369), (905, 369), (905, 514), (922, 514)])

    # Reporting fan-in
    poly_arrow(c, [(1092, 514), (1112, 514), (1112, 338)])
    poly_arrow(c, [(1092, 419), (1126, 419), (1126, 338)])
    poly_arrow(c, [(886, 314), (1034, 314)])
    poly_arrow(c, [(886, 205), (978, 205), (978, 299), (1034, 299)])
    poly_arrow(c, [(202, 217), (680, 217), (680, 112), (1000, 112), (1000, 282), (1034, 282)])
    arrow(c, 1114, 260, 1114, 209)

    # Compact legend
    legend_y = 76
    items = [(INPUT, "Inputs / QC"), (WMS, "WMS"), (S16, "16S"),
             (ITS, "ITS"), (INTEGRATION, "Integration"), (OUTPUT_FILL, "Outputs")]
    x = 42
    for fill, text in items:
        c.setFillColor(fill)
        c.setStrokeColor(LINE)
        c.roundRect(x, legend_y, 14, 14, 3, fill=1, stroke=1)
        label(c, x + 20, legend_y + 2, text, size=8.5)
        x += 115

    c.setFillColor(MUTED)
    c.setFont(REGULAR, 7.5)
    c.drawRightString(PAGE_W - 42, 28, "Derived from main.nf and modules/local/*.nf")
    c.save()
    print(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help=f"PDF destination (default: {DEFAULT_OUTPUT.relative_to(ROOT)})",
    )
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    render(output.resolve())


if __name__ == "__main__":
    main()

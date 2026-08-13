#!/usr/bin/env python3
"""Convert simple repository-generated SVG figures to vector PDF.

This intentionally supports the small SVG subset emitted by
`scripts/build_manuscript_figures.py`: rectangles, lines, circles, polylines,
and text. It avoids external converter dependencies in lightweight workspaces.
"""

from __future__ import annotations

import argparse
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path


SVG_NS = "{http://www.w3.org/2000/svg}"


def tag(element: ET.Element) -> str:
    return element.tag.removeprefix(SVG_NS)


def number(value: str | None, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    return float(str(value).replace("px", ""))


def color(value: str | None) -> tuple[float, float, float] | None:
    if not value or value == "none":
        return None
    value = value.strip()
    if re.fullmatch(r"#[0-9a-fA-F]{6}", value):
        return (
            int(value[1:3], 16) / 255,
            int(value[3:5], 16) / 255,
            int(value[5:7], 16) / 255,
        )
    if value == "white":
        return (1, 1, 1)
    return (0, 0, 0)


def pdf_string(value: str) -> str:
    return value.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


class PDF:
    def __init__(self, width: float, height: float):
        self.width = width
        self.height = height
        self.commands: list[str] = []

    def y(self, value: float) -> float:
        return self.height - value

    def set_fill(self, value: str | None) -> None:
        rgb = color(value)
        if rgb is not None:
            self.commands.append(f"{rgb[0]:.4f} {rgb[1]:.4f} {rgb[2]:.4f} rg")

    def set_stroke(self, value: str | None) -> None:
        rgb = color(value)
        if rgb is not None:
            self.commands.append(f"{rgb[0]:.4f} {rgb[1]:.4f} {rgb[2]:.4f} RG")

    def rect(self, element: ET.Element) -> None:
        x = number(element.get("x"))
        y = number(element.get("y"))
        w = number(element.get("width"))
        h = number(element.get("height"))
        fill = color(element.get("fill"))
        stroke = color(element.get("stroke"))
        sw = number(element.get("stroke-width"), 1.0)
        if fill is not None:
            self.set_fill(element.get("fill"))
            self.commands.append(f"{x:.3f} {self.y(y + h):.3f} {w:.3f} {h:.3f} re f")
        if stroke is not None and sw > 0:
            self.set_stroke(element.get("stroke"))
            self.commands.append(f"{sw:.3f} w {x:.3f} {self.y(y + h):.3f} {w:.3f} {h:.3f} re S")

    def line(self, element: ET.Element) -> None:
        x1 = number(element.get("x1"))
        y1 = number(element.get("y1"))
        x2 = number(element.get("x2"))
        y2 = number(element.get("y2"))
        self.set_stroke(element.get("stroke"))
        self.commands.append(f"{number(element.get('stroke-width'), 1.0):.3f} w")
        self.commands.append(f"{x1:.3f} {self.y(y1):.3f} m {x2:.3f} {self.y(y2):.3f} l S")
        if element.get("marker-end"):
            self.arrowhead(x1, y1, x2, y2, element.get("stroke"))

    def arrowhead(self, x1: float, y1: float, x2: float, y2: float, stroke: str | None) -> None:
        angle = math.atan2(y2 - y1, x2 - x1)
        length = 14
        spread = math.radians(26)
        pts = [(x2, y2)]
        for sign in (1, -1):
            theta = angle + math.pi + sign * spread
            pts.append((x2 + length * math.cos(theta), y2 + length * math.sin(theta)))
        self.set_fill(stroke)
        self.commands.append(
            f"{pts[0][0]:.3f} {self.y(pts[0][1]):.3f} m "
            f"{pts[1][0]:.3f} {self.y(pts[1][1]):.3f} l "
            f"{pts[2][0]:.3f} {self.y(pts[2][1]):.3f} l h f"
        )

    def circle(self, element: ET.Element) -> None:
        cx = number(element.get("cx"))
        cy = number(element.get("cy"))
        r = number(element.get("r"))
        k = 0.5522847498
        y = self.y(cy)
        self.set_fill(element.get("fill"))
        self.commands.append(
            f"{cx + r:.3f} {y:.3f} m "
            f"{cx + r:.3f} {y + k*r:.3f} {cx + k*r:.3f} {y + r:.3f} {cx:.3f} {y + r:.3f} c "
            f"{cx - k*r:.3f} {y + r:.3f} {cx - r:.3f} {y + k*r:.3f} {cx - r:.3f} {y:.3f} c "
            f"{cx - r:.3f} {y - k*r:.3f} {cx - k*r:.3f} {y - r:.3f} {cx:.3f} {y - r:.3f} c "
            f"{cx + k*r:.3f} {y - r:.3f} {cx + r:.3f} {y - k*r:.3f} {cx + r:.3f} {y:.3f} c f"
        )
        if color(element.get("stroke")) is not None:
            self.set_stroke(element.get("stroke"))
            self.commands.append(f"{number(element.get('stroke-width'), 1.0):.3f} w")

    def polyline(self, element: ET.Element) -> None:
        points = []
        for pair in (element.get("points") or "").split():
            x, y = pair.split(",", 1)
            points.append((float(x), float(y)))
        if not points:
            return
        self.set_stroke(element.get("stroke"))
        self.commands.append(f"{number(element.get('stroke-width'), 1.0):.3f} w")
        first, *rest = points
        parts = [f"{first[0]:.3f} {self.y(first[1]):.3f} m"]
        parts += [f"{x:.3f} {self.y(y):.3f} l" for x, y in rest]
        self.commands.append(" ".join(parts) + " S")

    def text(self, element: ET.Element) -> None:
        value = "".join(element.itertext())
        size = number(element.get("font-size"), 12)
        weight = element.get("font-weight", "400")
        font = "F2" if weight in {"600", "700", "bold"} else "F1"
        x = number(element.get("x"))
        y = number(element.get("y"))
        anchor = element.get("text-anchor", "start")
        if anchor != "start":
            width = len(value) * size * 0.52
            x -= width / 2 if anchor == "middle" else width
        self.set_fill(element.get("fill"))
        self.commands.append(f"BT /{font} {size:.3f} Tf {x:.3f} {self.y(y):.3f} Td ({pdf_string(value)}) Tj ET")

    def write(self, output: Path) -> None:
        stream = "\n".join(self.commands).encode("latin-1", "replace")
        objects: list[bytes] = [
            b"<< /Type /Catalog /Pages 2 0 R >>",
            b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
            (
                f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {self.width:.3f} {self.height:.3f}] "
                f"/Resources << /Font << /F1 4 0 R /F2 5 0 R >> >> /Contents 6 0 R >>"
            ).encode(),
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>",
            b"<< /Length " + str(len(stream)).encode() + b" >>\nstream\n" + stream + b"\nendstream",
        ]
        data = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
        offsets = [0]
        for index, obj in enumerate(objects, start=1):
            offsets.append(len(data))
            data.extend(f"{index} 0 obj\n".encode())
            data.extend(obj)
            data.extend(b"\nendobj\n")
        xref = len(data)
        data.extend(f"xref\n0 {len(objects)+1}\n0000000000 65535 f \n".encode())
        for offset in offsets[1:]:
            data.extend(f"{offset:010d} 00000 n \n".encode())
        data.extend(
            f"trailer << /Size {len(objects)+1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n".encode()
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(data)


def convert(input_path: Path, output_path: Path) -> None:
    root = ET.parse(input_path).getroot()
    pdf = PDF(number(root.get("width")), number(root.get("height")))
    for element in root.iter():
        name = tag(element)
        if name == "rect":
            pdf.rect(element)
        elif name == "line":
            pdf.line(element)
        elif name == "circle":
            pdf.circle(element)
        elif name == "polyline":
            pdf.polyline(element)
        elif name == "text":
            pdf.text(element)
    pdf.write(output_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    convert(args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

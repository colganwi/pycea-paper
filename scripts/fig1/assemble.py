"""Assemble Figure 1 from panel files -> manuscript/figures/figure1.{svg,pdf}.

Layout (183 mm wide):
  row 1: A ecosystem schematic | B supported data types (plots/fig1/datatypes.png, user-made: molecular recording on
         the left, live-cell tracking on the right, each labeled above its half). A is 10% taller than the
         height both images would share to fill the row, and B takes the remaining width. The letters and B labels
         sit on one line above B; A's drawn content starts level with the top of the letters.
  row 2: C example workflow, four steps across the full width; tool logos below
         the plots (Pycea logo centered under steps 3 and 4, flanked by bars)
"""

import struct
import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import FIGURES, PLOTS, SvgFigure

W, H = 183, 118
M = 3  # outer margin (mm)
LETTER, GAP = 4.5, 2  # letter gutter, gap between panels (mm)
TITLE_H = 5  # mm above the row 1 images: panel letters and the B labels sit on this line
A_SCALE = 1.1  # A is drawn this much taller than the shared height that would fill the row; B takes the remaining width
fig = SvgFigure(W, H)


def aspect(png):
    width, height = struct.unpack(">II", Path(png).read_bytes()[16:24])  # PNG IHDR
    return width / height


def content_top(png):
    """Top of the drawn content (fraction of the height) of a PNG."""
    rgba = np.asarray(Image.open(png).convert("RGBA"))
    drawn = ((rgba[:, :, 3] > 0) & (rgba[:, :, :3].min(axis=2) < 250)).any(axis=1)
    return np.argmax(drawn) / rgba.shape[0]


def half_centers(png):
    """Horizontal centers (fractions of the width) of the content left and right of the widest empty column gap."""
    rgba = np.asarray(Image.open(png).convert("RGBA"))
    drawn = ((rgba[:, :, 3] > 0) & (rgba[:, :, :3].min(axis=2) < 250)).any(axis=0)
    cols = np.where(drawn)[0]
    split = np.argmax(np.diff(cols))  # widest gap between drawn columns
    left, right = (cols[0], cols[split]), (cols[split + 1], cols[-1])
    return [(a + b + 1) / 2 / rgba.shape[1] for a, b in (left, right)]


# --- Row 1 ---
eco, types = PLOTS / "ecosystem.png", PLOTS / "fig1" / "datatypes.png"
row_w = W - 2 * M - 2 * LETTER - GAP
a_h = A_SCALE * row_w / (aspect(eco) + aspect(types))
y_img = M + TITLE_H
LETTER_TOP = M + 3.2 - 10 * 0.716 * 25.4 / 72  # cap top of the 10 pt bold letters (baseline M + 3.2)
fig.letter("A", M, M + 3.2)
a_top = LETTER_TOP - content_top(eco) * a_h  # A's drawn content starts level with the top of the letters
a_w, _ = fig.image(eco, M + LETTER, a_top, height=a_h)
xb = M + LETTER + a_w + GAP
fig.letter("B", xb, M + 3.2)
b_w, b_h = fig.image(types, xb + LETTER, y_img, width=row_w - a_w)
for frac, label in zip(half_centers(types), ["Molecular recording", "Live-cell tracking"], strict=True):
    fig.text(label, xb + LETTER + frac * b_w, M + 3.2, size=7, weight="bold", anchor="middle")
row1_h = max(a_top + a_h, y_img + b_h) - M
print(f"row 1: A {a_w:.1f} x {a_h:.1f} mm, B {b_w:.1f} x {b_h:.1f} mm; right edge {xb + LETTER + b_w:.1f} mm")

# --- Row 2: workflow ---
y2 = M + row1_h + 2
fig.letter("C", M, y2 + 4)
steps = [
    ("1. Lineage tree reconstruction", "c1_reconstruction.svg"),
    ("2. UMAP and clustering", "c2_umap.svg"),
    ("3. Lineage visualization", "c3_tree.svg"),
    ("4. Lineage analysis", "c4_linkage.svg"),
]
step_w = (W - 2 * M) / len(steps)
plot_top, plot_h = y2 + 7, 40
for i, (title, svg) in enumerate(steps):
    x = M + i * step_w
    fig.text(title, x + step_w / 2, y2 + 4.5, size=7, weight="bold", anchor="middle")
    pw, ph = fig.place(PLOTS / "fig1" / svg, x + 2, plot_top, height=plot_h)
    if pw > step_w - 4:  # keep wide panels inside their column
        fig.parts.pop()
        fig.place(PLOTS / "fig1" / svg, x + 2, plot_top, width=step_w - 4)
    if i < len(steps) - 1:
        fig.arrow(x + step_w - 3, plot_top + plot_h / 2, x + step_w + 3)

# --- Logos below the plots ---
logo_y = plot_top + plot_h + 1
logo_h = 10
centers = [M + (i + 0.5) * step_w for i in range(4)]
for name, cx, h in [("cassipopeia_logo.png", centers[0], 5), ("scanpy_logo.png", centers[1], 9)]:
    w, _ = fig.image(PLOTS / name, 0, -100, height=h)  # measure, then replace below
    fig.parts.pop()
    fig.image(PLOTS / name, cx - w / 2, logo_y + (logo_h - h) / 2, height=h)

# Pycea logo centered between steps 3 and 4, with bars on either side spanning each step
cx = M + 3 * step_w
w, _ = fig.image(PLOTS / "pycea_logo.png", 0, -100, height=logo_h)
fig.parts.pop()
fig.image(PLOTS / "pycea_logo.png", cx - w / 2, logo_y, height=logo_h)
bar_y = logo_y + logo_h / 2
pad = 2
fig.line(M + 2 * step_w + pad, bar_y, cx - w / 2 - pad, bar_y)
fig.line(cx + w / 2 + pad, bar_y, M + 4 * step_w - pad, bar_y)

fig.h = logo_y + logo_h + M  # fit page height to content
fig.save(FIGURES / "figure1")
print("wrote", FIGURES / "figure1.pdf")

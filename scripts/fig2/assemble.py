"""Assemble Figure 2 from panel files -> manuscript/figures/figure2.{svg,pdf}.

Layout (183 mm wide), work in progress:
  row 1: A TreeData schematic | B alignment options | C multiple trees
  Panel letters for B and C sit in a left gutter so the panels start at the top of the row.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import FIGURES, PLOTS, SvgFigure, visible_extent

W = 183
M = 3
GAP = 3  # horizontal gap between columns
fig = SvgFigure(W, 100)

# --- Row 1, panel B geometry (computed first so that A can be matched to its height) ---
options = [
    ("b1_leaves.svg", "alignment = “leaves”", "only leaves are observed"),
    ("b2_nodes.svg", "alignment = “nodes”", "leaves and internal nodes are observed"),
    ("b3_subset.svg", "alignment = “subset”", "a subset of leaves and internal nodes are observed"),
]
top = M + 3  # top of the B and C panels (their letters sit in the left gutter)
LETTER_GUTTER = 5  # mm reserved left of B and C panels for the panel letters
b_width = 42  # mm (width of every B and C panel, so scale and UMAP size are identical; A is scaled to the height of B)
LABEL_H, SPACING = 5.5, 1.0
heights = []
for svg, *_ in options:
    pw, ph = SvgFigure._size(PLOTS / "fig2" / svg)
    heights.append(b_width * ph / pw)
b_height = sum(LABEL_H + h for h in heights) + SPACING * (len(options) - 1)

# --- A: TreeData schematic (PNG), 20% smaller than the height of B ---
A_SCALE = 0.92  # fraction of the height of B (15% larger than the previous 0.8)
a_height = A_SCALE * b_height
fig.letter("A", M, top + 4)  # same baseline as the B and C letters
a_width, _ = fig.image(PLOTS / "treedata_schema.png", M, M + 4, height=a_height)

# --- Shift B and C right so the visible right edge of C matches the right edge of D (W - M) ---
c_right = max(visible_extent(PLOTS / "fig2" / svg)[1] for svg in ["c1_non_overlap.svg", "c2_overlap.svg"])
shift = max(0, (W - M) - (W - M - b_width + c_right * b_width))

# --- B ---
xb = M + a_width + GAP + shift
fig.letter("B", xb, top + 4)
xb_panels = xb + LETTER_GUTTER
y = top
for (svg, title, desc), h in zip(options, heights, strict=True):
    fig.text(title, xb_panels, y + 2.3, size=6.5, weight="bold")
    fig.text(desc, xb_panels, y + 4.8, size=5.5, color="#555555")
    fig.place(PLOTS / "fig2" / svg, xb_panels, y + LABEL_H, width=b_width)
    y += LABEL_H + h + SPACING

# --- C: multiple trees ---
xc_panels = W - M - b_width + shift  # visible right edge of C aligns with the right edge of row 2 (panel D)
xc = xc_panels - LETTER_GUTTER
fig.letter("C", xc, top + 4)
y = top
for svg, title, desc in [
    ("c1_non_overlap.svg", "allow_overlap = False", "multiple trees, disjoint observations"),
    ("c2_overlap.svg", "allow_overlap = True", "multiple trees, same observations"),
]:
    fig.text(title, xc_panels, y + 2.3, size=6.5, weight="bold")
    fig.text(desc, xc_panels, y + 4.8, size=5.5, color="#555555")
    _, h = fig.place(PLOTS / "fig2" / svg, xc_panels, y + LABEL_H, width=b_width)
    y += LABEL_H + h + SPACING

# --- Row 2: D, a tree plot built up with successive Pycea calls ---
y2 = max(M + 4 + a_height, y) + 3
fig.letter("D", M, y2 + 4)
steps = [
    ("d1_branches.svg", "pycea.pl.branches", "branches colored by length"),
    ("d2_nodes.svg", "pycea.pl.nodes", "Add nodes colored by cluster"),
    ("d3_star.svg", "pycea.pl.nodes", "Add a single highlighted node"),
    ("d4_annotation.svg", "pycea.pl.annotation", "Add gene expression at the leaves"),
]
# All four panels share one height (their widths differ because panel 4 includes a legend);
# panels are packed from the left, with a gap for the arrows.
d_x0 = M + LETTER_GUTTER
ARROW_GAP = 6
sizes = [SvgFigure._size(PLOTS / "fig2" / svg) for svg, *_ in steps]
aspect_sum = sum(pw / ph for pw, ph in sizes)
d_h = (W - M - d_x0 - ARROW_GAP * (len(steps) - 1)) / aspect_sum
x = d_x0
for i, ((svg, func, desc), (pw, ph)) in enumerate(zip(steps, sizes, strict=True)):
    w = d_h * pw / ph
    center = x + w / 2
    if svg == "d4_annotation.svg":  # center over the plot, excluding the legends to its right
        plot_w, _ = SvgFigure._size(PLOTS / "fig2" / "d4_annotation_plot_only.svg")
        center = x + w * (plot_w / pw) / 2
    fig.text(func, center, y2 + 3, size=7, weight="bold", anchor="middle", mono=True)
    fig.text(desc, center, y2 + 6, size=5.5, anchor="middle", color="#555555")
    fig.place(PLOTS / "fig2" / svg, x, y2 + 7.5, width=w)
    if i < len(steps) - 1:
        fig.arrow(x + w + 1, y2 + 7.5 + d_h / 2, x + w + ARROW_GAP - 1, width=0.6, head=2)
    x += w + ARROW_GAP
d_bottom = y2 + 7.5 + d_h

fig.h = d_bottom + M
fig.save(FIGURES / "figure2")
print("wrote", FIGURES / "figure2.pdf", f"(A: {a_width:.0f} x {a_height:.0f} mm; B: {b_height:.0f} mm tall)")

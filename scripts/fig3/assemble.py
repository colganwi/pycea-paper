"""Assemble Figure 3 from panel SVGs -> manuscript/figures/figure3.{svg,pdf}.

The panels (except the polar tree A, which has no text and is scaled to its box width) are drawn at their final size,
so they are placed 1:1 and keep their font sizes. Box sizes and row layout are in layout.py; letters sit on the
top-left corner of their panel (G has its own gutter).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import FIGURES, PLOTS, SvgFigure
from layout import GAP, LETTER, M, ROW_GAP, ROWS, SUB_GAP, W

P = PLOTS / "fig3"
PT = 25.4 / 72  # mm per SVG point (1:1)
fig = SvgFigure(W, 100)

y = M
for row in ROWS:
    x, row_h = M, 0
    for i, (letter, name, width, _, gutter) in enumerate(row):
        if i:
            x += GAP if letter else SUB_GAP
        if letter:
            fig.letter(letter, x, y + 3.2)
            x += LETTER if gutter else 0
        pw, ph = SvgFigure._size(P / f"{name}.svg")
        if name != "a_tree" and abs(pw * PT - width) > 0.05:
            print(f"warning: {name} is {pw * PT:.1f} mm wide, its box is {width} mm")
        _, h = fig.place(P / f"{name}.svg", x, y, width=width)
        if name == "a_tree":  # which ring is which
            fig.text("Inner ring: germ layer; outer ring: lineage", x + width / 2, y + h + 3.2, size=5.5, anchor="middle")
            h += 3.2
        x += width
        row_h = max(row_h, h)
    print(f"row at y={y:.1f}: right edge {x:.1f} mm (page edge {W - M}), height {row_h:.1f} mm")
    y += row_h + ROW_GAP

fig.h = y - ROW_GAP + M
fig.save(FIGURES / "figure3")
print(f"wrote {FIGURES / 'figure3.pdf'} ({W} x {fig.h:.1f} mm)")

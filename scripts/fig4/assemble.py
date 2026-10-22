"""Assemble Figure 4 from panel SVGs -> manuscript/figures/figure4.{svg,pdf}.

The panels are drawn at their final size (boxes in layout.py), so they are placed 1:1 and keep their font sizes.
Letters sit in a gutter left of their panel, or on the empty top-left corner of the polar trees.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import FIGURES, PLOTS, SvgFigure
from layout import GAP, LETTER, M, NUDGE, ROW_GAP, ROWS, W

P = PLOTS / "fig4"
PT = 25.4 / 72  # mm per SVG point (1:1)
fig = SvgFigure(W, 100)

y, bottom = M, M
for row in ROWS:
    x, row_h = M, 0
    for i, (letter, name, width, _, gutter) in enumerate(row):
        x += GAP if i else 0
        fig.letter(letter, x, y + 3.2)
        x += LETTER if gutter else 0
        pw, _ = SvgFigure._size(P / f"{name}.svg")
        if abs(pw * PT - width) > 0.05:
            print(f"warning: {name} is {pw * PT:.1f} mm wide, its box is {width} mm")
        dx, dy = NUDGE.get(name, (0, 0))
        _, h = fig.place(P / f"{name}.svg", x + dx, y + dy, width=width)
        x += width
        row_h = max(row_h, h)
        bottom = max(bottom, y + dy + h)
    print(f"row at y={y:.1f}: right edge {x:.1f} mm (page edge {W - M}), height {row_h:.1f} mm")
    y += row_h + ROW_GAP

fig.h = bottom + M
fig.save(FIGURES / "figure4")
print(f"wrote {FIGURES / 'figure4.pdf'} ({W} x {fig.h:.1f} mm)")

"""Assemble Figure S1 (benchmarks) from panel SVGs -> manuscript/figures/figureS1.{svg,pdf}. One row, three panels."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import FIGURES, PLOTS, SvgFigure

W, M, GUTTER = 183, 3, 5
fig = SvgFigure(W, 100)
col = (W - 2 * M) / 3
h_max = 0
for i, (letter, svg) in enumerate([("A", "a_treedata.svg"), ("B", "b_plot.svg"), ("C", "c_functions.svg")]):
    x = M + i * col
    fig.letter(letter, x, M + 4)
    _, h = fig.place(PLOTS / "figs1" / svg, x + GUTTER, M, width=col - GUTTER - 2)
    h_max = max(h_max, h)
fig.h = M + h_max + M
fig.save(FIGURES / "figureS1")
print("wrote", FIGURES / "figureS1.pdf")

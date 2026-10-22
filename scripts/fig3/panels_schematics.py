"""Figure 3 schematics (toy tree drawn rotated: root at the top, leaves down), the pairwise linkage heatmap and the
heritability waterfall. Light: no embryo data is loaded, run directly.

Needs data/fig3_palettes.json (panels_data.py), data/fig3_linkage_uns.pkl (prepare_linkage.py) and
data/fig3_moran.csv (prepare_moran.py). Writes SVGs to plots/fig3/, each exactly its box in layout.py:

d_schematic  ancestral linkage of each cell to target cell type B: LCA time with the nearest B cell, normalized
e_schematic  pairwise linkage: linkage to each cell type, averaged per cell type
e_heatmap    pairwise linkage between lineages (embryo data)
f_schematic  Moran's I of a heritable and a random gene
f_waterfall  ranked heritability (Moran's I) of the 5,000 highly variable genes

All toy values are computed with pycea on a simulated 12-cell tree (cell types A, B, C in leaf order); all tree
schematics have identical tree axes (TREE_W x TREE_H mm).
"""

import contextlib
import io
import json
import pickle
import sys
import warnings
from pathlib import Path

import cassiopeia as cas
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pycea as py
import treedata as td
from matplotlib.patches import FancyArrowPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import DATA, mm_axes, mm_figure, mm_text, save_panel_fixed, set_style
from config import colors
from layout import BOX, TITLE_X, TREE_H, TREE_W

warnings.filterwarnings("ignore")
SUBDIR = "fig3"
set_style()
plt.rcParams["savefig.dpi"] = 600
INK = "#444444"
LIGHT = "#E3E3E3"
TITLE, LABEL, SMALL = 6.5, 6, 5.5  # font sizes (pt)
LW = 0.7  # tree line width (pt)
DOT = 15  # leaf marker area (pt^2)
TYPES = dict(zip("ABC", colors[:3]))  # toy cell types, in the project palette
LINK_CMAP, TIME_CMAP = "RdBu_r", "Purples"
ROW_H, ROW_GAP = 2.6, 0.6  # annotation rows below the trees (mm)

palettes = json.load(open(DATA / "fig3_palettes.json"))


def quiet(f, *args, **kwargs):
    with contextlib.redirect_stderr(io.StringIO()):  # permutation progress bars
        return f(*args, **kwargs)


# ---------------- toy tree ----------------
toy = cas.sim.birth_death_process(num_extant=12, random_seed=0, alignment="leaves")
g = toy.obst["simulated"]
_fig, _ax = plt.subplots()
py.pl.branches(toy, depth_key="time", ax=_ax, legend=False)
POS = {node: (pos, depth) for (_, node), (depth, pos) in _ax._attrs["node_coords"].items()}  # node -> (x, y) rotated
plt.close(_fig)
ORDERED = sorted(toy.obs_names, key=lambda leaf: POS[leaf][0])
DX = np.diff([POS[leaf][0] for leaf in ORDERED]).mean()
XLIM = (POS[ORDERED[0]][0] - DX / 2, POS[ORDERED[-1]][0] + DX / 2)  # one cell of width DX per leaf
DMAX = max(depth for _, depth in POS.values())
ROOT = next(n for n, d in g.in_degree() if d == 0)
# Cell types in leaf order: B is one pure clade and A and C are intermixed on both sides of it, so the cells left of the
# B clade share their LCA with the nearest B cell high in the tree (room for the "LCA" label in D). The pairwise matrix
# then has clear structure under pycea's default null: A-C linked, B unlinked from both.
TYPE_OF = dict(zip(ORDERED, "ACAACBBBBCAC", strict=True))
GROUPS = {k: [leaf for leaf in ORDERED if TYPE_OF[leaf] == k] for k in TYPES}
toy.obs["type"] = pd.Categorical([TYPE_OF[c] for c in toy.obs_names], categories=list(TYPES))

# ---------------- toy ancestral linkage (pycea) ----------------
raw, link = {}, {}
for k in TYPES:
    quiet(py.tl.ancestral_linkage, toy, groupby="type", target=k, normalize=False, key_added=f"raw_{k}", depth_key="time")
    quiet(py.tl.ancestral_linkage, toy, groupby="type", target=k, key_added=f"to_{k}", test="permutation",
          n_permutations=2000, random_state=0, depth_key="time")  # fmt: skip
    raw[k], link[k] = toy.obs[f"raw_{k}_linkage"], toy.obs[f"to_{k}_linkage"]
quiet(py.tl.ancestral_linkage, toy, groupby="type", symmetrize=False, test="permutation", n_permutations=2000,
      random_state=0, depth_key="time")  # fmt: skip
pairwise = toy.uns["type_linkage"]  # source rows x target columns
# The schematic matrix is the mean of the drawn per-cell rows; it equals pycea's pairwise matrix up to permutation noise
mean_link = pd.DataFrame({t: [link[t][GROUPS[s]].mean() for s in TYPES] for t in TYPES}, index=list(TYPES))
assert np.allclose(mean_link.to_numpy(), pairwise.loc[list(TYPES), list(TYPES)].to_numpy(), atol=0.05), (mean_link, pairwise)
LINK_LIM = max(np.abs(v).max() for v in link.values())


# ---------------- drawing helpers ----------------
def tree_axes(fig, x, y, highlight=None):
    """Rotated toy tree in TREE_W x TREE_H mm axes, top-left at (x, y) mm; leaves sit on the bottom edge."""
    ax = mm_axes(fig, x, y, TREE_W, TREE_H)
    for parent, child in g.edges:
        (xp, yp), (xc, yc) = POS[parent], POS[child]
        ax.plot([xp, xc, xc], [yp, yp, yc], color=INK, lw=LW, solid_joinstyle="miter", solid_capstyle="butt")
    ax.plot([POS[ROOT][0]] * 2, [-0.04 * DMAX, 0], color=INK, lw=LW, solid_capstyle="butt")
    ax.set_xlim(*XLIM)
    ax.set_ylim(DMAX, -0.04 * DMAX)
    ax.set_axis_off()
    return ax


def draw_leaves(ax, color_of, size=DOT):
    for leaf in ORDERED:
        ax.scatter(*POS[leaf], s=size, color=color_of[leaf], zorder=4, linewidths=0.3, edgecolors=INK, clip_on=False)


def highlight_path(ax, a, b, color, lw=1.5):
    """Draw the path from node a up to its ancestor b."""
    node = a
    while node != b:
        parent = next(g.predecessors(node))
        (xp, yp), (xc, yc) = POS[parent], POS[node]
        ax.plot([xp, xc, xc], [yp, yp, yc], color=color, lw=lw, solid_capstyle="round", solid_joinstyle="round", zorder=3)
        node = parent


def row_axes(fig, x, y, width=TREE_W, height=ROW_H):
    ax = mm_axes(fig, x, y, width, height)
    ax.set_xlim(*XLIM)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    return ax


def cell_row(fig, x, y, colors_of, header=None):
    """One cell per leaf, aligned with the leaves of a tree drawn at x."""
    ax = row_axes(fig, x, y)
    for leaf in ORDERED:
        ax.add_patch(Rectangle((POS[leaf][0] - DX / 2, 0), DX, 1, color=colors_of[leaf], lw=0))
    if header:
        mm_text(fig, x - 0.8, y + ROW_H / 2, header, ha="right", va="center", fontsize=SMALL)
    return ax


def cmap_colors(values, cmap, vmin, vmax):
    norm = plt.Normalize(vmin, vmax)
    return {leaf: plt.get_cmap(cmap)(norm(values[leaf])) for leaf in ORDERED}


def type_dot(fig, x, y, k, size=DOT):
    """A colored dot for cell type k with its letter to the right; (x, y) mm is the dot center."""
    fw, fh = fig.get_size_inches() * 25.4
    fig.add_artist(plt.Line2D([x / fw], [1 - y / fh], marker="o", ls="", markersize=np.sqrt(size), color=TYPES[k],
                              markeredgecolor=INK, markeredgewidth=0.3, transform=fig.transFigure))  # fmt: skip
    mm_text(fig, x + 1.0, y, k, ha="left", va="center", fontsize=SMALL)


def colorbar_key(fig, x, y, width, cmap, left, right, header, header_x):
    """Small horizontal color key over the full colormap, `left`/`right` labels at its ends, header right-aligned at
    `header_x` (like the row headers)."""
    ax = mm_axes(fig, x, y, width, 1.3)
    ax.imshow(np.linspace(0, 1, 64)[None, :], cmap=cmap, aspect="auto")
    ax.set_axis_off()
    mm_text(fig, header_x, y + 0.65, header, ha="right", va="center", fontsize=SMALL)
    mm_text(fig, x - 0.6, y + 0.65, left, ha="right", va="center", fontsize=5, color=INK)
    mm_text(fig, x + width + 0.6, y + 0.65, right, ha="left", va="center", fontsize=5, color=INK)


def down_arrow(fig, x, y0, y1, lw=0.6):
    fw, fh = fig.get_size_inches() * 25.4
    fig.add_artist(FancyArrowPatch((x / fw, 1 - y0 / fh), (x / fw, 1 - y1 / fh), transform=fig.transFigure, arrowstyle="-|>",
                                   mutation_scale=5, color=INK, lw=lw, shrinkA=0, shrinkB=0))  # fmt: skip


type_colors = {leaf: TYPES[TYPE_OF[leaf]] for leaf in ORDERED}

# ---------------- D: ancestral linkage to target type B ----------------
w, h = BOX["d_schematic"]
fig = mm_figure(w, h)
tx = w - TREE_W - 0.2  # tree x (row headers to its left)
mm_text(fig, TITLE_X, 2.6, "Ancestral linkage", fontsize=TITLE, va="baseline")
for i, k in enumerate(TYPES):
    type_dot(fig, tx + 0.8 + 5 * i, 6.3, k)
mm_text(fig, tx + 16, 6.3, "target: B", fontsize=SMALL, va="center", color=INK)
ty = 9.0
ax = tree_axes(fig, tx, ty)
target = "B"
source = ORDERED[4]  # the C cell just left of the B clade: its LCA with the nearest B cell is high in the tree
assert TYPE_OF[source] == "C"
# nearest target cell: most recent LCA, ties to the cell closest in leaf order
nearest = max(GROUPS[target], key=lambda c: (g.nodes[nx.lowest_common_ancestor(g, source, c)]["time"], -abs(POS[c][0] - POS[source][0])))
lca = nx.lowest_common_ancestor(g, source, nearest)
highlight_path(ax, source, lca, colors[5])
highlight_path(ax, nearest, lca, colors[5])
draw_leaves(ax, type_colors)
ax.scatter(*POS[source], s=42, facecolor="none", edgecolor="black", lw=0.6, zorder=5, clip_on=False)
ax.scatter(*POS[lca], s=38, marker="*", color=colors[5], edgecolor="black", lw=0.35, zorder=6)
ax.annotate("LCA", POS[lca], xytext=(POS[lca][0] - 1.6 * DX, POS[lca][1] - 0.16 * DMAX), fontsize=SMALL, ha="right", va="center",
            arrowprops={"arrowstyle": "-", "lw": 0.4, "color": INK, "shrinkA": 0.5, "shrinkB": 2.5})  # fmt: skip
y = ty + TREE_H + 2.0
target_colors = {leaf: TYPES[target] if TYPE_OF[leaf] == target else LIGHT for leaf in ORDERED}
cell_row(fig, tx, y, target_colors, "Target")
cell_row(fig, tx, y + ROW_H + ROW_GAP, cmap_colors(raw[target], TIME_CMAP, 0, DMAX), "LCA time")
cell_row(fig, tx, y + 2 * (ROW_H + ROW_GAP), cmap_colors(link[target], LINK_CMAP, -LINK_LIM, LINK_LIM), "Linkage")
y += 3 * ROW_H + 2 * ROW_GAP + 4.8
colorbar_key(fig, tx + 7.0, y, 11, TIME_CMAP, "early", "recent", "LCA time", tx - 0.8)
colorbar_key(fig, tx + 7.0, y + 4.2, 11, LINK_CMAP, "\u2212", "+", "Linkage", tx - 0.8)
mm_text(fig, tx + TREE_W / 2 - 3, y + 10.8, "Linkage = LCA time − permuted mean", ha="center", va="center", fontsize=SMALL)
save_panel_fixed(fig, "d_schematic", SUBDIR)

# ---------------- E: pairwise linkage = mean linkage per cell type, for every target type ----------------
def row_group_label(fig, label, y0, y1, row_centers, tx):
    """Rotated `label` with a bracket spanning rows y0..y1 (mm), and the target letters A, B, C next to the rows."""
    fw, fh = fig.get_size_inches() * 25.4
    xb, tick = 3.1, 0.5  # bracket x and end-tick length (mm)
    for xs, ys in [([xb, xb], [y0, y1]), ([xb, xb + tick], [y0, y0]), ([xb, xb + tick], [y1, y1])]:
        fig.add_artist(plt.Line2D([x / fw for x in xs], [1 - y / fh for y in ys], color=INK, lw=0.4, transform=fig.transFigure))
    mm_text(fig, 1.5, (y0 + y1) / 2, label, rotation=90, ha="center", va="center", fontsize=SMALL)
    for k, yc in zip(TYPES, row_centers, strict=True):
        mm_text(fig, tx - 0.7, yc, k, ha="right", va="center", fontsize=SMALL)


w, h = BOX["e_schematic"]
fig = mm_figure(w, h)
tx = 5.6  # tree x: the rotated row labels are to its left, so the schematic keeps clear of the heatmap
mm_text(fig, TITLE_X, 2.6, "Pairwise linkage", fontsize=TITLE, va="baseline")
ty = 5.0
ax = tree_axes(fig, tx, ty)
draw_leaves(ax, type_colors)
y_rows = ty + TREE_H + 2.0
for i, k in enumerate(TYPES):
    cell_row(fig, tx, y_rows + i * (ROW_H + ROW_GAP), cmap_colors(link[k], LINK_CMAP, -LINK_LIM, LINK_LIM))
rows_h = 3 * ROW_H + 2 * ROW_GAP
row_group_label(fig, "Linkage to", y_rows, y_rows + rows_h, [y_rows + i * (ROW_H + ROW_GAP) + ROW_H / 2 for i in range(3)], tx)
# mean over the cells of each (source) cell type -> one matrix column per cell type
y_arrow = y_rows + rows_h + 1.0
down_arrow(fig, tx + TREE_W / 2, y_arrow, y_arrow + 4.6)
mm_text(fig, tx + TREE_W / 2 + 1.2, y_arrow + 2.3, "mean per cell type", ha="left", va="center", fontsize=SMALL, color=INK)
col_w = TREE_W / 3
col_center = {k: tx + (i + 0.5) * col_w for i, k in enumerate(TYPES)}
y_mat = y_arrow + 5.4
cell_h, cell_gap = 4.2, 0.5
mat_h = 3 * cell_h + 2 * cell_gap
norm = plt.Normalize(-LINK_LIM, LINK_LIM)
mat_ax = mm_axes(fig, tx, y_mat, TREE_W, mat_h)
mat_ax.set_xlim(0, TREE_W)
mat_ax.set_ylim(mat_h, 0)
mat_ax.set_axis_off()
for r, tgt in enumerate(TYPES):
    for c, src in enumerate(TYPES):
        mat_ax.add_patch(Rectangle((c * col_w + cell_gap / 2, r * (cell_h + cell_gap)), col_w - cell_gap, cell_h, lw=0,
                                   color=plt.get_cmap(LINK_CMAP)(norm(mean_link.loc[src, tgt]))))  # fmt: skip
row_group_label(fig, "Mean linkage to", y_mat, y_mat + mat_h, [y_mat + r * (cell_h + cell_gap) + cell_h / 2 for r in range(3)], tx)
y_src = y_mat + mat_h + 1.8
for k in TYPES:
    type_dot(fig, col_center[k] - 0.5, y_src, k, size=10)
mm_text(fig, tx + TREE_W / 2, y_src + 2.8, "source cell type", ha="center", va="center", fontsize=SMALL, color=INK)
save_panel_fixed(fig, "e_schematic", SUBDIR)

# ---------------- E heatmap: pairwise linkage between lineages (embryo) ----------------
w, h = BOX["e_heatmap"]
with open(DATA / "fig3_linkage_uns.pkl", "rb") as f:
    uns = pickle.load(f)
dummy = td.TreeData(obs=pd.DataFrame(index=["cell"]), uns=uns)
fig = mm_figure(w, h)
side, block = 35.0, 1.6  # heatmap side; width of the lineage color blocks
x0, y0 = 0.3 + block + 0.3, 0.4
ax = mm_axes(fig, x0, y0, side, side)
py.pl.ancestral_linkage(dummy, groupby="lineage", yticklabels=False, cbar=False, ax=ax)
order = [t.get_text() for t in ax.get_xticklabels()]  # rows/columns in pycea's clustered order
n = len(order)
ax.set_xticks([])
for spine in ax.spines.values():
    spine.set_visible(False)
left = mm_axes(fig, 0.3, y0, block, side)  # lineage colors of the rows
bottom = mm_axes(fig, x0, y0 + side + 0.3, side, block)  # ... and of the columns
for i, lin in enumerate(order):
    left.add_patch(Rectangle((0, i), 1, 1, color=palettes["lineage"][lin], lw=0))
    bottom.add_patch(Rectangle((i, 0), 1, 1, color=palettes["lineage"][lin], lw=0))
left.set_xlim(0, 1)
left.set_ylim(n, 0)
bottom.set_xlim(0, n)
bottom.set_ylim(0, 1)
bottom.set_xticks(np.arange(n) + 0.5)
bottom.set_xticklabels(order, rotation=90, fontsize=SMALL)
bottom.tick_params(length=0, pad=1.5)
bottom.set_yticks([])
for a in (left, bottom):
    for spine in a.spines.values():
        spine.set_visible(False)
left.set_axis_off()
cax = mm_axes(fig, x0 + side + 1.2, y0 + 6, 1.4, side - 12)
cbar = fig.colorbar(ax.collections[0], cax=cax)
cbar.set_label("Ancestral linkage", fontsize=LABEL, labelpad=2)
cbar.ax.tick_params(labelsize=SMALL, length=1.5, width=0.4, pad=1)
cbar.outline.set_linewidth(0.4)
save_panel_fixed(fig, "e_heatmap", SUBDIR)

# ---------------- F: Moran's I of a heritable and a random gene (toy values; Moran's I by pycea) ----------------
quiet(py.tl.tree_neighbors, toy, n_neighbors=3, depth_key="time", random_state=0)
rng = np.random.default_rng(2)  # heritable gene: Brownian motion along the tree
bm = {}
for node in nx.topological_sort(g):
    parent = next(g.predecessors(node), None)
    bm[node] = 0.0 if parent is None else bm[parent] + rng.normal(0, np.sqrt(g.nodes[node]["time"] - g.nodes[parent]["time"]))
heritable = pd.Series({leaf: bm[leaf] for leaf in ORDERED})
heritable = (heritable - heritable.min()) / (heritable.max() - heritable.min())
best = None
for seed in range(200):  # a random gene with Moran's I close to 0
    candidate = pd.Series(np.random.default_rng(seed).random(len(ORDERED)), index=ORDERED)
    toy.obs["random"] = candidate.reindex(toy.obs_names)
    toy.obs["heritable"] = heritable.reindex(toy.obs_names)
    I = quiet(py.tl.autocorr, toy, keys=["heritable", "random"], copy=True)["autocorr"]
    if best is None or abs(I["random"]) < abs(best[1]["random"]):
        best = (candidate, I)
random_gene, moran = best
w, h = BOX["f_schematic"]
fig = mm_figure(w, h)
tx = TITLE_X
genes = [("Heritable gene", heritable, moran["heritable"]), ("Random gene", random_gene, moran["random"])]
tree_tops = [5.6, 33.0]
for (title, values, I), ty in zip(genes, tree_tops, strict=True):
    mm_text(fig, tx, ty - 2.0, f"{title}: Moran's I = {round(I, 2) + 0.0:.2f}", fontsize=SMALL, va="baseline")
    ax = tree_axes(fig, tx, ty)
    draw_leaves(ax, {leaf: plt.get_cmap("viridis")(values[leaf]) for leaf in ORDERED}, size=DOT + 4)
cax = mm_axes(fig, tx + TREE_W + 3.0, tree_tops[0] + 8, 1.4, tree_tops[1] + TREE_H - tree_tops[0] - 16)
cb = fig.colorbar(plt.cm.ScalarMappable(cmap="viridis"), cax=cax, ticks=[])
cb.outline.set_linewidth(0.4)
cb.set_label("Expression", fontsize=SMALL, labelpad=2)
mm_text(fig, tx + TREE_W + 3.7, tree_tops[0] + 7.2, "high", ha="center", va="baseline", fontsize=5, color=INK)
mm_text(fig, tx + TREE_W + 3.7, tree_tops[1] + TREE_H - 8 + 2.3, "low", ha="center", va="baseline", fontsize=5, color=INK)
save_panel_fixed(fig, "f_schematic", SUBDIR)

# ---------------- F waterfall: ranked heritability of the HVGs ----------------
table = pd.read_csv(DATA / "fig3_moran.csv", index_col="gene")["autocorr"].sort_values(ascending=False)
w, h = BOX["f_waterfall"]
fig = mm_figure(w, h)
ax = mm_axes(fig, 9.5, 2.0, w - 13.0, h - 10.5)
ranks = np.arange(len(table))
ax.bar(ranks, table.to_numpy(), width=1.0, color="#B0B0B0", linewidth=0)
ax.bar(ranks[:4], table.to_numpy()[:4], width=40, color=colors[1], linewidth=0)
for i, gene in enumerate(table.index[:4]):
    ax.annotate(gene, (ranks[i] + 20, table.iloc[i]), xytext=(700, table.iloc[0] - 0.08 * (i + 0.2)), fontsize=SMALL, style="italic",
                va="center", arrowprops={"arrowstyle": "-", "lw": 0.35, "color": INK, "shrinkA": 1, "shrinkB": 0.5})  # fmt: skip
ax.set_xlim(-60, len(table) + 60)
ax.set_ylim(0, None)
ax.set_xlabel("Genes (ranked)", fontsize=LABEL, labelpad=1.5)
ax.set_ylabel("Heritability (Moran's I)", fontsize=LABEL, labelpad=2)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(length=2, width=0.5, pad=1.5)
save_panel_fixed(fig, "f_waterfall", SUBDIR)
print("toy Moran's I:", moran.round(3).to_dict(), flush=True)
print("done", flush=True)

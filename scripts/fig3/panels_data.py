"""Figure 3 panels from the full embryo (382k cells): the polar tree and the UMAPs.

The embryo is only loaded (Slurm) to redraw the polar tree (`--tree`) or to (re)build the per-cell cache (`--cache`, or
when the cache is missing); otherwise the UMAPs are drawn from cached tables in about a minute (`--only b c d g`
restricts the panels).

Needs data/fig3_obs.csv (prepare_linkage.py) and data/fig3_top_expr.csv (prepare_moran.py). Writes:

data/fig3_cells.csv.gz   per-cell UMAP coordinates, germ layer and lineage (cache)
data/fig3_palettes.json  germ layer / lineage colors (also used by panels_schematics.py)
plots/fig3/a_tree        polar tree, leaves annotated with germ layer and lineage (only with --tree; slow; 1200 dpi)
plots/fig3/b_umap_germ_layer, c_umap_lineage   UMAPs colored by germ layer / lineage, legend below
plots/fig3/d_umap        UMAP colored by neural crest ancestral linkage (scanpy, sort_order=False)
plots/fig3/g_umaps       2x2 UMAPs colored by the expression of the four most heritable genes

The UMAP panels are drawn at their final size (boxes in layout.py); all UMAP axes have the aspect ratio of the data, with
equal x and y scales.
"""

import json
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad
import pycea as py
import scanpy as sc
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import DATA, mm_axes, mm_figure, mm_text, save_panel_full, save_panel_fixed, set_style
from layout import BOX, UMAP_ASPECT, UMAP_W

warnings.filterwarnings("ignore")
SUBDIR = "fig3"
CELLS = DATA / "fig3_cells.csv.gz"
SEED = 0
LABEL, SMALL = 6, 5.5  # font sizes (pt)
set_style()
plt.rcParams["savefig.dpi"] = 600  # resolution of rasterized elements
ONLY = set(sys.argv[sys.argv.index("--only") + 1 :]) if "--only" in sys.argv else None  # e.g. --only g

if "--tree" in sys.argv or "--cache" in sys.argv or not CELLS.exists():
    tdata = py.datasets.colgan26(embryos="E9.5-R1")
    umap = tdata.obsm["X_umap"]
    cells = pd.DataFrame(
        {"umap1": umap[:, 0], "umap2": umap[:, 1], "germ_layer": tdata.obs["germ_layer"], "lineage": tdata.obs["lineage"]},
        index=tdata.obs_names,
    )
    cells.to_csv(CELLS)
    palettes = {"germ_layer": py.get.palette(tdata, key="germ_layer"), "lineage": py.get.palette(tdata, key="lineage")}
    with open(DATA / "fig3_palettes.json", "w") as f:
        json.dump(palettes, f)
    print("cached", len(cells), "cells; umap ranges", umap.min(0), umap.max(0), flush=True)
    # --- A: polar tree with germ layer and lineage rings, drawn at its final size; rasters at 1200 dpi ---
    if "--tree" in sys.argv:
        side = BOX["a_tree"][0] / 25.4
        fig, ax = plt.subplots(figsize=(side, side), subplot_kw={"projection": "polar"})
        py.pl.branches(tdata, depth_key="time", polar=True, linewidth=0.05, color="#555555", ax=ax)
        py.pl.annotation(tdata, keys=["germ_layer", "lineage"], width=0.08, gap=0.02, label=False, legend=False, ax=ax)
        for c in ax.collections:
            c.set_rasterized(True)
        plt.rcParams["savefig.dpi"] = 1200
        save_panel_full(fig, "a_tree", SUBDIR)
        plt.rcParams["savefig.dpi"] = 600
    del tdata


def want(name):
    return ONLY is None or name in ONLY


cells = pd.read_csv(CELLS, index_col=0)
cells["nc_linkage"] = pd.read_csv(DATA / "fig3_obs.csv", index_col=0)["nc_linkage"].reindex(cells.index)
expr = pd.read_csv(DATA / "fig3_top_expr.csv", index_col=0).reindex(cells.index)
genes = list(expr.columns)
palettes = json.load(open(DATA / "fig3_palettes.json"))
umap = cells[["umap1", "umap2"]].to_numpy()
lo, hi = umap.min(0), umap.max(0)
pad = 0.015 * (hi - lo)
XLIM, YLIM = (lo[0] - pad[0], hi[0] + pad[0]), (lo[1] - pad[1], hi[1] + pad[1])  # same aspect as UMAP_ASPECT


def umap_axes(fig, x, y, width):
    """UMAP axes `width` mm wide (height from the data aspect ratio), top-left at (x, y) mm."""
    ax = mm_axes(fig, x, y, width, width / UMAP_ASPECT)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_axis_off()
    return ax


def colorbar(fig, mappable, rect, label, orientation, ticks=None):
    cax = mm_axes(fig, *rect)
    cbar = fig.colorbar(mappable, cax=cax, orientation=orientation, ticks=ticks)
    cbar.set_label(label, fontsize=LABEL, labelpad=2)
    cbar.ax.tick_params(labelsize=SMALL, length=1.5, width=0.4, pad=1)
    cbar.outline.set_linewidth(0.4)
    return cbar


def category_umap(color, name, legend):
    """UMAP colored by a category; legend "below" (two columns) or "right" (one column)."""
    palette = palettes[color]
    w, h = BOX[name]
    fig = mm_figure(w, h)
    ax = umap_axes(fig, 0.3, 0.3, UMAP_W)
    values = cells[color].astype(str).map(palette).fillna("#D3D3D3").to_numpy()
    order = np.random.default_rng(SEED).permutation(len(values))
    ax.scatter(umap[order, 0], umap[order, 1], c=values[order], s=0.3, linewidths=0, rasterized=True)
    cats = [c for c in cells[color].value_counts().index if c in palette]
    handles = [Line2D([], [], marker="o", ls="", markersize=3.2, color=palette[k], label=k) for k in cats]
    where = {"below": {"loc": "upper left", "bbox_to_anchor": (0, -0.012), "ncol": 2},
             "right": {"loc": "center left", "bbox_to_anchor": (1.0, 0.5), "ncol": 1}}[legend]  # fmt: skip
    leg = ax.legend(handles=handles, frameon=False, fontsize=SMALL, handlelength=0.8, handletextpad=0.3, labelspacing=0.32,
                    columnspacing=0.8, borderaxespad=0, borderpad=0, **where)  # fmt: skip
    fig.canvas.draw()
    box = leg.get_window_extent().transformed(fig.transFigure.inverted())
    if box.x1 > 1 or box.y0 < 0:
        print(f"warning: {name} legend extends past the panel box ({box})", flush=True)
    save_panel_fixed(fig, name, SUBDIR)


# --- B, C ---
if want("b"):
    category_umap("germ_layer", "b_umap_germ_layer", "below")
if want("c"):
    category_umap("lineage", "c_umap_lineage", "right")

# --- D: neural crest ancestral linkage (py.tl.ancestral_linkage(tdata, groupby="lineage", target="Neural crest",
#     key_added="nc") in prepare_linkage.py), drawn by Scanpy in the original cell order (sort_order=False). Cells
#     without a value (not in a tree; light gray) come last in the data, so they are moved first to be drawn below ---
if want("d"):
    w, h = BOX["d_umap"]
    fig = mm_figure(w, h)
    umap_w = w - 0.5
    ax = umap_axes(fig, 0.25, 0.3, umap_w)
    first = np.argsort(cells["nc_linkage"].notna().to_numpy(), kind="stable")  # missing values first, order kept
    adata = ad.AnnData(obs=cells[["nc_linkage"]].iloc[first], obsm={"X_umap": umap[first]})
    sc.pl.umap(adata, color="nc_linkage", sort_order=False, vmin=-2, vmax=2, cmap="RdBu_r", s=0.3, frameon=False,
               title="", colorbar_loc=None, ax=ax, show=False)  # fmt: skip
    for c in ax.collections:
        c.set_rasterized(True)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xlabel("")
    ax.set_ylabel("")
    y_cbar = 0.3 + umap_w / UMAP_ASPECT + 2.0
    colorbar(fig, ax.collections[-1], (w / 2 - 14, y_cbar, 28, 1.4), "Neural crest ancestral linkage", "horizontal", ticks=[-2, -1, 0, 1, 2])
    save_panel_fixed(fig, "d_umap", SUBDIR)

# --- G: 2x2 UMAPs of the most heritable genes (names in the empty top-left corner), shared scale, colorbar right ---
if want("g"):
    w, h = BOX["g_umaps"]
    fig = mm_figure(w, h)
    gap_x, gap_y = 1.0, 1.2
    side_h = (h - gap_y) / 2
    side_w = side_h * UMAP_ASPECT
    vmax = max(np.nanpercentile(expr[gene], 99.5) for gene in genes)
    for k, gene in enumerate(genes):
        row, col = divmod(k, 2)
        x, y = col * (side_w + gap_x), row * (side_h + gap_y)
        ax = umap_axes(fig, x, y, side_w)
        values = expr[gene].to_numpy()
        order = np.argsort(np.nan_to_num(values))
        mappable = ax.scatter(umap[order, 0], umap[order, 1], c=values[order], s=0.15, cmap="viridis", vmin=0, vmax=vmax,
                              linewidths=0, rasterized=True)  # fmt: skip
        mm_text(fig, x + 0.3, y + 2.4, gene, ha="left", va="baseline", fontsize=LABEL, style="italic")
    colorbar(fig, mappable, (2 * side_w + gap_x + 1.5, h / 2 - 12, 1.4, 24), "Expression (log1p)", "vertical")
    save_panel_fixed(fig, "g_umaps", SUBDIR)
print("done", flush=True)

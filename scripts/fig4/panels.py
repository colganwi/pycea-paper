"""Figure 4 panels: growth dynamics of 4T1 tumors (Koblan et al. 2025; follows the growth-dynamics tutorial).

Needs data/fig4_obs.csv and data/fig4_gene_corr.csv from prepare.py. Writes SVGs to plots/fig4/, each drawn at the
final size of its box in layout.py (fonts are the printed sizes):

a_tree_characters  tree + character matrix, branches colored by clade
b_spatial_clade    spatial coordinates colored by clade
c_n_extant         number of extant cells per clade over time
d_tree_fitness     polar tree (clades) with a fitness ring
e_spatial_fitness  spatial coordinates colored by fitness
f_waterfall        correlation of gene expression with fitness
g_tree_cldn4       polar tree with fitness and Cldn4 rings
h_spatial_cldn4    spatial coordinates colored by Cldn4 expression
"""

import sys
import warnings
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import pycea as py
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import DATA, mm_axes, mm_figure, save_panel_fixed, set_style
from config import colors
from layout import BOX, SPATIAL_ASPECT, SPATIAL_W

warnings.filterwarnings("ignore")
SUBDIR = "fig4"
LABEL, SMALL = 6, 5.5  # font sizes (pt), as in Figure 3
set_style()
plt.rcParams.update({"font.size": LABEL, "axes.labelsize": LABEL, "xtick.labelsize": SMALL, "ytick.labelsize": SMALL,
                     "legend.fontsize": SMALL, "xtick.major.size": 2, "ytick.major.size": 2, "xtick.major.pad": 1.5,
                     "ytick.major.pad": 1.5, "axes.labelpad": 2})  # fmt: skip
plt.rcParams["savefig.dpi"] = 600  # resolution of rasterized elements (cells, leaf annotations) in the SVGs
edit_cmap = mcolors.ListedColormap(
    ["white", "lightgray", "#CD2626", "#E69F00", "#FFE600", "#009E73", "#83A4FF", "#1874CD", "#8E0496", "#DB65D2"]
)

tdata = py.datasets.koblan25()
# Keep only the three tumor sections that have lineage data (drops the section without lineage-traced cells)
lineage_sections = tdata.obs.loc[tdata.obs["tree"].notna(), "sample"].unique()
tdata = tdata[tdata.obs["sample"].isin(lineage_sections)].copy()
# Clades: cells sharing a common ancestor at day 5 (as in the growth-dynamics tutorial)
py.tl.clades(tdata, depth=5, depth_key="time", update=False)
tdata.obs["fitness"] = pd.read_csv(DATA / "fig4_obs.csv", index_col=0)["fitness"].reindex(tdata.obs_names)  # prepare.py
clade_palette = py.get.palette(tdata, key="clade", cmap="rainbow")
xy = tdata.obsm["spatial"]
pad = 0.01 * np.ptp(xy, axis=0)
XLIM, YLIM = (xy[:, 0].min() - pad[0], xy[:, 0].max() + pad[0]), (xy[:, 1].min() - pad[1], xy[:, 1].max() + pad[1])


def rasterize_since(ax, n_before):
    """Rasterize artists added after the first `n_before` collections (leaf annotations have ~10k cells each)."""
    for collection in ax.collections[n_before:]:
        collection.set_rasterized(True)
    for child in ax.get_children():
        if child.__class__.__name__ in ("QuadMesh", "AxesImage"):
            child.set_rasterized(True)


def spatial(color, name, label, mappable=None, ticks=True, **kwargs):
    """Spatial scatter of the three sections (equal x/y scale) with a slim colorbar labeled `label` at its right."""
    w, h = BOX[name]
    fig = mm_figure(w, h)
    strip_h = SPATIAL_W / SPATIAL_ASPECT  # the 1% padding on both axes keeps the data aspect ratio
    ax = mm_axes(fig, 0, (h - strip_h) / 2, SPATIAL_W, strip_h)
    sc.pl.spatial(
        tdata, color=color, spot_size=40, frameon=False, ax=ax, show=False, title="", colorbar_loc=None, **kwargs
    )
    inverted = ax.get_ylim()[0] > ax.get_ylim()[1]
    ax.set_xlim(*XLIM)
    ax.set_ylim(*(YLIM[::-1] if inverted else YLIM))
    ax.set_aspect("auto")
    for collection in ax.collections:
        collection.set_rasterized(True)  # keep the SVG small: many points
    if mappable is None:
        mappable = [c for c in ax.collections if c.get_array() is not None][-1]
    cax = mm_axes(fig, SPATIAL_W + 1.5, h / 2 - 12, 1.4, 24)
    cbar = fig.colorbar(mappable, cax=cax)
    cbar.set_label(label, fontsize=LABEL, labelpad=2)
    cbar.ax.tick_params(labelsize=SMALL, length=1.5, width=0.4, pad=1)
    cbar.outline.set_linewidth(0.4)
    if ticks:
        cbar.locator = plt.MaxNLocator(4, integer=True)
        cbar.update_ticks()
    save_panel_fixed(fig, name, SUBDIR)


def polar_axes(name):
    w, h = BOX[name]
    fig = mm_figure(w, h)
    return fig, mm_axes(fig, 0, (h - w) / 2, w, w, projection="polar")


def plot_axes(name):
    """Axes for the line/bar plots C and F (room for tick labels and axis labels on the left and below)."""
    w, h = BOX[name]
    fig = mm_figure(w, h)
    return fig, mm_axes(fig, 9.5, 1.5, w - 10.5, h - 9.5)


# --- A: tree + character matrix colored by clade ---
w, h = BOX["a_tree_characters"]
fig = mm_figure(w, h)
ax = mm_axes(fig, 0, 0.5, 64.0, h - 1.0)
py.pl.tree(tdata, depth_key="time", branch_color="clade", palette=clade_palette, legend=False, ax=ax)
n = len(ax.collections)
py.pl.annotation(tdata, keys="characters", cmap=edit_cmap, label=False, legend=False, ax=ax)
rasterize_since(ax, n)
# Legend for the character states (cmap: ND white, UE light gray, then the 8 edited states LM1-LM8)
states = ["LM%d" % i for i in range(1, 9)] + ["UE", "ND"]
state_colors = list(edit_cmap.colors[2:]) + [edit_cmap.colors[1], edit_cmap.colors[0]]
handles = [Patch(facecolor=c, edgecolor="black", linewidth=0.3, label=l) for l, c in zip(states, state_colors)]
ax.legend(
    handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=SMALL, handlelength=1,
    handleheight=1, labelspacing=0.35, borderaxespad=0,
)  # fmt: skip
save_panel_fixed(fig, "a_tree_characters", SUBDIR)

# --- B: spatial, clades ---
spatial(
    "clade", "b_spatial_clade", "Clades", ticks=True, legend_loc=None,
    mappable=plt.cm.ScalarMappable(norm=mcolors.Normalize(0, tdata.obs["clade"].nunique()), cmap="rainbow"),
)  # fmt: skip

# --- C: number of extant cells per clade ---
py.tl.n_extant(tdata, depth_key="time", groupby="clade", bins=20, dropna=True)
fig, ax = plot_axes("c_n_extant")
py.pl.n_extant(tdata, legend=False, ax=ax)
ax.set_xlabel("Time (days)")
ax.set_ylabel("Number of extant lineages")
save_panel_fixed(fig, "c_n_extant", SUBDIR)

# --- D: polar tree with clades and a fitness ring ---
fig, ax = polar_axes("d_tree_fitness")
py.pl.tree(tdata, depth_key="time", branch_color="clade", palette=clade_palette, polar=True, legend=False, ax=ax)
n = len(ax.collections)
py.pl.annotation(tdata, keys="fitness", width=0.2, cmap="magma", legend=False, ax=ax)
rasterize_since(ax, n)
save_panel_fixed(fig, "d_tree_fitness", SUBDIR)

# --- E: spatial, fitness ---
spatial("fitness", "e_spatial_fitness", "Fitness", cmap="magma")

# --- F: waterfall of gene-fitness correlations ---
corr = pd.read_csv(DATA / "fig4_gene_corr.csv", index_col="gene")["correlation"].sort_values(ascending=False)
fig, ax = plot_axes("f_waterfall")
ranks = np.arange(len(corr))
ax.bar(ranks, corr.to_numpy(), width=1.0, color="#B0B0B0", linewidth=0)
ax.bar(ranks[corr.index == "Cldn4"], corr[corr.index == "Cldn4"], width=3.0, color=colors[1])
ax.annotate("Cldn4", (ranks[corr.index == "Cldn4"][0], corr["Cldn4"]), xytext=(8, 0), textcoords="offset points", fontsize=SMALL)
ax.axhline(0, color="black", linewidth=0.4)
ax.set_xlabel("Genes (ranked)")
ax.set_ylabel("Correlation with fitness")
ax.spines[["top", "right"]].set_visible(False)
save_panel_fixed(fig, "f_waterfall", SUBDIR)

# --- G: polar tree with fitness and Cldn4 rings ---
fig, ax = polar_axes("g_tree_cldn4")
py.pl.tree(tdata, depth_key="time", branch_color="clade", palette=clade_palette, polar=True, legend=False, ax=ax)
n = len(ax.collections)
py.pl.annotation(tdata, keys="fitness", width=0.15, cmap="magma", legend=False, ax=ax)
py.pl.annotation(tdata, keys="Cldn4", width=0.15, vmax=20, legend=False, ax=ax)
rasterize_since(ax, n)
save_panel_fixed(fig, "g_tree_cldn4", SUBDIR)

# --- H: spatial, Cldn4 expression ---
spatial("Cldn4", "h_spatial_cldn4", "Cldn4 expression", cmap="viridis", vmax=20)

"""Figure 1C: example workflow panels (SVGs in plots/fig1/).

c1_reconstruction: tree + character matrix for all cells
c2_umap:           UMAP colored by leiden clusters (Scanpy)
c3_tree:           circular tree with leaves annotated by cluster (Pycea)
c4_linkage:        ancestral linkage heatmap (Pycea)
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pycea as py
import scanpy as sc
import treedata as td

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import CLUSTER_PALETTE, DATA, save_panel, set_style

SUBDIR = "fig1"
SEED = 0
set_style()

# Load simulated data (200 cells, leaf-aligned tree)
tdata = td.read_h5td(DATA / "fig1_sim.h5td")

# --- Scanpy: UMAP and clustering with leiden clusters ---
sc.pp.pca(tdata, random_state=SEED)
sc.pp.neighbors(tdata, random_state=SEED)
sc.tl.umap(tdata, random_state=SEED)
sc.tl.leiden(tdata, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
clusters = sorted(tdata.obs["leiden"].unique(), key=int)
palette = dict(zip(clusters, CLUSTER_PALETTE))
tdata.uns["leiden_colors"] = [palette[c] for c in clusters]

# --- c1: tree + character matrix for all cells ---
# pycea requires a single shared categorical dtype across the columns of an obsm array
tdata_chars = tdata.copy()
chars = tdata_chars.obsm["characters"].astype(str)
states = sorted(set(chars.values.ravel()) - {"*"}, key=int) + ["*"]
tdata_chars.obsm["characters"] = chars.astype(pd.CategoricalDtype(states))
state_colors = {s: c for s, c in zip(states[:-1], sc.pl.palettes.default_102, strict=False)} | {"*": "#E8E8E8"}
fig, ax = plt.subplots(figsize=(1.8, 1.6))
py.pl.tree(tdata_chars, depth_key="time", branch_linewidth=0.6, ax=ax)
py.pl.annotation(
    tdata_chars,
    keys="characters",
    slot="obsm",
    width=0.25,
    gap=0.02,
    label=False,
    palette=state_colors,
    ax=ax,
    legend=False,
)
ax.set_axis_off()
save_panel(fig, "c1_reconstruction", SUBDIR)

# --- c2: UMAP ---
fig, ax = plt.subplots(figsize=(1.6, 1.6))
sc.pl.umap(tdata, color="leiden", ax=ax, show=False, legend_loc=None, title="", size=30, frameon=False)
save_panel(fig, "c2_umap", SUBDIR)

# --- c3: circular tree annotated by cluster ---
fig, ax = plt.subplots(figsize=(1.8, 1.8), subplot_kw={"projection": "polar"})
py.pl.tree(tdata, polar=True, depth_key="time", branch_linewidth=0.4, ax=ax)
py.pl.annotation(tdata, keys="leiden", width=0.08, gap=0.03, label=False, legend=False, ax=ax)
ax.set_axis_off()
save_panel(fig, "c3_tree", SUBDIR)

# --- c4: ancestral linkage heatmap ---
py.pp.add_depth(tdata)
py.tl.ancestral_linkage(tdata, groupby="leiden", test="permutation", n_permutations=100, random_state=SEED)
fig, ax = plt.subplots(figsize=(1.8, 1.6))
py.pl.ancestral_linkage(tdata, groupby="leiden", ax=ax)
save_panel(fig, "c4_linkage", SUBDIR)

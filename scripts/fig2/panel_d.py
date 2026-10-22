"""Figure 2D: a tree plot built up with successive Pycea calls (SVGs in plots/fig2/).

Uses the same simulated data as Figure 1 (200 leaf-aligned cells), non-polar layout. Each panel adds one call:
d1_branches:   py.pl.branches   - branches colored by branch length
d2_nodes:      py.pl.nodes      - leaves and internal nodes colored by cluster
d3_star:       py.pl.nodes      - a single node marked with a large black star
d4_annotation: py.pl.annotation - expression of five genes at the leaves (RdBu_r); legends for all features
"""

import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pycea as py
import scanpy as sc
import treedata as td

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import CLUSTER_PALETTE, DATA, save_panel, save_panel_full, set_style

warnings.filterwarnings("ignore")
SEED = 0
# Panels are drawn 30% larger than their final size so that legends (fixed font sizes) shrink when assembled;
# line widths and marker sizes are scaled up by the same factor so the tree itself looks unchanged.
SCALE = 1.3
FIGSIZE = (1.75 * SCALE, 2.4 * SCALE)  # tall enough that the legends of the last panel stack in one column
LW = 0.6 * SCALE
NODE = 3.5 * SCALE**2  # marker sizes are areas
STAR = 60 * SCALE**2
set_style()

# --- Same embedding/clustering as Figure 1 so cluster labels and colors match ---
tdata = td.read_h5td(DATA / "fig1_sim.h5td")
sc.pp.pca(tdata, random_state=SEED)
sc.pp.neighbors(tdata, random_state=SEED)
sc.tl.umap(tdata, random_state=SEED)
sc.tl.leiden(tdata, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
clusters = sorted(tdata.obs["leiden"].unique(), key=int)
palette = dict(zip(clusters, CLUSTER_PALETTE))

# Internal-node clusters: most common cluster among descendant leaves
py.tl.ancestral_states(tdata, "leiden", method="mode")

# Branch lengths (the simulated tree stores node times only)
tree = tdata.obst["simulated"]
nx.set_edge_attributes(tree, {(u, v): tree.nodes[v]["time"] - tree.nodes[u]["time"] for u, v in tree.edges}, "length")

# Node to highlight: the largest subtree whose leaves all belong to one cluster (a cluster founder)
leaf_set = {n for n, d in tree.out_degree() if d == 0}
best, best_size = None, 0
for node in tree.nodes:
    if node in leaf_set:
        continue
    sub_leaves = (nx.descendants(tree, node)) & leaf_set
    if len({tree.nodes[leaf]["leiden"] for leaf in sub_leaves}) == 1 and len(sub_leaves) > best_size:
        best, best_size = node, len(sub_leaves)

# Five genes with the highest variance across cells
X = np.asarray(tdata.X)
genes = list(tdata.var_names[np.argsort(X.var(axis=0))[::-1][:5]])


def new_axes():
    fig, ax = plt.subplots(figsize=FIGSIZE)
    return fig, ax


# --- d1: branches colored by length ---
fig, ax = new_axes()
py.pl.branches(tdata, depth_key="time", color="length", cmap="viridis", linewidth=LW, legend=False, ax=ax)
save_panel(fig, "d1_branches", "fig2")

# --- d2: + nodes (leaves and internal) colored by cluster ---
fig, ax = new_axes()
py.pl.branches(tdata, depth_key="time", color="length", cmap="viridis", linewidth=LW, legend=False, ax=ax)
py.pl.nodes(tdata, nodes="all", color="leiden", slot="obst", palette=palette, size=NODE, linewidths=0, legend=False, ax=ax)
save_panel(fig, "d2_nodes", "fig2")

# --- d3: + a single node marked with a large black star ---
fig, ax = new_axes()
py.pl.branches(tdata, depth_key="time", color="length", cmap="viridis", linewidth=LW, legend=False, ax=ax)
py.pl.nodes(tdata, nodes="all", color="leiden", slot="obst", palette=palette, size=NODE, linewidths=0, legend=False, ax=ax)
py.pl.nodes(tdata, nodes=[best], color="black", style="*", size=STAR, legend=False, ax=ax)
save_panel(fig, "d3_star", "fig2")

# --- d4: + leaf annotation of gene expression; legends for all features ---
# Legends are stacked in the order the plotting functions are called, so annotation is drawn before the nodes to get
# the legend order length, expression, cluster. Also saved without legends (d4_annotation_plot_only.svg) so the
# assembly can center the title over the plot alone.
for legend, name in [(True, "d4_annotation"), (False, "d4_annotation_plot_only")]:
    fig, ax = new_axes()
    py.pl.branches(tdata, depth_key="time", color="length", cmap="viridis", linewidth=LW, legend=legend, ax=ax)
    py.pl.annotation(
        tdata, keys=genes, width=0.2, gap=0.03, label=False, cmap="RdBu_r", share_cmap=True, legend=legend, ax=ax
    )
    py.pl.nodes(tdata, nodes="all", color="leiden", slot="obst", palette=palette, size=NODE, linewidths=0, legend=legend, ax=ax)
    py.pl.nodes(tdata, nodes=[best], color="black", style="*", size=STAR, legend=False, ax=ax)
    save_panel_full(fig, name, "fig2")
print("highlighted node:", best, "genes:", genes)

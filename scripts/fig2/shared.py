"""Shared helpers for the Figure 2 panels: one embedding of the simulated data, sampling, tree pruning, and the
tree + UMAP renderer (thin gray lines connect tree nodes to their observations in the UMAP)."""

import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pycea as py
import scanpy as sc
import treedata as td
from matplotlib.patches import ConnectionPatch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import CLUSTER_PALETTE, DATA, save_panel, set_style

warnings.filterwarnings("ignore")
SUBDIR = "fig2"
SEED = 0
N_CELLS = 20  # leaves in the small trees
LINE_KW = dict(color="#A0A0A0", linewidth=0.3, zorder=-1)  # behind the axes (whose backgrounds are transparent)
set_style()

# --- Embed the node-aligned simulation (leaves + internal nodes) so all panels share one UMAP ---
full = td.read_h5td(DATA / "fig1_sim_nodes.h5td")
sc.pp.pca(full, random_state=SEED)
sc.pp.neighbors(full, random_state=SEED)
sc.tl.umap(full, random_state=SEED)
sc.tl.leiden(full, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
palette = dict(zip(sorted(full.obs["leiden"].unique(), key=int), CLUSTER_PALETTE))
tree = full.obst["simulated"]
root = next(n for n, d in tree.in_degree() if d == 0)
leaves = [n for n, d in tree.out_degree() if d == 0]

# Common square UMAP window (same scale in x and y) used by every panel
_xy = full.obsm["X_umap"]
_center, _half = (_xy.min(0) + _xy.max(0)) / 2, 0.55 * np.ptp(_xy, axis=0).max()
UMAP_XLIM = (_center[0] - _half, _center[0] + _half)
UMAP_YLIM = (_center[1] - _half, _center[1] + _half)
TREE_NODE_SCALE = 1.5  # tree node circles relative to the `node_size` passed to render
UMAP_POINT_SCALE = 3.0  # UMAP circles relative to `node_size`
UMAP_SIDE = 1.0  # inches; UMAP axes are square
TREE_WIDTH = 1.1  # inches
UMAP_LEFT = 1.3  # inches; left edge of the UMAP axes (leaves the middle gap for the gray lines)

# --- Sample leaves: one per cluster (as many as possible), remaining ones at random ---
rng = np.random.default_rng(SEED)
leaf_clusters = full.obs.loc[leaves, "leiden"]
sampled = [rng.choice(leaf_clusters.index[leaf_clusters == c]) for c in palette]
rest = [leaf for leaf in leaves if leaf not in sampled]
sampled = list(sampled) + list(rng.choice(rest, N_CELLS - len(sampled), replace=False))
sampled = sampled[:N_CELLS]


def prune(tree, keep_leaves, from_mrca=False):
    """Subtree spanning `keep_leaves`: keeps leaves, branch points and the root; collapses unifurcations.

    With `from_mrca`, the root is the most recent common ancestor of the leaves (so trees do not share nodes).
    """
    ancestors = set(keep_leaves)
    for leaf in keep_leaves:
        ancestors |= nx.ancestors(tree, leaf)
    sub = tree.subgraph(ancestors)
    root_ = next(n for n, d in sub.in_degree() if d == 0)
    if from_mrca:
        while sub.out_degree(root_) == 1:
            root_ = next(sub.successors(root_))
        sub = sub.subgraph(nx.descendants(sub, root_) | {root_})
    kept = [n for n in sub if sub.out_degree(n) != 1 or n == root_]
    pruned = nx.DiGraph()
    for node in kept:
        pruned.add_node(node, **tree.nodes[node])
        parent = next(sub.predecessors(node), None)
        while parent is not None and parent not in kept:
            parent = next(sub.predecessors(parent), None)
        if parent is not None:
            pruned.add_edge(parent, node)
    return pruned


def node_labels(pruned):
    """Cluster label for every node of `pruned`: leaves keep their own, internal nodes get the descendant mode."""
    leaf_nodes = [n for n, d in pruned.out_degree() if d == 0]
    obs = full.obs.loc[leaf_nodes, ["leiden"]].copy()
    obs.index = obs.index.astype(str)
    tdata = td.TreeData(obs=obs, obst={"tree": pruned}, alignment="leaves")
    py.tl.ancestral_states(tdata, "leiden", method="mode")
    return py.get.node_df(tdata)["leiden"].astype(str).rename(index=str)


def make_tdata(obs_names, trees, labels, alignment, **kwargs):
    obs = pd.DataFrame({"leiden": labels.loc[[str(n) for n in obs_names]].to_numpy()}, index=[str(n) for n in obs_names])
    return td.TreeData(
        obs=obs,
        obsm={"X_umap": full.obsm["X_umap"][full.obs_names.get_indexer(obs_names)]},
        obst=trees if isinstance(trees, dict) else {"tree": trees},
        alignment=alignment,
        **kwargs,
    )


def render(tdata, name, node_set, node_size, figsize, groups=None, labels=None):
    """Trees (left) + UMAP (right) with gray lines joining tree nodes to UMAP cells.

    `groups` is a list of lists of tree keys; each group is drawn on its own tree axes (stacked vertically).
    By default all trees share one axes. `labels` optionally names each group.
    """
    groups = groups or [None]
    fig = plt.figure(figsize=figsize)
    width, height = figsize
    # UMAP: identical square axes (same size, same limits for x and y) in every panel, vertically centered
    ax_u = fig.add_axes([UMAP_LEFT / width, (height - UMAP_SIDE) / 2 / height, UMAP_SIDE / width, UMAP_SIDE / height])
    ax_u.patch.set_alpha(0)
    # Trees: stacked in the remaining space on the left
    gap = 0.12 * height / len(groups) if len(groups) > 1 else 0
    slot = (height - gap * (len(groups) - 1)) / len(groups)
    tree_axes, coords = [], {}
    for i, keys in enumerate(groups):
        ax_t = fig.add_axes([0, (height - (i + 1) * slot - i * gap) / height, TREE_WIDTH / width, slot / height])
        ax_t.patch.set_alpha(0)
        keys = keys[0] if keys is not None and len(keys) == 1 else keys  # single tree: pass the key itself
        py.pl.branches(tdata, depth_key="time", linewidth=0.5, tree=keys, ax=ax_t)
        py.pl.nodes(
            tdata, nodes=node_set, color="leiden", slot="obs", palette=palette, size=node_size * TREE_NODE_SCALE, legend=False,
            linewidths=0, tree=keys, ax=ax_t,
        )  # fmt: skip
        ax_t.set_ylim(-0.1 * np.pi, 2.1 * np.pi)  # room for node circles at the first/last leaf
        tree_axes.append(ax_t)
        coords[i] = {str(node): xy for (_, node), xy in ax_t._attrs["node_coords"].items()}  # keyed by (tree, node)
        if labels:
            ax_t.text(0.0, 1.0, labels[i], transform=ax_t.transAxes, fontsize=6, fontweight="bold", va="top")
    umap = tdata.obsm["X_umap"]
    ax_u.scatter(umap[:, 0], umap[:, 1], c=[palette[c] for c in tdata.obs["leiden"]], s=node_size * UMAP_POINT_SCALE, linewidths=0)
    ax_u.set_xlim(*UMAP_XLIM)
    ax_u.set_ylim(*UMAP_YLIM)
    ax_u.set_axis_off()
    for i, cell in enumerate(tdata.obs_names):
        for ax_t, c in zip(tree_axes, coords.values(), strict=True):
            if cell in c:
                fig.add_artist(
                    ConnectionPatch(
                        xyA=c[cell], coordsA=ax_t.transData, xyB=tuple(umap[i]), coordsB=ax_u.transData, **LINE_KW
                    )
                )
    return save_panel(fig, name, SUBDIR)

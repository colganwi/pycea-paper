"""Figure 2B: three TreeData alignment options (SVGs in plots/fig2/).

Each panel shows a lineage tree sampled from the Figure 1 simulation next to the UMAP of the same observations,
with thin gray lines connecting tree nodes to their observations (cells) in the UMAP. Internal nodes are colored by
the mode cluster of their descendant leaves (pycea.tl.ancestral_states, method="mode").

b1_leaves: alignment="leaves" - 20 leaves; only leaves are observations
b2_nodes:  alignment="nodes"  - same 20-leaf tree; leaves and internal nodes are observations
b3_subset: alignment="subset" - 100-leaf tree; 30 randomly chosen nodes (leaves + internal) are observed
"""

from shared import *  # noqa: F403  (embedding, sampling, prune, node_labels, make_tdata, render)

N_TREE_LEAVES = 100  # leaves of the larger tree (b3)
N_OBSERVED = 30  # randomly observed nodes in b3

# --- b1: leaves ---
small_tree = prune(tree, sampled)
small_labels = node_labels(small_tree)
render(make_tdata(sampled, small_tree, small_labels, "leaves"), "b1_leaves", "leaves", node_size=4, figsize=(2.3, 0.9))

# --- b2: nodes (leaves + internal nodes are observations) ---
render(make_tdata(list(small_tree.nodes), small_tree, small_labels, "nodes"), "b2_nodes", "all", node_size=4, figsize=(2.3, 0.9))

# --- b3: subset (larger tree; a random subset of leaves + internal nodes is observed) ---
big_leaves = list(rng.choice(leaves, N_TREE_LEAVES, replace=False))
big_tree = prune(tree, big_leaves)
big_labels = node_labels(big_tree)
observed = list(rng.choice(list(big_tree.nodes), N_OBSERVED, replace=False))
render(make_tdata(observed, big_tree, big_labels, "subset"), "b3_subset", observed, node_size=3, figsize=(2.3, 1.2))

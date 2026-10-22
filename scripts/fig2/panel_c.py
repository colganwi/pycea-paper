"""Figure 2C: TreeData supports multiple trees (SVGs in plots/fig2/).

c1_non_overlap: allow_overlap=False - multiple (two) trees with disjoint sets of cells
c2_overlap:     allow_overlap=True  - multiple (two) trees (NJ and UPGMA, reconstructed with Cassiopeia from the same
                                      character matrix) for the same set of cells
"""

import cassiopeia as cas
from shared import *  # noqa: F403  (embedding, sampling, prune, make_tdata, render)

# --- c1: two non-overlapping trees (the two main clades of the lineage) ---
first_split = root
while tree.out_degree(first_split) == 1:  # skip the root unifurcation
    first_split = next(tree.successors(first_split))
clade_a, clade_b = (nx.descendants(tree, c) | {c} for c in tree.successors(first_split))
cells_a = [c for c in sampled if c in clade_a]
cells_b = [c for c in sampled if c in clade_b]
trees = {"tree1": prune(tree, cells_a, from_mrca=True), "tree2": prune(tree, cells_b, from_mrca=True)}
leaf_labels = full.obs["leiden"].astype(str).rename(index=str)
tdata = make_tdata(cells_a + cells_b, trees, leaf_labels, "leaves", allow_overlap=False)
render(tdata, "c1_non_overlap", "leaves", node_size=4, figsize=(2.3, 1.5))

# --- c2: two trees for the same cells (NJ and UPGMA from the character matrix) ---
chars = full.obsm["characters"].loc[sampled].astype(str)
tdata = make_tdata(sampled, {}, leaf_labels, "leaves", allow_overlap=True)
tdata.obsm["characters"] = chars.set_index(tdata.obs_names)
tdata.uns["unmodified_state"] = full.uns["unmodified_state"]
tdata.uns["missing_state"] = full.uns["missing_state"]
cas.solver.nj(tdata, characters_key="characters", key_added="NJ")
cas.solver.upgma(tdata, characters_key="characters", key_added="UPGMA")
for key in ["NJ", "UPGMA"]:
    py.pp.add_depth(tdata, tree=key, key_added="time")
render(tdata, "c2_overlap", "leaves", node_size=4, figsize=(2.3, 2.2), groups=[["NJ"], ["UPGMA"]], labels=["NJ", "UPGMA"])

"""Simulate the lineage tracing + expression dataset used in Figure 1C.

Outputs:
  data/fig1_sim.h5td        200 cells, leaf-aligned tree, expression + 20 characters (Figure 1C)
  data/fig1_sim_nodes.h5td  same simulation (same seed) with alignment="nodes", i.e. internal nodes are also
                            observations (used for the node-alignment schematic in Figure 2B)
"""

import sys
from pathlib import Path

import cassiopeia as cas
import networkx as nx
import treedata as td

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import DATA

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 1

def simulate(alignment):
    # Tree: 200 extant cells
    tdata = cas.sim.birth_death_process(num_extant=200, random_seed=SEED, alignment=alignment)

    # Lineage tracing: 20 independent characters
    cas.sim.stochastic_tracing(
        tdata,
        number_of_cassettes=20,
        size_of_cassette=1,
        mutation_rate=0.12,  # ~70% of sites edited at the leaves
        number_of_states=20,
        state_generating_distribution=lambda: 1.0,
        random_seed=SEED,
        tree_key="simulated",
    )

    # Expression: cells follow a hierarchical cell-state trajectory.
    # - Every cell in the (ultrametric) lineage tree sits at the end of its trajectory branch, so cell states are
    #   set by (i) which terminal branch a lineage commits to and (ii) latent_noise.
    # - Fate choices are made when a lineage passes a trajectory branch point, so branch points must occur late
    #   enough (>= ~0.35 of the experiment) that many lineages exist and choose independently.
    # - Graded split times (0.35 ... 1.0) give hierarchical structure: early splits -> distinct clusters, late
    #   splits -> neighboring sub-populations, which gives touching clusters and continuous bridges in the UMAP.
    # - High momentum (0.95) makes latent drift persistent along each branch (directional trajectories);
    #   latent_noise (0.8) sets the within-state spread relative to the between-state distances.
    split_times = [0.35, 0.5, 0.65, 0.8, 0.9, 1.0]  # time of each successive round of binary splits
    trajectory = nx.DiGraph()
    times = {"root": 0, "n": split_times[0]}
    trajectory.add_edge("root", "n")
    frontier = ["n"]
    for t in split_times[1:]:
        new = []
        for node in frontier:
            for i in range(2):
                child = f"{node}{i}"
                trajectory.add_edge(node, child)
                times[child] = t
                new.append(child)
        frontier = new
    nx.set_node_attributes(trajectory, times, "time")
    nx.set_edge_attributes(trajectory, {(u, v): 0.5 for u, v in trajectory.edges}, "prob")  # equal fate probabilities
    trajectory = td.TreeData(obst={"tree": trajectory}, alignment="nodes")
    trajectory = cas.sim.brownian_expression(trajectory, momentum=0.95, random_seed=SEED)
    tdata = cas.sim.trajectory_expression(
        tdata, trajectory, n_genes=500, latent_noise=0.8, prob_key="prob", random_seed=SEED
    )
    # uns["priors"] has non-string keys, which h5td cannot serialize; not needed downstream
    tdata.uns.pop("priors", None)
    return tdata


DATA.mkdir(exist_ok=True)
for alignment, name in [("leaves", "fig1_sim"), ("nodes", "fig1_sim_nodes")]:
    tdata = simulate(alignment)
    tdata.write_h5td(DATA / f"{name}.h5td")
    print(name, tdata)

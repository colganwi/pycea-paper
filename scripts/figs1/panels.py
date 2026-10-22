"""Figure S1 panels (SVGs in plots/figs1/): benchmarks, each a log-log line plot of time by tree size.

a_treedata:    TreeData vs AnnData disk write / read and slicing (data/figs1_treedata.csv)
b_plot:        tree plotting time, Pycea vs Cassiopeia 2.1.0 (data/figs1_pycea.csv, data/figs1_cassiopeia.csv)
c_functions:   runtime of key Pycea analysis functions (data/figs1_pycea.csv)
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import DATA, colors, save_panel_full, set_style

set_style()
SUBDIR = "figs1"
SIZE = (2.3, 1.9)
XLABEL, YLABEL = "Number of cells (leaves)", "Time (s)"


def setup(ax):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.spines[["top", "right"]].set_visible(False)


def line(ax, d, color, ls="-", label=None):
    d = d.sort_values("n_leaves")
    ax.plot(d["n_leaves"], d["seconds"], color=color, linestyle=ls, linewidth=1.1, marker="o", markersize=2.2, label=label)


# --- A: TreeData vs AnnData ---
df = pd.read_csv(DATA / "figs1_treedata.csv")
ops = {"write": "Write", "read": "Read", "slice": "Slice"}
fig, ax = plt.subplots(figsize=SIZE)
for (op, label), color in zip(ops.items(), colors, strict=False):
    for obj, ls in [("TreeData", "-"), ("AnnData", (0, (2, 1.5)))]:
        line(ax, df[(df["object"] == obj) & (df["operation"] == op)], color, ls)
setup(ax)
handles = [Line2D([], [], color=c, linewidth=1.1, label=label) for (op, label), c in zip(ops.items(), colors, strict=False)]
handles += [Line2D([], [], color="black", linestyle=ls, linewidth=1.1, label=obj) for obj, ls in [("TreeData", "-"), ("AnnData", (0, (2, 1.5)))]]
ax.legend(handles=handles, frameon=False, loc="upper left", fontsize=6, handlelength=1.8, borderaxespad=0.2)
save_panel_full(fig, "a_treedata", SUBDIR)

# --- B: tree plotting, Pycea vs Cassiopeia (2.1.0) ---
pyc = pd.read_csv(DATA / "figs1_pycea.csv")
cas_ = pd.read_csv(DATA / "figs1_cassiopeia.csv")
fig, ax = plt.subplots(figsize=SIZE)
line(ax, pyc[pyc["task"] == "plot"], colors[0], label="Pycea")
line(ax, cas_[cas_["task"] == "plot"], colors[1], label="Cassiopeia 2.1.0")
setup(ax)
ax.legend(frameon=False, loc="upper left", fontsize=6, handlelength=1.8, borderaxespad=0.2)
save_panel_full(fig, "b_plot", SUBDIR)

# --- C: key Pycea functions ---
names = {
    "add_depth": "pp.add_depth",
    "clades": "tl.clades",
    "ancestral_states": "tl.ancestral_states",
    "tree_neighbors": "tl.tree_neighbors",
    "ancestral_linkage": "tl.ancestral_linkage",
    "fitness": "tl.fitness",
    "expansion_test": "tl.expansion_test",
}
fig, ax = plt.subplots(figsize=SIZE)
for (task, label), color in zip(names.items(), colors, strict=False):
    line(ax, pyc[pyc["task"] == task], color, label=label)
setup(ax)
ax.legend(frameon=False, loc="upper left", fontsize=5.5, handlelength=1.5, borderaxespad=0.2, labelspacing=0.25)
save_panel_full(fig, "c_functions", SUBDIR)

"""Figure S1B/C benchmark (Pycea side): tree plotting time and runtime of key analysis functions by tree size.

Output: data/figs1_pycea.csv (task = "plot" or a function name). Cassiopeia plotting is timed separately
(benchmark_cassiopeia.py, old Cassiopeia environment).
"""

import sys
import warnings
from pathlib import Path

import cassiopeia as cas
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pycea as py

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_common import GIVE_UP, timeit
from common import DATA

warnings.filterwarnings("ignore")
SIZES = [100, 300, 1000, 3000, 10000, 30000, 100000]
rng = np.random.default_rng(0)


def plot(tdata):
    fig, ax = plt.subplots()
    py.pl.tree(tdata, depth_key="time", ax=ax)
    fig.canvas.draw()
    plt.close(fig)


rows = []
for n in SIZES:
    print(f"n={n}: simulating", flush=True)
    tdata = cas.sim.birth_death_process(num_extant=n, random_seed=0, alignment="leaves")
    tdata.obs["group"] = pd.Categorical(rng.choice(list("abcde"), n))
    tdata.obs["value"] = rng.normal(size=n)
    tasks = {
        "plot": lambda: plot(tdata),
        "add_depth": lambda: py.pp.add_depth(tdata),
        "clades": lambda: py.tl.clades(tdata, depth=3, depth_key="time"),
        "ancestral_states": lambda: py.tl.ancestral_states(tdata, "group", method="mode"),
        "tree_neighbors": lambda: py.tl.tree_neighbors(tdata, n_neighbors=10, depth_key="time"),
        "ancestral_linkage": lambda: py.tl.ancestral_linkage(tdata, groupby="group", depth_key="time"),
        "fitness": lambda: py.tl.fitness(tdata, depth_key="time", method="lbi"),
        "expansion_test": lambda: py.tl.expansion_test(tdata, copy=True),
    }
    slowest = 0
    for name, fn in tasks.items():
        sec = timeit(fn)
        slowest = max(slowest, sec)
        rows.append({"n_leaves": n, "task": name, "seconds": sec})
    print(pd.DataFrame(rows).query("n_leaves == @n").to_string(index=False), flush=True)
    pd.DataFrame(rows).to_csv(DATA / "figs1_pycea.csv", index=False)
    if slowest > GIVE_UP:
        break

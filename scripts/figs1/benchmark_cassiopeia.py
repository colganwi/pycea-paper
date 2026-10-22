"""Figure S1B benchmark (Cassiopeia side): tree plotting time with Cassiopeia 2.1.0 (before the 3.0 refactor).

Run in an environment with Cassiopeia < 3.0 (e.g. conda env `cas-testing`). Output: data/figs1_cassiopeia.csv.
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_common import GIVE_UP, timeit
from common import DATA

warnings.filterwarnings("ignore")
SIZES = [100, 300, 1000, 3000, 10000, 30000, 100000]
print("Cassiopeia", cas.__version__, flush=True)
assert cas.__version__.startswith("2."), "needs Cassiopeia < 3.0 for the old plotting functions"


def plot(tree):
    fig = cas.pl.plot_matplotlib(tree, orient="right")
    plt.gcf().canvas.draw()
    plt.close("all")


rows = []
for n in SIZES:
    np.random.seed(0)
    sim = cas.sim.BirthDeathFitnessSimulator(
        birth_waiting_distribution=lambda scale: np.random.exponential(scale), initial_birth_scale=1, num_extant=n
    )
    tree = sim.simulate_tree()
    sec = timeit(lambda: plot(tree))
    rows.append({"n_leaves": n, "task": "plot", "seconds": sec})
    print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(DATA / "figs1_cassiopeia.csv", index=False)
    if sec > GIVE_UP:
        break

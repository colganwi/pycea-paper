"""Figure 4 analysis (growth dynamics; follows the pycea growth-dynamics tutorial).

Loads the Koblan et al. 2025 4T1 tumor dataset, then computes clades (shared ancestor at day 5), SBD fitness, and the
correlation of each gene with fitness. Results are written to data/ so the panel scripts do not recompute them:

  data/fig4_obs.csv        per-cell clade and fitness
  data/fig4_gene_corr.csv  correlation of every gene with fitness
"""

import sys
import warnings
from pathlib import Path

import pycea as py
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import DATA

warnings.filterwarnings("ignore")

tdata = py.datasets.koblan25()
py.tl.clades(tdata, depth=5, depth_key="time", update=False)
py.tl.fitness(tdata, depth_key="time", method="sbd")

DATA.mkdir(exist_ok=True)
tdata.obs[["clade", "fitness"]].to_csv(DATA / "fig4_obs.csv")
expr = sc.get.obs_df(tdata, keys=list(tdata.var_names))
expr.corrwith(tdata.obs["fitness"]).sort_values(ascending=False).rename("correlation").to_csv(
    DATA / "fig4_gene_corr.csv", index_label="gene"
)
print(tdata.obs["clade"].nunique(), "clades;", tdata.obs["fitness"].notna().sum(), "cells with fitness")

"""Figure 3 analysis (ancestral linkage; follows the mouse-embryo tutorial) for embryo E9.5-R1.

Computes neural crest ancestral linkage per cell and pairwise lineage-level linkage, and saves the results:

  data/fig3_obs.csv          per-cell neural crest linkage ("nc_linkage")
  data/fig3_linkage_uns.pkl  the `uns` entries of the pairwise lineage linkage (restored by the panel scripts)
"""

import pickle
import sys
import warnings
from pathlib import Path

import pycea as py

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import DATA

warnings.filterwarnings("ignore")
tdata = py.datasets.colgan26(embryos="E9.5-R1")
print(tdata, flush=True)

py.tl.ancestral_linkage(tdata, groupby="lineage", target="Neural crest", key_added="nc")
print("nc linkage done", flush=True)
before = set(tdata.uns)
py.tl.ancestral_linkage(tdata, groupby="lineage")
new_keys = sorted(set(tdata.uns) - before)
print("pairwise linkage done; new uns keys:", new_keys, flush=True)

DATA.mkdir(exist_ok=True)
tdata.obs[["nc_linkage"]].to_csv(DATA / "fig3_obs.csv")
with open(DATA / "fig3_linkage_uns.pkl", "wb") as f:
    pickle.dump({k: tdata.uns[k] for k in new_keys}, f)

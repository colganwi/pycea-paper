"""Figure S1A benchmark: TreeData vs AnnData disk read/write and slicing time by tree size (20,000 genes).

Output: data/figs1_treedata.csv. Run via Slurm: sbatch scripts/figs1/benchmark.sbatch <python> scripts/figs1/benchmark_treedata.py
"""

import os
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import anndata as ad
import cassiopeia as cas
import numpy as np
import pandas as pd
import scipy.sparse as sp
import treedata as td

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_common import GIVE_UP, timeit
from common import DATA

warnings.filterwarnings("ignore")
SIZES = [100, 300, 1000, 3000, 10000, 30000, 100000]  # simulated expression: 20,000 genes (sparse)
N_GENES = 20000
DENSITY = 0.05  # fraction of nonzero entries (simulated sparse counts, as in scRNA-seq)
SLICE_FRACTION = 0.1  # slice keeps a random 10% of the cells
rng = np.random.default_rng(0)
tmp = Path(tempfile.mkdtemp(prefix=f"{os.environ.get('USER', 'user')}_pycea_bench_"))

rows = []
try:
    for n in SIZES:
        print(f"n={n}: simulating", flush=True)
        tree = cas.sim.birth_death_process(num_extant=n, random_seed=0, alignment="leaves")
        X = sp.random(n, N_GENES, density=DENSITY, format="csr", random_state=0, data_rvs=lambda k: rng.poisson(3, k) + 1).astype(np.float32)
        tdata = td.TreeData(X=X, obs=tree.obs.copy(), obst=dict(tree.obst), alignment="leaves")
        adata = ad.AnnData(X=X, obs=tree.obs.copy())
        idx = np.sort(rng.choice(n, int(n * SLICE_FRACTION), replace=False))
        td_path, ad_path = tmp / "t.h5td", tmp / "a.h5ad"
        ops = {
            ("TreeData", "write"): lambda: tdata.write_h5td(td_path),
            ("AnnData", "write"): lambda: adata.write_h5ad(ad_path),
        }
        slowest = 0
        for (obj, op), fn in ops.items():
            sec = timeit(fn)
            slowest = max(slowest, sec)
            rows.append({"n_leaves": n, "object": obj, "operation": op, "seconds": sec})
        ops = {
            ("TreeData", "read"): lambda: td.read_h5td(td_path),
            ("AnnData", "read"): lambda: ad.read_h5ad(ad_path),
            ("TreeData", "slice"): lambda: tdata[idx].copy(),
            ("AnnData", "slice"): lambda: adata[idx].copy(),
        }
        for (obj, op), fn in ops.items():
            sec = timeit(fn)
            slowest = max(slowest, sec)
            rows.append({"n_leaves": n, "object": obj, "operation": op, "seconds": sec})
        print(pd.DataFrame(rows).query("n_leaves == @n").to_string(index=False), flush=True)
        pd.DataFrame(rows).to_csv(DATA / "figs1_treedata.csv", index=False)
        if slowest > GIVE_UP:
            break
finally:
    shutil.rmtree(tmp, ignore_errors=True)

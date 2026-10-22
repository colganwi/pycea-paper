"""Figure 3 analysis (gene heritability) for embryo E9.5-R1: Moran's I of gene expression over the lineage tree.

Needs expression, which is not in the pycea dataset; it is read from the devmap-paper data (log1p expression, all
genes; the 5,000 highly variable genes in data/hvg.csv are tested). Writes:

  data/fig3_placeholder_expr.csv  expression of 4 random genes (placeholder for panels G/H until results are in)
  data/fig3_moran.csv             Moran's I (heritability) per gene: autocorr, pval_norm, var_norm (pycea.tl.autocorr)
  data/fig3_top_expr.csv          expression of the 4 genes with the highest Moran's I
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pycea as py
import scipy.sparse as sp
import treedata as td

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import DATA

warnings.filterwarnings("ignore")
DEVMAP = Path("/lab/solexa_weissman/wcolgan/devmap-paper/data")
N_TOP = 4

hvg = pd.read_csv(DEVMAP / "hvg.csv", header=None)[0].tolist()
tdata = td.read_h5td(DEVMAP / "embryos" / "E9.5-R1.h5td")
print(tdata, flush=True)
tdata = tdata[:, [g for g in hvg if g in tdata.var_names]].copy()
print("HVG subset:", tdata.shape, flush=True)

# The embryo file stores raw UMI counts: library-size normalize (using each cell's total counts over all genes,
# since only the HVGs are kept) to 10,000 counts per cell, then log1p
tdata.X = sp.diags(1e4 / tdata.obs["total_counts"].to_numpy()) @ sp.csr_matrix(tdata.X)
tdata.X.data = np.log1p(tdata.X.data)


def expression(genes):
    return pd.DataFrame(np.asarray(tdata[:, genes].X.todense()), index=tdata.obs_names, columns=genes)


# Placeholder: 4 random genes
rng = np.random.default_rng(0)
expression(list(rng.choice(tdata.var_names, N_TOP, replace=False))).to_csv(DATA / "fig3_placeholder_expr.csv")

# Moran's I over the tree neighbor graph (5 nearest neighbors by tree distance)
py.tl.tree_neighbors(tdata, n_neighbors=5, depth_key="time")
print("tree neighbors done", flush=True)
# Default `keys` tests every gene in tdata.X (kept sparse; the tdata was already subset to the HVGs)
result = py.tl.autocorr(tdata, connect_key="tree_connectivities", method="moran", copy=True)
print(result.head(), flush=True)
result.sort_values("autocorr", ascending=False).to_csv(DATA / "fig3_moran.csv", index_label="gene")
top = result.sort_values("autocorr", ascending=False).head(N_TOP).index.tolist()
print("top genes:", top, flush=True)
expression(top).to_csv(DATA / "fig3_top_expr.csv")

# Regenerate figures (uses the pycea-paper conda env; override with PY=...)
PY ?= PYTHONHASHSEED=0 /lab/solexa_weissman/wcolgan/tools/miniconda3/envs/pycea-paper/bin/python  # fixed hash seed: obs order in simulated data depends on set ordering

fig1: manuscript/figures/figure1.pdf

data/fig1_sim.h5td data/fig1_sim_nodes.h5td &: scripts/fig1/simulate.py
	$(PY) scripts/fig1/simulate.py

plots/fig1/c4_linkage.svg: scripts/fig1/panel_c.py data/fig1_sim.h5td scripts/common.py config.py
	$(PY) scripts/fig1/panel_c.py

manuscript/figures/figure1.pdf: scripts/fig1/assemble.py plots/fig1/c4_linkage.svg plots/ecosystem.png plots/fig1/datatypes.png plots/pycea_logo.png plots/scanpy_logo.png plots/cassipopeia_logo.png scripts/common.py
	$(PY) scripts/fig1/assemble.py

fig2: manuscript/figures/figure2.pdf

plots/fig2/b3_subset.svg: scripts/fig2/panel_b.py scripts/fig2/shared.py data/fig1_sim_nodes.h5td scripts/common.py config.py
	$(PY) scripts/fig2/panel_b.py

plots/fig2/d4_annotation.svg: scripts/fig2/panel_d.py data/fig1_sim.h5td scripts/common.py config.py
	$(PY) scripts/fig2/panel_d.py

plots/fig2/c2_overlap.svg: scripts/fig2/panel_c.py scripts/fig2/shared.py data/fig1_sim_nodes.h5td scripts/common.py config.py
	$(PY) scripts/fig2/panel_c.py

manuscript/figures/figure2.pdf: scripts/fig2/assemble.py plots/fig2/b3_subset.svg plots/fig2/c2_overlap.svg plots/fig2/d4_annotation.svg plots/treedata_schema.png scripts/common.py
	$(PY) scripts/fig2/assemble.py

# Figure 4: the analysis and panel scripts are heavy (146k cells) and run through Slurm; run them first:
#   sbatch scripts/fig4/prepare.sbatch   # -> data/fig4_obs.csv, data/fig4_gene_corr.csv
#   sbatch scripts/fig4/panels.sbatch    # -> plots/fig4/*.svg
# `make fig4` only assembles the panels into the figure.
fig4: manuscript/figures/figure4.pdf

manuscript/figures/figure4.pdf: scripts/fig4/assemble.py scripts/fig4/layout.py $(wildcard plots/fig4/*.svg) scripts/common.py
	$(PY) scripts/fig4/assemble.py

# Figure S1 (benchmarks): the benchmark scripts run through Slurm (see scripts/figs1/benchmark.sbatch for the three
# sbatch commands; Cassiopeia 2.1.0 timings use the `cas-testing` conda env). `make figs1` re-plots and assembles.
figs1: manuscript/figures/figureS1.pdf

plots/figs1/c_functions.svg: scripts/figs1/panels.py data/figs1_treedata.csv data/figs1_pycea.csv data/figs1_cassiopeia.csv scripts/common.py config.py
	$(PY) scripts/figs1/panels.py

manuscript/figures/figureS1.pdf: scripts/figs1/assemble.py plots/figs1/c_functions.svg scripts/common.py
	$(PY) scripts/figs1/assemble.py

# Figure 3 (mouse embryo E9.5-R1): analysis and panels run through Slurm (see scripts/fig3/prepare.sbatch for the
# commands: linkage, moran, then panels). `make fig3` only assembles the panels.
fig3: manuscript/figures/figure3.pdf

manuscript/figures/figure3.pdf: scripts/fig3/assemble.py scripts/fig3/layout.py $(wildcard plots/fig3/*.svg) scripts/common.py
	$(PY) scripts/fig3/assemble.py

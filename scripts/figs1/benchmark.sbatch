#!/bin/bash
#SBATCH --job-name=wcolgan_pycea_figs1
#SBATCH --partition=20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=0-08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Submit from the repository root with the python to use and the benchmark script:
#   ENVS=/lab/solexa_weissman/wcolgan/tools/miniconda3/envs
#   sbatch --job-name=wcolgan_pycea_figs1_td scripts/figs1/benchmark.sbatch $ENVS/pycea-paper/bin/python scripts/figs1/benchmark_treedata.py
#   sbatch --job-name=wcolgan_pycea_figs1_py scripts/figs1/benchmark.sbatch $ENVS/pycea-paper/bin/python scripts/figs1/benchmark_pycea.py
#   sbatch --job-name=wcolgan_pycea_figs1_cas scripts/figs1/benchmark.sbatch $ENVS/cas-testing/bin/python scripts/figs1/benchmark_cassiopeia.py
PYTHONHASHSEED=0 "$1" "$2"

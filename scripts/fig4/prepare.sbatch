#!/bin/bash
#SBATCH --job-name=wcolgan_pycea_fig4_prep
#SBATCH --partition=20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=0-02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Submit from the repository root: sbatch scripts/fig4/prepare.sbatch
PYTHONHASHSEED=0 /lab/solexa_weissman/wcolgan/tools/miniconda3/envs/pycea-paper/bin/python scripts/fig4/prepare.py

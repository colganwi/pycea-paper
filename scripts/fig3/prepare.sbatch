#!/bin/bash
#SBATCH --job-name=wcolgan_pycea_fig3
#SBATCH --partition=20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Submit from the repository root (the job name suffix identifies the analysis):
#   sbatch --job-name=wcolgan_pycea_fig3_linkage scripts/fig3/prepare.sbatch scripts/fig3/prepare_linkage.py
#   sbatch --job-name=wcolgan_pycea_fig3_moran   scripts/fig3/prepare.sbatch scripts/fig3/prepare_moran.py
#   sbatch --job-name=wcolgan_pycea_fig3_cache  --mem=48G scripts/fig3/prepare.sbatch scripts/fig3/panels_data.py --cache [--tree]   # after the two above
# The cache (data/fig3_cells.csv.gz) makes the UMAP panels quick (~1.5 min, 16G); redraw them on a compute node with
#   srun --job-name=wcolgan_pycea_fig3_umaps --partition=20 --cpus-per-task=2 --mem=16G --time=0:30:00 \
#     bash -c "PYTHONHASHSEED=0 <env python> scripts/fig3/panels_data.py [--only b c d g]"
# and run scripts/fig3/panels_schematics.py directly (seconds), then `make fig3`.
PYTHONHASHSEED=0 /lab/solexa_weissman/wcolgan/tools/miniconda3/envs/pycea-paper/bin/python "$@"

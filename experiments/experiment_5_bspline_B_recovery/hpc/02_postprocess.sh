#!/bin/bash
#SBATCH --job-name=exp5_pair200
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/postprocess/slurm-%j.out
#SBATCH --error=logs/postprocess/slurm-%j.err

set -euo pipefail
module --force purge
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="$HOME/packages-R4.1"
cd "$SLURM_SUBMIT_DIR"

Rscript code/03_combine_results.R
Rscript code/04_compare_with_gamma.R
Rscript code/05_compare_three_models.R

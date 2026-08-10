#!/bin/bash
#SBATCH --job-name=exp2_finish
#SBATCH --time=03:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/slurm-finish-%j.out
#SBATCH --error=logs/slurm-finish-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:?Submit this script from the Experiment 2 directory}"

module --force purge
module load StdEnv/2020
module load r/4.1.2

export R_LIBS_USER="$HOME/packages-R4.1"

Rscript code/03_Combine_MIF2_Array_Results.R
Rscript code/04_Plot_Best_Fit.R
Rscript code/05_regenerate_figures.R

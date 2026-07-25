#!/bin/bash
#SBATCH --job-name=exp2_finish
#SBATCH --time=03:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --chdir=/global/home/hpc6245/experiment2_large_scale_gamma_B_recovery_HPC
#SBATCH --output=/global/home/hpc6245/experiment2_large_scale_gamma_B_recovery_HPC/logs/slurm-finish-%j.out
#SBATCH --error=/global/home/hpc6245/experiment2_large_scale_gamma_B_recovery_HPC/logs/slurm-finish-%j.err

set -euo pipefail

module --force purge
module load StdEnv/2020
module load r/4.1.2

export R_LIBS_USER="$HOME/packages-R4.1"

Rscript code/03_Combine_MIF2_Array_Results.R
Rscript code/04_Plot_Best_Fit.R

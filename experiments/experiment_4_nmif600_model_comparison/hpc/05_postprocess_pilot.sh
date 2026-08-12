#!/bin/bash
#SBATCH --job-name=exp4_pilotpost
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/pilot/postprocess-%j.out
#SBATCH --error=logs/pilot/postprocess-%j.err

set -euo pipefail
module --force purge
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="$HOME/packages-R4.1"
cd "$SLURM_SUBMIT_DIR"

TASKS="1,50,100,150,200"

Rscript code/04_combine_results.R \
  gamma results_raw/gamma results/pilot_combined/gamma "$TASKS"

Rscript code/04_combine_results.R \
  constant results_raw/constant results/pilot_combined/constant "$TASKS"

Rscript code/05_compare_models.R \
  results/pilot_combined/gamma \
  results/pilot_combined/constant \
  results/pilot_comparison \
  figures/pilot \
  false

Rscript code/06_make_convergence_diagnostics.R \
  results/pilot_combined/gamma \
  results/pilot_combined/constant \
  figures/pilot

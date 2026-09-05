#!/bin/bash
# Independent production-parameter pilot for Experiment 5 task 1.
# This file is never invoked by submit_all.sh and must be submitted manually.

#SBATCH --job-name=exp5_task1pilot
#SBATCH --time=08:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/pilot/slurm-%j.out
#SBATCH --error=logs/pilot/slurm-%j.err

set -euo pipefail
if [[ "${EXP5_TASK1_PILOT_GUARD:-}" != "YES" ]]; then
  printf '%s\n' \
    'Task-1 pilot guard missing; submit with hpc/submit_task1_pilot.sh.' >&2
  exit 2
fi
module --force purge
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="$HOME/packages-R4.1"
cd "$SLURM_SUBMIT_DIR"

# Pin the production computation explicitly for this isolated pilot.
export EXP5_NMIF=600
export EXP5_NP_MIF=5000
export EXP5_N_START=10
export EXP5_NP_EVAL=50000
export EXP5_N_PF_EVALS=5

Rscript --vanilla -e 'stopifnot(requireNamespace("pomp", quietly = TRUE)); cat(R.version.string, "\npomp ", as.character(packageVersion("pomp")), "\n", sep = "")'

Rscript --vanilla code/02_fit_bspline_B.R \
  1 \
  ../experiment_4_nmif600_model_comparison/shared_data \
  results_pilot/bspline \
  results/paired_input_manifest.csv

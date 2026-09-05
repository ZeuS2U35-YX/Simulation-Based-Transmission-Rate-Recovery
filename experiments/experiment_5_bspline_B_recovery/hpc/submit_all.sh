#!/bin/bash
# Validate the promoted task-1 production pilot, submit tasks 2-200, then run
# one exact 200-row paired comparison.

set -euo pipefail
cd "$(dirname "$0")/.."

if [[ "${1:-}" == "--dry-run" ]]; then
  printf '%s\n' \
    'Experiment 5 submission plan (no commands executed):' \
    '  1. Validate the 200 paired inputs and promoted task 1.' \
    '  2. Submit B-spline tasks 2-200 as one Slurm array.' \
    '  3. Submit post-processing with an afterok dependency.' \
    'Run with EXP5_CONFIRM_SUBMIT=YES to submit.'
  exit 0
fi
if [[ $# -ne 0 ]]; then
  printf 'Usage: bash hpc/submit_all.sh [--dry-run]\n' >&2
  exit 2
fi
if [[ "${EXP5_CONFIRM_SUBMIT:-}" != "YES" ]]; then
  printf '%s\n' \
    'Submission stopped before validation, file creation, or sbatch.' \
    'Review with: bash hpc/submit_all.sh --dry-run' \
    'Submit with: EXP5_CONFIRM_SUBMIT=YES bash hpc/submit_all.sh' >&2
  exit 2
fi

mkdir -p logs/bspline logs/postprocess results/bspline \
  results/combined/bspline results/comparison results

module --force purge
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="$HOME/packages-R4.1"
command -v sbatch >/dev/null
Rscript --vanilla -e 'stopifnot(requireNamespace("pomp", quietly = TRUE)); cat(R.version.string, "\npomp ", as.character(packageVersion("pomp")), "\n", sep = "")'

# Validate the accepted Experiment 4 batch and refresh the paired manifest.
Rscript --vanilla code/01_validate_paired_inputs.R

# Refuse to submit the 199-task array until the independently run production
# pilot is present in the canonical result tree. Promotion is deliberately a
# separate, non-overwriting manual step documented in README.md.
Rscript --vanilla code/00_validate_completed_bspline_task.R \
  1 results/bspline/task_001 results/paired_input_manifest.csv

clean_job_id() { echo "${1%%;*}"; }

bspline_raw=$(sbatch --parsable --array=2-200%20 hpc/01_bspline_array.sh)
bspline_job=$(clean_job_id "$bspline_raw")
post_raw=$(sbatch --parsable --dependency=afterok:"$bspline_job" \
  hpc/02_postprocess.sh)
post_job=$(clean_job_id "$post_raw")

printf 'BSPLINE_JOB_ID=%s\nPOSTPROCESS_JOB_ID=%s\n' \
  "$bspline_job" "$post_job" > results/full_job_ids.env

printf '%s\n' \
  "Validated the promoted task-1 production pilot." \
  "Submitted one 199-task B-spline array (tasks 2-200): $bspline_job" \
  "Submitted one dependent paired postprocess job: $post_job" \
  "No observed-data or Gamma-noise jobs were submitted."

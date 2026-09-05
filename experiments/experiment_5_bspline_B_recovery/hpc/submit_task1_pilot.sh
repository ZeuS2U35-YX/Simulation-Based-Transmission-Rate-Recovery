#!/bin/bash
# Safely preview or submit the independent Experiment 5 task-1 pilot.

set -euo pipefail
cd "$(dirname "$0")/.."

if [[ "${1:-}" == "--dry-run" ]]; then
  printf '%s\n' \
    'Experiment 5 task-1 pilot plan (no commands executed):' \
    '  1. Check that sbatch is available.' \
    '  2. Create only the pilot log and result directories.' \
    '  3. Submit task 1 with production fitting settings.' \
    'Run with EXP5_CONFIRM_TASK1_PILOT=YES to submit.'
  exit 0
fi
if [[ $# -ne 0 ]]; then
  printf 'Usage: bash hpc/submit_task1_pilot.sh [--dry-run]\n' >&2
  exit 2
fi
if [[ "${EXP5_CONFIRM_TASK1_PILOT:-}" != "YES" ]]; then
  printf '%s\n' \
    'Submission stopped before directory creation or sbatch.' \
    'Review with: bash hpc/submit_task1_pilot.sh --dry-run' \
    'Submit with: EXP5_CONFIRM_TASK1_PILOT=YES bash hpc/submit_task1_pilot.sh' >&2
  exit 2
fi

command -v sbatch >/dev/null
mkdir -p logs/pilot results_pilot/bspline

pilot_raw=$(sbatch --parsable --export=ALL,EXP5_TASK1_PILOT_GUARD=YES \
  hpc/03_task1_production_pilot.sh)
pilot_job="${pilot_raw%%;*}"

printf 'TASK1_PILOT_JOB_ID=%s\n' "$pilot_job" \
  > results_pilot/task1_pilot_job_id.env
printf '%s\n' \
  "Submitted the Experiment 5 task-1 production pilot: $pilot_job" \
  'Job ID saved to results_pilot/task1_pilot_job_id.env.'

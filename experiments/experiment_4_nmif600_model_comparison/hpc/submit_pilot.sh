#!/bin/bash
# Submit five full-setting pilot tasks. These outputs are production-quality
# and will be skipped, not overwritten, when submit_all.sh is run later.

set -euo pipefail
cd "$(dirname "$0")/.."

if [[ "${1:-}" == "--dry-run" ]]; then
  printf '%s\n' \
    'Experiment 4 pilot submission plan (no commands executed):' \
    '  1. Create the fixed parameter list.' \
    '  2. Submit tasks 1, 50, 100, 150, and 200 for shared data and both models.' \
    '  3. Submit pilot post-processing after both model arrays succeed.' \
    'Run with EXP4_CONFIRM_PILOT_SUBMIT=YES to submit.'
  exit 0
fi
if [[ $# -ne 0 ]]; then
  printf 'Usage: bash hpc/submit_pilot.sh [--dry-run]\n' >&2
  exit 2
fi
if [[ "${EXP4_CONFIRM_PILOT_SUBMIT:-}" != "YES" ]]; then
  printf '%s\n' \
    'Submission stopped before file creation or sbatch.' \
    'Review with: bash hpc/submit_pilot.sh --dry-run' \
    'Submit with: EXP4_CONFIRM_PILOT_SUBMIT=YES bash hpc/submit_pilot.sh' >&2
  exit 2
fi

mkdir -p logs/data logs/gamma logs/constant logs/pilot \
  shared_data results_raw/gamma results_raw/constant \
  results/pilot_combined/gamma results/pilot_combined/constant \
  results/pilot_comparison figures/pilot results

module --force purge
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="$HOME/packages-R4.1"
command -v sbatch >/dev/null
Rscript --vanilla -e 'stopifnot(requireNamespace("pomp", quietly = TRUE)); cat(R.version.string, "\npomp ", as.character(packageVersion("pomp")), "\n", sep = "")'
Rscript code/00_create_paramlist.R results/paramlist.csv

clean_job_id() { echo "${1%%;*}"; }

data_raw=$(sbatch --parsable --array=1,50,100,150,200 hpc/01_generate_data_array.sh)
data_job=$(clean_job_id "$data_raw")

gamma_raw=$(sbatch --parsable --array=1,50,100,150,200 --dependency=afterok:"$data_job" hpc/02_gamma_array.sh)
gamma_job=$(clean_job_id "$gamma_raw")

constant_raw=$(sbatch --parsable --array=1,50,100,150,200 --dependency=afterok:"$data_job" hpc/03_constant_array.sh)
constant_job=$(clean_job_id "$constant_raw")

post_raw=$(sbatch --parsable --dependency=afterok:"$gamma_job":"$constant_job" hpc/05_postprocess_pilot.sh)
post_job=$(clean_job_id "$post_raw")

cat > results/pilot_job_ids.env <<EOF
DATA_JOB_ID=$data_job
GAMMA_JOB_ID=$gamma_job
CONSTANT_JOB_ID=$constant_job
POST_JOB_ID=$post_job
EOF

printf 'Pilot submitted.\nData job: %s\nGamma job: %s\nConstant job: %s\nPilot postprocess: %s\nJob IDs saved to results/pilot_job_ids.env\n' \
  "$data_job" "$gamma_job" "$constant_job" "$post_job"

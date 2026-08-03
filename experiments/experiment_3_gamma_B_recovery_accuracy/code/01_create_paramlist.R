# ============================================================
# Create paramlist.csv for the 200-task simulation study
#
# Each array task receives:
#   - one seed for generating its accepted simulated epidemic data set;
#   - separate, reproducible random-number streams for mif2,
#     likelihood evaluation, and final filtering.
#
# Usage:
# Rscript code/01_create_paramlist.R results/paramlist.csv
#
# Optional:
# Rscript code/01_create_paramlist.R another_filename.csv
# ============================================================

options(
  stringsAsFactors = FALSE
)

args <- commandArgs(
  trailingOnly = TRUE
)

if (length(args) >= 1) {
  output_file <- args[[1]]
} else {
  output_file <- file.path("results", "paramlist.csv")
}

output_parent <- dirname(output_file)

if (!identical(output_parent, ".")) {
  dir.create(
    output_parent,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

n_tasks <- 200L

task_id <- seq_len(
  n_tasks
)

paramlist <- data.frame(
  task_id = task_id,

  # Task 1 uses 1001, task 2 uses 1002, ..., task 200 uses 1200.
  simulation_seed = 1000L + task_id,

  # Separate random-number streams for fitting and evaluation.
  mif_seed_base = 20260728L + 1000000L * task_id,

  evaluation_seed_base =
    20260800L + 1000000L * task_id,

  final_pf_seed =
    20260900L + 1000000L * task_id,

  stringsAsFactors = FALSE
)

write.csv(
  paramlist,
  output_file,
  row.names = FALSE
)

cat(
  "Created ",
  output_file,
  " with ",
  nrow(paramlist),
  " task-specific seed rows.\n",
  sep = ""
)
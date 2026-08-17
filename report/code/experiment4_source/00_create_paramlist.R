# ============================================================
# Create the model-specific seed table for Experiment 4
# ============================================================

options(stringsAsFactors = FALSE)
source("config/experiment_config.R")

args <- commandArgs(trailingOnly = TRUE)
output_file <- if (length(args) >= 1L) args[[1]] else file.path("results", "paramlist.csv")

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

task_id <- seq_len(experiment_config$n_tasks)

constant_mif_seed_base <- 20260628L + 1000000L * task_id

paramlist <- data.frame(
  task_id = task_id,
  simulation_seed = 1000L + task_id,

  # These reproduce the organized Gamma-model experiment's seed design.
  gamma_mif_seed_base = 20260728L + 1000000L * task_id,
  gamma_evaluation_seed_base = 20260800L + 1000000L * task_id,
  gamma_final_pf_seed = 20260900L + 1000000L * task_id,

  # These reproduce the earlier constant-B experiment's seed design.
  constant_mif_seed_base = constant_mif_seed_base,
  constant_evaluation_seed_base = constant_mif_seed_base + 1000L,
  constant_final_pf_seed = constant_mif_seed_base + 900000L,

  stringsAsFactors = FALSE
)

write.csv(paramlist, output_file, row.names = FALSE)
cat("Created ", output_file, " with ", nrow(paramlist), " rows.\n", sep = "")

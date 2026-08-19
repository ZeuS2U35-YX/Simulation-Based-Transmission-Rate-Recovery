# ============================================================
# Read-only validation of one completed Experiment 5 task
#
# Usage:
# Rscript code/00_validate_completed_bspline_task.R \
#   <task_id> <task_directory> [input_manifest]
# ============================================================

options(stringsAsFactors = FALSE)

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
)
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "path_helpers.R"))

experiment_directory <- get_experiment_directory()
source(file.path(experiment_directory, "config", "experiment_config.R"))
source(file.path(experiment_directory, "code", "io_helpers.R"))
source(file.path(experiment_directory, "code", "model_components.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript code/00_validate_completed_bspline_task.R ",
    "<task_id> <task_directory> [input_manifest]"
  )
}

task_id <- suppressWarnings(as.integer(args[[1]]))
if (is.na(task_id) || task_id < 1L || task_id > experiment_config$n_tasks) {
  stop("task_id must be between 1 and ", experiment_config$n_tasks, ".")
}
task_dir <- args[[2]]
manifest_file <- if (length(args) >= 3L) args[[3]] else {
  file.path(experiment_directory, "results", "paired_input_manifest.csv")
}

required_files <- c(
  "COMPLETE", "start_selection_audit.csv", "likelihood_evaluations.csv",
  "best_fit_summary.csv", "B_path.csv", "best_mif2.rds", "run_config.csv"
)
missing <- required_files[!file.exists(file.path(task_dir, required_files))]
if (length(missing) > 0L) {
  stop("Task directory is incomplete. Missing: ", paste(missing, collapse = ", "))
}

manifest <- read_csv(manifest_file)
require_exact_task_ids(
  manifest$task_id, experiment_config$n_tasks, "Paired input manifest"
)
manifest_row <- manifest[manifest$task_id == task_id, , drop = FALSE]

audit <- read_csv(file.path(task_dir, "start_selection_audit.csv"))
evaluations <- read_csv(file.path(task_dir, "likelihood_evaluations.csv"))
best <- read_csv(file.path(task_dir, "best_fit_summary.csv"))
path <- read_csv(file.path(task_dir, "B_path.csv"))
run_config <- read_csv(file.path(task_dir, "run_config.csv"))

if (nrow(audit) != experiment_config$n_start) {
  stop("Expected ", experiment_config$n_start, " internal-start rows.")
}
if (nrow(evaluations) !=
    experiment_config$n_start * experiment_config$n_pf_evals) {
  stop("Likelihood-evaluation row count is not the production count.")
}
if (nrow(best) != 1L || best$status[[1]] != "success") {
  stop("Task must contain exactly one successful best-fit row.")
}

expected_times <- seq(
  experiment_config$observation_dt,
  experiment_config$n_weeks,
  by = experiment_config$observation_dt
)
if (nrow(path) != 70L || is.unsorted(path$week, strictly = TRUE) ||
    anyDuplicated(path$week) ||
    max(abs(path$week - expected_times)) > 1e-12) {
  stop("B_path.csv does not contain the configured 70 observation times.")
}
if (!all(path$path_semantics == "deterministic_selected_bspline_trajectory")) {
  stop("B_path.csv has unexpected path semantics.")
}
expected_true <- true_B_at_times(path$week)
legacy_true <- ifelse(
  path$week < experiment_config$truth[["t_switch"]],
  experiment_config$truth[["Beta_high"]],
  experiment_config$truth[["Beta_low"]]
)
uses_current_truth <- max(abs(path$B_true - expected_true)) <= 1e-12
uses_legacy_week5_truth <- max(abs(path$B_true - legacy_true)) <= 1e-12
if (!uses_current_truth && !uses_legacy_week5_truth) {
  stop("B_path.csv does not match the current or legacy true B(t) encoding.")
}

task_tables <- list(audit = audit, evaluations = evaluations, best = best, path = path)
for (table_name in names(task_tables)) {
  x <- task_tables[[table_name]]
  if (!all(x$task_id == task_id) || !all(x$model == "bspline_B")) {
    stop(table_name, " contains a task_id or model mismatch.")
  }
  if (!all(x$simulation_seed == manifest_row$simulation_seed[[1]]) ||
      !all(x$observed_data_md5 == manifest_row$observed_data_md5[[1]])) {
    stop(table_name, " does not match the paired-input manifest.")
  }
}

if (!all(audit$Nmif == experiment_config$Nmif) ||
    best$Nmif[[1]] != experiment_config$Nmif) {
  stop("Task does not use the configured Nmif.")
}
if (best$selection_rule[[1]] != "maximum_independently_evaluated_logLik") {
  stop("Best-fit row has the wrong selection rule.")
}

eligible <- audit[
  as.logical(audit$eligible_for_selection) &
    is.finite(audit$independently_evaluated_logLik),
  ,
  drop = FALSE
]
if (nrow(eligible) == 0L) stop("Task has no eligible internal start.")
expected_best <- eligible[
  which.max(eligible$independently_evaluated_logLik), , drop = FALSE
]
if (best$best_run[[1]] != expected_best$run[[1]] ||
    !isTRUE(all.equal(
      as.numeric(best$logLik[[1]]),
      as.numeric(expected_best$independently_evaluated_logLik[[1]]),
      tolerance = 1e-12
    ))) {
  stop("Best-fit row is not the maximum independently evaluated candidate.")
}

evaluation_counts <- table(evaluations$run)
if (length(evaluation_counts) != experiment_config$n_start ||
    any(evaluation_counts != experiment_config$n_pf_evals) ||
    anyDuplicated(evaluations$evaluation_seed)) {
  stop("Task lacks distinct complete PF-evaluation slots per start.")
}

rds_files <- list.files(task_dir, pattern = "\\.rds$", full.names = FALSE)
if (!identical(rds_files, "best_mif2.rds")) {
  stop("Task must retain only best_mif2.rds as its fitted object.")
}

required_settings <- c(
  Nmif = experiment_config$Nmif,
  Np_mif = experiment_config$Np_mif,
  n_internal_starts = experiment_config$n_start,
  Np_eval = experiment_config$Np_eval,
  n_pf_evals_per_start = experiment_config$n_pf_evals
)
if (!all(c("setting", "value") %in% names(run_config))) {
  stop("run_config.csv lacks setting/value columns.")
}
for (setting in names(required_settings)) {
  actual <- run_config$value[run_config$setting == setting]
  actual_numeric <- suppressWarnings(as.numeric(actual))
  if (length(actual) != 1L ||
      length(actual_numeric) != 1L || !is.finite(actual_numeric) ||
      actual_numeric != required_settings[[setting]]) {
    stop("run_config.csv has the wrong production value for ", setting, ".")
  }
}

cat(
  "Validated completed production task ", task_id,
  " at ", normalizePath(task_dir, mustWork = TRUE), ".\n",
  if (uses_legacy_week5_truth) {
    paste0(
      "The raw task uses the legacy B(5)=2 label; combination will ",
      "normalize B(5) to 4 without changing the fitted B-spline.\n"
    )
  } else {
    ""
  },
  sep = ""
)

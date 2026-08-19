# ============================================================
# Validate the single set of 200 paired inputs for Experiment 5
#
# This script never simulates or copies observed data. It reads the
# accepted Experiment 4 data in place and writes only a 200-row manifest.
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

source_experiment_directory <- file.path(
  dirname(experiment_directory), experiment_config$source_experiment_id
)
args <- commandArgs(trailingOnly = TRUE)
shared_data_root <- if (length(args) >= 1L) args[[1]] else {
  file.path(source_experiment_directory, "shared_data")
}
gamma_best_file <- if (length(args) >= 2L) args[[2]] else {
  file.path(
    source_experiment_directory, "results", "combined", "gamma",
    "combined_best_fit_summary.csv"
  )
}
output_manifest <- if (length(args) >= 3L) args[[3]] else {
  file.path(experiment_directory, "results", "paired_input_manifest.csv")
}

task_directories <- list.dirs(shared_data_root, recursive = FALSE, full.names = FALSE)
task_directories <- sort(task_directories[grepl("^task_[0-9]{3}$", task_directories)])
expected_directories <- sprintf("task_%03d", seq_len(experiment_config$n_tasks))
if (!identical(task_directories, expected_directories)) {
  stop(
    "The selected Experiment 4 shared-data root must contain exactly task_001 ",
    "through task_200, with no missing or additional task directories."
  )
}

gamma_best <- read_csv(gamma_best_file)
required_gamma_columns <- c(
  "task_id", "simulation_seed", "observed_data_md5", "model",
  "Nmif", "logLik", "status"
)
missing_gamma_columns <- setdiff(required_gamma_columns, names(gamma_best))
if (length(missing_gamma_columns) > 0L) {
  stop(
    "Gamma best-fit table is missing: ",
    paste(missing_gamma_columns, collapse = ", ")
  )
}
require_exact_task_ids(
  gamma_best$task_id, experiment_config$n_tasks, "Gamma best-fit table"
)
if (!all(gamma_best$model == "gamma_noise") ||
    !all(gamma_best$status == "success")) {
  stop("All 200 Gamma-noise reference rows must have status='success'.")
}
if (!all(gamma_best$Nmif == experiment_config$Nmif)) {
  stop(
    "All 200 Gamma-noise reference rows must use Nmif=",
    experiment_config$Nmif, "."
  )
}

validate_task <- function(task_id) {
  shared <- read_exp4_shared_task(shared_data_root, task_id)
  metadata <- shared$metadata
  gamma_row <- gamma_best[gamma_best$task_id == task_id, , drop = FALSE]

  simulation_seed <- as.integer(metadata$simulation_seed[[1]])
  simulation_attempt <- as.integer(metadata$simulation_attempt[[1]])
  if (simulation_seed != as.integer(gamma_row$simulation_seed[[1]])) {
    stop("Simulation-seed mismatch between shared data and Gamma task ", task_id)
  }
  if (!identical(
    shared$observed_data_md5,
    as.character(gamma_row$observed_data_md5[[1]])
  )) {
    stop("Observed-data checksum mismatch between shared data and Gamma task ", task_id)
  }

  data.frame(
    task_id = task_id,
    simulation_seed = simulation_seed,
    simulation_attempt = simulation_attempt,
    observed_data_md5 = shared$observed_data_md5,
    simulated_data_md5 = shared$simulated_data_md5,
    acceptance_threshold = shared$acceptance_threshold,
    recomputed_max_H = shared$recomputed_max_H,
    accepted = as.logical(metadata$accepted[[1]]),
    source_experiment = experiment_config$source_experiment_id,
    source_task_directory = sprintf("task_%03d", task_id),
    gamma_model = as.character(gamma_row$model[[1]]),
    gamma_Nmif = as.integer(gamma_row$Nmif[[1]]),
    gamma_logLik = as.numeric(gamma_row$logLik[[1]]),
    stringsAsFactors = FALSE
  )
}

validation_cores <- read_positive_integer_env("EXP5_VALIDATION_CORES", 1L)
if (validation_cores > 1L && .Platform$OS.type != "windows") {
  manifest_rows <- parallel::mclapply(
    seq_len(experiment_config$n_tasks),
    validate_task,
    mc.cores = min(validation_cores, experiment_config$n_tasks)
  )
} else {
  manifest_rows <- lapply(
    seq_len(experiment_config$n_tasks),
    validate_task
  )
}

manifest <- do.call(rbind, manifest_rows)
require_exact_task_ids(
  manifest$task_id, experiment_config$n_tasks, "Paired input manifest"
)
if (nrow(manifest) != 200L || anyDuplicated(manifest$observed_data_md5)) {
  stop(
    "The manifest must contain exactly 200 distinct accepted observed-data ",
    "checksums."
  )
}
if (!all(manifest$acceptance_threshold == 20) ||
    !all(manifest$recomputed_max_H > manifest$acceptance_threshold)) {
  stop(
    "Every manifest row must independently satisfy max(H) > 20 with ",
    "acceptance_threshold=20."
  )
}

dir.create(dirname(output_manifest), recursive = TRUE, showWarnings = FALSE)
write.csv(manifest, output_manifest, row.names = FALSE)
cat(
  "Validated one shared set of exactly 200 accepted Experiment 4 data sets.\n",
  "No observed data were generated or copied.\n",
  "Manifest: ", output_manifest, "\n",
  sep = ""
)

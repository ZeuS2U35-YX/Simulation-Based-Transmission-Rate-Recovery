#!/usr/bin/env Rscript

# Recompute the week-8 sensitivity comparison from tracked Gamma-noise and
# B-spline trajectories. This analysis truncates the common 10-week grid at
# week 8; it does not refit either model or modify the canonical full-window
# results.

options(stringsAsFactors = FALSE, digits = 17)

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Could not determine the script path.")
  }
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

script_path <- get_script_path()
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)

output_dir <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = FALSE)
} else {
  file.path(
    repo_root,
    "experiments", "experiment_5_bspline_B_recovery", "results",
    "comparison_three_models"
  )
}

cutoff_week <- if (length(args) >= 2L) as.numeric(args[[2]]) else 8
if (length(cutoff_week) != 1L || !is.finite(cutoff_week) ||
    cutoff_week <= 0 || cutoff_week > 10) {
  stop("cutoff_week must be one finite number in (0, 10].")
}

read_required_csv <- function(path) {
  if (!file.exists(path)) stop("Required file does not exist: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

manifest_path <- file.path(
  repo_root,
  "experiments", "experiment_5_bspline_B_recovery", "results",
  "paired_input_manifest.csv"
)
gamma_path <- file.path(
  repo_root,
  "experiments", "experiment_4_nmif600_model_comparison", "results",
  "combined", "gamma", "combined_B_paths.csv"
)
bspline_path <- file.path(
  repo_root,
  "experiments", "experiment_5_bspline_B_recovery", "results",
  "combined", "bspline", "combined_B_paths.csv"
)

manifest <- read_required_csv(manifest_path)
gamma <- read_required_csv(gamma_path)
bspline <- read_required_csv(bspline_path)

manifest_required <- c(
  "task_id", "simulation_seed", "observed_data_md5", "acceptance_threshold",
  "recomputed_max_H", "accepted"
)
missing_manifest <- setdiff(manifest_required, names(manifest))
if (length(missing_manifest) > 0L) {
  stop(
    "The paired manifest is missing: ",
    paste(missing_manifest, collapse = ", ")
  )
}
expected_ids <- seq_len(200L)
actual_ids <- sort(as.integer(manifest$task_id))
if (nrow(manifest) != 200L || anyDuplicated(actual_ids) ||
    !identical(actual_ids, expected_ids)) {
  stop("The paired input manifest must contain task IDs 1 through 200 once each.")
}
if (!all(manifest$accepted) || !all(manifest$acceptance_threshold == 20) ||
    !all(manifest$recomputed_max_H > manifest$acceptance_threshold)) {
  stop("The paired manifest does not satisfy the accepted-outbreak contract.")
}

validate_paths <- function(
  paths, model_name, required_semantics, required_seed_column = NULL
) {
  required <- c(
    "task_id", "simulation_seed", "observed_data_md5", "model", "week",
    "B_estimate", "B_true", "path_semantics"
  )
  missing <- setdiff(required, names(paths))
  if (length(missing) > 0L) {
    stop(model_name, " paths are missing: ", paste(missing, collapse = ", "))
  }
  if (nrow(paths) != 14000L || !all(paths$model == model_name) ||
      !all(paths$path_semantics == required_semantics) ||
      any(!is.finite(paths$week)) || any(!is.finite(paths$B_estimate))) {
    stop(model_name, " paths fail the row, model, semantics, or finite-value contract.")
  }
  if (!is.null(required_seed_column) &&
      !required_seed_column %in% names(paths)) {
    stop(model_name, " paths are missing ", required_seed_column, ".")
  }

  split_paths <- split(paths, paths$task_id)
  split_ids <- sort(as.integer(names(split_paths)))
  if (!identical(split_ids, expected_ids)) {
    stop(model_name, " paths must contain task IDs 1 through 200.")
  }

  expected_times <- seq(1 / 7, 10, by = 1 / 7)
  checked <- lapply(names(split_paths), function(task_name) {
    task_id <- as.integer(task_name)
    x <- split_paths[[task_name]]
    x <- x[order(x$week), , drop = FALSE]
    manifest_row <- manifest[manifest$task_id == task_id, , drop = FALSE]
    if (nrow(x) != 70L || anyDuplicated(x$week) ||
        max(abs(x$week - expected_times)) > 1e-12) {
      stop(model_name, " task ", task_id, " does not use the 70-time grid.")
    }
    if (length(unique(x$simulation_seed)) != 1L ||
        x$simulation_seed[[1]] != manifest_row$simulation_seed[[1]] ||
        length(unique(x$observed_data_md5)) != 1L ||
        x$observed_data_md5[[1]] != manifest_row$observed_data_md5[[1]]) {
      stop(model_name, " task ", task_id, " does not match the paired manifest.")
    }
    if (!is.null(required_seed_column) &&
        length(unique(x[[required_seed_column]])) != 1L) {
      stop(model_name, " task ", task_id, " does not have one stable path seed.")
    }
    expected_truth <- ifelse(x$week <= 5, 4, 2)
    if (any(!is.finite(x$B_true)) ||
        max(abs(x$B_true - expected_truth)) > 1e-12) {
      stop(model_name, " task ", task_id, " has an unexpected truth convention.")
    }
    x
  })

  out <- do.call(rbind, checked)
  rownames(out) <- NULL
  out
}

gamma <- validate_paths(
  gamma,
  "gamma_noise",
  "particle_filtering_mean",
  "filter_seed"
)
bspline <- validate_paths(
  bspline,
  "bspline_B",
  "deterministic_selected_bspline_trajectory"
)

path_key <- function(x) {
  x <- x[order(x$task_id, x$week), , drop = FALSE]
  paste(x$task_id, sprintf("%.14f", x$week), sep = "|")
}
if (!identical(path_key(gamma), path_key(bspline))) {
  stop("Gamma-noise and B-spline paths do not use the same task/time pairs.")
}

calculate_truncated_rmse <- function(paths) {
  paths <- paths[paths$week <= cutoff_week + 1e-12, , drop = FALSE]
  split_paths <- split(paths, paths$task_id)
  values <- vapply(split_paths, function(x) {
    sqrt(mean((x$B_estimate - x$B_true)^2))
  }, numeric(1))
  counts <- vapply(split_paths, nrow, integer(1))
  if (length(values) != 200L || length(unique(counts)) != 1L) {
    stop("The truncated comparison does not contain one common grid for 200 tasks.")
  }
  list(values = values, n_times = unique(counts))
}

gamma_rmse <- calculate_truncated_rmse(gamma)
bspline_rmse <- calculate_truncated_rmse(bspline)
if (!identical(names(gamma_rmse$values), names(bspline_rmse$values)) ||
    gamma_rmse$n_times != bspline_rmse$n_times) {
  stop("The two truncated RMSE vectors are not paired.")
}

delta <- gamma_rmse$values - bspline_rmse$values
winner <- ifelse(
  delta < 0,
  "gamma_noise",
  ifelse(delta > 0, "bspline_B", "tie")
)

task_metrics <- data.frame(
  task_id = as.integer(names(gamma_rmse$values)),
  cutoff_week = cutoff_week,
  n_observation_times = gamma_rmse$n_times,
  gamma_RMSE = as.numeric(gamma_rmse$values),
  bspline_RMSE = as.numeric(bspline_rmse$values),
  delta_RMSE_gamma_minus_bspline = as.numeric(delta),
  winner_lower_RMSE = winner,
  stringsAsFactors = FALSE
)
task_metrics <- task_metrics[order(task_metrics$task_id), , drop = FALSE]

summary <- data.frame(
  cutoff_week = cutoff_week,
  n_tasks = nrow(task_metrics),
  n_observation_times_per_task = gamma_rmse$n_times,
  gamma_mean_RMSE = mean(task_metrics$gamma_RMSE),
  bspline_mean_RMSE = mean(task_metrics$bspline_RMSE),
  gamma_lower_RMSE_count = sum(winner == "gamma_noise"),
  bspline_lower_RMSE_count = sum(winner == "bspline_B"),
  tie_count = sum(winner == "tie"),
  truth_convention = "B(t)=4 through week 5; B(t)=2 after week 5",
  scope = "paired latent-path recovery truncated at week 8",
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  task_metrics,
  file.path(output_dir, "week8_sensitivity_task_metrics.csv"),
  row.names = FALSE
)
write.csv(
  summary,
  file.path(output_dir, "week8_sensitivity_summary.csv"),
  row.names = FALSE
)

cat(
  sprintf(
    paste0(
      "Week-8 sensitivity complete: mean RMSE %.6f (Gamma-noise) vs ",
      "%.6f (B-spline); Gamma-noise lower in %d/%d paired tasks.\n"
    ),
    summary$gamma_mean_RMSE,
    summary$bspline_mean_RMSE,
    summary$gamma_lower_RMSE_count,
    summary$n_tasks
  )
)

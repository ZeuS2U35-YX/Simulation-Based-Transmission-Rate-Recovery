# ============================================================
# Fair path-level comparison of Gamma-noise, B-spline, and
# constant-B recovery on the same 200 accepted simulations.
#
# This script does not run MIF2 or particle filtering. It reads
# existing combined trajectories and recomputes every recovery
# metric after applying the common B(5) = 4 truth definition.
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
gamma_dir <- if (length(args) >= 1L) args[[1]] else {
  file.path(source_experiment_directory, "results", "combined", "gamma")
}
bspline_dir <- if (length(args) >= 2L) args[[2]] else {
  file.path(experiment_directory, "results", "combined", "bspline")
}
constant_dir <- if (length(args) >= 3L) args[[3]] else {
  file.path(source_experiment_directory, "results", "combined", "constant")
}
output_dir <- if (length(args) >= 4L) args[[4]] else {
  file.path(experiment_directory, "results", "comparison_three_models")
}
manifest_file <- if (length(args) >= 5L) args[[5]] else {
  file.path(experiment_directory, "results", "paired_input_manifest.csv")
}

manifest <- read_csv(manifest_file)
require_exact_task_ids(
  manifest$task_id, experiment_config$n_tasks, "Paired input manifest"
)
if (nrow(manifest) != 200L) {
  stop("Paired input manifest must contain exactly 200 rows.")
}

read_model_outputs <- function(directory) {
  list(
    best = read_csv(file.path(directory, "combined_best_fit_summary.csv")),
    paths = read_csv(file.path(directory, "combined_B_paths.csv"))
  )
}

gamma <- read_model_outputs(gamma_dir)
bspline <- read_model_outputs(bspline_dir)
constant <- read_model_outputs(constant_dir)

model_contract <- data.frame(
  model = c("gamma_noise", "bspline_B", "constant_B"),
  display_model = c("Gamma-noise", "B-spline", "Constant-B"),
  suffix = c("gamma", "bspline", "constant"),
  path_semantics = c(
    "particle_filtering_mean",
    "deterministic_selected_bspline_trajectory",
    "repeated_static_estimate"
  ),
  stringsAsFactors = FALSE
)

validate_best <- function(best, model_name) {
  required <- c(
    "task_id", "simulation_seed", "observed_data_md5", "model", "status"
  )
  missing <- setdiff(required, names(best))
  if (length(missing) > 0L) {
    stop(model_name, " best summaries are missing: ", paste(missing, collapse = ", "))
  }
  require_exact_task_ids(
    best$task_id, experiment_config$n_tasks,
    paste(model_name, "best summaries")
  )
  if (nrow(best) != 200L || !all(best$model == model_name) ||
      !all(best$status == "success")) {
    stop(
      model_name,
      " must have exactly 200 successful best-fit rows with the expected model label."
    )
  }
  manifest_index <- match(best$task_id, manifest$task_id)
  if (!all(best$simulation_seed == manifest$simulation_seed[manifest_index]) ||
      !all(best$observed_data_md5 ==
        manifest$observed_data_md5[manifest_index])) {
    stop(model_name, " best summaries do not match the paired-input manifest.")
  }
  invisible(TRUE)
}

validate_best(gamma$best, "gamma_noise")
validate_best(bspline$best, "bspline_B")
validate_best(constant$best, "constant_B")

if (!all(c("filter_seed", "path_semantics") %in% names(gamma$paths)) ||
    !all(gamma$paths$path_semantics == "particle_filtering_mean")) {
  stop(
    "Gamma main-analysis paths must all be particle filtering means; ",
    "sampled latent trajectories are not permitted."
  )
}
if (!"path_semantics" %in% names(bspline$paths) ||
    !all(
      bspline$paths$path_semantics ==
        "deterministic_selected_bspline_trajectory"
    )) {
  stop("B-spline paths do not have deterministic selected-trajectory semantics.")
}

expected_times <- seq(
  experiment_config$observation_dt,
  experiment_config$n_weeks,
  by = experiment_config$observation_dt
)
if (length(expected_times) != 70L) {
  stop("The configured observation grid must contain exactly 70 times.")
}

validate_paths_and_calculate_metrics <- function(
  paths, model_name, display_model, required_semantics = NULL,
  require_static = FALSE
) {
  required <- c(
    "task_id", "simulation_seed", "observed_data_md5", "model",
    "week", "B_estimate", "B_true"
  )
  missing <- setdiff(required, names(paths))
  if (length(missing) > 0L) {
    stop(model_name, " paths are missing: ", paste(missing, collapse = ", "))
  }
  if (nrow(paths) != 14000L || !all(paths$model == model_name)) {
    stop(model_name, " must have exactly 14,000 path rows with the expected model label.")
  }
  if (!is.null(required_semantics) &&
      (!"path_semantics" %in% names(paths) ||
       !all(paths$path_semantics == required_semantics))) {
    stop(model_name, " paths do not have the required path semantics.")
  }

  split_paths <- split(paths, paths$task_id)
  if (!identical(sort(as.integer(names(split_paths))), seq_len(200L))) {
    stop(model_name, " paths must contain exactly task IDs 1 through 200.")
  }

  normalized_list <- vector("list", length(split_paths))
  metrics_list <- vector("list", length(split_paths))
  names(normalized_list) <- names(split_paths)
  names(metrics_list) <- names(split_paths)

  for (task_name in names(split_paths)) {
    x <- split_paths[[task_name]]
    task_id <- as.integer(task_name)
    x <- x[order(x$week), , drop = FALSE]
    manifest_row <- manifest[manifest$task_id == task_id, , drop = FALSE]

    if (nrow(x) != 70L || is.unsorted(x$week, strictly = TRUE) ||
        anyDuplicated(x$week) ||
        max(abs(x$week - expected_times)) > 1e-12) {
      stop(
        model_name, " task ", task_id,
        " does not contain exactly the 70 configured observation times."
      )
    }
    if (length(unique(x$simulation_seed)) != 1L ||
        as.integer(x$simulation_seed[[1]]) !=
          as.integer(manifest_row$simulation_seed[[1]]) ||
        length(unique(x$observed_data_md5)) != 1L ||
        as.character(x$observed_data_md5[[1]]) !=
          as.character(manifest_row$observed_data_md5[[1]])) {
      stop(
        model_name, " task ", task_id,
        " path provenance does not match the paired-input manifest."
      )
    }
    if (!all(is.finite(x$B_estimate))) {
      stop(model_name, " task ", task_id, " has non-finite B estimates.")
    }
    if (require_static && diff(range(x$B_estimate)) > 1e-12) {
      stop(
        model_name, " task ", task_id,
        " is not a repeated static estimate across the 70 observation times."
      )
    }

    # Apply one common truth definition in memory. The input combined files
    # are never modified by this script.
    x$B_true <- ifelse(
      x$week <= experiment_config$truth[["t_switch"]],
      experiment_config$truth[["Beta_high"]],
      experiment_config$truth[["Beta_low"]]
    )
    error <- x$B_estimate - x$B_true
    through <- x$week <= experiment_config$truth[["t_switch"]]
    after <- x$week > experiment_config$truth[["t_switch"]]
    if (sum(through) != 35L || sum(after) != 35L) {
      stop(model_name, " task ", task_id, " does not split into 35 + 35 times.")
    }

    metrics_list[[task_name]] <- data.frame(
      task_id = task_id,
      simulation_seed = as.integer(x$simulation_seed[[1]]),
      observed_data_md5 = as.character(x$observed_data_md5[[1]]),
      model = model_name,
      display_model = display_model,
      n_observation_times = nrow(x),
      RSS = sum(error^2),
      RMSE = sqrt(mean(error^2)),
      mean_error = mean(error),
      AOB = abs(mean(error)),
      mean_error_through_5 = mean(error[through]),
      mean_error_after_5 = mean(error[after]),
      stringsAsFactors = FALSE
    )
    normalized_list[[task_name]] <- x
  }

  metrics <- do.call(rbind, metrics_list)
  normalized_paths <- do.call(rbind, normalized_list)
  rownames(metrics) <- NULL
  rownames(normalized_paths) <- NULL
  require_exact_task_ids(
    metrics$task_id, experiment_config$n_tasks,
    paste(model_name, "task-level recovery metrics")
  )
  list(metrics = metrics, normalized_paths = normalized_paths)
}

gamma_checked <- validate_paths_and_calculate_metrics(
  gamma$paths, "gamma_noise", "Gamma-noise",
  required_semantics = "particle_filtering_mean"
)
bspline_checked <- validate_paths_and_calculate_metrics(
  bspline$paths, "bspline_B", "B-spline",
  required_semantics = "deterministic_selected_bspline_trajectory"
)
constant_checked <- validate_paths_and_calculate_metrics(
  constant$paths, "constant_B", "Constant-B", require_static = TRUE
)

path_key <- function(paths) {
  paths <- paths[order(paths$task_id, paths$week), , drop = FALSE]
  paste(paths$task_id, sprintf("%.14f", paths$week), sep = "|")
}
gamma_key <- path_key(gamma_checked$normalized_paths)
if (!identical(gamma_key, path_key(bspline_checked$normalized_paths)) ||
    !identical(gamma_key, path_key(constant_checked$normalized_paths))) {
  stop("The three models do not use exactly the same task/time pairs.")
}

metrics_long <- rbind(
  gamma_checked$metrics,
  bspline_checked$metrics,
  constant_checked$metrics
)
metrics_long$model <- factor(
  metrics_long$model,
  levels = model_contract$model
)
metrics_long <- metrics_long[order(metrics_long$model, metrics_long$task_id), ]
metrics_long$model <- as.character(metrics_long$model)
rownames(metrics_long) <- NULL

metric_columns <- c(
  "RSS", "RMSE", "mean_error", "AOB",
  "mean_error_through_5", "mean_error_after_5"
)
paired_wide <- manifest[, c("task_id", "simulation_seed", "observed_data_md5")]
for (i in seq_len(nrow(model_contract))) {
  model_name <- model_contract$model[[i]]
  suffix <- model_contract$suffix[[i]]
  part <- metrics_long[
    metrics_long$model == model_name,
    c("task_id", metric_columns),
    drop = FALSE
  ]
  names(part)[names(part) != "task_id"] <- paste0(
    names(part)[names(part) != "task_id"], "_", suffix
  )
  paired_wide <- merge(paired_wide, part, by = "task_id", sort = TRUE)
}
require_exact_task_ids(
  paired_wide$task_id, experiment_config$n_tasks,
  "Three-model paired metric table"
)
if (nrow(paired_wide) != 200L) {
  stop("The paired three-model metric table must contain exactly 200 rows.")
}

summary_rows <- lapply(split(metrics_long, metrics_long$model), function(x) {
  data.frame(
    model = x$model[[1]],
    display_model = x$display_model[[1]],
    n_tasks = nrow(x),
    mean_RSS = mean(x$RSS),
    mean_RMSE = mean(x$RMSE),
    median_RMSE = median(x$RMSE),
    RMSE_Q1 = as.numeric(quantile(x$RMSE, 0.25, names = FALSE)),
    RMSE_Q3 = as.numeric(quantile(x$RMSE, 0.75, names = FALSE)),
    RMSE_IQR = IQR(x$RMSE),
    mean_AOB = mean(x$AOB),
    mean_error_overall = mean(x$mean_error),
    mean_error_through_5 = mean(x$mean_error_through_5),
    mean_error_after_5 = mean(x$mean_error_after_5),
    stringsAsFactors = FALSE
  )
})
overall_summary <- do.call(rbind, summary_rows)
overall_summary <- merge(
  model_contract[, c("model", "display_model")],
  overall_summary[, setdiff(names(overall_summary), "display_model")],
  by = "model", sort = FALSE
)
overall_summary <- overall_summary[
  match(model_contract$model, overall_summary$model), , drop = FALSE
]
rownames(overall_summary) <- NULL

pair_contract <- data.frame(
  comparison = c(
    "Gamma-noise vs Constant-B",
    "B-spline vs Constant-B",
    "B-spline vs Gamma-noise"
  ),
  model_a = c("gamma_noise", "bspline_B", "bspline_B"),
  model_b = c("constant_B", "constant_B", "gamma_noise"),
  suffix_a = c("gamma", "bspline", "bspline"),
  suffix_b = c("constant", "constant", "gamma"),
  stringsAsFactors = FALSE
)

paired_difference_list <- vector("list", nrow(pair_contract))
pairwise_summary_list <- vector("list", nrow(pair_contract))
for (i in seq_len(nrow(pair_contract))) {
  contract_row <- pair_contract[i, , drop = FALSE]
  rmse_a <- paired_wide[[paste0("RMSE_", contract_row$suffix_a)]]
  rmse_b <- paired_wide[[paste0("RMSE_", contract_row$suffix_b)]]
  difference <- rmse_a - rmse_b
  tolerance <- sqrt(.Machine$double.eps)
  winner <- ifelse(
    abs(difference) <= tolerance, "tie",
    ifelse(difference < 0, contract_row$model_a, contract_row$model_b)
  )
  paired_difference_list[[i]] <- data.frame(
    task_id = paired_wide$task_id,
    simulation_seed = paired_wide$simulation_seed,
    observed_data_md5 = paired_wide$observed_data_md5,
    comparison = contract_row$comparison,
    model_a = contract_row$model_a,
    model_b = contract_row$model_b,
    RMSE_model_a = rmse_a,
    RMSE_model_b = rmse_b,
    paired_RMSE_difference_model_a_minus_model_b = difference,
    winner_lower_RMSE = winner,
    stringsAsFactors = FALSE
  )
  pairwise_summary_list[[i]] <- data.frame(
    comparison = contract_row$comparison,
    model_a = contract_row$model_a,
    model_b = contract_row$model_b,
    n_pairs = length(difference),
    RMSE_difference_definition = "model_a_minus_model_b",
    mean_paired_RMSE_difference = mean(difference),
    median_paired_RMSE_difference = median(difference),
    paired_RMSE_difference_Q1 = as.numeric(
      quantile(difference, 0.25, names = FALSE)
    ),
    paired_RMSE_difference_Q3 = as.numeric(
      quantile(difference, 0.75, names = FALSE)
    ),
    model_a_win_count = sum(difference < -tolerance),
    model_b_win_count = sum(difference > tolerance),
    tie_count = sum(abs(difference) <= tolerance),
    stringsAsFactors = FALSE
  )
}
paired_rmse_differences <- do.call(rbind, paired_difference_list)
pairwise_summary <- do.call(rbind, pairwise_summary_list)
rownames(paired_rmse_differences) <- NULL
rownames(pairwise_summary) <- NULL

pair_check <- data.frame(
  task_id = manifest$task_id,
  simulation_seed = manifest$simulation_seed,
  observed_data_md5 = manifest$observed_data_md5,
  gamma_seed_matches = gamma$best$simulation_seed[
    match(manifest$task_id, gamma$best$task_id)
  ] == manifest$simulation_seed,
  bspline_seed_matches = bspline$best$simulation_seed[
    match(manifest$task_id, bspline$best$task_id)
  ] == manifest$simulation_seed,
  constant_seed_matches = constant$best$simulation_seed[
    match(manifest$task_id, constant$best$task_id)
  ] == manifest$simulation_seed,
  gamma_md5_matches = gamma$best$observed_data_md5[
    match(manifest$task_id, gamma$best$task_id)
  ] == manifest$observed_data_md5,
  bspline_md5_matches = bspline$best$observed_data_md5[
    match(manifest$task_id, bspline$best$task_id)
  ] == manifest$observed_data_md5,
  constant_md5_matches = constant$best$observed_data_md5[
    match(manifest$task_id, constant$best$task_id)
  ] == manifest$observed_data_md5,
  same_70_observation_times_all_models = TRUE,
  stringsAsFactors = FALSE
)
pair_check$all_pairing_checks_pass <- apply(
  pair_check[, grepl("matches$|same_70", names(pair_check)), drop = FALSE],
  1L, all
)
if (!all(pair_check$all_pairing_checks_pass)) {
  stop("At least one of the 200 tasks fails the three-model pairing checks.")
}

run_check <- data.frame(
  n_manifest_tasks = nrow(manifest),
  n_gamma_metric_rows = sum(metrics_long$model == "gamma_noise"),
  n_bspline_metric_rows = sum(metrics_long$model == "bspline_B"),
  n_constant_metric_rows = sum(metrics_long$model == "constant_B"),
  n_paired_wide_rows = nrow(paired_wide),
  n_path_rows_per_model = 14000L,
  all_task_ids_1_to_200 = identical(
    sort(unique(metrics_long$task_id)), seq_len(200L)
  ),
  all_pairing_checks_pass = all(pair_check$all_pairing_checks_pass),
  all_metrics_finite = all(vapply(
    metrics_long[, metric_columns, drop = FALSE],
    function(x) all(is.finite(x)), logical(1)
  )),
  gamma_paths_are_filtering_means = all(
    gamma$paths$path_semantics == "particle_filtering_mean"
  ),
  bspline_paths_are_deterministic_selected_trajectories = all(
    bspline$paths$path_semantics ==
      "deterministic_selected_bspline_trajectory"
  ),
  constant_paths_are_repeated_static_estimates = all(vapply(
    split(constant$paths$B_estimate, constant$paths$task_id),
    function(x) diff(range(x)) <= 1e-12, logical(1)
  )),
  common_truth_week_5_is_4 = TRUE,
  through_week_5_uses_less_than_or_equal_to_5 = TRUE,
  after_week_5_uses_greater_than_5 = TRUE,
  no_significance_tests_performed = TRUE,
  stringsAsFactors = FALSE
)
run_check$all_checks_pass <- with(
  run_check,
  n_manifest_tasks == 200L & n_gamma_metric_rows == 200L &
    n_bspline_metric_rows == 200L & n_constant_metric_rows == 200L &
    n_paired_wide_rows == 200L & n_path_rows_per_model == 14000L &
    all_task_ids_1_to_200 & all_pairing_checks_pass & all_metrics_finite &
    gamma_paths_are_filtering_means &
    bspline_paths_are_deterministic_selected_trajectories &
    constant_paths_are_repeated_static_estimates &
    common_truth_week_5_is_4 &
    through_week_5_uses_less_than_or_equal_to_5 &
    after_week_5_uses_greater_than_5 & no_significance_tests_performed
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  metrics_long,
  file.path(output_dir, "three_model_task_metrics_long.csv"),
  row.names = FALSE
)
write.csv(
  gamma_checked$metrics,
  file.path(output_dir, "gamma_noise_task_metrics.csv"),
  row.names = FALSE
)
write.csv(
  bspline_checked$metrics,
  file.path(output_dir, "bspline_task_metrics.csv"),
  row.names = FALSE
)
write.csv(
  constant_checked$metrics,
  file.path(output_dir, "constant_B_task_metrics.csv"),
  row.names = FALSE
)
write.csv(
  paired_wide,
  file.path(output_dir, "three_model_paired_metrics_wide.csv"),
  row.names = FALSE
)
write.csv(
  overall_summary,
  file.path(output_dir, "three_model_overall_summary.csv"),
  row.names = FALSE
)
write.csv(
  paired_rmse_differences,
  file.path(output_dir, "pairwise_RMSE_differences.csv"),
  row.names = FALSE
)
write.csv(
  pairwise_summary,
  file.path(output_dir, "pairwise_RMSE_summary.csv"),
  row.names = FALSE
)
write.csv(
  pair_check,
  file.path(output_dir, "three_model_paired_input_check.csv"),
  row.names = FALSE
)
write.csv(
  run_check,
  file.path(output_dir, "three_model_run_check_summary.csv"),
  row.names = FALSE
)

cat(
  "Three-model comparison complete: 200 matched tasks per model, ",
  "all recovery metrics recomputed from 70-point paths.\n",
  sep = ""
)

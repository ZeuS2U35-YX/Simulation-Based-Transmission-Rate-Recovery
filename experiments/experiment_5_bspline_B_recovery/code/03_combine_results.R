# ============================================================
# Validate and combine exactly 200 final B-spline task results
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

args <- commandArgs(trailingOnly = TRUE)
raw_root <- if (length(args) >= 1L) args[[1]] else {
  file.path(experiment_directory, "results", "bspline")
}
combined_dir <- if (length(args) >= 2L) args[[2]] else {
  file.path(experiment_directory, "results", "combined", "bspline")
}
manifest_file <- if (length(args) >= 3L) args[[3]] else {
  file.path(experiment_directory, "results", "paired_input_manifest.csv")
}

manifest <- read_csv(manifest_file)
require_exact_task_ids(
  manifest$task_id, experiment_config$n_tasks, "Paired input manifest"
)
if (nrow(manifest) != 200L) {
  stop("Paired input manifest must contain exactly 200 rows.")
}

required_files <- c(
  "COMPLETE", "start_selection_audit.csv", "likelihood_evaluations.csv",
  "best_fit_summary.csv", "B_path.csv", "best_mif2.rds", "run_config.csv"
)
audit_list <- list()
evaluation_list <- list()
best_list <- list()
path_list <- list()
config_list <- list()
trace_list <- list()
problems <- character(0)

for (task_id in seq_len(experiment_config$n_tasks)) {
  task_dir <- file.path(raw_root, sprintf("task_%03d", task_id))
  missing <- required_files[!file.exists(file.path(task_dir, required_files))]
  if (length(missing) > 0L) {
    problems <- c(
      problems,
      paste0(
        "task_", sprintf("%03d", task_id), " missing: ",
        paste(missing, collapse = ", ")
      )
    )
    next
  }

  audit <- read_csv(file.path(task_dir, "start_selection_audit.csv"))
  evaluations <- read_csv(file.path(task_dir, "likelihood_evaluations.csv"))
  best <- read_csv(file.path(task_dir, "best_fit_summary.csv"))
  path <- read_csv(file.path(task_dir, "B_path.csv"))
  run_config <- read_csv(file.path(task_dir, "run_config.csv"))
  manifest_row <- manifest[manifest$task_id == task_id, , drop = FALSE]

  if (nrow(audit) != experiment_config$n_start) {
    problems <- c(
      problems,
      paste0(
        "task_", sprintf("%03d", task_id), " has ", nrow(audit),
        " internal-start rows; expected ", experiment_config$n_start
      )
    )
  }
  if (nrow(evaluations) !=
      experiment_config$n_start * experiment_config$n_pf_evals) {
    problems <- c(
      problems,
      paste0(
        "task_", sprintf("%03d", task_id),
        " has an unexpected likelihood-evaluation row count"
      )
    )
  }
  if (nrow(best) != 1L) {
    problems <- c(
      problems,
      paste0(
        "task_", sprintf("%03d", task_id),
        " must have exactly one final best-fit row"
      )
    )
    next
  }
  if (nrow(path) != 70L ||
      is.unsorted(path$week, strictly = TRUE) || anyDuplicated(path$week)) {
    problems <- c(
      problems,
      paste0(
        "task_", sprintf("%03d", task_id),
        " must have exactly 70 ordered B-path rows"
      )
    )
  }
  expected_times <- seq(
    experiment_config$observation_dt,
    experiment_config$n_weeks,
    by = experiment_config$observation_dt
  )
  expected_true <- ifelse(
    path$week <= experiment_config$truth[["t_switch"]],
    experiment_config$truth[["Beta_high"]],
    experiment_config$truth[["Beta_low"]]
  )
  legacy_true <- ifelse(
    path$week < experiment_config$truth[["t_switch"]],
    experiment_config$truth[["Beta_high"]],
    experiment_config$truth[["Beta_low"]]
  )
  uses_current_truth <- nrow(path) == 70L &&
    max(abs(path$B_true - expected_true)) <= 1e-12
  uses_legacy_week5_truth <- nrow(path) == 70L &&
    max(abs(path$B_true - legacy_true)) <= 1e-12
  if (nrow(path) == 70L &&
      (max(abs(path$week - expected_times)) > 1e-12 ||
       !("path_semantics" %in% names(path)) ||
       !all(
         path$path_semantics == "deterministic_selected_bspline_trajectory"
       ) ||
       (!uses_current_truth && !uses_legacy_week5_truth))) {
    problems <- c(
      problems,
      paste0(
        "task_", sprintf("%03d", task_id),
        " has incorrect B-path times, semantics, or true B(t); B(5) must be 4"
      )
    )
  }
  if (nrow(path) == 70L &&
      (uses_current_truth || uses_legacy_week5_truth)) {
    path$B_true <- expected_true
    path_error <- path$B_estimate - path$B_true
    best$RSS <- sum(path_error^2)
    best$RMSE <- sqrt(mean(path_error^2))
    best$mean_error <- mean(path_error)
    best$AOB <- abs(best$mean_error)
    best$mean_error_through_5 <- mean(
      path_error[
        path$week <= experiment_config$truth[["t_switch"]]
      ]
    )
    best$mean_error_after_5 <- mean(
      path_error[
        path$week > experiment_config$truth[["t_switch"]]
      ]
    )
    # Retain the historical column name for backward compatibility, but do
    # not trust the value stored in the raw task summary.
    best$B_rmse <- best$RMSE
  }

  task_tables <- list(audit = audit, evaluations = evaluations, best = best, path = path)
  for (table_name in names(task_tables)) {
    x <- task_tables[[table_name]]
    if (!all(x$task_id == task_id)) {
      problems <- c(
        problems,
        paste0(
          "task_", sprintf("%03d", task_id), " ", table_name,
          " contains a task_id mismatch"
        )
      )
    }
    if (!all(x$model == "bspline_B")) {
      problems <- c(
        problems,
        paste0(
          "task_", sprintf("%03d", task_id), " ", table_name,
          " contains a model mismatch"
        )
      )
    }
    if (!all(x$simulation_seed == manifest_row$simulation_seed[[1]]) ||
        !all(x$observed_data_md5 == manifest_row$observed_data_md5[[1]])) {
      problems <- c(
        problems,
        paste0(
          "task_", sprintf("%03d", task_id), " ", table_name,
          " does not match the accepted-data manifest"
        )
      )
    }
  }

  if (!all(audit$Nmif == experiment_config$Nmif) ||
      best$Nmif[[1]] != experiment_config$Nmif) {
    problems <- c(
      problems,
      paste0(
        "task_", sprintf("%03d", task_id),
        " is not the configured Nmif=", experiment_config$Nmif, " fit"
      )
    )
  }
  if (best$status[[1]] != "success" ||
      best$selection_rule[[1]] !=
        "maximum_independently_evaluated_logLik") {
    problems <- c(
      problems,
      paste0(
        "task_", sprintf("%03d", task_id),
        " does not contain one successfully selected final fit"
      )
    )
  }

  eligible <- audit[
    as.logical(audit$eligible_for_selection) &
      is.finite(audit$independently_evaluated_logLik),
    ,
    drop = FALSE
  ]
  if (nrow(eligible) == 0L) {
    problems <- c(
      problems,
      paste0("task_", sprintf("%03d", task_id), " has no eligible start")
    )
  } else {
    expected_best <- eligible[which.max(
      eligible$independently_evaluated_logLik
    ), , drop = FALSE]
    same_run <- best$best_run[[1]] == expected_best$run[[1]]
    same_logLik <- isTRUE(all.equal(
      as.numeric(best$logLik[[1]]),
      as.numeric(expected_best$independently_evaluated_logLik[[1]]),
      tolerance = 1e-12
    ))
    if (!same_run || !same_logLik) {
      problems <- c(
        problems,
        paste0(
          "task_", sprintf("%03d", task_id),
          " final fit is not the maximum independently evaluated candidate"
        )
      )
    }
  }

  evaluation_counts <- table(evaluations$run)
  if (length(evaluation_counts) != experiment_config$n_start ||
      any(evaluation_counts != experiment_config$n_pf_evals) ||
      anyDuplicated(evaluations$evaluation_seed)) {
    problems <- c(
      problems,
      paste0(
        "task_", sprintf("%03d", task_id),
        " does not have distinct complete PF-evaluation slots per start"
      )
    )
  }

  rds_files <- list.files(task_dir, pattern = "\\.rds$", full.names = FALSE)
  if (!identical(rds_files, "best_mif2.rds")) {
    problems <- c(
      problems,
      paste0(
        "task_", sprintf("%03d", task_id),
        " must retain only best_mif2.rds as its fitted object"
      )
    )
  }

  audit_list[[as.character(task_id)]] <- audit
  evaluation_list[[as.character(task_id)]] <- evaluations
  best_list[[as.character(task_id)]] <- best
  path_list[[as.character(task_id)]] <- path
  run_config$task_id <- task_id
  config_list[[as.character(task_id)]] <- run_config

  trace_file <- file.path(task_dir, "selected_mif2_trace.csv")
  if (file.exists(trace_file)) {
    trace_list[[as.character(task_id)]] <- read_csv(trace_file)
  }
}

dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(problems, file.path(combined_dir, "combine_problems.txt"))
if (length(problems) > 0L) {
  stop(
    "Combination found ", length(problems),
    " problem(s). See ", file.path(combined_dir, "combine_problems.txt"), "."
  )
}

combined_audit <- do.call(rbind, audit_list)
combined_evaluations <- do.call(rbind, evaluation_list)
combined_best <- do.call(rbind, best_list)
combined_path <- do.call(rbind, path_list)
combined_config <- do.call(rbind, config_list)
require_exact_task_ids(
  combined_best$task_id, experiment_config$n_tasks,
  "Combined final B-spline results"
)
if (nrow(combined_best) != 200L) {
  stop("Combined final B-spline result must contain exactly 200 rows.")
}
expected_path_rows <- experiment_config$n_tasks * 70L
path_rows_per_task <- table(combined_path$task_id)
if (nrow(combined_path) != expected_path_rows ||
    length(path_rows_per_task) != experiment_config$n_tasks ||
    any(path_rows_per_task != 70L)) {
  stop(
    "Combined B-spline paths must contain exactly 14,000 rows: ",
    "70 observation times for each of 200 tasks."
  )
}

write.csv(
  combined_audit,
  file.path(combined_dir, "combined_start_selection_audit.csv"),
  row.names = FALSE
)
write.csv(
  combined_evaluations,
  file.path(combined_dir, "combined_likelihood_evaluations.csv"),
  row.names = FALSE
)
write.csv(
  combined_best,
  file.path(combined_dir, "combined_best_fit_summary.csv"),
  row.names = FALSE
)
write.csv(
  combined_path,
  file.path(combined_dir, "combined_B_paths.csv"),
  row.names = FALSE
)
write.csv(
  combined_config,
  file.path(combined_dir, "combined_run_config.csv"),
  row.names = FALSE
)
if (length(trace_list) > 0L) {
  write.csv(
    do.call(rbind, trace_list),
    file.path(combined_dir, "combined_selected_mif2_traces.csv"),
    row.names = FALSE
  )
}

manifest_index <- match(combined_path$task_id, manifest$task_id)
metric_columns <- c(
  "RSS", "RMSE", "mean_error", "AOB",
  "mean_error_through_5", "mean_error_after_5"
)
run_check <- data.frame(
  model = "bspline_B",
  n_accepted_observed_data = nrow(manifest),
  n_expected_tasks = experiment_config$n_tasks,
  n_final_bspline_fits = nrow(combined_best),
  n_combined_B_path_rows = nrow(combined_path),
  n_expected_B_path_rows = expected_path_rows,
  n_tasks_with_70_observation_times = sum(path_rows_per_task == 70L),
  n_unique_task_ids = length(unique(combined_best$task_id)),
  n_unique_task_seed_pairs = nrow(unique(
    combined_best[, c("task_id", "simulation_seed")]
  )),
  n_unique_task_checksum_pairs = nrow(unique(
    combined_best[, c("task_id", "observed_data_md5")]
  )),
  all_nmif_configured = all(combined_best$Nmif == experiment_config$Nmif),
  all_status_success = all(combined_best$status == "success"),
  all_model_bspline_B = all(combined_best$model == "bspline_B") &&
    all(combined_path$model == "bspline_B"),
  all_path_semantics_expected = all(
    combined_path$path_semantics ==
      "deterministic_selected_bspline_trajectory"
  ),
  all_simulation_seeds_match_manifest = all(
    combined_path$simulation_seed == manifest$simulation_seed[manifest_index]
  ),
  all_observed_data_md5_match_manifest = all(
    combined_path$observed_data_md5 ==
      manifest$observed_data_md5[manifest_index]
  ),
  all_truth_uses_week_5_B4 = all(
    combined_path$B_true == ifelse(
      combined_path$week <= experiment_config$truth[["t_switch"]],
      experiment_config$truth[["Beta_high"]],
      experiment_config$truth[["Beta_low"]]
    )
  ),
  all_requested_metrics_recomputed = all(
    metric_columns %in% names(combined_best)
  ) && all(vapply(
    combined_best[, metric_columns, drop = FALSE],
    function(x) all(is.finite(x)),
    logical(1)
  )),
  exactly_one_final_fit_per_task = TRUE,
  stringsAsFactors = FALSE
)
run_check$all_checks_pass <- with(
  run_check,
  n_accepted_observed_data == 200L &
    n_expected_tasks == 200L &
    n_final_bspline_fits == 200L &
    n_combined_B_path_rows == 14000L &
    n_expected_B_path_rows == 14000L &
    n_tasks_with_70_observation_times == 200L &
    n_unique_task_ids == 200L &
    n_unique_task_seed_pairs == 200L &
    n_unique_task_checksum_pairs == 200L &
    all_nmif_configured & all_status_success & all_model_bspline_B &
    all_path_semantics_expected & all_simulation_seeds_match_manifest &
    all_observed_data_md5_match_manifest & all_truth_uses_week_5_B4 &
    all_requested_metrics_recomputed & exactly_one_final_fit_per_task
)
write.csv(
  run_check, file.path(combined_dir, "run_check_summary.csv"),
  row.names = FALSE
)

cat(
  "Combined exactly 200 tasks into exactly 200 final B-spline fits.\n",
  "Internal-start audit rows are retained separately and are not ",
  "simulation replicates or final fits.\n",
  sep = ""
)

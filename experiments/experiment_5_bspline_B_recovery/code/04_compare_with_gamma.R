# ============================================================
# Build the final 200-row paired Gamma-noise/B-spline comparison
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
output_dir <- if (length(args) >= 3L) args[[3]] else {
  file.path(experiment_directory, "results", "comparison")
}
manifest_file <- if (length(args) >= 4L) args[[4]] else {
  file.path(experiment_directory, "results", "paired_input_manifest.csv")
}

manifest <- read_csv(manifest_file)
gamma_best <- read_csv(file.path(gamma_dir, "combined_best_fit_summary.csv"))
gamma_paths <- read_csv(file.path(gamma_dir, "combined_B_paths.csv"))
bspline_best <- read_csv(
  file.path(bspline_dir, "combined_best_fit_summary.csv")
)
bspline_paths <- read_csv(file.path(bspline_dir, "combined_B_paths.csv"))

if (!all(c("filter_seed", "path_semantics") %in% names(gamma_paths)) ||
    !all(gamma_paths$path_semantics == "particle_filtering_mean")) {
  stop(
    "Gamma main-analysis paths must be particle filtering means, not ",
    "sampled latent trajectories."
  )
}
if (!"path_semantics" %in% names(bspline_paths) ||
    !all(
      bspline_paths$path_semantics ==
        "deterministic_selected_bspline_trajectory"
    )) {
  stop("B-spline paths have unexpected trajectory semantics.")
}

require_exact_task_ids(
  manifest$task_id, experiment_config$n_tasks, "Paired input manifest"
)
require_exact_task_ids(
  gamma_best$task_id, experiment_config$n_tasks, "Gamma final fits"
)
require_exact_task_ids(
  bspline_best$task_id, experiment_config$n_tasks, "B-spline final fits"
)
if (nrow(manifest) != 200L || nrow(gamma_best) != 200L ||
    nrow(bspline_best) != 200L) {
  stop(
    "The manifest, Gamma final fits, and B-spline final fits must each ",
    "contain exactly 200 rows."
  )
}
if (!all(gamma_best$model == "gamma_noise") ||
    !all(gamma_best$status == "success")) {
  stop("All 200 Gamma reference fits must be successful gamma_noise fits.")
}
if (!all(bspline_best$model == "bspline_B") ||
    !all(bspline_best$status == "success") ||
    !all(
      bspline_best$selection_rule ==
        "maximum_independently_evaluated_logLik"
    )) {
  stop(
    "All 200 B-spline rows must be successfully selected by maximum ",
    "independently evaluated log likelihood."
  )
}

pair_check <- merge(
  manifest[, c("task_id", "simulation_seed", "observed_data_md5")],
  gamma_best[, c(
    "task_id", "simulation_seed", "observed_data_md5", "Nmif"
  )],
  by = "task_id",
  suffixes = c("_manifest", "_gamma"),
  sort = TRUE
)
pair_check <- merge(
  pair_check,
  bspline_best[, c(
    "task_id", "simulation_seed", "observed_data_md5", "Nmif"
  )],
  by = "task_id",
  suffixes = c("_gamma", "_bspline"),
  sort = TRUE
)
names(pair_check)[names(pair_check) == "simulation_seed"] <-
  "simulation_seed_bspline"
names(pair_check)[names(pair_check) == "observed_data_md5"] <-
  "observed_data_md5_bspline"

pair_check$same_simulation_seed <-
  pair_check$simulation_seed_manifest == pair_check$simulation_seed_gamma &
  pair_check$simulation_seed_manifest == pair_check$simulation_seed_bspline
pair_check$same_observed_data_checksum <-
  pair_check$observed_data_md5_manifest ==
    pair_check$observed_data_md5_gamma &
  pair_check$observed_data_md5_manifest ==
    pair_check$observed_data_md5_bspline
pair_check$both_nmif_configured <-
  pair_check$Nmif_gamma == experiment_config$Nmif &
  pair_check$Nmif_bspline == experiment_config$Nmif

if (nrow(pair_check) != 200L ||
    !all(pair_check$same_simulation_seed) ||
    !all(pair_check$same_observed_data_checksum) ||
    !all(pair_check$both_nmif_configured)) {
  stop(
    "The two models do not form exactly 200 pairs with identical task_id, ",
    "simulation seed, observed-data checksum, and configured Nmif."
  )
}

calculate_metrics <- function(paths, model_name) {
  required <- c(
    "task_id", "simulation_seed", "observed_data_md5",
    "week", "B_estimate", "B_true"
  )
  missing <- setdiff(required, names(paths))
  if (length(missing) > 0L) {
    stop(model_name, " paths are missing: ", paste(missing, collapse = ", "))
  }

  split_paths <- split(paths, paths$task_id)
  if (!identical(sort(as.integer(names(split_paths))), seq_len(200L))) {
    stop(model_name, " paths do not contain exactly tasks 1 through 200.")
  }
  out <- do.call(rbind, lapply(split_paths, function(x) {
    task_id <- as.integer(x$task_id[[1]])
    manifest_row <- manifest[manifest$task_id == task_id, , drop = FALSE]
    expected_times <- seq(
      experiment_config$observation_dt,
      experiment_config$n_weeks,
      by = experiment_config$observation_dt
    )
    if (nrow(x) != 70L || is.unsorted(x$week, strictly = TRUE) ||
        anyDuplicated(x$week) ||
        max(abs(x$week - expected_times)) > 1e-12) {
      stop(
        model_name, " task ", task_id,
        " does not contain exactly 70 ordered observation times."
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
    expected_true <- ifelse(
      x$week <= experiment_config$truth[["t_switch"]],
      experiment_config$truth[["Beta_high"]],
      experiment_config$truth[["Beta_low"]]
    )
    if (max(abs(x$B_true - expected_true)) > 1e-12) {
      stop(
        model_name, " task ", task_id,
        " does not encode B(5)=4 and the configured true B(t)."
      )
    }
    error <- x$B_estimate - x$B_true
    data.frame(
      task_id = task_id,
      RSS = sum(error^2),
      RMSE = sqrt(mean(error^2)),
      mean_error = mean(error),
      AOB = abs(mean(error)),
      mean_error_through_5 = mean(
        error[x$week <= experiment_config$truth[["t_switch"]]]
      ),
      mean_error_after_5 = mean(
        error[x$week > experiment_config$truth[["t_switch"]]]
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  require_exact_task_ids(out$task_id, 200L, paste(model_name, "metrics"))
  out
}

gamma_metrics <- calculate_metrics(gamma_paths, "Gamma-noise")
bspline_metrics <- calculate_metrics(bspline_paths, "B-spline")
names(gamma_metrics)[names(gamma_metrics) != "task_id"] <- paste0(
  names(gamma_metrics)[names(gamma_metrics) != "task_id"], "_gamma"
)
names(bspline_metrics)[names(bspline_metrics) != "task_id"] <- paste0(
  names(bspline_metrics)[names(bspline_metrics) != "task_id"], "_bspline"
)

paired <- manifest[, c(
  "task_id", "simulation_seed", "observed_data_md5"
)]
paired <- merge(paired, gamma_metrics, by = "task_id", sort = TRUE)
paired <- merge(paired, bspline_metrics, by = "task_id", sort = TRUE)
paired <- merge(
  paired,
  gamma_best[, c("task_id", "logLik", "logLik_se", "status")],
  by = "task_id",
  sort = TRUE
)
names(paired)[names(paired) %in% c("logLik", "logLik_se", "status")] <-
  paste0(
    names(paired)[names(paired) %in% c("logLik", "logLik_se", "status")],
    "_gamma"
  )
paired <- merge(
  paired,
  bspline_best[, c(
    "task_id", "logLik", "logLik_se", "best_run",
    "n_internal_starts", "selection_rule", "status"
  )],
  by = "task_id",
  sort = TRUE
)
names(paired)[names(paired) %in% c(
  "logLik", "logLik_se", "best_run", "n_internal_starts",
  "selection_rule", "status"
)] <- paste0(
  names(paired)[names(paired) %in% c(
    "logLik", "logLik_se", "best_run", "n_internal_starts",
    "selection_rule", "status"
  )],
  "_bspline"
)

paired$delta_RMSE_bspline_minus_gamma <-
  paired$RMSE_bspline - paired$RMSE_gamma
paired$delta_RSS_bspline_minus_gamma <-
  paired$RSS_bspline - paired$RSS_gamma
paired$delta_AOB_bspline_minus_gamma <-
  paired$AOB_bspline - paired$AOB_gamma
paired$delta_logLik_bspline_minus_gamma <-
  paired$logLik_bspline - paired$logLik_gamma
paired$bspline_lower_RMSE <- paired$RMSE_bspline < paired$RMSE_gamma
paired$bspline_lower_RSS <- paired$RSS_bspline < paired$RSS_gamma
paired$bspline_lower_AOB <- paired$AOB_bspline < paired$AOB_gamma

require_exact_task_ids(
  paired$task_id, experiment_config$n_tasks, "Final paired comparison"
)
if (nrow(paired) != 200L) {
  stop("Final paired comparison must contain exactly 200 rows.")
}

overall <- data.frame(
  quantity = c(
    "mean RMSE", "mean RSS", "mean AOB",
    "mean independently evaluated log likelihood",
    "B-spline win proportion by RMSE",
    "B-spline win proportion by RSS",
    "B-spline win proportion by AOB"
  ),
  gamma = c(
    mean(paired$RMSE_gamma),
    mean(paired$RSS_gamma),
    mean(paired$AOB_gamma),
    mean(paired$logLik_gamma),
    NA_real_, NA_real_, NA_real_
  ),
  bspline = c(
    mean(paired$RMSE_bspline),
    mean(paired$RSS_bspline),
    mean(paired$AOB_bspline),
    mean(paired$logLik_bspline),
    mean(paired$bspline_lower_RMSE),
    mean(paired$bspline_lower_RSS),
    mean(paired$bspline_lower_AOB)
  ),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  pair_check,
  file.path(output_dir, "paired_input_check.csv"),
  row.names = FALSE
)
write.csv(
  paired,
  file.path(output_dir, "paired_gamma_bspline_comparison.csv"),
  row.names = FALSE
)
write.csv(
  overall,
  file.path(output_dir, "overall_gamma_bspline_summary.csv"),
  row.names = FALSE
)

cat(
  "Final paired comparison contains exactly 200 rows.\n",
  "Each row contains one existing Gamma-noise fit and one selected ",
  "B-spline fit for the same accepted observed data.\n",
  sep = ""
)

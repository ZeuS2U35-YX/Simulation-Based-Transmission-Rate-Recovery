# ============================================================
# Reconstruct Experiment 4 Gamma-noise recovery paths
# as sampled latent trajectories
#
# This script does not rerun MIF2 or the five independent likelihood
# evaluations. For each saved best fit, it runs the final particle filter
# with filter.traj=TRUE and extracts one ancestry-preserving B trajectory.
# Only the 70 observation-time values are written to the combined path file;
# t0 is retained in the provenance table for validation.
#
# Usage from the Experiment 4 root:
#   Rscript code/09_regenerate_sampled_B_trajectories.R \
#     <combined_B_output.csv> <provenance_output.csv> [workers] [task_spec]
# ============================================================

options(stringsAsFactors = FALSE, digits = 17)

suppressPackageStartupMessages(library(pomp))

source(file.path("config", "experiment_config.R"))
source(file.path("code", "model_components.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript code/09_regenerate_sampled_B_trajectories.R ",
    "<combined_B_output.csv> <provenance_output.csv> [workers] [task_spec]"
  )
}

output_path <- args[[1]]
provenance_path <- args[[2]]
workers <- if (length(args) >= 3L) as.integer(args[[3]]) else 1L
task_spec <- if (length(args) >= 4L) args[[4]] else "1:200"
if (!is.finite(workers) || workers < 1L) stop("workers must be a positive integer.")

parse_task_spec <- function(x) {
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    ends <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    return(seq.int(ends[[1]], ends[[2]]))
  }
  out <- suppressWarnings(as.integer(strsplit(x, ",", fixed = TRUE)[[1]]))
  if (any(is.na(out))) stop("Invalid task_spec: ", x)
  sort(unique(out))
}

task_ids <- parse_task_spec(task_spec)
if (any(task_ids < 1L | task_ids > experiment_config$n_tasks)) {
  stop("task_spec contains IDs outside 1:", experiment_config$n_tasks, ".")
}

best_path <- file.path(
  "results", "combined", "gamma", "combined_best_fit_summary.csv"
)
paramlist_path <- file.path("results", "paramlist.csv")
best <- read.csv(best_path, check.names = FALSE)
paramlist <- read.csv(paramlist_path, check.names = FALSE)

expected_times <- seq(
  from = experiment_config$observation_interval,
  to = experiment_config$n_weeks,
  by = experiment_config$observation_interval
)
if (length(expected_times) != 70L) stop("Configured observation grid is not length 70.")

recover_one <- function(task_id) {
  tryCatch({
    best_row <- best[best$task_id == task_id, , drop = FALSE]
    seed_row <- paramlist[paramlist$task_id == task_id, , drop = FALSE]
    if (nrow(best_row) != 1L || nrow(seed_row) != 1L) {
      stop("expected one best-fit row and one seed row")
    }
    if (!identical(as.character(best_row$status[[1]]), "success")) {
      stop("saved best fit is not successful")
    }

    observed_path <- file.path(
      "shared_data", sprintf("task_%03d", task_id), "observed_data.csv"
    )
    observed <- read.csv(observed_path, check.names = FALSE)
    if (!all(c("week", "reports") %in% names(observed))) {
      stop("observed data lack week or reports")
    }
    observed <- observed[, c("week", "reports"), drop = FALSE]
    if (nrow(observed) != 70L ||
        max(abs(observed$week - expected_times)) > 1e-12) {
      stop("observed data do not use the configured 70-time grid")
    }

    observed_md5 <- unname(tools::md5sum(observed_path))
    if (!identical(observed_md5, as.character(best_row$observed_data_md5[[1]]))) {
      stop("observed-data checksum differs from the best-fit record")
    }

    model <- make_gamma_model(observed, experiment_config)
    theta <- gamma_baseline_parameters(experiment_config)
    theta[["B0"]] <- best_row$B0_hat[[1]]
    theta[["sigma_beta"]] <- best_row$sigma_beta_hat[[1]]
    trajectory_seed <- as.integer(seed_row$gamma_final_pf_seed[[1]])

    set.seed(trajectory_seed)
    pf <- pfilter(
      model,
      params = theta,
      Np = experiment_config$Np_final,
      filter.traj = TRUE
    )
    sampled <- filter_traj(pf, vars = "B", format = "data.frame")

    if (!all(c("time", "value") %in% names(sampled))) {
      stop("filter_traj returned an unexpected schema")
    }
    if ("name" %in% names(sampled) &&
        !identical(unique(as.character(sampled$name)), "B")) {
      stop("filter_traj returned a state other than B")
    }
    expected_with_t0 <- c(0, expected_times)
    if (nrow(sampled) != 71L || any(!is.finite(sampled$time)) ||
        any(!is.finite(sampled$value)) ||
        max(abs(sampled$time - expected_with_t0)) > 1e-12 ||
        is.unsorted(sampled$time, strictly = TRUE) || anyDuplicated(sampled$time)) {
      stop("sampled trajectory failed the 71-time trajectory check")
    }
    if (abs(sampled$value[[1]] - best_row$B0_hat[[1]]) > 1e-4) {
      stop("sampled trajectory does not begin at the fitted B0")
    }

    final_logLik <- as.numeric(logLik(pf))
    if (!is.finite(final_logLik) ||
        abs(final_logLik - best_row$final_pf_logLik[[1]]) > 1e-8) {
      stop("final particle-filter log likelihood was not reproduced")
    }

    observation_path <- sampled[-1L, , drop = FALSE]
    path <- data.frame(
      task_id = task_id,
      simulation_seed = as.integer(best_row$simulation_seed[[1]]),
      simulation_attempt = as.integer(best_row$simulation_attempt[[1]]),
      observed_data_md5 = observed_md5,
      model = "gamma_noise",
      Nmif = as.integer(best_row$Nmif[[1]]),
      week = as.numeric(observation_path$time),
      B_estimate = as.numeric(observation_path$value),
      B_true = true_B_at_times(observation_path$time, experiment_config),
      trajectory_seed = trajectory_seed,
      path_semantics = "ancestry_preserving_sampled_latent_trajectory",
      stringsAsFactors = FALSE
    )

    provenance <- data.frame(
      task_id = task_id,
      best_run = as.integer(best_row$best_run[[1]]),
      B0_hat = best_row$B0_hat[[1]],
      sigma_beta_hat = best_row$sigma_beta_hat[[1]],
      trajectory_seed = trajectory_seed,
      Np = experiment_config$Np_final,
      final_pf_logLik = final_logLik,
      t0_B = sampled$value[[1]],
      n_trajectory_times_including_t0 = nrow(sampled),
      n_metric_times = nrow(path),
      path_semantics = "ancestry_preserving_sampled_latent_trajectory",
      parameter_uncertainty_integrated = FALSE,
      filtering_mean_used_for_metrics = FALSE,
      stringsAsFactors = FALSE
    )

    list(ok = TRUE, path = path, provenance = provenance, error = NA_character_)
  }, error = function(e) {
    list(
      ok = FALSE,
      path = NULL,
      provenance = NULL,
      error = paste0("task ", task_id, ": ", conditionMessage(e))
    )
  })
}

if (.Platform$OS.type == "unix" && workers > 1L) {
  recovered <- parallel::mclapply(
    task_ids,
    recover_one,
    mc.cores = workers,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
} else {
  recovered <- lapply(task_ids, recover_one)
}

errors <- vapply(recovered, function(x) x$error, character(1))
errors <- errors[!is.na(errors)]
if (length(errors) > 0L) {
  stop(
    "Sampled-trajectory reconstruction failed:\n",
    paste(errors, collapse = "\n")
  )
}

combined_path <- do.call(rbind, lapply(recovered, `[[`, "path"))
combined_provenance <- do.call(rbind, lapply(recovered, `[[`, "provenance"))
combined_path <- combined_path[order(combined_path$task_id, combined_path$week), ]
combined_provenance <- combined_provenance[order(combined_provenance$task_id), ]
rownames(combined_path) <- NULL
rownames(combined_provenance) <- NULL

if (nrow(combined_path) != 70L * length(task_ids) ||
    any(table(combined_path$task_id) != 70L)) {
  stop("Combined output does not contain exactly 70 rows per requested task.")
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(provenance_path), recursive = TRUE, showWarnings = FALSE)
write.csv(combined_path, output_path, row.names = FALSE)
write.csv(combined_provenance, provenance_path, row.names = FALSE)

cat(
  "Reconstructed ", length(task_ids),
  " sampled Gamma-noise B trajectories.\n",
  "Combined observation-time paths: ", output_path, "\n",
  "Trajectory provenance: ", provenance_path, "\n",
  sep = ""
)

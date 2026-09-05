#!/usr/bin/env Rscript

# Fast release gate for tracked, post-processed evidence. This script does not
# run simulation, MIF2, particle filtering, or Slurm jobs.

options(stringsAsFactors = FALSE, digits = 17, warn = 1)

fail <- function(...) stop(..., call. = FALSE)
expect <- function(condition, ...) {
  if (!isTRUE(condition)) fail(...)
  invisible(TRUE)
}

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) fail("Could not determine the script path.")
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

script_path <- get_script_path()
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

read_required_csv <- function(path) {
  if (!file.exists(path)) fail("Required file does not exist: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

relative <- function(...) file.path(repo_root, ...)

# Syntax gate for every project R file, including the plotting profile.
r_files <- list.files(
  repo_root,
  pattern = "[.]R$",
  recursive = TRUE,
  full.names = TRUE
)
r_files <- r_files[!grepl("[/\\]renv[/\\]library[/\\]", r_files)]
for (path in r_files) invisible(parse(file = path))
profile <- relative("experiments", "experiment_5_bspline_B_recovery", ".Rprofile")
if (file.exists(profile)) invisible(parse(file = profile))

model_component_path <- relative(
  "experiments", "experiment_4_nmif600_model_comparison", "code",
  "model_components.R"
)
model_component_source <- readLines(model_component_path, warn = FALSE)
expect(
  sum(grepl(
    "if[[:space:]]*\\(t[[:space:]]*<[[:space:]]*t_switch\\)",
    model_component_source
  )) == 1L,
  "The data-generating process must retain the t < t_switch step."
)
expect(
  any(grepl(
    "times[[:space:]]*<=[[:space:]]*config\\$true_parameters",
    model_component_source
  )),
  "The reporting helper must retain the normalized B(5)=4 convention."
)

config_env <- new.env(parent = baseenv())
sys.source(
  relative(
    "experiments", "experiment_5_bspline_B_recovery", "config",
    "experiment_config.R"
  ),
  envir = config_env
)
config <- config_env$experiment_config
expect(config$n_tasks == 200L, "Experiment 5 must define exactly 200 tasks.")

manifest_path <- relative(
  "experiments", "experiment_5_bspline_B_recovery", "results",
  "paired_input_manifest.csv"
)
manifest <- read_required_csv(manifest_path)
manifest_required <- c(
  "task_id", "simulation_seed", "observed_data_md5", "simulated_data_md5",
  "acceptance_threshold", "recomputed_max_H", "accepted"
)
missing_manifest <- setdiff(manifest_required, names(manifest))
expect(
  length(missing_manifest) == 0L,
  "Paired manifest is missing: ", paste(missing_manifest, collapse = ", ")
)
manifest_ids <- sort(as.integer(manifest$task_id))
expect(
  nrow(manifest) == 200L && !anyDuplicated(manifest_ids) &&
    identical(manifest_ids, seq_len(200L)),
  "Paired manifest must contain task IDs 1 through 200 exactly once."
)
expect(
  all(manifest$accepted) && all(manifest$acceptance_threshold == 20) &&
    all(manifest$recomputed_max_H > manifest$acceptance_threshold),
  "Paired manifest does not satisfy the accepted-outbreak contract."
)
md5_pattern <- "^[0-9a-f]{32}$"
expect(
  all(grepl(md5_pattern, manifest$observed_data_md5)) &&
    all(grepl(md5_pattern, manifest$simulated_data_md5)),
  "Paired manifest contains malformed MD5 values."
)

path_specs <- list(
  gamma_noise = list(
    path = relative(
      "experiments", "experiment_4_nmif600_model_comparison", "results",
      "combined", "gamma", "combined_B_paths.csv"
    ),
    semantics = "particle_filtering_mean",
    seed_column = "filter_seed",
    require_static = FALSE
  ),
  bspline_B = list(
    path = relative(
      "experiments", "experiment_5_bspline_B_recovery", "results",
      "combined", "bspline", "combined_B_paths.csv"
    ),
    semantics = "deterministic_selected_bspline_trajectory",
    seed_column = NULL,
    require_static = FALSE
  ),
  constant_B = list(
    path = relative(
      "experiments", "experiment_4_nmif600_model_comparison", "results",
      "combined", "constant", "combined_B_paths.csv"
    ),
    semantics = NULL,
    seed_column = NULL,
    require_static = TRUE
  )
)

expected_times <- seq(
  config$observation_dt,
  config$n_weeks,
  by = config$observation_dt
)
expect(length(expected_times) == 70L, "Expected a 70-time observation grid.")

validate_path_table <- function(model_name, spec) {
  paths <- read_required_csv(spec$path)
  required <- c(
    "task_id", "simulation_seed", "observed_data_md5", "model", "week",
    "B_estimate", "B_true"
  )
  missing <- setdiff(required, names(paths))
  expect(
    length(missing) == 0L,
    model_name, " paths are missing: ", paste(missing, collapse = ", ")
  )
  expect(
    nrow(paths) == 14000L && all(paths$model == model_name),
    model_name, " must contain 14,000 rows with the expected model label."
  )
  expect(
    all(is.finite(paths$week)) && all(is.finite(paths$B_estimate)),
    model_name, " contains non-finite path values."
  )
  if (!is.null(spec$semantics)) {
    expect(
      "path_semantics" %in% names(paths) &&
        all(paths$path_semantics == spec$semantics),
      model_name, " does not have the required path semantics."
    )
  }
  if (!is.null(spec$seed_column)) {
    expect(
      spec$seed_column %in% names(paths),
      model_name, " is missing ", spec$seed_column, "."
    )
  }

  paths <- paths[order(paths$task_id, paths$week), , drop = FALSE]
  task_counts <- table(paths$task_id)
  expect(
    identical(as.integer(names(task_counts)), seq_len(200L)) &&
      all(task_counts == 70L),
    model_name, " must contain 70 rows for each task ID 1 through 200."
  )
  expect(
    max(abs(paths$week - rep(expected_times, times = 200L))) <= 1e-12,
    model_name, " does not use the configured observation grid."
  )

  manifest_index <- match(paths$task_id, manifest$task_id)
  expect(
    all(paths$simulation_seed == manifest$simulation_seed[manifest_index]) &&
      all(paths$observed_data_md5 == manifest$observed_data_md5[manifest_index]),
    model_name, " path provenance does not match the paired manifest."
  )

  expected_truth <- ifelse(paths$week <= config$truth[["t_switch"]],
    config$truth[["Beta_high"]], config$truth[["Beta_low"]]
  )
  expect(
    max(abs(paths$B_true - expected_truth)) <= 1e-12,
    model_name, " does not use the normalized B(5)=4 reporting truth."
  )

  if (!is.null(spec$seed_column)) {
    seed_counts <- vapply(
      split(paths[[spec$seed_column]], paths$task_id),
      function(x) length(unique(x)),
      integer(1)
    )
    expect(
      all(seed_counts == 1L),
      model_name, " does not have one stable path seed per task."
    )
  }
  if (spec$require_static) {
    ranges <- vapply(
      split(paths$B_estimate, paths$task_id),
      function(x) diff(range(x)),
      numeric(1)
    )
    expect(
      all(ranges <= 1e-12),
      model_name, " is not a repeated static estimate within each task."
    )
  }
  paths
}

validated_paths <- Map(validate_path_table, names(path_specs), path_specs)
path_key <- function(paths) {
  paste(paths$task_id, sprintf("%.14f", paths$week), sep = "|")
}
reference_key <- path_key(validated_paths$gamma_noise)
expect(
  identical(reference_key, path_key(validated_paths$bspline_B)) &&
    identical(reference_key, path_key(validated_paths$constant_B)),
  "The three path tables do not use identical task/time pairs."
)

compare_csv <- function(expected_path, actual_path, tolerance = 1e-12) {
  expected <- read_required_csv(expected_path)
  actual <- read_required_csv(actual_path)
  expect(
    identical(names(actual), names(expected)) && nrow(actual) == nrow(expected),
    "CSV schema or row-count mismatch: ", basename(expected_path)
  )
  for (column in names(expected)) {
    left <- expected[[column]]
    right <- actual[[column]]
    expect(
      identical(is.na(left), is.na(right)),
      "NA pattern mismatch in ", basename(expected_path), ":", column
    )
    keep <- !is.na(left)
    if (is.numeric(left) || is.integer(left)) {
      delta <- abs(as.numeric(left[keep]) - as.numeric(right[keep]))
      scale <- pmax(1, abs(as.numeric(left[keep])))
      expect(
        length(delta) == 0L || all(delta <= tolerance * scale),
        "Numeric mismatch in ", basename(expected_path), ":", column
      )
    } else {
      expect(
        identical(as.character(left[keep]), as.character(right[keep])),
        "Value mismatch in ", basename(expected_path), ":", column
      )
    }
  }
  invisible(TRUE)
}

run_r_script <- function(script, args) {
  rscript <- Sys.which("Rscript")
  if (!nzchar(rscript)) fail("Rscript is not available.")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(repo_root)
  output <- system2(
    rscript,
    args = c("--vanilla", shQuote(script), vapply(args, shQuote, character(1))),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    fail(
      "Command failed: ", script, "\n",
      paste(output, collapse = "\n")
    )
  }
  if (length(output) > 0L) cat(paste(output, collapse = "\n"), "\n")
  invisible(output)
}

exp5 <- relative("experiments", "experiment_5_bspline_B_recovery")
comparison_dir <- file.path(exp5, "results", "comparison_three_models")
comparison_tmp <- tempfile("three_model_validation_")
dir.create(comparison_tmp, recursive = TRUE)
on.exit(unlink(comparison_tmp, recursive = TRUE, force = TRUE), add = TRUE)

run_r_script(
  file.path(
    "experiments", "experiment_5_bspline_B_recovery", "code",
    "05_compare_three_models.R"
  ),
  c(
    file.path(
      "experiments", "experiment_4_nmif600_model_comparison", "results",
      "combined", "gamma"
    ),
    file.path(
      "experiments", "experiment_5_bspline_B_recovery", "results",
      "combined", "bspline"
    ),
    file.path(
      "experiments", "experiment_4_nmif600_model_comparison", "results",
      "combined", "constant"
    ),
    comparison_tmp,
    file.path(
      "experiments", "experiment_5_bspline_B_recovery", "results",
      "paired_input_manifest.csv"
    )
  )
)

expected_comparison_files <- sort(list.files(
  comparison_dir,
  pattern = "[.]csv$",
  full.names = FALSE
))
generated_comparison_files <- sort(list.files(
  comparison_tmp,
  pattern = "[.]csv$",
  full.names = FALSE
))
expected_core_files <- setdiff(
  expected_comparison_files,
  c("week8_sensitivity_task_metrics.csv", "week8_sensitivity_summary.csv")
)
expect(
  identical(expected_core_files, generated_comparison_files),
  "The recomputed three-model output set differs from the tracked output set."
)
for (name in expected_core_files) {
  compare_csv(file.path(comparison_dir, name), file.path(comparison_tmp, name))
}

overall <- read_required_csv(file.path(comparison_tmp, "three_model_overall_summary.csv"))
headline <- c(
  gamma_noise = 0.521755255106023,
  bspline_B = 0.712764112697815,
  constant_B = 1.18583153801711
)
for (model in names(headline)) {
  actual <- overall$mean_RMSE[overall$model == model]
  expect(
    length(actual) == 1L && abs(actual - headline[[model]]) <= 1e-12,
    "Headline mean RMSE changed for ", model, "."
  )
}

pairwise <- read_required_csv(file.path(comparison_tmp, "pairwise_RMSE_summary.csv"))
expected_wins <- data.frame(
  model_a = c("gamma_noise", "bspline_B", "bspline_B"),
  model_b = c("constant_B", "constant_B", "gamma_noise"),
  model_a_win_count = c(200L, 188L, 38L),
  model_b_win_count = c(0L, 12L, 162L),
  tie_count = c(0L, 0L, 0L),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(expected_wins))) {
  row <- pairwise[
    pairwise$model_a == expected_wins$model_a[[i]] &
      pairwise$model_b == expected_wins$model_b[[i]],
    ,
    drop = FALSE
  ]
  expect(
    nrow(row) == 1L && row$n_pairs == 200L &&
      row$model_a_win_count == expected_wins$model_a_win_count[[i]] &&
      row$model_b_win_count == expected_wins$model_b_win_count[[i]] &&
      row$tie_count == expected_wins$tie_count[[i]],
    "Paired RMSE win counts changed for ",
    expected_wins$model_a[[i]], " versus ", expected_wins$model_b[[i]], "."
  )
}

week8_tmp <- tempfile("week8_validation_")
dir.create(week8_tmp, recursive = TRUE)
on.exit(unlink(week8_tmp, recursive = TRUE, force = TRUE), add = TRUE)
run_r_script(
  file.path("scripts", "compute_week8_sensitivity.R"),
  c(week8_tmp, "8")
)
for (name in c(
  "week8_sensitivity_task_metrics.csv",
  "week8_sensitivity_summary.csv"
)) {
  compare_csv(file.path(comparison_dir, name), file.path(week8_tmp, name))
}

week8 <- read_required_csv(file.path(week8_tmp, "week8_sensitivity_summary.csv"))
expect(
  nrow(week8) == 1L &&
    abs(week8$gamma_mean_RMSE - 0.550096374886285) <= 1e-12 &&
    abs(week8$bspline_mean_RMSE - 0.580798337692794) <= 1e-12 &&
    week8$gamma_lower_RMSE_count == 114L &&
    week8$bspline_lower_RMSE_count == 86L &&
    week8$tie_count == 0L,
  "Week-8 sensitivity results changed unexpectedly."
)

cat(
  "Release validation passed: syntax, manifest, path semantics, pairing, ",
  "full-window results, and week-8 sensitivity all match tracked evidence.\n",
  sep = ""
)

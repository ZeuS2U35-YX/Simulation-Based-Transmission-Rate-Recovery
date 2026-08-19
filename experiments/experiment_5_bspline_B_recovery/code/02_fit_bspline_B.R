# ============================================================
# Fit one B-spline model to one accepted Experiment 4 data set
#
# Usage:
# Rscript code/02_fit_bspline_B.R \
#   <task_id> [shared_data_root] [output_root] [input_manifest]
#
# Multiple starts are internal optimization attempts for this one task.
# Only the candidate with the largest independently evaluated likelihood
# becomes the task's final fit.
# ============================================================

library(pomp)
options(stringsAsFactors = FALSE)

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
)
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "path_helpers.R"))

experiment_directory <- get_experiment_directory()
source(file.path(experiment_directory, "config", "experiment_config.R"))
source(file.path(experiment_directory, "code", "model_components.R"))
source(file.path(experiment_directory, "code", "io_helpers.R"))

source_experiment_directory <- file.path(
  dirname(experiment_directory), experiment_config$source_experiment_id
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop(
    "Usage: Rscript code/02_fit_bspline_B.R ",
    "<task_id> [shared_data_root] [output_root] [input_manifest]"
  )
}

task_id <- suppressWarnings(as.integer(args[[1]]))
if (is.na(task_id) || task_id < 1L ||
    task_id > experiment_config$n_tasks) {
  stop("task_id must be between 1 and ", experiment_config$n_tasks, ".")
}
shared_data_root <- if (length(args) >= 2L) args[[2]] else {
  file.path(source_experiment_directory, "shared_data")
}
output_root <- if (length(args) >= 3L) args[[3]] else {
  file.path(experiment_directory, "results", "bspline")
}
manifest_file <- if (length(args) >= 4L) args[[4]] else {
  file.path(experiment_directory, "results", "paired_input_manifest.csv")
}

fit_config <- get_fit_config(experiment_config)
manifest_row <- read_manifest_task(
  manifest_file, task_id, experiment_config$n_tasks
)
shared <- read_exp4_shared_task(shared_data_root, task_id)

simulation_seed <- as.integer(shared$metadata$simulation_seed[[1]])
simulation_attempt <- as.integer(shared$metadata$simulation_attempt[[1]])
if (simulation_seed != as.integer(manifest_row$simulation_seed[[1]]) ||
    simulation_attempt != as.integer(manifest_row$simulation_attempt[[1]]) ||
    !identical(
      shared$observed_data_md5,
      as.character(manifest_row$observed_data_md5[[1]])
    ) ||
    !identical(
      shared$simulated_data_md5,
      as.character(manifest_row$simulated_data_md5[[1]])
    ) ||
    !identical(
      shared$acceptance_threshold,
      as.numeric(manifest_row$acceptance_threshold[[1]])
    ) ||
    !identical(
      as.numeric(shared$recomputed_max_H),
      as.numeric(manifest_row$recomputed_max_H[[1]])
    ) ||
    !isTRUE(as.logical(manifest_row$accepted[[1]]))) {
  stop(
    "Task ", task_id,
    " does not match the validated paired-input manifest."
  )
}

atomic <- start_atomic_task(output_root, task_id)
if (atomic$skip) quit(save = "no", status = 0L)
on.exit({
  if (dir.exists(atomic$temp_dir)) {
    unlink(atomic$temp_dir, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)

observed_data <- shared$observed_data
bspline_model <- make_bspline_model(observed_data, experiment_config)
start_seed <- task_start_seed(experiment_config$seeds$starts, task_id)
start_values <- make_start_values(
  n_start = fit_config$n_start,
  seed = start_seed,
  config = experiment_config
)
mif_rw_sd <- make_mif_rw_sd(experiment_config)

cat(
  "Experiment = ", experiment_config$experiment_id,
  "\nModel = deterministic six-parameter B-spline",
  "\nTask ID = ", task_id,
  "\nSimulation seed = ", simulation_seed,
  "\nObserved-data MD5 = ", shared$observed_data_md5,
  "\nNmif = ", fit_config$Nmif,
  "\nNp_mif = ", fit_config$Np_mif,
  "\nInternal starts = ", nrow(start_values),
  "\nNp_eval = ", fit_config$Np_eval,
  "\nIndependent PF evaluations per start = ", fit_config$n_pf_evals,
  "\n", sep = ""
)

mif_objects <- vector("list", nrow(start_values))
selection_rows <- vector("list", nrow(start_values))
evaluation_rows <- vector("list", nrow(start_values))

for (s in seq_len(nrow(start_values))) {
  start_coefficients <- as.numeric(start_values[s, coefficient_names])
  names(start_coefficients) <- coefficient_names
  theta_start <- make_parameter_vector(start_coefficients, experiment_config)
  mif_seed <- task_start_seed(
    experiment_config$seeds$mif2, task_id, start_values$start_id[[s]]
  )
  evaluation_seed_base <- task_start_seed(
    experiment_config$seeds$evaluation,
    task_id,
    start_values$start_id[[s]]
  )

  cat(
    "Running internal start ", s, "/", nrow(start_values),
    " (", start_values$start_source[[s]], ")\n", sep = ""
  )

  set.seed(mif_seed)
  mif_now <- tryCatch(
    mif2(
      bspline_model,
      params = theta_start,
      Np = fit_config$Np_mif,
      Nmif = fit_config$Nmif,
      rw.sd = mif_rw_sd,
      cooling.type = experiment_config$cooling_type,
      cooling.fraction.50 = experiment_config$cooling_fraction_50
    ),
    error = function(e) structure(
      list(message = conditionMessage(e)), class = "mif2_error"
    )
  )

  base_row <- data.frame(
    run = s,
    start_id = start_values$start_id[[s]],
    start_source = start_values$start_source[[s]],
    mif_seed = mif_seed,
    evaluation_seed_base = evaluation_seed_base,
    as.list(setNames(
      start_coefficients, paste0("start_", coefficient_names)
    )),
    check.names = FALSE
  )

  if (inherits(mif_now, "mif2_error")) {
    selection_rows[[s]] <- cbind(
      base_row,
      as.data.frame(as.list(setNames(
        rep(NA_real_, length(coefficient_names)),
        paste0(coefficient_names, "_hat")
      ))),
      independently_evaluated_logLik = NA_real_,
      logLik_se = NA_real_,
      n_successful_pf_evals = 0L,
      eligible_for_selection = FALSE,
      status = "mif2_failed",
      error_message = mif_now$message
    )
    evaluation_rows[[s]] <- data.frame(
      run = s,
      evaluation_rep = seq_len(fit_config$n_pf_evals),
      evaluation_seed = evaluation_seed_base + seq_len(fit_config$n_pf_evals),
      logLik = NA_real_,
      status = "not_run_mif2_failed",
      stringsAsFactors = FALSE
    )
    next
  }

  fitted_parameters <- tryCatch(coef(mif_now), error = function(e) NULL)
  valid_coefficients <- !is.null(fitted_parameters) &&
    all(is.finite(fitted_parameters[coefficient_names]))
  if (!valid_coefficients) {
    selection_rows[[s]] <- cbind(
      base_row,
      as.data.frame(as.list(setNames(
        rep(NA_real_, length(coefficient_names)),
        paste0(coefficient_names, "_hat")
      ))),
      independently_evaluated_logLik = NA_real_,
      logLik_se = NA_real_,
      n_successful_pf_evals = 0L,
      eligible_for_selection = FALSE,
      status = "invalid_mif2_estimate",
      error_message = NA_character_
    )
    evaluation_rows[[s]] <- data.frame(
      run = s,
      evaluation_rep = seq_len(fit_config$n_pf_evals),
      evaluation_seed = evaluation_seed_base + seq_len(fit_config$n_pf_evals),
      logLik = NA_real_,
      status = "not_run_invalid_estimate",
      stringsAsFactors = FALSE
    )
    next
  }

  likelihood <- evaluate_parameter_vector(
    model = bspline_model,
    params = fitted_parameters,
    Np = fit_config$Np_eval,
    n_evals = fit_config$n_pf_evals,
    seed_base = evaluation_seed_base
  )
  evaluation_rows[[s]] <- data.frame(
    run = s,
    evaluation_rep = seq_len(fit_config$n_pf_evals),
    evaluation_seed = evaluation_seed_base + seq_len(fit_config$n_pf_evals),
    logLik = likelihood$replicates,
    status = ifelse(is.finite(likelihood$replicates), "success", "failed"),
    stringsAsFactors = FALSE
  )

  complete_evaluation <- likelihood$n_successful == fit_config$n_pf_evals &&
    is.finite(likelihood$logLik)
  fitted_coefficients <- fitted_parameters[coefficient_names]
  selection_rows[[s]] <- cbind(
    base_row,
    as.data.frame(as.list(setNames(
      as.numeric(fitted_coefficients), paste0(coefficient_names, "_hat")
    ))),
    independently_evaluated_logLik = likelihood$logLik,
    logLik_se = likelihood$logLik_se,
    n_successful_pf_evals = likelihood$n_successful,
    eligible_for_selection = complete_evaluation,
    status = if (complete_evaluation) "success" else "incomplete_evaluation",
    error_message = NA_character_
  )
  mif_objects[[s]] <- mif_now

  cat(
    "Finished internal start ", s,
    ": independently evaluated logLik=",
    round(likelihood$logLik, 4),
    ", successful evaluations=", likelihood$n_successful, "/",
    fit_config$n_pf_evals, "\n", sep = ""
  )
}

selection_audit <- do.call(rbind, selection_rows)
likelihood_evaluations <- do.call(rbind, evaluation_rows)
eligible <- which(
  selection_audit$eligible_for_selection &
    is.finite(selection_audit$independently_evaluated_logLik)
)
if (length(eligible) == 0L) {
  write.csv(
    selection_audit,
    file.path(atomic$temp_dir, "start_selection_audit.csv"),
    row.names = FALSE
  )
  write.csv(
    likelihood_evaluations,
    file.path(atomic$temp_dir, "likelihood_evaluations.csv"),
    row.names = FALSE
  )
  stop(
    "No start completed all independent likelihood evaluations for task ",
    task_id, "."
  )
}

selected_index <- eligible[[which.max(
  selection_audit$independently_evaluated_logLik[eligible]
)]]
selected_run <- selection_audit$run[[selected_index]]
best_mif <- mif_objects[[selected_run]]
best_parameters <- coef(best_mif)
best_coefficients <- best_parameters[coefficient_names]

B_path <- data.frame(
  week = observed_data$week,
  B_estimate = reconstruct_B(
    best_coefficients, observed_data$week, experiment_config
  ),
  B_true = true_B_at_times(observed_data$week, experiment_config),
  path_semantics = "deterministic_selected_bspline_trajectory",
  stringsAsFactors = FALSE
)
B_rmse <- sqrt(mean((B_path$B_estimate - B_path$B_true)^2))

best_summary <- data.frame(
  best_run = selected_run,
  start_id = selection_audit$start_id[[selected_index]],
  start_source = selection_audit$start_source[[selected_index]],
  as.list(setNames(
    as.numeric(best_coefficients), paste0(coefficient_names, "_hat")
  )),
  logLik = selection_audit$independently_evaluated_logLik[[selected_index]],
  logLik_se = selection_audit$logLik_se[[selected_index]],
  n_successful_pf_evals =
    selection_audit$n_successful_pf_evals[[selected_index]],
  B_rmse = B_rmse,
  n_internal_starts = nrow(start_values),
  selection_rule = "maximum_independently_evaluated_logLik",
  status = "success",
  check.names = FALSE
)

# Retain only the evidence needed to audit start selection. Alternate fitted
# coefficient vectors are deliberately not persisted as fit results.
selection_audit$selected <- selection_audit$run == selected_run
selection_audit <- selection_audit[, setdiff(
  names(selection_audit), paste0(coefficient_names, "_hat")
), drop = FALSE]

common <- list(
  task_id = task_id,
  simulation_seed = simulation_seed,
  simulation_attempt = simulation_attempt,
  observed_data_md5 = shared$observed_data_md5,
  model = "bspline_B",
  Nmif = fit_config$Nmif
)
add_common <- function(x) {
  for (name in rev(names(common))) {
    x <- cbind(setNames(data.frame(common[[name]]), name), x)
  }
  x
}

selection_audit <- add_common(selection_audit)
likelihood_evaluations <- add_common(likelihood_evaluations)
best_summary <- add_common(best_summary)
B_path <- add_common(B_path)

write.csv(
  selection_audit,
  file.path(atomic$temp_dir, "start_selection_audit.csv"),
  row.names = FALSE
)
write.csv(
  likelihood_evaluations,
  file.path(atomic$temp_dir, "likelihood_evaluations.csv"),
  row.names = FALSE
)
write.csv(
  best_summary,
  file.path(atomic$temp_dir, "best_fit_summary.csv"),
  row.names = FALSE
)
write.csv(
  B_path,
  file.path(atomic$temp_dir, "B_path.csv"),
  row.names = FALSE
)

# This is the only fitted MIF2 object retained for the task.
saveRDS(best_mif, file.path(atomic$temp_dir, "best_mif2.rds"))

if (task_id %in% experiment_config$diagnostic_task_ids) {
  selected_trace <- safe_trace_data(
    best_mif,
    task_id = task_id,
    run = selected_run,
    start_values = as.list(
      selection_audit[selected_index, paste0("start_", coefficient_names)]
    )
  )
  if (!is.null(selected_trace)) {
    selected_trace$simulation_seed <- simulation_seed
    selected_trace$observed_data_md5 <- shared$observed_data_md5
    selected_trace$Nmif <- fit_config$Nmif
    write.csv(
      selected_trace,
      file.path(atomic$temp_dir, "selected_mif2_trace.csv"),
      row.names = FALSE
    )
  }
}

write_run_config(
  file.path(atomic$temp_dir, "run_config.csv"),
  list(
    experiment_id = experiment_config$experiment_id,
    source_experiment_id = experiment_config$source_experiment_id,
    model = "bspline_B",
    task_id = task_id,
    simulation_seed = simulation_seed,
    observed_data_md5 = shared$observed_data_md5,
    Nmif = fit_config$Nmif,
    Np_mif = fit_config$Np_mif,
    n_internal_starts = nrow(start_values),
    Np_eval = fit_config$Np_eval,
    n_pf_evals_per_start = fit_config$n_pf_evals,
    selected_run = selected_run,
    selection_rule = "maximum_independently_evaluated_logLik",
    start_seed = start_seed,
    mif_seed_formula = "base + 1000000*task_id + 1000*start_id",
    evaluation_seed_formula =
      "base + 1000000*task_id + 1000*start_id + evaluation_rep"
  )
)

commit_atomic_task(atomic$temp_dir, atomic$final_dir)
cat(
  "Selected task ", task_id, " start ", selected_run,
  " with independently evaluated logLik=", round(best_summary$logLik[[1]], 4),
  ".\nSaved one final B-spline fit to ", atomic$final_dir, "\n",
  sep = ""
)

# ============================================================
# Fit the Gamma-noise model to one shared Experiment 4 data set
#
# Usage:
# Rscript code/02_run_gamma_task.R \
#   <task_id> <shared_data_root> <output_root> <paramlist_file>
# ============================================================

library(pomp)
options(stringsAsFactors = FALSE)

source("config/experiment_config.R")
source("code/model_components.R")
source("code/io_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript code/02_run_gamma_task.R <task_id> <shared_data_root> <output_root> <paramlist_file>")
}

task_id <- suppressWarnings(as.integer(args[[1]]))
if (is.na(task_id) || task_id < 1L || task_id > experiment_config$n_tasks) {
  stop("task_id must be between 1 and ", experiment_config$n_tasks, ".")
}
shared_data_root <- if (length(args) >= 2L) args[[2]] else "shared_data"
output_root <- if (length(args) >= 3L) args[[3]] else file.path("results_raw", "gamma")
paramlist_file <- if (length(args) >= 4L) args[[4]] else file.path("results", "paramlist.csv")

seed_row <- read_paramlist_task(
  paramlist_file,
  task_id,
  c(
    "task_id", "simulation_seed", "gamma_mif_seed_base",
    "gamma_evaluation_seed_base", "gamma_final_pf_seed"
  )
)

simulation_seed <- as.integer(seed_row$simulation_seed[[1]])
mif_seed_base <- as.integer(seed_row$gamma_mif_seed_base[[1]])
evaluation_seed_base <- as.integer(seed_row$gamma_evaluation_seed_base[[1]])
final_pf_seed <- as.integer(seed_row$gamma_final_pf_seed[[1]])

shared <- read_shared_data(shared_data_root, task_id)
if (as.integer(shared$metadata$simulation_seed[[1]]) != simulation_seed) {
  stop("Simulation seed in shared data does not match paramlist for task ", task_id, ".")
}

atomic <- start_atomic_task(output_root, task_id)
if (atomic$skip) quit(save = "no", status = 0L)
on.exit({
  if (dir.exists(atomic$temp_dir)) unlink(atomic$temp_dir, recursive = TRUE, force = TRUE)
}, add = TRUE)

config <- experiment_config
observed_data <- shared$observed_data
gamma_model <- make_gamma_model(observed_data, config)
theta_gamma <- gamma_baseline_parameters(config)
start_values <- config$gamma_start_values
mif_rw_sd <- rw_sd(
  B0 = ivp(config$gamma_rw_sd_B0_ivp),
  sigma_beta = config$gamma_rw_sd_sigma_beta
)
save_traces <- task_id %in% config$diagnostic_task_ids

cat(
  "Experiment = ", config$experiment_id,
  "\nModel = Gamma-noise model",
  "\nTask ID = ", task_id,
  "\nObserved-data MD5 = ", shared$observed_data_md5,
  "\nNmif = ", config$Nmif,
  "\nNp_mif = ", config$Np_mif,
  "\nNp_eval = ", config$Np_eval,
  "\nn_pf_evals = ", config$n_pf_evals,
  "\nNp_final = ", config$Np_final,
  "\nSave convergence traces = ", save_traces,
  "\n", sep = ""
)

mif_rows <- vector("list", nrow(start_values))
eval_rows <- vector("list", nrow(start_values))
trace_rows <- vector("list", nrow(start_values))

for (s in seq_len(nrow(start_values))) {
  theta_start <- theta_gamma
  theta_start[["B0"]] <- start_values$B0[[s]]
  theta_start[["sigma_beta"]] <- start_values$sigma_beta[[s]]

  cat(
    "Running Gamma-noise model start ", s, "/", nrow(start_values),
    ": B0=", theta_start[["B0"]],
    ", sigma_beta=", theta_start[["sigma_beta"]], "\n", sep = ""
  )

  set.seed(mif_seed_base)
  mif_now <- tryCatch(
    mif2(
      gamma_model,
      params = theta_start,
      Np = config$Np_mif,
      Nmif = config$Nmif,
      rw.sd = mif_rw_sd,
      cooling.type = config$cooling_type,
      cooling.fraction.50 = config$cooling_fraction_50
    ),
    error = function(e) structure(
      list(message = conditionMessage(e)),
      class = "mif2_error"
    )
  )

  if (inherits(mif_now, "mif2_error")) {
    mif_rows[[s]] <- data.frame(
      run = s,
      start_B0 = theta_start[["B0"]],
      start_sigma_beta = theta_start[["sigma_beta"]],
      B0_hat = NA_real_, sigma_beta_hat = NA_real_, Beta_hat = NA_real_,
      logLik = NA_real_, logLik_se = NA_real_,
      n_successful_pf_evals = 0L,
      status = "mif2_failed",
      error_message = mif_now$message,
      stringsAsFactors = FALSE
    )
    eval_rows[[s]] <- data.frame(
      run = s, evaluation_rep = seq_len(config$n_pf_evals),
      logLik = NA_real_, status = "not_run_mif2_failed",
      stringsAsFactors = FALSE
    )
    next
  }

  if (save_traces) {
    trace_rows[[s]] <- safe_trace_data(
      mif_now, task_id, "gamma_noise", s,
      list(B0 = theta_start[["B0"]], sigma_beta = theta_start[["sigma_beta"]])
    )
  }

  mif_coef <- tryCatch(coef(mif_now), error = function(e) NULL)
  valid_coef <- !is.null(mif_coef) &&
    is.finite(mif_coef[["B0"]]) && is.finite(mif_coef[["sigma_beta"]]) &&
    mif_coef[["B0"]] > 0 && mif_coef[["sigma_beta"]] > 0

  if (!valid_coef) {
    mif_rows[[s]] <- data.frame(
      run = s,
      start_B0 = theta_start[["B0"]],
      start_sigma_beta = theta_start[["sigma_beta"]],
      B0_hat = NA_real_, sigma_beta_hat = NA_real_, Beta_hat = NA_real_,
      logLik = NA_real_, logLik_se = NA_real_,
      n_successful_pf_evals = 0L,
      status = "invalid_mif2_estimate",
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    eval_rows[[s]] <- data.frame(
      run = s, evaluation_rep = seq_len(config$n_pf_evals),
      logLik = NA_real_, status = "not_run_invalid_estimate",
      stringsAsFactors = FALSE
    )
    next
  }

  theta_eval <- theta_gamma
  theta_eval[["B0"]] <- unname(mif_coef[["B0"]])
  theta_eval[["sigma_beta"]] <- unname(mif_coef[["sigma_beta"]])

  eval_logLik <- rep(NA_real_, config$n_pf_evals)
  eval_status <- rep("failed", config$n_pf_evals)

  for (j in seq_len(config$n_pf_evals)) {
    set.seed(evaluation_seed_base + j)
    pf_eval <- tryCatch(
      pfilter(gamma_model, params = theta_eval, Np = config$Np_eval),
      error = function(e) NULL
    )
    if (!is.null(pf_eval)) {
      eval_logLik[[j]] <- as.numeric(logLik(pf_eval))
      if (is.finite(eval_logLik[[j]])) eval_status[[j]] <- "success"
    }
  }

  eval_rows[[s]] <- data.frame(
    run = s,
    evaluation_rep = seq_len(config$n_pf_evals),
    logLik = eval_logLik,
    status = eval_status,
    stringsAsFactors = FALSE
  )

  finite_eval <- eval_logLik[is.finite(eval_logLik)]
  if (length(finite_eval) == 0L) {
    mif_rows[[s]] <- data.frame(
      run = s,
      start_B0 = theta_start[["B0"]],
      start_sigma_beta = theta_start[["sigma_beta"]],
      B0_hat = theta_eval[["B0"]],
      sigma_beta_hat = theta_eval[["sigma_beta"]],
      Beta_hat = NA_real_,
      logLik = NA_real_, logLik_se = NA_real_,
      n_successful_pf_evals = 0L,
      status = "all_evaluation_pfilters_failed",
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }

  likelihood_summary <- pomp::logmeanexp(finite_eval, se = TRUE)
  mif_rows[[s]] <- data.frame(
    run = s,
    start_B0 = theta_start[["B0"]],
    start_sigma_beta = theta_start[["sigma_beta"]],
    B0_hat = theta_eval[["B0"]],
    sigma_beta_hat = theta_eval[["sigma_beta"]],
    Beta_hat = NA_real_,
    logLik = as.numeric(likelihood_summary[[1]]),
    logLik_se = as.numeric(likelihood_summary[[2]]),
    n_successful_pf_evals = length(finite_eval),
    status = "success",
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )

  cat(
    "Finished start ", s,
    ": B0_hat=", round(theta_eval[["B0"]], 5),
    ", sigma_beta_hat=", round(theta_eval[["sigma_beta"]], 5),
    ", independent logLik=", round(mif_rows[[s]]$logLik[[1]], 4), "\n",
    sep = ""
  )
}

mif_results <- do.call(rbind, mif_rows)
evaluation_results <- do.call(rbind, eval_rows)
valid <- mif_results[is.finite(mif_results$logLik), , drop = FALSE]

empty_path <- data.frame(
  week = numeric(0), B_estimate = numeric(0), B_true = numeric(0)
)

if (nrow(valid) == 0L) {
  best_summary <- data.frame(
    best_run = NA_integer_, start_B0 = NA_real_, start_sigma_beta = NA_real_,
    start_Beta = NA_real_, B0_hat = NA_real_, sigma_beta_hat = NA_real_,
    Beta_hat = NA_real_, logLik = NA_real_, logLik_se = NA_real_,
    n_successful_pf_evals = 0L, final_pf_logLik = NA_real_,
    fit_success = FALSE, final_pf_success = FALSE,
    status = "all_mif2_starts_failed", stringsAsFactors = FALSE
  )
  B_path <- empty_path
} else {
  best <- valid[which.max(valid$logLik), , drop = FALSE]
  theta_best <- theta_gamma
  theta_best[["B0"]] <- best$B0_hat[[1]]
  theta_best[["sigma_beta"]] <- best$sigma_beta_hat[[1]]

  set.seed(final_pf_seed)
  pf_best <- tryCatch(
    pfilter(
      gamma_model,
      params = theta_best,
      Np = config$Np_final,
      filter.mean = TRUE
    ),
    error = function(e) structure(
      list(message = conditionMessage(e)), class = "final_pf_error"
    )
  )

  if (inherits(pf_best, "final_pf_error")) {
    B_path <- empty_path
    final_pf_logLik <- NA_real_
    final_pf_success <- FALSE
    final_status <- "final_pf_failed"
  } else {
    fm <- filter_mean(pf_best)
    B_path <- data.frame(
      week = as.numeric(time(pf_best)),
      B_estimate = as.numeric(fm["B", ])
    )
    B_path$B_true <- true_B_at_times(B_path$week, config)
    final_pf_logLik <- as.numeric(logLik(pf_best))
    final_pf_success <- TRUE
    final_status <- "success"
  }

  best_summary <- data.frame(
    best_run = best$run[[1]],
    start_B0 = best$start_B0[[1]],
    start_sigma_beta = best$start_sigma_beta[[1]],
    start_Beta = NA_real_,
    B0_hat = best$B0_hat[[1]],
    sigma_beta_hat = best$sigma_beta_hat[[1]],
    Beta_hat = NA_real_,
    logLik = best$logLik[[1]],
    logLik_se = best$logLik_se[[1]],
    n_successful_pf_evals = best$n_successful_pf_evals[[1]],
    final_pf_logLik = final_pf_logLik,
    fit_success = TRUE,
    final_pf_success = final_pf_success,
    status = final_status,
    stringsAsFactors = FALSE
  )
}

common <- list(
  task_id = task_id,
  simulation_seed = simulation_seed,
  simulation_attempt = as.integer(shared$metadata$simulation_attempt[[1]]),
  observed_data_md5 = shared$observed_data_md5,
  model = "gamma_noise",
  Nmif = config$Nmif
)

add_common <- function(x) {
  if (nrow(x) == 0L) {
    common_empty <- as.data.frame(
      lapply(common, function(value) value[FALSE]),
      stringsAsFactors = FALSE
    )
    return(cbind(common_empty, x))
  }

  for (name in rev(names(common))) {
    x <- cbind(setNames(data.frame(common[[name]]), name), x)
  }
  x
}

mif_results <- add_common(mif_results)
evaluation_results <- add_common(evaluation_results)
best_summary <- add_common(best_summary)
B_path <- add_common(B_path)

write.csv(mif_results, file.path(atomic$temp_dir, "mif2_results.csv"), row.names = FALSE)
write.csv(evaluation_results, file.path(atomic$temp_dir, "evaluation_logliks.csv"), row.names = FALSE)
write.csv(best_summary, file.path(atomic$temp_dir, "best_fit_summary.csv"), row.names = FALSE)
write.csv(B_path, file.path(atomic$temp_dir, "B_path.csv"), row.names = FALSE)

trace_non_null <- Filter(Negate(is.null), trace_rows)
if (length(trace_non_null) > 0L) {
  traces_all <- do.call(rbind, trace_non_null)
  traces_all$simulation_seed <- simulation_seed
  traces_all$observed_data_md5 <- shared$observed_data_md5
  traces_all$Nmif <- config$Nmif
  write.csv(traces_all, file.path(atomic$temp_dir, "mif2_traces.csv"), row.names = FALSE)
}

write_run_config(
  file.path(atomic$temp_dir, "run_config.csv"),
  list(
    experiment_id = config$experiment_id,
    model = "gamma_noise",
    task_id = task_id,
    simulation_seed = simulation_seed,
    observed_data_md5 = shared$observed_data_md5,
    Nmif = config$Nmif,
    Np_mif = config$Np_mif,
    Np_eval = config$Np_eval,
    n_pf_evals = config$n_pf_evals,
    Np_final = config$Np_final,
    starts_B0 = paste(unique(start_values$B0), collapse = ","),
    starts_sigma_beta = paste(unique(start_values$sigma_beta), collapse = ","),
    rw_sd_B0_ivp = config$gamma_rw_sd_B0_ivp,
    rw_sd_sigma_beta = config$gamma_rw_sd_sigma_beta,
    cooling_type = config$cooling_type,
    cooling_fraction_50 = config$cooling_fraction_50,
    mif_seed_base = mif_seed_base,
    evaluation_seed_base = evaluation_seed_base,
    final_pf_seed = final_pf_seed
  )
)

commit_atomic_task(atomic$temp_dir, atomic$final_dir)
cat("Gamma-noise model task status: ", best_summary$status[[1]], "\nSaved to ", atomic$final_dir, "\n", sep = "")

# ============================================================
# Experiment 5 configuration
#
# One paired B-spline analysis on the 200 accepted observed
# data sets already created and accepted by Experiment 4.
# No simulation is performed by Experiment 5.
# ============================================================

experiment_config <- list(
  experiment_id = "experiment_5_bspline_B_recovery",
  source_experiment_id = "experiment_4_nmif600_model_comparison",
  n_tasks = 200L,
  diagnostic_task_ids = c(1L, 50L, 100L, 150L, 200L),

  n_weeks = 10,
  observation_dt = 1 / 7,
  process_dt = 1 / 30,
  initial_infectious = 10,
  truth = c(
    Beta_high = 4,
    Beta_low = 2,
    t_switch = 5
  ),
  fixed_parameters = c(
    mu_IR = 3,
    N = 10000,
    rho = 0.5,
    k = 10
  ),

  spline = list(
    nbasis = 6L,
    degree = 3L,
    range = c(0, 10),
    rw_sd = 0.05,
    start_B_range = c(1.5, 5.5),
    start_deviation_sd = 0.35,
    coefficient_limits = log(c(0.5, 8))
  ),

  # All starts belong to the same simulation task. Seeds are separated by
  # task and start so every MIF2 run and every likelihood evaluation is an
  # independent optimization/evaluation attempt, not another replicate.
  seeds = list(
    starts = 20261100L,
    mif2 = 20261200L,
    evaluation = 20261300L
  ),

  Nmif = 600L,
  Np_mif = 5000L,
  n_start = 10L,
  Np_eval = 50000L,
  n_pf_evals = 5L,
  cooling_type = "geometric",
  cooling_fraction_50 = 0.5
)

read_positive_integer_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)

  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed <= 0L) {
    stop(name, " must be a positive integer when it is set.")
  }
  parsed
}

get_fit_config <- function(config = experiment_config) {
  list(
    Nmif = read_positive_integer_env("EXP5_NMIF", config$Nmif),
    Np_mif = read_positive_integer_env("EXP5_NP_MIF", config$Np_mif),
    n_start = read_positive_integer_env("EXP5_N_START", config$n_start),
    Np_eval = read_positive_integer_env("EXP5_NP_EVAL", config$Np_eval),
    n_pf_evals = read_positive_integer_env(
      "EXP5_N_PF_EVALS", config$n_pf_evals
    )
  )
}

task_start_seed <- function(base_seed, task_id, start_id = 0L) {
  value <- as.double(base_seed) + 1000000 * as.double(task_id) +
    1000 * as.double(start_id)
  if (!is.finite(value) || value > .Machine$integer.max) {
    stop("Derived task/start seed exceeds the R integer range.")
  }
  as.integer(value)
}

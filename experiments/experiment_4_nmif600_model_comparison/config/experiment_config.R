# ============================================================
# Experiment 4 configuration
#
# Scientific question:
# Gamma-noise vs. constant-B model comparison,
# both fitted with Nmif = 600, on exactly the same 200 accepted
# simulated epidemic data sets.
# ============================================================

experiment_config <- list(
  experiment_id = "experiment_4_nmif600_model_comparison",
  n_tasks = 200L,
  diagnostic_task_ids = c(1L, 50L, 100L, 150L, 200L),

  n_weeks = 10,
  observation_interval = 1 / 7,
  process_delta_t = 1 / 30,

  true_parameters = c(
    Beta_high = 4,
    Beta_low = 2,
    t_switch = 5,
    mu_IR = 3,
    N = 10000,
    rho = 0.5,
    k = 10
  ),

  acceptance_threshold = 20,
  max_simulation_attempts = 10000L,

  Nmif = 600L,
  Np_mif = 5000L,
  Np_eval = 50000L,
  n_pf_evals = 5L,
  Np_final = 50000L,

  cooling_type = "geometric",
  cooling_fraction_50 = 0.5,

  gamma_start_values = expand.grid(
    B0 = c(2, 4, 6),
    sigma_beta = c(0.10, 0.30, 0.45)
  ),
  gamma_rw_sd_B0_ivp = 0.20,
  gamma_rw_sd_sigma_beta = 0.05,

  constant_start_values = data.frame(
    Beta = c(1, 2, 3, 4, 5, 6)
  ),
  constant_rw_sd_Beta = 0.05
)

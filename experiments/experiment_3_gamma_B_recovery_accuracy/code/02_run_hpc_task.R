# ============================================================
# Simulation study: 200 accepted simulation replicates fitted with MIF2
# HPC array-job version
#
# One Slurm array task:
#   1. generates one accepted simulated epidemic data set;
#   2. fits the Gamma-transition model from 9 mif2 starts;
#   3. evaluates every fitted parameter vector with repeated pfilters;
#   4. selects the best start;
#   5. saves numerical CSV outputs only.
#
# Usage:
# Rscript code/02_run_hpc_task.R <task_id> <output_folder> <paramlist_file>
#
# Example:
# Rscript code/02_run_hpc_task.R 1 Results results/paramlist.csv
# ============================================================

library(pomp)

options(
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 1. Read array-task arguments and task-specific seeds
# ------------------------------------------------------------

args <- commandArgs(
  trailingOnly = TRUE
)

if (length(args) < 1) {
  stop(
    paste0(
      "Usage: Rscript code/02_run_hpc_task.R ",
      "<task_id> <output_folder> <paramlist_file>"
    )
  )
}

task_id <- suppressWarnings(
  as.integer(args[[1]])
)

if (is.na(task_id) || task_id < 1) {
  stop("task_id must be a positive integer.")
}

if (length(args) >= 2) {
  output_folder <- args[[2]]
} else {
  output_folder <- "Results"
}

if (length(args) >= 3) {
  paramlist_file <- args[[3]]
} else {
  paramlist_file <- "paramlist.csv"
}

if (!file.exists(paramlist_file)) {
  stop(
    "Could not find paramlist file: ",
    paramlist_file
  )
}

paramlist <- read.csv(
  paramlist_file,
  check.names = FALSE
)

required_paramlist_columns <- c(
  "task_id",
  "simulation_seed",
  "mif_seed_base",
  "evaluation_seed_base",
  "final_pf_seed"
)

missing_paramlist_columns <- setdiff(
  required_paramlist_columns,
  names(paramlist)
)

if (length(missing_paramlist_columns) > 0) {
  stop(
    "paramlist.csv is missing column(s): ",
    paste(missing_paramlist_columns, collapse = ", ")
  )
}

task_row <- paramlist[
  paramlist$task_id == task_id,
  ,
  drop = FALSE
]

if (nrow(task_row) != 1) {
  stop(
    "Expected exactly one row for task_id = ",
    task_id,
    " in ",
    paramlist_file,
    "."
  )
}

simulation_seed <- as.integer(
  task_row$simulation_seed[[1]]
)

mif_seed_base <- as.integer(
  task_row$mif_seed_base[[1]]
)

evaluation_seed_base <- as.integer(
  task_row$evaluation_seed_base[[1]]
)

final_pf_seed <- as.integer(
  task_row$final_pf_seed[[1]]
)

all_seeds <- c(
  simulation_seed,
  mif_seed_base,
  evaluation_seed_base,
  final_pf_seed
)

if (any(!is.finite(all_seeds))) {
  stop("At least one task-specific seed is invalid.")
}

task_output_folder <- file.path(
  output_folder,
  sprintf("task_%03d", task_id)
)

dir.create(
  task_output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "Task ID = ", task_id,
  "\nOutput folder = ", task_output_folder,
  "\nSimulation seed = ", simulation_seed,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 2. Define the time grid for the simulated data
# ------------------------------------------------------------

n_weeks <- 10

template <- data.frame(
  week = seq(
    from = 1 / 7,
    to = n_weeks,
    by = 1 / 7
  ),
  reports = 0
)

# ------------------------------------------------------------
# 2.1 Define the data-generating process model
# ------------------------------------------------------------

sir_step <- Csnippet("
  double Beta_now;
  double dN_SI;
  double dN_IR;
  double p_SI;
  double p_IR;

  // Piecewise transmission rate
  if (t < t_switch) {
    Beta_now = Beta_high;
  } else {
    Beta_now = Beta_low;
  }

  // Transition probabilities
  p_SI = 1.0 - exp(-Beta_now * I / N * dt);
  p_IR = 1.0 - exp(-mu_IR * dt);

  // New infections and recoveries
  dN_SI = rbinom(S, p_SI);
  dN_IR = rbinom(I, p_IR);

  // Update hidden states
  S = S - dN_SI;
  I = I + dN_SI - dN_IR;
  R = R + dN_IR;

  // Accumulate new infections
  H = H + dN_SI;
")

sir_rinit_piecewise <- Csnippet("
  S = N - 10;
  I = 10;
  R = 0;
  H = 0;
")

sir_rmeas <- Csnippet("
  reports = rnbinom_mu(k, rho * H);
")

sir_dmeas <- Csnippet("
  lik = dnbinom_mu(
    reports,
    k,
    rho * H,
    give_log
  );
")

simple_SIR <- pomp(
  data = template,
  times = "week",
  t0 = 0,

  rinit = sir_rinit_piecewise,

  rprocess = euler(
    sir_step,
    delta.t = 1 / 30
  ),

  rmeasure = sir_rmeas,
  dmeasure = sir_dmeas,

  accumvars = "H",

  statenames = c(
    "S",
    "I",
    "R",
    "H"
  ),

  paramnames = c(
    "Beta_high",
    "Beta_low",
    "t_switch",
    "mu_IR",
    "N",
    "rho",
    "k"
  )
)

# ------------------------------------------------------------
# 3. Set the true data-generating parameters
# ------------------------------------------------------------

theta <- c(
  Beta_high = 4,
  Beta_low = 2,
  t_switch = 5,
  mu_IR = 3,
  N = 10000,
  rho = 0.5,
  k = 10
)

# ------------------------------------------------------------
# 4. Generate one accepted fake epidemic data set
# ------------------------------------------------------------
#
# Keep generating fake epidemics until:
#
#   max(H) > 20
#
# This means each task contributes one accepted epidemic,
# not merely one attempted simulation.

acceptance_threshold <- 20

max_simulation_attempts <- 10000L

set.seed(
  simulation_seed
)

accepted <- FALSE

simulation_attempt <- 0L

simulated_data <- NULL

while (
  !accepted &&
    simulation_attempt < max_simulation_attempts
) {
  simulation_attempt <- simulation_attempt + 1L

  simulated_data <- simulate(
    simple_SIR,
    params = theta,
    nsim = 1,
    format = "data.frame",
    include.data = FALSE
  )

  accepted <- is.finite(
    max(simulated_data$H)
  ) &&
    max(simulated_data$H) > acceptance_threshold
}

if (!accepted) {
  stop(
    "Task ",
    task_id,
    " did not generate an accepted simulation after ",
    max_simulation_attempts,
    " attempts."
  )
}

observed_data <- simulated_data[
  ,
  c("week", "reports"),
  drop = FALSE
]

cat(
  "Accepted simulation after ",
  simulation_attempt,
  " attempt(s). Maximum H = ",
  max(simulated_data$H),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 5. Define the Gamma-transition filtering model
# ------------------------------------------------------------

sir_step_gamma <- Csnippet("
  double shape_B;
  double scale_B;

  double dN_SI;
  double dN_IR;

  double p_SI;
  double p_IR;

  // Update the hidden transmission rate B
  if (sigma_beta > 0.0) {

    shape_B =
      1.0 /
      (sigma_beta * sigma_beta * dt);

    scale_B =
      B *
      sigma_beta *
      sigma_beta *
      dt;

    B = rgamma(
      shape_B,
      scale_B
    );
  }

  // Transition probabilities
  p_SI = 1.0 - exp(-B * I / N * dt);
  p_IR = 1.0 - exp(-mu_IR * dt);

  // New infections and recoveries
  dN_SI = rbinom(S, p_SI);
  dN_IR = rbinom(I, p_IR);

  // Update SIR states
  S = S - dN_SI;
  I = I + dN_SI - dN_IR;
  R = R + dN_IR;

  // Accumulate new infections
  H = H + dN_SI;
")

sir_rinit_gamma <- Csnippet("
  S = N - 10;
  I = 10;
  R = 0;
  H = 0;

  // Initial hidden transmission rate
  B = B0;
")

theta_gamma <- c(
  B0 = 4,
  sigma_beta = 0.2,
  mu_IR = 3,
  N = 10000,
  rho = 0.5,
  k = 10
)

pf_gamma_model <- pomp(
  data = observed_data,
  times = "week",
  t0 = 0,

  rinit = sir_rinit_gamma,

  rprocess = euler(
    sir_step_gamma,
    delta.t = 1 / 30
  ),

  rmeasure = sir_rmeas,
  dmeasure = sir_dmeas,

  accumvars = "H",

  statenames = c(
    "S",
    "I",
    "R",
    "H",
    "B"
  ),

  paramnames = c(
    "B0",
    "sigma_beta",
    "mu_IR",
    "N",
    "rho",
    "k"
  ),

  partrans = parameter_trans(
    log = c(
      "B0",
      "sigma_beta"
    )
  )
)

# ------------------------------------------------------------
# 6. Set up the formal mif2 search
# ------------------------------------------------------------

start_values <- expand.grid(
  B0 = c(2, 4, 6),
  sigma_beta = c(0.10, 0.30, 0.45)
)

# Number of iterated-filtering iterations
Nmif <- 100

# Number of particles inside mif2
Np_mif <- 5000

# Number of particles used to evaluate each fitted parameter vector
Np_eval <- 50000

# Number of repeated pfilter evaluations for each mif2 result
n_pf_evals <- 5

# Number of particles for the final filtered B path
Np_final <- 50000

mif_rw_sd <- rw_sd(
  B0 = ivp(0.20),
  sigma_beta = 0.05
)

# ------------------------------------------------------------
# 7. Fit one fake data set from all starting values
# ------------------------------------------------------------

fit_one_dataset <- function(
  pf_gamma_model,
  theta_gamma,
  theta_true,
  start_values,
  mif_rw_sd,
  Nmif,
  Np_mif,
  Np_eval,
  n_pf_evals,
  Np_final,
  mif_seed_base,
  evaluation_seed_base,
  final_pf_seed
) {

  n_starts <- nrow(
    start_values
  )

  mif_result_list <- vector(
    mode = "list",
    length = n_starts
  )

  for (s in seq_len(n_starts)) {

    theta_start <- theta_gamma

    theta_start[["B0"]] <- start_values$B0[[s]]

    theta_start[["sigma_beta"]] <-
      start_values$sigma_beta[[s]]

    cat(
      "\nRunning mif2 start ",
      s,
      " of ",
      n_starts,
      ": B0 = ",
      theta_start[["B0"]],
      ", sigma_beta = ",
      theta_start[["sigma_beta"]],
      "\n",
      sep = ""
    )

    # Reuse the same task-specific MIF2 seed across starting values.
    # This keeps the seed setting fixed for the multi-start comparison.
    set.seed(
      mif_seed_base
    )

    mif_now <- tryCatch(
      mif2(
        pf_gamma_model,
        params = theta_start,
        Np = Np_mif,
        Nmif = Nmif,
        rw.sd = mif_rw_sd,
        cooling.type = "geometric",
        cooling.fraction.50 = 0.5
      ),
      error = function(e) {
        structure(
          list(
            message = conditionMessage(e)
          ),
          class = "mif2_error"
        )
      }
    )

    if (inherits(mif_now, "mif2_error")) {
      mif_result_list[[s]] <- data.frame(
        run = s,
        start_B0 = theta_start[["B0"]],
        start_sigma_beta = theta_start[["sigma_beta"]],
        B0_hat = NA_real_,
        sigma_beta_hat = NA_real_,
        logLik = NA_real_,
        logLik_se = NA_real_,
        n_successful_pf_evals = 0L,
        status = "mif2_failed",
        error_message = mif_now$message,
        stringsAsFactors = FALSE
      )

      next
    }

    mif_coef_now <- tryCatch(
      coef(
        mif_now
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      is.null(mif_coef_now) ||
        !is.finite(mif_coef_now[["B0"]]) ||
        !is.finite(mif_coef_now[["sigma_beta"]]) ||
        mif_coef_now[["B0"]] <= 0 ||
        mif_coef_now[["sigma_beta"]] <= 0
    ) {
      mif_result_list[[s]] <- data.frame(
        run = s,
        start_B0 = theta_start[["B0"]],
        start_sigma_beta = theta_start[["sigma_beta"]],
        B0_hat = NA_real_,
        sigma_beta_hat = NA_real_,
        logLik = NA_real_,
        logLik_se = NA_real_,
        n_successful_pf_evals = 0L,
        status = "invalid_mif2_estimate",
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )

      next
    }

    theta_eval_now <- theta_gamma

    theta_eval_now[["B0"]] <-
      unname(
        mif_coef_now[["B0"]]
      )

    theta_eval_now[["sigma_beta"]] <-
      unname(
        mif_coef_now[["sigma_beta"]]
      )

    eval_logLik <- rep(
      NA_real_,
      n_pf_evals
    )

    for (j in seq_len(n_pf_evals)) {

      # Reuse the same evaluation-seed sequence across starting values.
      set.seed(
        evaluation_seed_base + j
      )

      pf_eval <- tryCatch(
        pfilter(
          pf_gamma_model,
          params = theta_eval_now,
          Np = Np_eval
        ),
        error = function(e) {
          message(
            "Evaluation pfilter failed for start ",
            s,
            ", repetition ",
            j,
            ": ",
            conditionMessage(e)
          )
          NULL
        }
      )

      if (!is.null(pf_eval)) {
        eval_logLik[[j]] <- as.numeric(
          logLik(pf_eval)
        )
      }
    }

    eval_logLik <- eval_logLik[
      is.finite(eval_logLik)
    ]

    if (length(eval_logLik) == 0) {
      mif_result_list[[s]] <- data.frame(
        run = s,
        start_B0 = theta_start[["B0"]],
        start_sigma_beta = theta_start[["sigma_beta"]],
        B0_hat = theta_eval_now[["B0"]],
        sigma_beta_hat = theta_eval_now[["sigma_beta"]],
        logLik = NA_real_,
        logLik_se = NA_real_,
        n_successful_pf_evals = 0L,
        status = "all_evaluation_pfilters_failed",
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )

      next
    }

    likelihood_summary <- pomp::logmeanexp(
      eval_logLik,
      se = TRUE
    )

    estimated_logLik <- as.numeric(
      likelihood_summary[[1]]
    )

    estimated_logLik_se <- as.numeric(
      likelihood_summary[[2]]
    )

    mif_result_list[[s]] <- data.frame(
      run = s,
      start_B0 = theta_start[["B0"]],
      start_sigma_beta = theta_start[["sigma_beta"]],
      B0_hat = theta_eval_now[["B0"]],
      sigma_beta_hat = theta_eval_now[["sigma_beta"]],
      logLik = estimated_logLik,
      logLik_se = estimated_logLik_se,
      n_successful_pf_evals = length(eval_logLik),
      status = "success",
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )

    cat(
      "Finished start ",
      s,
      ": B0_hat = ",
      round(theta_eval_now[["B0"]], 4),
      ", sigma_beta_hat = ",
      round(theta_eval_now[["sigma_beta"]], 4),
      ", logLik = ",
      round(estimated_logLik, 4),
      ", logLik_se = ",
      round(estimated_logLik_se, 4),
      "\n",
      sep = ""
    )
  }

  mif_results <- do.call(
    rbind,
    mif_result_list
  )

  rownames(
    mif_results
  ) <- NULL

  valid_results <- mif_results[
    is.finite(mif_results$logLik),
    ,
    drop = FALSE
  ]

  empty_B_path <- data.frame(
    week = numeric(0),
    B_filtered_mean = numeric(0),
    B_true = numeric(0)
  )

  if (nrow(valid_results) == 0) {
    best_summary <- data.frame(
      best_run = NA_integer_,
      start_B0 = NA_real_,
      start_sigma_beta = NA_real_,
      B0_hat = NA_real_,
      sigma_beta_hat = NA_real_,
      logLik = NA_real_,
      logLik_se = NA_real_,
      n_successful_pf_evals = 0L,
      final_pf_logLik = NA_real_,
      fit_success = FALSE,
      final_pf_success = FALSE,
      status = "all_mif2_starts_failed",
      stringsAsFactors = FALSE
    )

    return(
      list(
        mif_results = mif_results,
        best_summary = best_summary,
        B_estimate = empty_B_path
      )
    )
  }

  best_row_index <- which.max(
    valid_results$logLik
  )

  best_fit <- valid_results[
    best_row_index,
    ,
    drop = FALSE
  ]

  theta_best <- theta_gamma

  theta_best[["B0"]] <- best_fit$B0_hat[[1]]

  theta_best[["sigma_beta"]] <-
    best_fit$sigma_beta_hat[[1]]

  set.seed(
    final_pf_seed
  )

  pf_best <- tryCatch(
    pfilter(
      pf_gamma_model,
      params = theta_best,
      Np = Np_final,
      filter.mean = TRUE
    ),
    error = function(e) {
      structure(
        list(
          message = conditionMessage(e)
        ),
        class = "final_pf_error"
      )
    }
  )

  if (inherits(pf_best, "final_pf_error")) {
    best_summary <- data.frame(
      best_run = best_fit$run[[1]],
      start_B0 = best_fit$start_B0[[1]],
      start_sigma_beta = best_fit$start_sigma_beta[[1]],
      B0_hat = best_fit$B0_hat[[1]],
      sigma_beta_hat = best_fit$sigma_beta_hat[[1]],
      logLik = best_fit$logLik[[1]],
      logLik_se = best_fit$logLik_se[[1]],
      n_successful_pf_evals =
        best_fit$n_successful_pf_evals[[1]],
      final_pf_logLik = NA_real_,
      fit_success = TRUE,
      final_pf_success = FALSE,
      status = "final_pf_failed",
      stringsAsFactors = FALSE
    )

    return(
      list(
        mif_results = mif_results,
        best_summary = best_summary,
        B_estimate = empty_B_path
      )
    )
  }

  fm <- filter_mean(
    pf_best
  )

  B_estimate <- data.frame(
    week = as.numeric(
      time(pf_best)
    ),
    B_filtered_mean = as.numeric(
      fm["B", ]
    )
  )

  B_estimate$B_true <- ifelse(
    B_estimate$week < theta_true[["t_switch"]],
    theta_true[["Beta_high"]],
    theta_true[["Beta_low"]]
  )

  best_summary <- data.frame(
    best_run = best_fit$run[[1]],
    start_B0 = best_fit$start_B0[[1]],
    start_sigma_beta = best_fit$start_sigma_beta[[1]],
    B0_hat = best_fit$B0_hat[[1]],
    sigma_beta_hat = best_fit$sigma_beta_hat[[1]],
    logLik = best_fit$logLik[[1]],
    logLik_se = best_fit$logLik_se[[1]],
    n_successful_pf_evals =
      best_fit$n_successful_pf_evals[[1]],
    final_pf_logLik = as.numeric(
      logLik(pf_best)
    ),
    fit_success = TRUE,
    final_pf_success = TRUE,
    status = "success",
    stringsAsFactors = FALSE
  )

  list(
    mif_results = mif_results,
    best_summary = best_summary,
    B_estimate = B_estimate
  )
}

fit_output <- fit_one_dataset(
  pf_gamma_model = pf_gamma_model,
  theta_gamma = theta_gamma,
  theta_true = theta,
  start_values = start_values,
  mif_rw_sd = mif_rw_sd,
  Nmif = Nmif,
  Np_mif = Np_mif,
  Np_eval = Np_eval,
  n_pf_evals = n_pf_evals,
  Np_final = Np_final,
  mif_seed_base = mif_seed_base,
  evaluation_seed_base = evaluation_seed_base,
  final_pf_seed = final_pf_seed
)

# ------------------------------------------------------------
# 8. Add identifiers and save task-level outputs
# ------------------------------------------------------------

mif_results <- fit_output$mif_results

mif_results$task_id <- task_id

mif_results$simulation_seed <- simulation_seed

mif_results <- mif_results[
  ,
  c(
    "task_id",
    "simulation_seed",
    "run",
    "start_B0",
    "start_sigma_beta",
    "B0_hat",
    "sigma_beta_hat",
    "logLik",
    "logLik_se",
    "n_successful_pf_evals",
    "status",
    "error_message"
  ),
  drop = FALSE
]

best_fit_summary <- fit_output$best_summary

best_fit_summary$task_id <- task_id

best_fit_summary$simulation_seed <- simulation_seed

best_fit_summary <- best_fit_summary[
  ,
  c(
    "task_id",
    "simulation_seed",
    "best_run",
    "start_B0",
    "start_sigma_beta",
    "B0_hat",
    "sigma_beta_hat",
    "logLik",
    "logLik_se",
    "n_successful_pf_evals",
    "final_pf_logLik",
    "fit_success",
    "final_pf_success",
    "status"
  ),
  drop = FALSE
]

B_estimate <- fit_output$B_estimate

B_estimate$task_id <- task_id

B_estimate$simulation_seed <- simulation_seed

B_estimate <- B_estimate[
  ,
  c(
    "task_id",
    "simulation_seed",
    "week",
    "B_filtered_mean",
    "B_true"
  ),
  drop = FALSE
]

write.csv(
  mif_results,
  file.path(
    task_output_folder,
    "mif2_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  best_fit_summary,
  file.path(
    task_output_folder,
    "best_fit_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  B_estimate,
  file.path(
    task_output_folder,
    "filtered_B_path.csv"
  ),
  row.names = FALSE
)

write.csv(
  simulated_data,
  file.path(
    task_output_folder,
    "simulated_data.csv"
  ),
  row.names = FALSE
)

cat(
  "\nSaved task outputs to: ",
  task_output_folder,
  "\n",
  sep = ""
)

cat(
  "Task status: ",
  best_fit_summary$status[[1]],
  "\n",
  sep = ""
)
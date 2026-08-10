# ============================================================
# Run one MIF2 search for one starting point
#
# Usage:
#   Rscript code/02_Run_MIF2_Array.R <task_id>
#
# task_id must be 1,...,9.
#
# All tasks use the same fixed dataset.
# Each task corresponds to one starting point.
# ============================================================

library(pomp)

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. Read task ID
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    paste(
      "Usage:",
      "Rscript code/02_Run_MIF2_Array.R <task_id>"
    )
  )
}

task_id <- suppressWarnings(
  as.integer(args[[1]])
)

if (
  is.na(task_id) ||
  task_id < 1 ||
  task_id > 9
) {
  stop("task_id must be an integer from 1 to 9.")
}

# ------------------------------------------------------------
# 2. Read the fixed dataset
# ------------------------------------------------------------

data_file <- "data/fixed_piecewise_B_dataset.csv"

if (!file.exists(data_file)) {
  stop(
    paste(
      "Fixed dataset not found:",
      data_file,
      "\nRun code/01_Generate_Fixed_Piecewise_Data.R first."
    )
  )
}

sim1 <- read.csv(
  data_file,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "week",
  "reports"
)

if (!all(required_columns %in% names(sim1))) {
  stop(
    "The fixed dataset must contain columns named week and reports."
  )
}

sim1_observed <- sim1[
  ,
  c(
    "week",
    "reports"
  ),
  drop = FALSE
]

# ------------------------------------------------------------
# 3. Measurement model
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 4. Gamma-noise process model
# ------------------------------------------------------------

sir_step_gamma <- Csnippet("

  double shape_B;
  double scale_B;

  double dN_SI;
  double dN_IR;

  double p_SI;
  double p_IR;

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

  p_SI = 1.0 - exp(-B * I / N * dt);
  p_IR = 1.0 - exp(-mu_IR * dt);

  dN_SI = rbinom(S, p_SI);
  dN_IR = rbinom(I, p_IR);

  S = S - dN_SI;
  I = I + dN_SI - dN_IR;
  R = R + dN_IR;

  H = H + dN_SI;
")

sir_rinit_gamma <- Csnippet("
  S = N - 10;
  I = 10;
  R = 0;
  H = 0;
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
  data = sim1_observed,
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
# 5. Starting points
# ------------------------------------------------------------

start_values <- expand.grid(
  B0 = c(
    2,
    4,
    6
  ),
  sigma_beta = c(
    0.10,
    0.30,
    0.45
  )
)

theta_start <- theta_gamma

theta_start[["B0"]] <-
  start_values$B0[task_id]

theta_start[["sigma_beta"]] <-
  start_values$sigma_beta[task_id]

# ------------------------------------------------------------
# 6. Formal MIF2 settings
# ------------------------------------------------------------

Nmif <- 100
Np_mif <- 50000

Np_eval <- 50000
n_pf_evals <- 5

mif_seed <- 20260628

evaluation_seeds <-
  20260800 +
  seq_len(n_pf_evals)

mif_rw_sd <- rw_sd(
  B0 = ivp(0.20),
  sigma_beta = 0.05
)

cat(
  "Task ID = ",
  task_id,
  "\nStarting B0 = ",
  theta_start[["B0"]],
  "\nStarting sigma_beta = ",
  theta_start[["sigma_beta"]],
  "\nNmif = ",
  Nmif,
  "\nNp_mif = ",
  Np_mif,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Run one MIF2 search
# ------------------------------------------------------------

set.seed(mif_seed)

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
    message(
      "MIF2 failed for task ",
      task_id,
      ": ",
      conditionMessage(e)
    )
    NULL
  }
)

task_output_folder <- file.path(
  "results",
  "array_output",
  sprintf(
    "task_%03d",
    task_id
  )
)

dir.create(
  task_output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

if (is.null(mif_now)) {

  failed_result <- data.frame(
    task_id = task_id,
    mif_seed = mif_seed,
    start_B0 = theta_start[["B0"]],
    start_sigma_beta = theta_start[["sigma_beta"]],
    B0_hat = NA_real_,
    sigma_beta_hat = NA_real_,
    logLik = NA_real_,
    logLik_se = NA_real_,
    successful_pf_evals = 0
  )

  write.csv(
    failed_result,
    file.path(
      task_output_folder,
      "mif2_result.csv"
    ),
    row.names = FALSE
  )

  stop("MIF2 failed. Failure result was saved.")
}

# ------------------------------------------------------------
# 8. Extract fitted parameters
# ------------------------------------------------------------

mif_coef_now <- coef(mif_now)

theta_eval_now <- theta_gamma

theta_eval_now[["B0"]] <-
  unname(
    mif_coef_now[["B0"]]
  )

theta_eval_now[["sigma_beta"]] <-
  unname(
    mif_coef_now[["sigma_beta"]]
  )

# ------------------------------------------------------------
# 9. Repeated pfilter likelihood evaluation
# ------------------------------------------------------------

eval_logLik <- rep(
  NA_real_,
  n_pf_evals
)

for (j in seq_len(n_pf_evals)) {

  set.seed(
    evaluation_seeds[j]
  )

  pf_eval <- tryCatch(
    pfilter(
      pf_gamma_model,
      params = theta_eval_now,
      Np = Np_eval
    ),
    error = function(e) {
      message(
        "Evaluation pfilter failed for task ",
        task_id,
        ", repetition ",
        j,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (!is.null(pf_eval)) {
    eval_logLik[j] <- as.numeric(
      logLik(pf_eval)
    )
  }
}

successful_logLik <- eval_logLik[
  is.finite(eval_logLik)
]

if (length(successful_logLik) > 0) {

  likelihood_summary <- pomp::logmeanexp(
    successful_logLik,
    se = TRUE
  )

  estimated_logLik <- as.numeric(
    likelihood_summary[1]
  )

  estimated_logLik_se <- as.numeric(
    likelihood_summary[2]
  )

} else {

  estimated_logLik <- NA_real_
  estimated_logLik_se <- NA_real_
}

# ------------------------------------------------------------
# 10. Save task output
# ------------------------------------------------------------

result_now <- data.frame(
  task_id = task_id,
  mif_seed = mif_seed,
  start_B0 = theta_start[["B0"]],
  start_sigma_beta = theta_start[["sigma_beta"]],
  B0_hat = theta_eval_now[["B0"]],
  sigma_beta_hat = theta_eval_now[["sigma_beta"]],
  logLik = estimated_logLik,
  logLik_se = estimated_logLik_se,
  successful_pf_evals = length(successful_logLik)
)

write.csv(
  result_now,
  file.path(
    task_output_folder,
    "mif2_result.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    repetition = seq_len(n_pf_evals),
    evaluation_seed = evaluation_seeds,
    logLik = eval_logLik
  ),
  file.path(
    task_output_folder,
    "pfilter_evaluations.csv"
  ),
  row.names = FALSE
)

saveRDS(
  mif_now,
  file.path(
    task_output_folder,
    "mif2_object.rds"
  )
)

cat(
  "\nTask completed successfully.\n",
  "Estimated B0 = ",
  round(theta_eval_now[["B0"]], 6),
  "\nEstimated sigma_beta = ",
  round(theta_eval_now[["sigma_beta"]], 6),
  "\nEstimated logLik = ",
  round(estimated_logLik, 6),
  "\n",
  sep = ""
)

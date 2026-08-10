# ============================================================
# Recreate and save the global-best MIF2 run from Experiment 3
#
# This script:
#   1. reads the combined 200 x 9 MIF2 results;
#   2. identifies the row with the largest finite evaluated log-likelihood;
#   3. regenerates that task's accepted simulated epidemic using the
#      original simulation seed;
#   4. reruns only the corresponding MIF2 start using the original seed;
#   5. saves the reconstructed MIF2 object, a summary CSV, the regenerated
#      simulated data, and a diagnostic PDF.
#
# Usage:
# Rscript code/05_recreate_global_best_mif2.R \
#   results/combined/combined_mif2_results.csv \
#   results/paramlist.csv \
#   results/recreated_mif2 \
#   figures
#
# Important:
# Exact bit-for-bit reproduction may require the same R version,
# pomp version, compiler, and computing environment as the original run.
# ============================================================

library(pomp)

options(
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 1. Read command-line arguments
# ------------------------------------------------------------

args <- commandArgs(
  trailingOnly = TRUE
)

results_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(
    "results",
    "combined",
    "combined_mif2_results.csv"
  )
}

paramlist_file <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path("results", "paramlist.csv")
}

output_folder <- if (length(args) >= 3) {
  args[[3]]
} else {
  file.path("results", "recreated_mif2")
}

figures_folder <- if (length(args) >= 4) {
  args[[4]]
} else {
  "figures"
}

if (!file.exists(results_file)) {
  stop(
    "Could not find combined MIF2 results file: ",
    results_file
  )
}

if (!file.exists(paramlist_file)) {
  stop(
    "Could not find paramlist file: ",
    paramlist_file
  )
}

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Identify the global-best row among all 200 x 9 runs
# ------------------------------------------------------------

all_results <- read.csv(
  results_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_result_columns <- c(
  "task_id",
  "simulation_seed",
  "run",
  "start_B0",
  "start_sigma_beta",
  "B0_hat",
  "sigma_beta_hat",
  "logLik"
)

missing_result_columns <- setdiff(
  required_result_columns,
  names(all_results)
)

if (length(missing_result_columns) > 0) {
  stop(
    "The combined results file is missing column(s): ",
    paste(
      missing_result_columns,
      collapse = ", "
    )
  )
}

valid_results <- all_results[
  is.finite(all_results$logLik),
  ,
  drop = FALSE
]

if (nrow(valid_results) == 0) {
  stop(
    "No finite log-likelihood values were found in ",
    results_file,
    "."
  )
}

global_best <- valid_results[
  which.max(valid_results$logLik),
  ,
  drop = FALSE
]

task_id <- as.integer(
  global_best$task_id[[1]]
)

simulation_seed_from_results <- as.integer(
  global_best$simulation_seed[[1]]
)

best_run <- as.integer(
  global_best$run[[1]]
)

start_B0 <- as.numeric(
  global_best$start_B0[[1]]
)

start_sigma_beta <- as.numeric(
  global_best$start_sigma_beta[[1]]
)

original_B0_hat <- as.numeric(
  global_best$B0_hat[[1]]
)

original_sigma_beta_hat <- as.numeric(
  global_best$sigma_beta_hat[[1]]
)

original_logLik <- as.numeric(
  global_best$logLik[[1]]
)

cat(
  "\nGlobal-best stored result\n",
  "Task ID: ", task_id, "\n",
  "Run: ", best_run, "\n",
  "Starting B0: ", start_B0, "\n",
  "Starting sigma_beta: ", start_sigma_beta, "\n",
  "Stored B0_hat: ", original_B0_hat, "\n",
  "Stored sigma_beta_hat: ", original_sigma_beta_hat, "\n",
  "Stored evaluated logLik: ", original_logLik, "\n",
  sep = ""
)

# Save the selected original row immediately.
write.csv(
  global_best,
  file.path(
    output_folder,
    "selected_original_result.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Read the original task-specific seeds
# ------------------------------------------------------------

paramlist <- read.csv(
  paramlist_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_param_columns <- c(
  "task_id",
  "simulation_seed",
  "mif_seed_base",
  "evaluation_seed_base",
  "final_pf_seed"
)

missing_param_columns <- setdiff(
  required_param_columns,
  names(paramlist)
)

if (length(missing_param_columns) > 0) {
  stop(
    "The paramlist file is missing column(s): ",
    paste(
      missing_param_columns,
      collapse = ", "
    )
  )
}

task_row <- paramlist[
  paramlist$task_id == task_id,
  ,
  drop = FALSE
]

if (nrow(task_row) != 1) {
  stop(
    "Expected exactly one paramlist row for task_id = ",
    task_id,
    "."
  )
}

simulation_seed <- as.integer(
  task_row$simulation_seed[[1]]
)

mif_seed_base <- as.integer(
  task_row$mif_seed_base[[1]]
)

if (
  is.finite(simulation_seed_from_results) &&
    simulation_seed != simulation_seed_from_results
) {
  stop(
    "Simulation-seed mismatch. Results file gives ",
    simulation_seed_from_results,
    ", while paramlist gives ",
    simulation_seed,
    "."
  )
}

# ------------------------------------------------------------
# 4. Define the original data-generating model
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

sir_step <- Csnippet("
  double Beta_now;
  double dN_SI;
  double dN_IR;
  double p_SI;
  double p_IR;

  if (t < t_switch) {
    Beta_now = Beta_high;
  } else {
    Beta_now = Beta_low;
  }

  p_SI = 1.0 - exp(-Beta_now * I / N * dt);
  p_IR = 1.0 - exp(-mu_IR * dt);

  dN_SI = rbinom(S, p_SI);
  dN_IR = rbinom(I, p_IR);

  S = S - dN_SI;
  I = I + dN_SI - dN_IR;
  R = R + dN_IR;

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

data_generating_model <- pomp(
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

theta_true <- c(
  Beta_high = 4,
  Beta_low = 2,
  t_switch = 5,
  mu_IR = 3,
  N = 10000,
  rho = 0.5,
  k = 10
)

# ------------------------------------------------------------
# 5. Regenerate the accepted simulated dataset
# ------------------------------------------------------------

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
    data_generating_model,
    params = theta_true,
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
    "Failed to regenerate an accepted epidemic after ",
    max_simulation_attempts,
    " attempts."
  )
}

observed_data <- simulated_data[
  ,
  c(
    "week",
    "reports"
  ),
  drop = FALSE
]

write.csv(
  simulated_data,
  file.path(
    output_folder,
    "regenerated_simulated_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  observed_data,
  file.path(
    output_folder,
    "regenerated_observed_data.csv"
  ),
  row.names = FALSE
)

cat(
  "\nRegenerated accepted simulation after ",
  simulation_attempt,
  " attempt(s). Maximum H = ",
  max(simulated_data$H),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 6. Define the original Gamma-noise fitting model
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
# 7. Recreate only the selected MIF2 run
# ------------------------------------------------------------

theta_start <- theta_gamma

theta_start[["B0"]] <- start_B0

theta_start[["sigma_beta"]] <- start_sigma_beta

Nmif <- 100

Np_mif <- 5000

mif_rw_sd <- rw_sd(
  B0 = ivp(0.20),
  sigma_beta = 0.05
)

set.seed(
  mif_seed_base
)

recreated_mif2 <- mif2(
  pf_gamma_model,

  params = theta_start,

  Np = Np_mif,

  Nmif = Nmif,

  rw.sd = mif_rw_sd,

  cooling.type = "geometric",

  cooling.fraction.50 = 0.5
)

recreated_coef <- coef(
  recreated_mif2
)

# ------------------------------------------------------------
# 8. Save the recreated object and diagnostics
# ------------------------------------------------------------

rds_file <- file.path(
  output_folder,
  sprintf(
    "global_best_task_%03d_run_%02d_mif2.rds",
    task_id,
    best_run
  )
)

saveRDS(
  recreated_mif2,
  rds_file
)

diagnostic_pdf <- file.path(
  figures_folder,
  sprintf(
    "mif2_diagnostic_task_%03d_run_%02d.pdf",
    task_id,
    best_run
  )
)

pdf(
  diagnostic_pdf,
  width = 8,
  height = 8
)

plot(
  recreated_mif2
)

dev.off()

comparison_summary <- data.frame(
  task_id = task_id,
  run = best_run,
  simulation_seed = simulation_seed,
  mif_seed_base = mif_seed_base,
  simulation_attempt = simulation_attempt,
  max_H = max(simulated_data$H),
  start_B0 = start_B0,
  start_sigma_beta = start_sigma_beta,
  original_B0_hat = original_B0_hat,
  recreated_B0_hat = unname(
    recreated_coef[["B0"]]
  ),
  B0_difference =
    unname(recreated_coef[["B0"]]) -
    original_B0_hat,
  original_sigma_beta_hat =
    original_sigma_beta_hat,
  recreated_sigma_beta_hat =
    unname(recreated_coef[["sigma_beta"]]),
  sigma_beta_difference =
    unname(recreated_coef[["sigma_beta"]]) -
    original_sigma_beta_hat,
  original_evaluated_logLik =
    original_logLik,
  rds_file = basename(rds_file),
  diagnostic_pdf = diagnostic_pdf,
  stringsAsFactors = FALSE
)

write.csv(
  comparison_summary,
  file.path(
    output_folder,
    "recreation_summary.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_folder,
    "sessionInfo.txt"
  )
)

cat(
  "\nRecreation finished.\n",
  "Saved MIF2 object: ", rds_file, "\n",
  "Saved diagnostic PDF: ", diagnostic_pdf, "\n",
  "Saved comparison summary: ",
  file.path(output_folder, "recreation_summary.csv"),
  "\n",
  sep = ""
)

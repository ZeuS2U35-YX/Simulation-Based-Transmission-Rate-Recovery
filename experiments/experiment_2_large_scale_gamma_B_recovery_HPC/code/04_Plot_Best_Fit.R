# ============================================================
# Run the final particle filter and save compact path summaries
# for the selected best fit.
#
# Inputs:
#   data/fixed_piecewise_B_dataset.csv
#   results/best_fit.csv
#   results/best_mif2_object.rds
#   results/combined_mif2_results.csv
#
# Outputs:
#   results/final_filtered_B_path.csv
#   results/final_filtered_infectious_path.csv
#   results/starting_value_summary.csv
#
# Figures are generated separately by code/05_regenerate_figures.R.
# ============================================================

library(pomp)

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. Read inputs
# ------------------------------------------------------------

data_file <- "data/fixed_piecewise_B_dataset.csv"
best_fit_file <- "results/best_fit.csv"
best_object_file <- "results/best_mif2_object.rds"
combined_results_file <- "results/combined_mif2_results.csv"

required_files <- c(
  data_file,
  best_fit_file,
  best_object_file,
  combined_results_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required files:",
      paste(missing_files, collapse = ", ")
    )
  )
}

sim1 <- read.csv(data_file)
best_fit <- read.csv(best_fit_file)
combined_results <- read.csv(combined_results_file)
mif_best <- readRDS(best_object_file)

sim1_observed <- sim1[, c("week", "reports"), drop = FALSE]

# ------------------------------------------------------------
# 2. Rebuild Gamma-noise model
# ------------------------------------------------------------

sir_rmeas <- Csnippet("reports = rnbinom_mu(k, rho * H);")

sir_dmeas <- Csnippet("lik = dnbinom_mu(reports, k, rho * H, give_log);")

sir_step_gamma <- Csnippet("
  double shape_B;
  double scale_B;
  double dN_SI;
  double dN_IR;
  double p_SI;
  double p_IR;

  if (sigma_beta > 0.0) {
    shape_B = 1.0 / (sigma_beta * sigma_beta * dt);
    scale_B = B * sigma_beta * sigma_beta * dt;
    B = rgamma(shape_B, scale_B);
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
  rprocess = euler(sir_step_gamma, delta.t = 1 / 30),
  rmeasure = sir_rmeas,
  dmeasure = sir_dmeas,
  accumvars = "H",
  statenames = c("S", "I", "R", "H", "B"),
  paramnames = c("B0", "sigma_beta", "mu_IR", "N", "rho", "k"),
  partrans = parameter_trans(log = c("B0", "sigma_beta"))
)

# ------------------------------------------------------------
# 3. Construct best parameter vector
# ------------------------------------------------------------

mif_coef_best <- coef(mif_best)

theta_best <- theta_gamma
theta_best[["B0"]] <- unname(mif_coef_best[["B0"]])
theta_best[["sigma_beta"]] <- unname(mif_coef_best[["sigma_beta"]])

# ------------------------------------------------------------
# 4. Final particle filter
# ------------------------------------------------------------

Np_final <- 50000
final_filter_seed <- 999

set.seed(final_filter_seed)

pf_best <- pfilter(
  pf_gamma_model,
  params = theta_best,
  Np = Np_final,
  filter.mean = TRUE
)

final_logLik <- as.numeric(logLik(pf_best))
fm <- filter_mean(pf_best)

# ------------------------------------------------------------
# 5. Extract paths and save CSV outputs
# ------------------------------------------------------------

theta_true <- c(Beta_high = 4, Beta_low = 2, t_switch = 5)

B_estimate <- data.frame(
  week = time(pf_best),
  B_true = ifelse(
    time(pf_best) < theta_true[["t_switch"]],
    theta_true[["Beta_high"]],
    theta_true[["Beta_low"]]
  ),
  B_filtered_mean = as.numeric(fm["B", ])
)

infectious_estimate <- data.frame(
  week = time(pf_best),
  gamma_infectious = as.numeric(fm["I", ])
)

true_infectious <- sim1[, c("week", "I"), drop = FALSE]
names(true_infectious)[2] <- "true_infectious"

infectious_estimate <- merge(
  infectious_estimate,
  true_infectious,
  by = "week",
  all.x = TRUE,
  sort = FALSE
)

dir.create("results", recursive = TRUE, showWarnings = FALSE)

write.csv(B_estimate, "results/final_filtered_B_path.csv", row.names = FALSE)
write.csv(infectious_estimate, "results/final_filtered_infectious_path.csv", row.names = FALSE)

combined_results$is_selected_best <- combined_results$task_id == best_fit$task_id[[1]]
combined_results <- combined_results[order(-combined_results$logLik), ]
write.csv(combined_results, "results/starting_value_summary.csv", row.names = FALSE)

cat(
  "Best task ID = ", best_fit$task_id[[1]],
  "\nEstimated B0 = ", theta_best[["B0"]],
  "\nEstimated sigma_beta = ", theta_best[["sigma_beta"]],
  "\nFinal pfilter logLik = ", final_logLik,
  "\n",
  sep = ""
)

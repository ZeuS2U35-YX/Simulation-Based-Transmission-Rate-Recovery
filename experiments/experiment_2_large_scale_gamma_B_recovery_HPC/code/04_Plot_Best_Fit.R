# ============================================================
# Run the final particle filter and plot the selected best fit
#
# Inputs:
#   data/fixed_piecewise_B_dataset.csv
#   results/best_fit.csv
#   results/best_mif2_object.rds
#
# Outputs:
#   results/final_filtered_B_path.csv
#   results/final_filtered_infectious_path.csv
#   figures/best_B_path.png
#   figures/best_infectious_path.png
#   figures/best_mif2_trace.png
# ============================================================

library(pomp)
library(ggplot2)

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. Read inputs
# ------------------------------------------------------------

data_file <- "data/fixed_piecewise_B_dataset.csv"
best_fit_file <- "results/best_fit.csv"
best_object_file <- "results/best_mif2_object.rds"

required_files <- c(
  data_file,
  best_fit_file,
  best_object_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required files:",
      paste(
        missing_files,
        collapse = ", "
      )
    )
  )
}

sim1 <- read.csv(
  data_file,
  stringsAsFactors = FALSE
)

best_fit <- read.csv(
  best_fit_file,
  stringsAsFactors = FALSE
)

mif_best <- readRDS(
  best_object_file
)

sim1_observed <- sim1[
  ,
  c(
    "week",
    "reports"
  ),
  drop = FALSE
]

# ------------------------------------------------------------
# 2. Rebuild Gamma-transition model
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
# 3. Construct best complete parameter vector
# ------------------------------------------------------------

mif_coef_best <- coef(
  mif_best
)

theta_best <- theta_gamma

theta_best[["B0"]] <-
  unname(
    mif_coef_best[["B0"]]
  )

theta_best[["sigma_beta"]] <-
  unname(
    mif_coef_best[["sigma_beta"]]
  )

# ------------------------------------------------------------
# 4. Final particle filter
# ------------------------------------------------------------

Np_final <- 50000
final_filter_seed <- 999

set.seed(
  final_filter_seed
)

pf_best <- pfilter(
  pf_gamma_model,
  params = theta_best,
  Np = Np_final,
  filter.mean = TRUE
)

final_logLik <- as.numeric(
  logLik(pf_best)
)

fm <- filter_mean(
  pf_best
)

# ------------------------------------------------------------
# 5. Extract paths
# ------------------------------------------------------------

theta_true <- c(
  Beta_high = 4,
  Beta_low = 2,
  t_switch = 5
)

B_estimate <- data.frame(
  week = time(pf_best),
  B_true = ifelse(
    time(pf_best) < theta_true[["t_switch"]],
    theta_true[["Beta_high"]],
    theta_true[["Beta_low"]]
  ),
  B_filtered_mean = as.numeric(
    fm["B", ]
  )
)

infectious_estimate <- data.frame(
  week = time(pf_best),
  gamma_infectious = as.numeric(
    fm["I", ]
  )
)

true_infectious <- sim1[
  ,
  c(
    "week",
    "I"
  ),
  drop = FALSE
]

names(true_infectious)[2] <- "true_infectious"

infectious_estimate <- merge(
  infectious_estimate,
  true_infectious,
  by = "week",
  all.x = TRUE,
  sort = FALSE
)

dir.create(
  "results",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  B_estimate,
  "results/final_filtered_B_path.csv",
  row.names = FALSE
)

write.csv(
  infectious_estimate,
  "results/final_filtered_infectious_path.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Plot B path
# ------------------------------------------------------------

dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)

B_path_plot <- ggplot(
  B_estimate,
  aes(x = week)
) +
  geom_step(
    aes(
      y = B_true,
      color = "True B(t)"
    ),
    linewidth = 1,
    direction = "hv"
  ) +
  geom_line(
    aes(
      y = B_filtered_mean,
      color = "Gamma-filtered B(t)"
    ),
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = theta_true[["t_switch"]],
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = c(
      "True B(t)" = "#00BFC4",
      "Gamma-filtered B(t)" = "#F8766D"
    ),
    breaks = c(
      "True B(t)",
      "Gamma-filtered B(t)"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(
      0,
      6,
      by = 1
    )
  ) +
  theme_bw(
    base_size = 14
  ) +
  labs(
    x = "Week",
    y = expression(B(t)),
    color = NULL
  ) +
  theme(
    legend.position = "top"
  )

ggsave(
  filename = "figures/best_B_path.png",
  plot = B_path_plot,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 7. Plot infectious path
# ------------------------------------------------------------

infectious_path_plot <- ggplot(
  infectious_estimate,
  aes(x = week)
) +
  geom_line(
    aes(
      y = true_infectious,
      color = "True I(t)"
    ),
    linewidth = 1
  ) +
  geom_line(
    aes(
      y = gamma_infectious,
      color = "Gamma-filtered I(t)"
    ),
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = theta_true[["t_switch"]],
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = c(
      "True I(t)" = "#00BFC4",
      "Gamma-filtered I(t)" = "#E69F00"
    ),
    breaks = c(
      "True I(t)",
      "Gamma-filtered I(t)"
    )
  ) +
  theme_bw(
    base_size = 14
  ) +
  labs(
    x = "Week",
    y = "Number infected and infectious",
    color = NULL
  ) +
  theme(
    legend.position = "top"
  )

ggsave(
  filename = "figures/best_infectious_path.png",
  plot = infectious_path_plot,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 8. Plot MIF2 traces
# ------------------------------------------------------------

png(
  filename = "figures/best_mif2_trace.png",
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  mif_best,
  pars = c(
    "B0",
    "sigma_beta"
  )
)

dev.off()

cat(
  "Best task ID = ",
  best_fit$task_id[[1]],
  "\nEstimated B0 = ",
  theta_best[["B0"]],
  "\nEstimated sigma_beta = ",
  theta_best[["sigma_beta"]],
  "\nFinal pfilter logLik = ",
  final_logLik,
  "\n",
  sep = ""
)

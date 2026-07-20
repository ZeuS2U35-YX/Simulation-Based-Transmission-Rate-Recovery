# ============================================================
# Estimating epidemic transmission rates using mif2
# True transmission rate is constant: B(t) = 4
# ============================================================


# ------------------------------------------------------------
# 1. Load packages and set options
# ------------------------------------------------------------

library(tidyverse)
library(pomp)


options(
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. Define the time grid for the simulated data
# ------------------------------------------------------------

n_weeks <- 10

template <- tibble(
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

  // Piecewise transmission-rate structure
  //
  // In this experiment, Beta_high and Beta_low
  // are both equal to 4, so B(t) is constant.
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


# ------------------------------------------------------------
# 2.2 Define the initial states for the data-generating model
# ------------------------------------------------------------

sir_rinit_piecewise <- Csnippet("
  S = N - 10;
  I = 10;
  R = 0;
  H = 0;
")


# ------------------------------------------------------------
# 2.3 Define the measurement model
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
# 2.4 Build the data-generating pomp object
# ------------------------------------------------------------

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
# 3. Set the true parameters
# ------------------------------------------------------------

theta <- c(
  Beta_high = 4,
  Beta_low = 4,
  t_switch = 5,
  mu_IR = 3,
  N = 10000,
  rho = 0.5,
  k = 10
)


# ------------------------------------------------------------
# 4. Simulate one epidemic data set
# ------------------------------------------------------------

set.seed(20260527)

sim1 <- simulate(
  simple_SIR,
  params = theta,
  nsim = 1,
  format = "data.frame",
  include.data = FALSE
)


# ------------------------------------------------------------
# 4.1 Keep only the observed reports
# ------------------------------------------------------------

sim1_observed <- sim1 |>
  select(
    week,
    reports
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


# ------------------------------------------------------------
# 5.1 Define initial states for the Gamma-transition model
# ------------------------------------------------------------

sir_rinit_gamma <- Csnippet("
  S = N - 10;
  I = 10;
  R = 0;
  H = 0;

  // Initial hidden transmission rate
  B = B0;
")


# ------------------------------------------------------------
# 5.2 Define the baseline filtering parameter vector
# ------------------------------------------------------------

theta_gamma <- c(
  B0 = 4,
  sigma_beta = 0.2,
  mu_IR = 3,
  N = 10000,
  rho = 0.5,
  k = 10
)


# ------------------------------------------------------------
# 6. Build the Gamma-transition pomp object
# ------------------------------------------------------------

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
# 6.1 Sanity check:
# run one ordinary particle filter before fitting
# ------------------------------------------------------------

set.seed(20260628)

pf_sanity_check <- pfilter(
  pf_gamma_model,
  params = theta_gamma,
  Np = 5000
)


# ------------------------------------------------------------
# 7. Set up the mif2 search
# ------------------------------------------------------------

# These are different starting values for mif2.
# They are not a grid search.
# Each mif2 run searches continuously from its own start.

start_values <- expand.grid(
  B0 = c(2, 4, 6),
  sigma_beta = c(0.10, 0.30, 0.45)
)


# Number of iterated-filtering iterations

Nmif <- 100


# Number of particles inside mif2

Np_mif <- 1000


# Number of particles used to evaluate each fitted parameter vector

Np_eval <- 50000


# Number of repeated pfilter evaluations for each mif2 result

n_pf_evals <- 5


# Number of particles for the final filtered paths

Np_final <- 50000


# Random-walk perturbation sizes used only during mif2.
#
# B0 is an initial-value parameter:
# it determines B only at time 0.
#
# sigma_beta is a regular time-constant parameter.
#
# Because both parameters are log-transformed,
# these perturbations occur on the estimation scale.

mif_rw_sd <- rw_sd(
  B0 = ivp(0.20),
  sigma_beta = 0.05
)


# ------------------------------------------------------------
# 8. Run mif2 from multiple starting values
# ------------------------------------------------------------

mif_fits <- vector(
  mode = "list",
  length = nrow(start_values)
)


mif_result_list <- vector(
  mode = "list",
  length = nrow(start_values)
)


for (i in seq_len(nrow(start_values))) {

  # Start with the full baseline parameter vector

  theta_start <- theta_gamma


  # Replace the two starting values for this mif2 run

  theta_start[["B0"]] <-
    start_values$B0[i]

  theta_start[["sigma_beta"]] <-
    start_values$sigma_beta[i]


  cat(
    "\nRunning mif2 from start ",
    i,
    " of ",
    nrow(start_values),
    ": B0 = ",
    theta_start[["B0"]],
    ", sigma_beta = ",
    theta_start[["sigma_beta"]],
    "\n",
    sep = ""
  )


  # All nine starting values use the same mif2 seed.
  # The seed is reset before every mif2 run.

  set.seed(
    20260628
  )


  # Run mif2

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
        "\nmif2 failed for start ",
        i,
        "."
      )

      message(
        "Exact mif2 error: ",
        conditionMessage(e)
      )

      NULL
    }
  )


  # If mif2 failed, record missing results

  if (is.null(mif_now)) {

    mif_result_list[[i]] <- tibble(
      run = i,

      start_B0 =
        theta_start[["B0"]],

      start_sigma_beta =
        theta_start[["sigma_beta"]],

      B0_hat = NA_real_,

      sigma_beta_hat = NA_real_,

      logLik = NA_real_,

      logLik_se = NA_real_
    )

    next
  }


  # Store successful mif2 result

  mif_fits[[i]] <- mif_now


  # Extract the fitted parameter values

  mif_coef_now <- coef(
    mif_now
  )


  # Build a complete parameter vector for pfilter evaluation.
  #
  # Keep mu_IR, N, rho, and k fixed.
  # Replace only B0 and sigma_beta.

  theta_eval_now <- theta_gamma

  theta_eval_now[["B0"]] <-
    unname(
      mif_coef_now[["B0"]]
    )

  theta_eval_now[["sigma_beta"]] <-
    unname(
      mif_coef_now[["sigma_beta"]]
    )


  # ----------------------------------------------------------
  # Evaluate this fitted parameter vector using repeated
  # ordinary particle filters
  # ----------------------------------------------------------

  eval_logLik <- rep(
    NA_real_,
    n_pf_evals
  )


  for (j in seq_len(n_pf_evals)) {

    # Every mif2 result uses the same evaluation seed
    # for the same repetition j.
    #
    # The five seeds are:
    # 20260801, 20260802, ..., 20260805.

    set.seed(
      20260800 + j
    )


    pf_eval <- tryCatch(
      pfilter(
        pf_gamma_model,

        params = theta_eval_now,

        Np = Np_eval
      ),

      error = function(e) {

        message(
          "\nEvaluation pfilter failed for mif2 run ",
          i,
          ", repetition ",
          j,
          "."
        )

        message(
          "Exact pfilter error: ",
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


  # Keep only successful pfilter evaluations

  eval_logLik <- eval_logLik[
    is.finite(eval_logLik)
  ]


  # Combine repeated likelihood estimates

  if (length(eval_logLik) > 0) {

    likelihood_summary <- pomp::logmeanexp(
      eval_logLik,
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


  # Store the result from this starting point

  mif_result_list[[i]] <- tibble(
    run = i,

    start_B0 =
      theta_start[["B0"]],

    start_sigma_beta =
      theta_start[["sigma_beta"]],

    B0_hat =
      theta_eval_now[["B0"]],

    sigma_beta_hat =
      theta_eval_now[["sigma_beta"]],

    logLik =
      estimated_logLik,

    logLik_se =
      estimated_logLik_se
  )
}


# ------------------------------------------------------------
# 8.1 Combine all mif2 results
# ------------------------------------------------------------

mif_results <- bind_rows(
  mif_result_list
)


# ------------------------------------------------------------
# 9. Select the best mif2 result
# ------------------------------------------------------------

best_fit <- mif_results |>
  filter(
    is.finite(logLik)
  ) |>
  slice_max(
    order_by = logLik,
    n = 1,
    with_ties = FALSE
  )


if (nrow(best_fit) == 0) {

  stop(
    "All mif2 runs or all pfilter evaluations failed."
  )
}


print(
  best_fit
)


# ------------------------------------------------------------
# 10. Build the best parameter vector
# ------------------------------------------------------------

best_run_id <- best_fit$run[[1]]

mif_best <- mif_fits[[best_run_id]]


mif_coef_best <- coef(
  mif_best
)


# Keep fixed parameters from theta_gamma
# and use the selected estimates for B0 and sigma_beta

theta_best <- theta_gamma

theta_best[["B0"]] <-
  unname(
    mif_coef_best[["B0"]]
  )

theta_best[["sigma_beta"]] <-
  unname(
    mif_coef_best[["sigma_beta"]]
  )


print(
  theta_best
)


cat(
  "\nBest result:",
  " run =", best_run_id,
  "| B0 =", round(theta_best[["B0"]], 4),
  "| sigma_beta =", round(theta_best[["sigma_beta"]], 4),
  "| logLik =", round(best_fit$logLik[[1]], 4),
  "| MCSE =", round(best_fit$logLik_se[[1]], 4),
  "\n"
)



# ------------------------------------------------------------
# 10.1 Save the best mif2 object
# ------------------------------------------------------------

# Output folders
mif2_output_dir <- "data/experiment_01"
mif2_figure_dir <- "figures/experiment_01"


# Create the folders if they do not already exist
dir.create(
  mif2_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  mif2_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 10.2 Plot the best mif2 object
# ------------------------------------------------------------

# ------------------------------------------------------------
# 10.2 Plot the best MIF2 object
# ------------------------------------------------------------

plot(
  mif_best,
  pars = c(
    "B0",
    "sigma_beta"
  )
)

# ------------------------------------------------------------
# 11. Run a final particle filter for the filtered paths
# ------------------------------------------------------------

set.seed(999)

pf_best <- pfilter(
  pf_gamma_model,

  params = theta_best,

  Np = Np_final,

  filter.mean = TRUE
)


# ------------------------------------------------------------
# 12. Extract the filtered transmission-rate path
# ------------------------------------------------------------

fm <- filter_mean(
  pf_best
)


B_estimate <- tibble(
  week = time(
    pf_best
  ),

  B_filtered_mean = as.numeric(
    fm["B", ]
  )
)


# Add the true B path.
#
# Beta_high and Beta_low are both 4,
# so the true transmission rate remains constant.

B_estimate <- B_estimate |>
  mutate(
    B_true = if_else(
      week < theta[["t_switch"]],
      theta[["Beta_high"]],
      theta[["Beta_low"]]
    )
  )


# ------------------------------------------------------------
# 12.1 Extract the true and filtered infectious paths
# ------------------------------------------------------------

infectious_estimate <- tibble(
  week = time(
    pf_best
  ),

  gamma_infectious = as.numeric(
    fm["I", ]
  )
) |>
  left_join(
    sim1 |>
      select(
        week,
        true_infectious = I
      ),
    by = "week"
  )


# ------------------------------------------------------------
# 13. Plot the true and filtered transmission-rate paths
# ------------------------------------------------------------

B_path_plot <- ggplot(
  B_estimate,
  aes(
    x = week
  )
) +

  # True B(t) = 4
  geom_step(
    aes(
      y = B_true,
      color = "True B(t)"
    ),
    linewidth = 1,
    direction = "hv"
  ) +

  # Gamma-filtered B(t)
  geom_line(
    aes(
      y = B_filtered_mean,
      color = "Gamma-filtered B(t)"
    ),
    linewidth = 0.7
  ) +

  # Reference time t_switch = 5
  geom_vline(
    xintercept = theta[["t_switch"]],
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
      from = 0,
      to = 6,
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


print(
  B_path_plot
)


# ------------------------------------------------------------
# 13.1 Plot the true and Gamma-filtered infectious paths
# ------------------------------------------------------------

infectious_path_plot <- ggplot(
  infectious_estimate,
  aes(
    x = week
  )
) +

  # True I(t)
  geom_line(
    aes(
      y = true_infectious,
      color = "True I(t)"
    ),
    linewidth = 1
  ) +

  # Gamma-filtered I(t)
  geom_line(
    aes(
      y = gamma_infectious,
      color = "Gamma-filtered I(t)"
    ),
    linewidth = 0.6
  ) +

  # Reference time t_switch = 5
  geom_vline(
    xintercept = theta[["t_switch"]],
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


print(
  infectious_path_plot
)
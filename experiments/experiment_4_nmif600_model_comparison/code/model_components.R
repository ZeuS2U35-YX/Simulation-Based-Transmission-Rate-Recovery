# ============================================================
# Shared POMP model components for Experiment 4
# ============================================================

make_observation_template <- function(config) {
  data.frame(
    week = seq(
      from = config$observation_interval,
      to = config$n_weeks,
      by = config$observation_interval
    ),
    reports = 0
  )
}

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

make_data_generating_model <- function(config) {
  template <- make_observation_template(config)

  sir_step_piecewise <- Csnippet("
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

  pomp(
    data = template,
    times = "week",
    t0 = 0,
    rinit = sir_rinit_piecewise,
    rprocess = euler(
      sir_step_piecewise,
      delta.t = config$process_delta_t
    ),
    rmeasure = sir_rmeas,
    dmeasure = sir_dmeas,
    accumvars = "H",
    statenames = c("S", "I", "R", "H"),
    paramnames = c(
      "Beta_high", "Beta_low", "t_switch",
      "mu_IR", "N", "rho", "k"
    )
  )
}

make_gamma_model <- function(observed_data, config) {
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

  pomp(
    data = observed_data,
    times = "week",
    t0 = 0,
    rinit = sir_rinit_gamma,
    rprocess = euler(
      sir_step_gamma,
      delta.t = config$process_delta_t
    ),
    rmeasure = sir_rmeas,
    dmeasure = sir_dmeas,
    accumvars = "H",
    statenames = c("S", "I", "R", "H", "B"),
    paramnames = c("B0", "sigma_beta", "mu_IR", "N", "rho", "k"),
    partrans = parameter_trans(
      log = c("B0", "sigma_beta")
    )
  )
}

make_constant_model <- function(observed_data, config) {
  sir_step_constant <- Csnippet("
    double dN_SI;
    double dN_IR;
    double p_SI;
    double p_IR;

    p_SI = 1.0 - exp(-Beta * I / N * dt);
    p_IR = 1.0 - exp(-mu_IR * dt);

    dN_SI = rbinom(S, p_SI);
    dN_IR = rbinom(I, p_IR);

    S = S - dN_SI;
    I = I + dN_SI - dN_IR;
    R = R + dN_IR;
    H = H + dN_SI;
  ")

  sir_rinit_constant <- Csnippet("
    S = N - 10;
    I = 10;
    R = 0;
    H = 0;
  ")

  pomp(
    data = observed_data,
    times = "week",
    t0 = 0,
    rinit = sir_rinit_constant,
    rprocess = euler(
      sir_step_constant,
      delta.t = config$process_delta_t
    ),
    rmeasure = sir_rmeas,
    dmeasure = sir_dmeas,
    accumvars = "H",
    statenames = c("S", "I", "R", "H"),
    paramnames = c("Beta", "mu_IR", "N", "rho", "k"),
    partrans = parameter_trans(
      log = "Beta"
    )
  )
}

gamma_baseline_parameters <- function(config) {
  c(
    B0 = 4,
    sigma_beta = 0.2,
    mu_IR = config$true_parameters[["mu_IR"]],
    N = config$true_parameters[["N"]],
    rho = config$true_parameters[["rho"]],
    k = config$true_parameters[["k"]]
  )
}

constant_baseline_parameters <- function(config) {
  c(
    Beta = 4,
    mu_IR = config$true_parameters[["mu_IR"]],
    N = config$true_parameters[["N"]],
    rho = config$true_parameters[["rho"]],
    k = config$true_parameters[["k"]]
  )
}

true_B_at_times <- function(times, config) {
  ifelse(
    times < config$true_parameters[["t_switch"]],
    config$true_parameters[["Beta_high"]],
    config$true_parameters[["Beta_low"]]
  )
}

# Observation-time recovery estimates are endpoint states. In the implemented
# update-before-events Gamma process, B(t_n) is the value that drove events in
# the final Euler substep ending at t_n. The matching truth is therefore the
# left limit B_true(t_n-), not the right-continuous path value B_true(t_n).
# For the present step change, the endpoint target remains Beta_high at the
# switch itself and becomes Beta_low only after the switch.
true_B_driver_at_endpoints <- function(times, config) {
  ifelse(
    times <= config$true_parameters[["t_switch"]],
    config$true_parameters[["Beta_high"]],
    config$true_parameters[["Beta_low"]]
  )
}

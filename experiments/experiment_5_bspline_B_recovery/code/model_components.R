coefficient_names <- paste0("b", seq_len(experiment_config$spline$nbasis))

make_observation_template <- function(config = experiment_config) {
  data.frame(
    week = seq(
      from = config$observation_dt,
      to = config$n_weeks,
      by = config$observation_dt
    ),
    reports = 0
  )
}

true_B_at_times <- function(week, config = experiment_config) {
  ifelse(
    week <= unname(config$truth[["t_switch"]]),
    unname(config$truth[["Beta_high"]]),
    unname(config$truth[["Beta_low"]])
  )
}

make_truth_model <- function(config = experiment_config) {
  truth_step <- pomp::Csnippet("
    double Beta_now;
    double dN_SI;
    double dN_IR;
    double p_SI;
    double p_IR;

    if (t <= t_switch) {
      Beta_now = Beta_high;
    } else {
      Beta_now = Beta_low;
    }

    p_SI = 1.0 - exp(-Beta_now * I / N * dt);
    p_IR = 1.0 - exp(-mu_IR * dt);

    dN_SI = rbinom(S, p_SI);
    dN_IR = rbinom(I, p_IR);

    S -= dN_SI;
    I += dN_SI - dN_IR;
    R += dN_IR;
    H += dN_SI;
  ")

  truth_rinit <- pomp::Csnippet("
    S = N - initial_infectious;
    I = initial_infectious;
    R = 0;
    H = 0;
  ")

  measurement_model <- make_measurement_model()

  pomp::pomp(
    data = make_observation_template(config),
    times = "week",
    t0 = 0,
    rinit = truth_rinit,
    rprocess = pomp::euler(truth_step, delta.t = config$process_dt),
    rmeasure = measurement_model$rmeasure,
    dmeasure = measurement_model$dmeasure,
    accumvars = "H",
    statenames = c("S", "I", "R", "H"),
    paramnames = c(
      "Beta_high", "Beta_low", "t_switch", "initial_infectious",
      "mu_IR", "N", "rho", "k"
    )
  )
}

make_measurement_model <- function() {
  list(
    rmeasure = pomp::Csnippet("
      reports = rnbinom_mu(k, rho * H);
    "),
    dmeasure = pomp::Csnippet("
      lik = dnbinom_mu(reports, k, rho * H, give_log);
    ")
  )
}

make_truth_parameters <- function(config = experiment_config) {
  c(
    config$truth,
    initial_infectious = config$initial_infectious,
    config$fixed_parameters
  )
}

make_spline_basis <- function(week, config = experiment_config) {
  pomp::bspline_basis(
    x = week,
    nbasis = config$spline$nbasis,
    degree = config$spline$degree,
    names = "xi%d",
    rg = config$spline$range
  )
}

make_spline_covariates <- function(config = experiment_config) {
  spline_week <- seq(
    from = config$spline$range[[1]],
    to = config$spline$range[[2]],
    by = config$process_dt
  )

  if (tail(spline_week, 1) < config$spline$range[[2]]) {
    spline_week <- c(spline_week, config$spline$range[[2]])
  }

  basis <- make_spline_basis(spline_week, config)
  covariate_columns <- c(
    list(week = spline_week),
    as.list(as.data.frame(basis)),
    list(times = "week")
  )

  do.call(pomp::covariate_table, covariate_columns)
}

make_bspline_model <- function(observed_data, config = experiment_config) {
  required_columns <- c("week", "reports")

  if (!all(required_columns %in% names(observed_data))) {
    stop("observed_data must contain columns: week and reports.")
  }

  spline_step <- pomp::Csnippet("
    double log_B;
    double B_now;
    double dN_SI;
    double dN_IR;
    double p_SI;
    double p_IR;

    log_B =
      b1 * xi1 + b2 * xi2 + b3 * xi3 +
      b4 * xi4 + b5 * xi5 + b6 * xi6;
    B_now = exp(log_B);

    p_SI = 1.0 - exp(-B_now * I / N * dt);
    p_IR = 1.0 - exp(-mu_IR * dt);

    dN_SI = rbinom(S, p_SI);
    dN_IR = rbinom(I, p_IR);

    S -= dN_SI;
    I += dN_SI - dN_IR;
    R += dN_IR;
    H += dN_SI;
  ")

  spline_rinit <- pomp::Csnippet("
    S = N - initial_infectious;
    I = initial_infectious;
    R = 0;
    H = 0;
  ")

  measurement_model <- make_measurement_model()

  pomp::pomp(
    data = observed_data[, c("week", "reports")],
    times = "week",
    t0 = 0,
    rinit = spline_rinit,
    rprocess = pomp::euler(spline_step, delta.t = config$process_dt),
    rmeasure = measurement_model$rmeasure,
    dmeasure = measurement_model$dmeasure,
    accumvars = "H",
    covar = make_spline_covariates(config),
    statenames = c("S", "I", "R", "H"),
    paramnames = c(
      coefficient_names, "initial_infectious", "mu_IR", "N", "rho", "k"
    )
  )
}

make_parameter_vector <- function(
  spline_coefficients,
  config = experiment_config
) {
  if (length(spline_coefficients) != config$spline$nbasis) {
    stop("Exactly six spline coefficients are required.")
  }

  names(spline_coefficients) <- coefficient_names

  c(
    spline_coefficients,
    initial_infectious = config$initial_infectious,
    config$fixed_parameters
  )
}

make_start_values <- function(
  n_start,
  seed,
  config = experiment_config
) {
  stopifnot(n_start >= 1L)
  set.seed(seed)

  starts <- matrix(
    NA_real_,
    nrow = n_start,
    ncol = config$spline$nbasis,
    dimnames = list(NULL, coefficient_names)
  )

  starts[1, ] <- log(3)

  if (n_start >= 2L) {
    starts[2, ] <- log(4)
  }

  if (n_start >= 3L) {
    starts[3, ] <- seq(log(4), log(2), length.out = config$spline$nbasis)
  }

  if (n_start >= 4L) {
    for (i in 4:n_start) {
      baseline <- runif(
        1,
        min = log(config$spline$start_B_range[[1]]),
        max = log(config$spline$start_B_range[[2]])
      )
      starts[i, ] <- baseline + rnorm(
        config$spline$nbasis,
        mean = 0,
        sd = config$spline$start_deviation_sd
      )
    }
  }

  starts[] <- pmax(
    config$spline$coefficient_limits[[1]],
    pmin(config$spline$coefficient_limits[[2]], starts)
  )

  data.frame(
    start_id = seq_len(n_start),
    start_source = c(
      "constant_B3",
      if (n_start >= 2L) "constant_B4",
      if (n_start >= 3L) "declining_B4_to_B2",
      if (n_start >= 4L) rep("random", n_start - 3L)
    ),
    starts,
    check.names = FALSE
  )
}

make_mif_rw_sd <- function(config = experiment_config) {
  pomp::rw_sd(
    b1 = config$spline$rw_sd,
    b2 = config$spline$rw_sd,
    b3 = config$spline$rw_sd,
    b4 = config$spline$rw_sd,
    b5 = config$spline$rw_sd,
    b6 = config$spline$rw_sd
  )
}

evaluate_parameter_vector <- function(
  model,
  params,
  Np,
  n_evals,
  seed_base
) {
  log_likelihoods <- rep(NA_real_, n_evals)

  for (j in seq_len(n_evals)) {
    set.seed(seed_base + j)

    filter_result <- tryCatch(
      pomp::pfilter(model, params = params, Np = Np),
      error = function(e) {
        message(
          "pfilter evaluation failed for repetition ", j, ": ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(filter_result)) {
      log_likelihoods[[j]] <- as.numeric(logLik(filter_result))
    }
  }

  successful <- log_likelihoods[is.finite(log_likelihoods)]

  if (length(successful) == 0L) {
    return(list(
      logLik = NA_real_,
      logLik_se = NA_real_,
      n_successful = 0L,
      replicates = log_likelihoods
    ))
  }

  summary <- pomp::logmeanexp(successful, se = TRUE)

  list(
    logLik = as.numeric(summary[[1]]),
    logLik_se = as.numeric(summary[[2]]),
    n_successful = length(successful),
    replicates = log_likelihoods
  )
}

reconstruct_B <- function(
  spline_coefficients,
  week,
  config = experiment_config
) {
  basis <- make_spline_basis(week, config)
  coefficients <- spline_coefficients[coefficient_names]

  if (any(!is.finite(coefficients))) {
    stop("All six spline coefficients must be finite.")
  }

  as.numeric(exp(basis %*% coefficients))
}

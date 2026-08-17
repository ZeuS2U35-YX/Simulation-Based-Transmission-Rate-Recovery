# Multi-start initialization for IF2

# Nine starting-value combinations for the Gamma-noise model
gamma_start_values <- expand.grid(
  B0 = c(2, 4, 6),
  sigma_beta = c(0.10, 0.30, 0.45)
)

# Six starting values for the constant-B model
constant_start_values <- data.frame(
  Beta = c(1, 2, 3, 4, 5, 6)
)

gamma_rw_sd <- rw_sd(
  B0 = ivp(0.20),
  sigma_beta = 0.05
)

constant_rw_sd <- rw_sd(
  Beta = 0.05
)

# Each Gamma-noise combination initializes one IF2 run
for (s in seq_len(nrow(gamma_start_values))) {
  theta_start <- theta_gamma
  theta_start[["B0"]] <- gamma_start_values$B0[[s]]
  theta_start[["sigma_beta"]] <- gamma_start_values$sigma_beta[[s]]

  mif_now <- mif2(
    gamma_model,
    params = theta_start,
    Np = config$Np_mif,
    Nmif = config$Nmif,
    rw.sd = gamma_rw_sd,
    cooling.type = config$cooling_type,
    cooling.fraction.50 = config$cooling_fraction_50
  )
}

# Each constant-B value initializes one IF2 run
for (s in seq_len(nrow(constant_start_values))) {
  theta_start <- theta_constant
  theta_start[["Beta"]] <- constant_start_values$Beta[[s]]

  mif_now <- mif2(
    constant_model,
    params = theta_start,
    Np = config$Np_mif,
    Nmif = config$Nmif,
    rw.sd = constant_rw_sd,
    cooling.type = config$cooling_type,
    cooling.fraction.50 = config$cooling_fraction_50
  )
}

# For each completed IF2 run, evaluate the fitted parameters with five
# independent particle filters. logmeanexp averages likelihood estimates on the
# likelihood scale and returns the result on the log scale; it is not an
# arithmetic mean of log likelihoods and it is unrelated to a filtering mean.
evaluation_loglik <- replicate(
  config$n_pf_evals,
  as.numeric(logLik(pfilter(
    model,
    params = fitted_params,
    Np = config$Np_eval
  )))
)
aggregated_loglik <- as.numeric(
  pomp::logmeanexp(evaluation_loglik)
)

# The stored logLik column contains aggregated_loglik for each candidate run.
best <- valid[which.max(valid$logLik), , drop = FALSE]

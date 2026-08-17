# Parameters used to generate the simulated data

theta <- c(
  Beta_high = 4,
  Beta_low = 2,
  t_switch = 5,
  mu_IR = 3,
  N = 10000,
  rho = 0.5,
  k = 10
)


# Baseline parameters for the Gamma-noise POMP model

theta_gamma <- c(
  B0 = 4,
  sigma_beta = 0.2,
  mu_IR = 3,
  N = 10000,
  rho = 0.5,
  k = 10
)


# Baseline parameters for the constant-B filtering model

theta_constant <- c(
  Beta = 4,
  mu_IR = 3,
  N = 10000,
  rho = 0.5,
  k = 10
)

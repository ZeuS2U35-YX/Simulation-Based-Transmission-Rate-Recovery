# Filtering latent process: Gamma-transition transmission rate

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

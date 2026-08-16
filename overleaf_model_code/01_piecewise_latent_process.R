# Data-generating latent process: piecewise transmission rate

sir_step <- Csnippet("
  double Beta_now;
  double dN_SI;
  double dN_IR;
  double p_SI;
  double p_IR;

  // Piecewise transmission rate
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

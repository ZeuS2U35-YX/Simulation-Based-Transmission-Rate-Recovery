# Filtering latent process: constant transmission rate

sir_step_constant <- Csnippet("
  double dN_SI;
  double dN_IR;
  double p_SI;
  double p_IR;

  // Beta is a fixed parameter and is not updated by the process
  p_SI = 1.0 - exp(-Beta * I / N * dt);
  p_IR = 1.0 - exp(-mu_IR * dt);

  // New infections and recoveries
  dN_SI = rbinom(S, p_SI);
  dN_IR = rbinom(I, p_IR);

  // Update SIR states and incidence accumulator
  S = S - dN_SI;
  I = I + dN_SI - dN_IR;
  R = R + dN_IR;
  H = H + dN_SI;
")

# Initial states for the Gamma-noise POMP model

sir_rinit_gamma <- Csnippet("
  S = N - 10;
  I = 10;
  R = 0;
  H = 0;

  // Initial hidden transmission rate
  B = B0;
")

# Negative-binomial measurement model

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

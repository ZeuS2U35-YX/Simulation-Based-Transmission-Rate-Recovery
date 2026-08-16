# Construction of the Gamma-transition filtering POMP model

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

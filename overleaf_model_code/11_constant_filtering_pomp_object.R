# Construction of the constant-B filtering POMP model

constant_model <- pomp(
  data = observed_data,
  times = "week",
  t0 = 0,

  rinit = sir_rinit_constant,

  rprocess = euler(
    sir_step_constant,
    delta.t = 1 / 30
  ),

  rmeasure = sir_rmeas,
  dmeasure = sir_dmeas,

  accumvars = "H",

  statenames = c(
    "S",
    "I",
    "R",
    "H"
  ),

  paramnames = c(
    "Beta",
    "mu_IR",
    "N",
    "rho",
    "k"
  ),

  partrans = parameter_trans(
    log = "Beta"
  )
)

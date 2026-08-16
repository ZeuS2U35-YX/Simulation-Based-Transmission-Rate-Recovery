# Construction of the data-generating POMP model

simple_SIR <- pomp(
  data = template,
  times = "week",
  t0 = 0,

  rinit = sir_rinit_piecewise,

  rprocess = euler(
    sir_step,
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
    "Beta_high",
    "Beta_low",
    "t_switch",
    "mu_IR",
    "N",
    "rho",
    "k"
  )
)

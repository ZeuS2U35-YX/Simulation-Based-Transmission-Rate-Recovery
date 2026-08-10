# ============================================================
# Generate one fixed epidemic dataset
#
# True transmission rate:
#   B(t) = 4 for t < 5
#   B(t) = 2 for t >= 5
#
# Output:
#   data/fixed_piecewise_B_dataset.csv
# ============================================================

library(pomp)

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. Observation-time grid
# ------------------------------------------------------------

n_weeks <- 10

template <- data.frame(
  week = seq(
    from = 1 / 7,
    to = n_weeks,
    by = 1 / 7
  ),
  reports = 0
)

# ------------------------------------------------------------
# 2. Data-generating process
# ------------------------------------------------------------

sir_step <- Csnippet("

  double Beta_now;
  double dN_SI;
  double dN_IR;
  double p_SI;
  double p_IR;

  if (t < t_switch) {
    Beta_now = Beta_high;
  } else {
    Beta_now = Beta_low;
  }

  p_SI = 1.0 - exp(-Beta_now * I / N * dt);
  p_IR = 1.0 - exp(-mu_IR * dt);

  dN_SI = rbinom(S, p_SI);
  dN_IR = rbinom(I, p_IR);

  S = S - dN_SI;
  I = I + dN_SI - dN_IR;
  R = R + dN_IR;

  H = H + dN_SI;
")

sir_rinit_piecewise <- Csnippet("
  S = N - 10;
  I = 10;
  R = 0;
  H = 0;
")

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

# ------------------------------------------------------------
# 3. True parameters
# ------------------------------------------------------------

theta <- c(
  Beta_high = 4,
  Beta_low = 2,
  t_switch = 5,
  mu_IR = 3,
  N = 10000,
  rho = 0.5,
  k = 10
)

# ------------------------------------------------------------
# 4. Simulate exactly one fixed dataset
# ------------------------------------------------------------

simulation_seed <- 20260527

set.seed(simulation_seed)

sim1 <- simulate(
  simple_SIR,
  params = theta,
  nsim = 1,
  format = "data.frame",
  include.data = FALSE
)

# ------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------

dir.create(
  "data",
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- "data/fixed_piecewise_B_dataset.csv"

write.csv(
  sim1,
  output_file,
  row.names = FALSE
)

cat(
  "Saved fixed dataset to: ",
  output_file,
  "\nSimulation seed: ",
  simulation_seed,
  "\n",
  sep = ""
)

print(head(sim1))

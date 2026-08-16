# Experiment 2 design

## Primary purpose

Repeat the Gamma-noise partially observed Markov process (POMP) recovery
workflow using 50,000 particles on an HPC cluster. The term `large_scale` in
the folder name refers to the particle count, not to the number of simulated
data sets.

The experiment also examines starting-value sensitivity by fitting one fixed
simulated data set from nine combinations of initial values for `B0` and
`sigma_beta`.

## Scope

This is a one-data-set, multi-start, large-particle experiment. It examines
whether the tested starts reach a similar independently evaluated
high-likelihood region.

It is not a replicated recovery-accuracy experiment and does not estimate
bias or RMSE across independently simulated data sets.

This experiment is retained as supporting evidence for the high-particle
implementation and starting-value behavior. Experiment 4 is the canonical
computational analysis and supersedes this experiment for final quantitative
Gamma-noise model claims.

## Data-generating model

- One fixed simulated epidemic data set.
- Observation period: 10 weeks.
- Observation interval: `1/7` week.
- Euler process step: `1/30` week.
- Initial state: `S = 9990`, `I = 10`, `R = 0`, `H = 0`.
- True transmission rate:
  - $B(t)=4$ before week 5;
  - $B(t)=2$ from week 5 onward.
- Measurement model:
  $Y_n\mid H_n\sim\operatorname{NegBin}(\text{mean}=\rho H_n,\text{size}=k)$.
- Fixed values: `mu_IR = 3`, `N = 10000`, `rho = 0.5`, `k = 10`.
- Simulation seed: `20260527`.

## Fitted model

Under the fitted Gamma-noise model, the one-step transition distribution for
the positive latent state `B` is Gamma. MIF2 estimates only:

- `B0`;
- `sigma_beta`.

The remaining parameters and the initial epidemic state are fixed at their
data-generating values. Both estimated parameters are log-transformed through
`parameter_trans(log = c("B0", "sigma_beta"))`, so MIF2 perturbations occur
on the transformed estimation scale. `B0` is an initial-value parameter that
determines the value of the latent transmission-rate process at time zero,
whereas `sigma_beta` is a regular time-constant parameter.

## Multi-start search design

- One fixed data set shared by all tasks.
- Nine starting points:
  - `B0` in `{2, 4, 6}`;
  - `sigma_beta` in `{0.10, 0.30, 0.45}`.
- One MIF2 search per starting point.
- One Slurm array task per starting point.
- Shared MIF2 seed across starting points: `20260628`.

## Large-particle settings

- `Nmif = 100`.
- `Np_mif = 50000`.
- Five likelihood evaluations per fitted parameter vector.
- `Np_eval = 50000`.
- `Np_final = 50000`.
- Evaluation seeds: `20260801` to `20260805`.
- Final filtering seed: `999`.

## Selection and final filtering

For each fitted parameter vector, five particle-filter likelihood estimates
are combined using `pomp::logmeanexp`. The candidate with the largest combined
evaluated log-likelihood is selected. A final particle filter with
`filter.mean = TRUE` then produces filtered transmission-rate and infectious
paths.

The selected candidate has the largest combined evaluated log likelihood among
the nine tested fits. This rule does not establish a unique global optimum.

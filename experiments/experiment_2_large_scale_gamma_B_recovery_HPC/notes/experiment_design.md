# Experiment 2 design

## Purpose

Assess whether MIF2 runs initialized from nine different starting points
converge to a similar high-likelihood region when fitting one fixed
piecewise-transmission epidemic dataset.

## Data-generating transmission rate

- \(B(t)=4\) before week 5.
- \(B(t)=2\) from week 5 onward.

## Search design

- One fixed simulated dataset.
- Nine starting points.
- One MIF2 search per starting point.
- Nine SLURM array tasks.
- `Nmif = 100`.
- `Np_mif = 50000`.
- Five likelihood evaluations per fitted parameter vector.
- `Np_eval = 50000`.

## Selection

The fitted parameter vector with the highest replicated particle-filter
log-likelihood estimate is selected. A final particle filter is then run
to produce filtered transmission-rate and infectious paths.

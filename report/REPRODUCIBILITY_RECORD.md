# Reproducibility record

This directory contains the complete LaTeX report and a report-local snapshot
of the Experiment 4 R source. The executable repository workflow remains
canonical:

- `../experiments/experiment_4_nmif600_model_comparison/code/`
- `../experiments/experiment_4_nmif600_model_comparison/config/`
- `../experiments/experiment_4_nmif600_model_comparison/shared_data/`
- `../experiments/experiment_4_nmif600_model_comparison/results/`

## Estimator definition

For every starting-value endpoint, five independent particle filters estimate
the likelihood. `pomp::logmeanexp` combines those likelihood estimates on the
likelihood scale, and the endpoint with the largest finite result is selected.
This operation does not average latent trajectories.

At the selected Gamma-noise parameters, a final 50,000-particle filter runs
with ancestry storage and `filter_traj()` extracts one internally consistent
sampled latent `B(t)` history. Its 70 observation-time values are used directly
for RSS, RMSE, signed mean error, and absolute overall bias. Time zero is shown
only in selected trajectory figures. No filtering mean enters the Experiment 4
recovery metrics.

## Verified retained summaries

- Mean RMSE: 0.653386 for Gamma-noise and 1.198924 for constant-B.
- Lower Gamma-noise RMSE: 195/200 tasks (97.5%).
- Mean absolute task-level mean error: 0.228762 and 0.606471.
- Lower Gamma-noise absolute mean error: 177/200 tasks (88.5%).
- Mean independent log likelihood: -191.083 for Gamma-noise and -209.228 for
  constant-B; the paired difference favors Gamma-noise in 199/200 tasks.

The values above were checked against
`../experiments/experiment_4_nmif600_model_comparison/results/comparison/`.

## Reproduction boundary

Task-specific seeds and data checksums are retained, but the original HPC run
did not preserve a complete software lockfile. Exact Monte Carlo reproduction
may depend on the R, `pomp`, compiler, and RNG versions. Preserve the supplied
seeds and accepted-data checksums when rerunning the workflow.

# 2026 Summer Epidemic Project

This repository develops and evaluates partially observed Markov process models for recovering a time-varying epidemic transmission rate.

This repository state is a **canonical computational results milestone**, not a final release. The report is currently in progress, and no licensing decision has yet been made.

## Experiment map

- **Experiment 4** (`experiments/experiment_4_nmif600_model_comparison/`) is the canonical final analysis and the primary source of quantitative evidence for the report. It is a **Gamma-noise vs. constant-B model comparison** using `Nmif = 600` on 200 shared accepted simulated outbreaks.
- **Experiment 3** is an earlier `Nmif = 100` recovery-accuracy study. It is retained as supporting evidence but is superseded by Experiment 4 for final Gamma-noise model numerical claims.
- **Experiments 1 and 2** are developmental studies of the fitting workflow, starting-value sensitivity, and high-particle implementation.

Experiments 1-3 should therefore not be presented as the repository's final quantitative analysis.

The constant-B model in Experiment 4 is a deliberately restricted comparator: it cannot represent temporal variation in `B(t)`. Traces from selected diagnostic tasks 1, 50, 100, 150, and 200 provide empirical support for `Nmif = 600`, but do not establish convergence for all 200 fitted data sets.

## Interpretation boundary

The repeated-simulation performance results are conditional on accepted informative outbreaks satisfying `max(H) > 20`. They do not estimate unconditional performance over all simulated epidemic trajectories.

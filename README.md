# A Simulation-Based Approach to Recovering Time-Varying Epidemic Transmission Rates

This repository develops and evaluates partially observed Markov process (POMP) models for recovering a time-varying epidemic transmission rate, `B(t)`, from noisy case-report data. Four simulation experiments progress from workflow development and computational diagnostics to a paired model comparison. The canonical analysis is Experiment 4, which compares a Gamma-noise model with latent time-varying `B(t)` against a deliberately restricted constant-B model on the same 200 accepted simulated outbreaks.

This repository is a **computational-results milestone**, not a final archival release. The accompanying report remains in development. Code is distributed under the [MIT License](LICENSE).

## Scientific question and evidence chain

The central question is whether a Gamma-noise POMP model can recover a prescribed change in transmission rate from partially observed epidemic data, and whether it does so more accurately than a model that assumes transmission is constant. Model parameters are estimated by iterated filtering (MIF2).

The experiments form a progression rather than four co-equal analyses:

| Experiment | Role in the evidence chain | Question addressed | Use in the report |
| --- | --- | --- | --- |
| [Experiment 1](experiments/experiment_1_gamma_B_recovery/) | Workflow development | Can the initial single-data-set workflows recover constant and piecewise transmission patterns from multiple starting values? | Developmental illustration |
| [Experiment 2](experiments/experiment_2_large_scale_gamma_B_recovery_HPC/) | High-particle computational diagnostic | On one fixed data set, do nine starting points reach similar independently evaluated likelihood regions with 50,000 particles? | Supporting implementation evidence |
| [Experiment 3](experiments/experiment_3_gamma_B_recovery_accuracy/) | Earlier repeated recovery study | How accurately do observation-time filtering means from the `Nmif = 100` Gamma-noise workflow recover prescribed `B(t)` values across 200 accepted outbreaks? | Supporting study; superseded by Experiment 4 for final numerical claims |
| [Experiment 4](experiments/experiment_4_nmif600_model_comparison/) | Canonical paired comparison | How does the `Nmif = 600` Gamma-noise model compare with the constant-B model on the same 200 accepted outbreaks? | Primary quantitative evidence |

Experiment 4 is therefore the source for final numerical comparisons. Experiments 1-3 document workflow development, particle-count and starting-value diagnostics, and supporting recovery behavior.

## Interpretation boundaries

- All findings arise from simulated stochastic SIR epidemics under the stated process and measurement models. They are not evidence of performance on real surveillance data.
- Experiments 3 and 4 retain only outbreaks satisfying `max(H) > 20`. Their recovery results are conditional on this informative-outbreak acceptance rule, not unconditional performance over all attempted simulations.
- Gamma-noise residual sum of squares (RSS), root mean squared error (RMSE), and signed-error summaries use per-task observation-time filtering means. Experiment 4 evaluates the constant-B model by repeating its fitted static estimate across the same 70 observation times.
- Selected-task curves serve different purposes. Experiment 3 Figure 01 and Experiment 4 Figure 08 show prespecified ancestry-preserving, finite-particle smoothing-trajectory approximations. Experiment 4 Figure 01 instead shows the task-1 Gamma filtering mean. None of these single-task displays is an across-task performance summary.
- Independent particle-filter log likelihoods are descriptive fitting diagnostics. They are not complexity-adjusted model-selection criteria.
- Experiment 4 convergence traces for tasks 1, 50, 100, 150, and 200 support the chosen computational settings for those tasks but do not establish convergence for all 200 fitted data sets.

## Repository guide

```text
experiments/
├── experiment_1_gamma_B_recovery/
├── experiment_2_large_scale_gamma_B_recovery_HPC/
├── experiment_3_gamma_B_recovery_accuracy/
└── experiment_4_nmif600_model_comparison/
README.md
SOFTWARE.md
LICENSE
```

Each experiment README documents its scientific purpose, model settings, execution order, retained outputs, figure provenance, limitations, and reproduction commands. Raw task-level HPC products and Slurm logs are excluded from Git when compact validated summaries provide the required evidence.

For the primary analysis, start with the [Experiment 4 README](experiments/experiment_4_nmif600_model_comparison/README.md). Its canonical evidence is stored in:

- [`results/combined/gamma/`](experiments/experiment_4_nmif600_model_comparison/results/combined/gamma/);
- [`results/combined/constant/`](experiments/experiment_4_nmif600_model_comparison/results/combined/constant/);
- [`results/comparison/`](experiments/experiment_4_nmif600_model_comparison/results/comparison/);
- [`results/selected_trajectory/`](experiments/experiment_4_nmif600_model_comparison/results/selected_trajectory/);
- [`figures/comparison/`](experiments/experiment_4_nmif600_model_comparison/figures/comparison/);
- [`figures/convergence/`](experiments/experiment_4_nmif600_model_comparison/figures/convergence/).

The root-level [trajectory audit](EXP3_EXP4_B_TRAJECTORY_AUDIT.md) and [repair record](EXP3_EXP4_B_TRAJECTORY_REPAIR.md) are retained as historical evidence of the distinction among filtering means, across-task averages, and ancestry-preserving trajectories. They are not substitutes for the current experiment READMEs.

## Lightweight reproduction

Run each command from the corresponding experiment directory. These commands rebuild figures and compact summaries from retained results unless explicitly stated otherwise; they do not rerun the full MIF2 studies.

```bash
# Experiment 1
Rscript code/04_regenerate_figures.R

# Experiment 2
Rscript code/05_regenerate_figures.R

# Experiment 3
Rscript code/04_analyze_results.R
```

For Experiment 4:

```bash
# Selected-task B(t) and infectious-state figures.
# These commands rerun only final 50,000-particle filters at saved parameters.
Rscript code/07_generate_task1_comparison_figures.R
Rscript code/08_generate_task117_comparison_figures.R

# Paired 200-task summaries and selected-task artifact validation.
Rscript code/05_compare_models.R \
  results/combined/gamma \
  results/combined/constant \
  results/comparison \
  figures/comparison

# Convergence diagnostics for the prespecified diagnostic tasks.
Rscript code/06_make_convergence_diagnostics.R \
  results/combined/gamma \
  results/combined/constant \
  figures/convergence
```

The complete workflows are computationally expensive and require Slurm. Consult the experiment-specific README before submitting an array job. Recorded and currently tested software environments are summarized in [SOFTWARE.md](SOFTWARE.md).

## Release status

The current milestone is suitable for supervisor review, report development, and independent inspection of the retained computational evidence. Exact historical package environments are only partly documented, so bit-for-bit reproduction of every original HPC particle-filter realization is not claimed.

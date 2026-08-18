# 2026 Summer Epidemic Project

This repository develops and evaluates partially observed Markov process models for recovering a time-varying epidemic transmission rate from noisy case-report data. The canonical final analysis is a Gamma-noise versus constant-B comparison, where the constant-`B` model is deliberately restricted.

This repository state is a **canonical computational results milestone**, not the final repository release. The report is currently in progress. The repository is distributed under the [MIT License](LICENSE).

## Scientific hierarchy

The experiments document a progression rather than four co-equal final analyses:

| Experiment | Scientific role | Main question | Status in the report |
| --- | --- | --- | --- |
| [Experiment 1](experiments/experiment_1_gamma_B_recovery/) | Developmental and exploratory | Do the initial one-dataset fitting scenarios recover constant and piecewise transmission paths from several starts? | Workflow development only |
| [Experiment 2](experiments/experiment_2_large_scale_gamma_B_recovery_HPC/) | Supporting high-particle diagnostic | On one fixed data set, do nine starts reach a similar high-likelihood region with 50,000 particles? | Supporting implementation evidence |
| [Experiment 3](experiments/experiment_3_gamma_B_recovery_accuracy/) | Earlier repeated recovery study | How accurately do `Nmif = 100` observation-time filtering means recover prescribed values across 200 accepted outbreaks? | Supporting study; superseded by Experiment 4 for final Gamma-noise numerical claims |
| [Experiment 4](experiments/experiment_4_nmif600_model_comparison/) | Canonical final computational analysis | How does the `Nmif = 600` Gamma-noise model compare with the constant-`B` model on the same 200 accepted outbreaks? | Primary quantitative evidence for the report |

Experiment 4 is therefore the canonical source for final numerical results. Experiments 1–3 remain useful for development history, diagnostics, and supporting interpretation, but should not be presented as equally important final analyses.

## Interpretation boundaries

- The constant-`B` model in Experiment 4 cannot represent temporal variation in `B(t)`; it is intentionally used as a restricted comparator.
- Recovery and comparison conclusions in Experiments 3 and 4 are conditional on accepted informative outbreaks satisfying `max(H) > 20`. They do not estimate unconditional performance over every attempted epidemic trajectory.
- Experiment 3 retains its earlier filtering-mean summaries. In Experiment 4, each accepted replicate contributes one ancestry-preserving sampled latent Gamma trajectory from the final plug-in particle filter; the constant-model estimate is repeated over the same 70 observation times. These paired replicate-level objects are the inputs to the RSS, RMSE, signed mean error, and absolute overall bias tables.
- Experiment 4 Figures 01 and 08 show the exact sampled trajectories used for the task-1 and task-117 recovery metrics. Each is a single finite-particle approximation conditional on one complete case series, not a filtering mean, uncertainty interval, or across-task average.
- Traces from selected Experiment 4 diagnostic tasks 1, 50, 100, 150, and 200 provide empirical support for `Nmif = 600`, but do not establish convergence for all 200 fitted data sets. These tasks are selected diagnostics, not a random sample.

## Repository layout

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

Each experiment README documents its scientific purpose, settings, execution order, compact retained results, figure provenance, limitations, and full or lightweight reproduction commands. Raw task-level HPC outputs and Slurm logs are intentionally kept out of Git when compact combined evidence is available.

The canonical Experiment 4 evidence is under:

- [`results/combined/gamma/`](experiments/experiment_4_nmif600_model_comparison/results/combined/gamma/)
- [`results/combined/constant/`](experiments/experiment_4_nmif600_model_comparison/results/combined/constant/)
- [`results/comparison/`](experiments/experiment_4_nmif600_model_comparison/results/comparison/)
- [`results/selected_trajectory/`](experiments/experiment_4_nmif600_model_comparison/results/selected_trajectory/)
- [`figures/comparison/`](experiments/experiment_4_nmif600_model_comparison/figures/comparison/)
- [`figures/convergence/`](experiments/experiment_4_nmif600_model_comparison/figures/convergence/)

The current sampled-trajectory results have mean RMSE 0.653 for the Gamma-noise model and 1.199 for the constant-`B` comparator. Gamma has lower paired RMSE in 195 of 200 accepted replicates (97.5%) and lower absolute replicate-level mean error in 177 of 200 (88.5%). The root-level trajectory audit and repair record document earlier stages of the analysis; their filtering-mean values are historical rather than active Experiment 4 results.

## Lightweight reproduction entry points

Run each command from the corresponding experiment directory. These commands rebuild figures and compact summaries from retained results; they do not rerun MIF2 or the full HPC studies.

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
# Reconstruct the 200 canonical sampled Gamma trajectories from saved fits.
# This reruns only the final 50,000-particle filters.
Rscript code/09_regenerate_sampled_B_trajectories.R \
  results/combined/gamma/combined_B_paths.csv \
  results/combined/gamma/sampled_B_trajectory_provenance.csv \
  4 \
  1:200

# Rebuild selected-task figures from the canonical sampled paths.
# These commands do not run another particle filter.
Rscript code/07_generate_task1_comparison_figures.R
Rscript code/08_generate_task117_comparison_figures.R

Rscript code/05_compare_models.R \
  results/combined/gamma \
  results/combined/constant \
  results/comparison \
  figures/comparison

Rscript code/06_make_convergence_diagnostics.R \
  results/combined/gamma \
  results/combined/constant \
  figures/convergence
```

The complete HPC workflows are computationally expensive and require Slurm. See each experiment README before submitting them. Software and environment evidence is summarized in [SOFTWARE.md](SOFTWARE.md).

## Release status

The computational milestone is suitable for supervisor review and continued report drafting. It is not a final archival release: the report remains in progress and exact historical package environments are only partly documented. The repository is licensed under the [MIT License](LICENSE).

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
- Experiment 3 retains its earlier filtering-mean summaries. Experiment 4 likewise uses the particle filtering mean of the latent Gamma state at each of the 70 observation times as its primary recovery estimate. The constant-model estimate is repeated over the same times. These paired point-estimator trajectories are the inputs to the RSS, RMSE, signed mean error, and absolute overall bias tables.
- For the prespecified Task 1 illustration, the report shows particle filtering-mean trajectories for both `B(t)` and `I(t)`, followed by a separate fixed-seed joint forward simulation at the fitted parameters. Neither panel changes the primary metrics: only the observation-time `B(t)` filtering means enter the Gamma-noise recovery calculations. Historical ancestry-path artifacts for tasks 1 and 117 remain in the repository for provenance but are not active report figures.
- Endpoint recovery is aligned to the rate that drove the final Euler substep ending at each observation time. The week-5 endpoint therefore has target `B = 4`, although the right-continuous prescribed path has already changed to `B(5) = 2`.
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

The current endpoint-aligned filtering-mean results have mean RMSE 0.522 for the Gamma-noise model and 1.186 for the constant-`B` comparator. Gamma has lower paired RMSE in all 200 accepted replicates and lower absolute replicate-level mean error in 180 of 200 (90.0%). Historical unaligned or sampled-trajectory summaries are not the active Experiment 4 results.

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

The filtering-mean and Task 1 trajectory-figure commands have one additional
input requirement. `code/10_regenerate_filtering_mean_B_paths.R` and
`code/11_generate_task1_filtering_mean_forward_figures.R` read the exact case series from
`shared_data/task_###/observed_data.csv`. The full Slurm workflow creates these
files during its shared-data stage, and the packaged completed-experiment
archive includes them. They are not committed in this compact repository
milestone. To rerun the final particle filters, first restore `shared_data/`
from that archive or run the full shared-data generation stage. Without those
files, use the retained canonical filtering means in
`results/combined/gamma/combined_B_filtering_means.csv` and skip particle-filter
reconstruction; the aggregate comparison commands operate on retained combined
results.

```bash
# Reconstruct the 200 primary Gamma filtering-mean trajectories from saved fits.
# Requires shared_data/task_###/observed_data.csv for tasks 1-200.
# This reruns only the final 50,000-particle filters, not MIF2.
Rscript code/10_regenerate_filtering_mean_B_paths.R \
  results/combined/gamma/combined_B_filtering_means.csv \
  results/combined/gamma/filtering_mean_provenance.csv \
  4 \
  1:200

# Build the active Task 1 manuscript figures: filtering means for B(t) and I(t),
# plus one fixed-seed joint forward simulation at the fitted parameters.
# Requires the exact Task 1 shared data.
Rscript code/11_generate_task1_filtering_mean_forward_figures.R shared_data

# Optional legacy provenance: reconstruct ancestry paths no longer used in the report.
Rscript code/09_regenerate_sampled_B_trajectories.R \
  results/combined/gamma/combined_B_paths.csv \
  results/combined/gamma/sampled_B_trajectory_provenance.csv \
  4 \
  1:200

# Optional: rebuild the legacy repository-only ancestry-path figures.
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

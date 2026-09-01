# Simulation-Based Transmission-Rate Recovery

[![Archived release v1.0.0 DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22127917.svg)](https://doi.org/10.5281/zenodo.22127917)

This repository contains simulation, fitting, recovery, validation, and
figure-generation workflows for recovering a time-varying epidemic
transmission rate, \(B(t)\), from partially observed epidemic data.

The study has **one shared data-generating stochastic SIR process** and
**three fitted transmission-rate models**:

1. a Gamma-noise model with a latent stochastic \(B(t)\), summarized for the
   primary analysis by observation-time particle filtering means;
2. a deterministic non-periodic cubic B-spline model for \(\log B(t)\), using
   six B-spline basis coefficients;
3. a constant-\(B\) model using one fitted transmission-rate value over the
   complete epidemic.

The six B-spline coefficients are components of one fitted B-spline model;
they are not six separate models. The stochastic SIR data-generating process
is also not counted as a fourth fitted model.

> **Release status.** The current public release is
> [v1.0.0](https://github.com/ZeuS2U35-YX/Simulation-Based-Transmission-Rate-Recovery/releases/tag/v1.0.0).
> A single-download full-replication package is being prepared for a future
> `v1.1.0` release, but `v1.1.0` has not been published. No DOI has been
> assigned to `v1.1.0`.

## Study design

One stochastic SIR data-generating process is used to produce the shared
simulation replicates. Each accepted replicate contains 70 observation times
over 10 weeks. Reported cases follow a negative-binomial measurement model.
Replicates are retained when \(\max(H)>20\), so the reported recovery
performance is conditional on this informative-outbreak acceptance rule.

The same 200 accepted simulated case series are used for all three fitted
models. Pairing is checked using task identifiers, simulation seeds, and
observed-data checksums.

### Week-5 convention and provenance

The transmission-rate step occurs at week 5, but two endpoint conventions are
present in the retained record and must be distinguished. Historical
step-change workflows and raw truth labels underlying Experiments 1–4 used
\(B(t)=4\) for \(t<5\) and \(B(t)=2\) for \(t\geq5\), so the exact week-5
point was labelled \(B(5)=2\). This statement does not apply to Experiment
1's separate constant-\(B=4\) scenario.

The current post-processing code and final Experiment 5 reporting normalize
the common comparison truth to \(B(5)=4\), using \(B(t)=4\) for \(t\leq5\)
and \(B(t)=2\) for \(t>5\), and recompute the reported recovery metrics under
that convention without changing the historical observed data, fitted
coefficients, or fitted model objects. Historical raw files are retained as
provenance; their truth labels or metrics must not be mixed with the normalized
final tables.

| Fitted model | Representation of \(B(t)\) | Primary path used for comparison |
| --- | --- | --- |
| Gamma-noise | Positive latent stochastic process | Observation-time particle filtering mean |
| B-spline | Deterministic non-periodic cubic B-spline for \(\log B(t)\) with six basis coefficients | Selected deterministic B-spline trajectory |
| Constant-\(B\) | One fitted scalar over the complete epidemic | Fitted scalar repeated at all observation times |

All three models are fitted within a partially observed Markov process
framework using multi-start iterated filtering. Starting values, particles,
Euler substeps, repeated likelihood evaluations, and observation times are
computational or within-replicate components. The independent analysis unit is
the accepted simulation replicate.

## Primary analysis

Experiment 4 and Experiment 5 have distinct roles:

- **Experiment 4 is the two-model computational foundation.** It generates the
  200 accepted shared data sets and fits the Gamma-noise and constant-\(B\)
  models.
- **Experiment 5 is the final three-model analysis.** It fits the deterministic
  six-basis cubic B-spline model to the same accepted Experiment 4 data and
  combines all three fitted models in one paired comparison.

The final design is therefore:

```text
one shared stochastic SIR data-generating process
  -> 200 accepted simulated epidemic replicates
  -> 200 Gamma-noise fits
  -> 200 deterministic six-basis cubic B-spline fits
  -> 200 constant-B fits
  -> exactly 200 paired three-model comparison rows
```

Under the normalized Experiment 5 reporting convention, the final comparison
gives the following mean replicate-level RMSE values:

| Model | Mean RMSE |
| --- | ---: |
| Gamma-noise | 0.5218 |
| Deterministic six-basis cubic B-spline | 0.7128 |
| Constant-\(B\) | 1.1858 |

Gamma-noise has lower RMSE than B-spline in 162 of 200 paired replicates.
B-spline has lower RMSE than constant-\(B\) in 188 of 200 paired replicates.
Gamma-noise has lower RMSE than constant-\(B\) in all 200 paired replicates.

These are latent transmission-path recovery results for the configured
single-step-change simulation. They are not predictive-performance
comparisons, do not establish a universal ranking of the three models, and
should not be generalized directly to real surveillance data.

## Experiment hierarchy

| Experiment | Role | Status |
| --- | --- | --- |
| [Experiment 1](experiments/experiment_1_gamma_B_recovery/) | Developmental Gamma-noise fitting scenarios and starting-value diagnostics | Supporting |
| [Experiment 2](experiments/experiment_2_large_scale_gamma_B_recovery_HPC/) | High-particle, single-data-set Gamma-noise diagnostic | Supporting |
| [Experiment 3](experiments/experiment_3_gamma_B_recovery_accuracy/) | Earlier repeated Gamma-noise recovery study | Supporting; superseded for final numerical claims |
| [Experiment 4](experiments/experiment_4_nmif600_model_comparison/) | Shared accepted data plus Gamma-noise and constant-\(B\) fits for 200 paired replicates | Primary computational foundation; two fitted models |
| [Experiment 5](experiments/experiment_5_bspline_B_recovery/) | B-spline extension and paired Gamma-noise/B-spline/constant-\(B\) comparison | Final primary analysis; three fitted models |

Experiments 1–3 document workflow development and supporting computational
evidence. Experiment 4 supplies the shared data and two-model foundation.
Experiment 5 supplies the B-spline extension and the final manuscript-facing
three-model results.

## Download options

### Current public release

The current public release is
[v1.0.0](https://github.com/ZeuS2U35-YX/Simulation-Based-Transmission-Rate-Recovery/releases/tag/v1.0.0).
Its exact archived record is available from Zenodo at
[https://doi.org/10.5281/zenodo.22127917](https://doi.org/10.5281/zenodo.22127917).

Version `v1.0.0` contains the tracked code, retained aggregate results, figure
source data, plotting scripts, software documentation, seeds, and pairing
information needed for the documented lightweight figure reproduction. It
does not contain a custom full-replication ZIP with every Git-ignored shared
input and raw task-level HPC output.

### Planned full-replication package

The planned future release asset is:

```text
Simulation-Based-Transmission-Rate-Recovery-v1.1.0-full-replication.zip
```

It is intended to provide one download containing one top-level folder while
preserving the internal experiment subdirectories and relative paths. The
package is also intended to include the shared accepted data, required raw HPC
outputs, retained combined results, figures, `MANIFEST.csv`, and
`SHA256SUMS`.

`MANIFEST.csv` uses the columns:

```text
path,size_bytes,sha256,role,provenance
```

This package has **not yet been published**. It will appear on the
[GitHub Releases page](https://github.com/ZeuS2U35-YX/Simulation-Based-Transmission-Rate-Recovery/releases)
only after validation and release approval.

GitHub's automatically generated “Source code” ZIP is a snapshot of
Git-tracked files. It should not be confused with the planned full-replication
package.

The package entry guide is [README_FIRST.md](README_FIRST.md).

## Choose a reproduction level

| Level | Purpose | Input | Computational scope |
| --- | --- | --- | --- |
| Quick reproduction | Regenerate the final three-model figures and figure source-data files | Git-tracked post-processed tables | No simulation, MIF2, particle filtering, or Slurm |
| Full reproduction | Regenerate accepted data and refit the Gamma-noise, B-spline, and constant-\(B\) models | Full-replication package plus an appropriate HPC environment | Slurm arrays, MIF2, repeated particle filtering, validation, and post-processing |

### Quick reproduction of the final figures

Requirements:

- R 4.5.2;
- Internet access for the first package restore;
- Git or a GitHub source-code download.

Run:

```bash
git clone https://github.com/ZeuS2U35-YX/Simulation-Based-Transmission-Rate-Recovery.git
cd Simulation-Based-Transmission-Rate-Recovery/experiments/experiment_5_bspline_B_recovery

Rscript -e 'renv::restore(prompt = FALSE)'
Rscript -e 'renv::status()'
Rscript code/06_plot_three_model_comparison.R
```

The command regenerates the two final comparison figures in PDF, SVG, PNG, and
600-dpi TIFF formats, together with their source-data CSV files under:

```text
figures/comparison_three_models/
```

For an additional data check, run:

```bash
Rscript -e 'x <- read.csv("figures/comparison_three_models/source_data/figure_2_mean_error_histograms_source_data.csv"); stopifnot(nrow(x) == 1800L, length(unique(x$task_id)) == 200L); cat("Reproduction checks passed.\n")'
```

See
[Experiment 5 figure-reproduction instructions](experiments/experiment_5_bspline_B_recovery/REPRODUCE_FIGURES.md)
for details.

### Full reproduction

The complete workflow is computationally expensive and is not a laptop-scale
figure-regeneration task. It requires a Slurm environment, an R compiler
toolchain compatible with `pomp`, and careful validation of all task-level
inputs and outputs.

The required order is:

1. run Experiment 4 to generate the accepted shared data and fit the
   Gamma-noise and constant-\(B\) models;
2. validate the paired Experiment 4 inputs;
3. run the Experiment 5 production pilot and B-spline task array;
4. combine and validate all 200 B-spline tasks;
5. construct the paired three-model tables and figures.

Follow the experiment-specific instructions rather than launching scripts from
the repository root:

- [Experiment 4 full workflow](experiments/experiment_4_nmif600_model_comparison/)
- [Experiment 5 full workflow](experiments/experiment_5_bspline_B_recovery/)
- [Software and computing environment](SOFTWARE.md)

The historical full HPC package environment was not preserved as a complete
lockfile or container image. Exact bit-for-bit reproduction of every historical
particle-filter realization is therefore not claimed. The Experiment 5
`renv.lock` is intentionally scoped to the final plotting workflow.

## Repository guide

```text
Simulation-Based-Transmission-Rate-Recovery/
├── README.md
├── README_FIRST.md
├── SOFTWARE.md
├── CITATION.cff
├── LICENSE
├── shared_code/
└── experiments/
    ├── experiment_1_gamma_B_recovery/
    ├── experiment_2_large_scale_gamma_B_recovery_HPC/
    ├── experiment_3_gamma_B_recovery_accuracy/
    ├── experiment_4_nmif600_model_comparison/
    └── experiment_5_bspline_B_recovery/
```

The final tracked three-model evidence is under Experiment 5:

- `results/combined/bspline/`;
- `results/comparison_three_models/`;
- `figures/comparison_three_models/`;
- `results/paired_input_manifest.csv`.

The figure directory includes vector and raster exports, source-data CSV files,
and visual quality-control notes.

## Interpretation boundaries

- All findings arise from simulated stochastic SIR epidemics under the stated
  process and measurement models.
- Recovery is evaluated conditionally on accepted outbreaks satisfying
  `max(H) > 20`.
- Historical step-change raw artifacts underlying Experiments 1–4 use
  `B(5) = 2`; current normalized post-processing and final Experiment 5
  reporting use `B(5) = 4`. Results from the two conventions must not be
  mixed.
- Gamma-noise recovery metrics use observation-time particle filtering means.
  Selected-task sampled trajectories are illustrative artifacts and are not
  the primary metric inputs.
- The B-spline path is deterministic conditional on its selected fitted
  coefficient vector.
- The constant-\(B\) estimate is repeated across observation times for the
  path-recovery comparison.
- Independent particle-filter likelihoods are descriptive fitting diagnostics,
  not complexity-adjusted model-selection criteria.
- No significance tests are used to turn the paired simulation results into a
  general model-ranking claim.

## Versioning, citation, and license

The current public release is `v1.0.0`. Cite its version-specific Zenodo record
when using that archived release:

- DOI: [10.5281/zenodo.22127917](https://doi.org/10.5281/zenodo.22127917)
- GitHub release:
  [v1.0.0](https://github.com/ZeuS2U35-YX/Simulation-Based-Transmission-Rate-Recovery/releases/tag/v1.0.0)

The `main` branch may contain documentation or packaging work performed after
the `v1.0.0` tag. A future `v1.1.0` release should be cited only after it is
actually published and assigned its own release metadata. No `v1.1.0` DOI is
claimed here.

Citation metadata are provided in [CITATION.cff](CITATION.cff). Code is
distributed under the [MIT License](LICENSE).

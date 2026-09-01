# README FIRST: Full Replication Package

## Release status

This file is the entry guide for the full-replication package being prepared
for:

```text
Simulation-Based-Transmission-Rate-Recovery-v1.1.0-full-replication.zip
```

As of the preparation of this guide, **v1.1.0 has not been published**. The
current public release remains
[v1.0.0](https://github.com/ZeuS2U35-YX/Simulation-Based-Transmission-Rate-Recovery/releases/tag/v1.0.0),
archived at
[https://doi.org/10.5281/zenodo.22127917](https://doi.org/10.5281/zenodo.22127917).

No DOI has been assigned to `v1.1.0`. If this file is being read from the
`main` branch or from a candidate archive, treat the package as pre-release
material until a `v1.1.0` GitHub release is publicly available.

## What this project compares

The study uses **one shared data-generating stochastic SIR process** to produce
200 accepted simulated epidemic replicates. It then fits **three
transmission-rate models** to every accepted replicate:

1. **Gamma-noise:** a latent stochastic \(B(t)\), summarized for the primary
   analysis by observation-time particle filtering means;
2. **B-spline:** a deterministic non-periodic cubic B-spline model for
   \(\log B(t)\), using six basis coefficients;
3. **Constant-\(B\):** one fitted transmission-rate value over the complete
   epidemic.

The six B-spline coefficients define one B-spline model. They are not six
models. The stochastic SIR data-generating process is not a fourth fitted
model.

### Week-5 convention and provenance

The retained files contain two explicit endpoint conventions at the week-5
change point. Historical step-change workflows and raw truth labels underlying
Experiments 1–4 used \(B(t)=4\) for \(t<5\) and \(B(t)=2\) for
\(t\geq5\), so the exact week-5 point was labelled \(B(5)=2\). This does not
describe Experiment 1's separate constant-\(B=4\) scenario.

The current post-processing code and final Experiment 5 reporting normalize
the shared comparison truth to \(B(5)=4\): \(B(t)=4\) for \(t\leq5\) and
\(B(t)=2\) for \(t>5\). They recompute the final metrics without changing
the historical observed data, fitted B-spline coefficients, or fitted model
objects. Historical raw products remain in the package as provenance. Do not
combine their truth labels or metrics with normalized Experiment 5 results.

## Experiment 4 and Experiment 5

The two primary experiment directories serve different purposes:

- `experiment_4_nmif600_model_comparison` is the **two-model computational
  foundation**. It generates the shared accepted data and fits the Gamma-noise
  and constant-\(B\) models.
- `experiment_5_bspline_B_recovery` is the **final three-model analysis**. It
  fits the B-spline model to the same accepted data and combines Gamma-noise,
  B-spline, and constant-\(B\) results.

Experiment 5 does not generate another set of 200 observations.

```text
200 accepted Experiment 4 data sets
  -> 200 Gamma-noise fits
  -> 200 B-spline fits
  -> 200 constant-B fits
  -> 200 paired three-model comparison rows
```

## Start here

Most readers should use the quick reproduction route. Use the full route only
if you intend to inspect or rerun the complete HPC fitting workflow.

| Route | What it does | Typical requirement |
| --- | --- | --- |
| Quick reproduction | Regenerates final three-model figures from retained post-processed tables | R 4.5.2; no Slurm |
| Full reproduction | Regenerates accepted simulations and refits all three models | Slurm HPC, `pomp`, compiler toolchain, substantial compute time |

## 1. Verify the package

A released archive should be accompanied by a published SHA-256 checksum.
Verify the downloaded ZIP before extracting it.

After extraction, the top-level directory should contain:

```text
README_FIRST.md
README.md
SOFTWARE.md
CITATION.cff
LICENSE
MANIFEST.csv
PACKAGE_METADATA.txt
SHA256SUMS
experiments/
scripts/
shared_code/
```

Verify the extracted files from the package root.

On macOS:

```bash
shasum -a 256 -c SHA256SUMS
```

On Linux:

```bash
sha256sum -c SHA256SUMS
```

`MANIFEST.csv` is the package inventory and uses exactly these columns:

```text
path,size_bytes,sha256,role,provenance
```

The fields identify the relative file path, file size in bytes, SHA-256 hash,
package role, and provenance classification for each regular payload file.
The manifest excludes itself and `SHA256SUMS` to avoid self-referential hashes;
`SHA256SUMS` covers `MANIFEST.csv` and every other regular file except itself.

If checksum verification fails, stop. Do not use the failed archive for
reproduction or citation.

## 2. Quick reproduction: regenerate the final figures

This route uses retained, post-processed three-model tables. It does not
simulate epidemics, rerun MIF2, rerun particle filters, or submit Slurm jobs.

### Requirements

- R 4.5.2;
- Internet access for the first `renv` restore;
- a standard R installation able to install the locked plotting packages.

From the extracted package root, run:

```bash
cd experiments/experiment_5_bspline_B_recovery

Rscript -e 'renv::restore(prompt = FALSE)'
Rscript -e 'renv::status()'
Rscript code/06_plot_three_model_comparison.R
```

Expected output directory:

```text
figures/comparison_three_models/
```

Expected principal figure files:

```text
01_three_model_RMSE_comparison.pdf
01_three_model_RMSE_comparison.svg
01_three_model_RMSE_comparison.png
01_three_model_RMSE_comparison.tiff
02_three_model_mean_error_comparison.pdf
02_three_model_mean_error_comparison.svg
02_three_model_mean_error_comparison.png
02_three_model_mean_error_comparison.tiff
```

The script also regenerates figure source-data CSV files under:

```text
figures/comparison_three_models/source_data/
```

Run the additional source-data check:

```bash
Rscript -e 'x <- read.csv("figures/comparison_three_models/source_data/figure_2_mean_error_histograms_source_data.csv"); stopifnot(nrow(x) == 1800L, length(unique(x$task_id)) == 200L); cat("Reproduction checks passed.\n")'
```

A successful figure reproduction validates the retained plotting inputs and
plot-generation workflow. It does not independently repeat the original model
fitting.

See `experiments/experiment_5_bspline_B_recovery/REPRODUCE_FIGURES.md` for the
complete plotting notes.

## 3. Inspect the retained final evidence

The final three-model comparison files are located under:

```text
experiments/experiment_5_bspline_B_recovery/results/comparison_three_models/
```

Important files include:

```text
three_model_overall_summary.csv
pairwise_RMSE_summary.csv
three_model_paired_input_check.csv
three_model_paired_metrics_wide.csv
three_model_run_check_summary.csv
three_model_task_metrics_long.csv
```

The validated combined B-spline outputs are under:

```text
experiments/experiment_5_bspline_B_recovery/results/combined/bspline/
```

The final run checks require:

- exactly 200 manifest tasks;
- exactly 200 metric rows for each fitted model;
- exactly 200 paired comparison rows;
- exactly 70 observation times per model and task;
- matching task IDs, seeds, and observed-data checksums;
- Gamma paths labelled as particle filtering means;
- B-spline paths labelled as selected deterministic trajectories;
- constant-\(B\) paths labelled as repeated static estimates;
- normalized final comparison truth with \(B(5)=4\);
- `week <= 5` for the through-week-5 period;
- `week > 5` for the after-week-5 period.

## 4. Full reproduction: rerun the fitted models

The full route is a computational workflow, not a figure-only command. Do not
start it until you have read both experiment READMEs and confirmed that the
required R and Slurm environments are available.

The historical Experiment 2–4 batch scripts recorded:

```bash
module load StdEnv/2020
module load r/4.1.2
```

However, the complete historical package environment for the original HPC
fits was not preserved as a lockfile or container. The Experiment 5
`renv.lock` applies only to final figure generation. Bit-for-bit reproduction
of every historical particle-filter realization is not claimed.

### Important preservation rule

Keep the downloaded archive unchanged as a reference copy. Perform a full
rerun in a separate working copy. Do not overwrite, delete, or silently replace
the archived shared data or completed raw task directories.

### Step A: run Experiment 4

Experiment 4 creates the accepted shared data and fits the Gamma-noise and
constant-\(B\) models.

```bash
cd experiments/experiment_4_nmif600_model_comparison
```

First run the documented full-setting pilot:

```bash
bash hpc/submit_pilot.sh
```

Review the pilot outputs and resource records. Then submit the complete
200-task workflow:

```bash
bash hpc/submit_all.sh
```

The workflow submits, in dependency order:

1. the shared-data generation array;
2. the Gamma-noise and constant-\(B\) fitting arrays;
3. combination, validation, paired two-model analysis, and diagnostics.

Follow the full instructions in:

```text
experiments/experiment_4_nmif600_model_comparison/README.md
```

### Step B: validate the Experiment 5 inputs

Experiment 5 must use the accepted Experiment 4 data directly. It must not
simulate a new set of observations.

```bash
cd ../experiment_5_bspline_B_recovery
Rscript code/01_validate_paired_inputs.R
```

The validation requires exactly 200 accepted tasks and checks task IDs,
simulation seeds, acceptance metadata, and data checksums.

### Step C: run the B-spline production workflow

Experiment 5 uses a production task-1 pilot followed by the formal array for
tasks 2–200.

Submit the task-1 production pilot:

```bash
mkdir -p logs/pilot results_pilot/bspline
sbatch hpc/03_task1_production_pilot.sh
```

After it completes, validate and promote task 1 exactly as documented in the
Experiment 5 README. Do not overwrite an existing formal task directory.

Then submit the formal B-spline array and dependent post-processing job:

```bash
bash hpc/submit_all.sh
```

The intended production total is:

```text
one validated and promoted task-1 pilot
+ 199 formal array tasks
= 200 B-spline fits
```

Follow the authoritative instructions in:

```text
experiments/experiment_5_bspline_B_recovery/README.md
```

### Step D: rebuild the final three-model analysis

After all 200 newly fitted B-spline tasks pass validation, combine the new
formal result tree rather than the bundled historical archive:

```bash
Rscript code/03_combine_results.R \
  results/bspline \
  results/combined/bspline \
  results/paired_input_manifest.csv

Rscript code/05_compare_three_models.R
Rscript code/06_plot_three_model_comparison.R
```

To reprocess the unchanged B-spline task outputs supplied in the full bundle
instead of a new fit, use `results_raw/bspline` as the first argument in a
separate archival-reprocessing run.

The final pairing checks must pass before the summary tables or figures are
treated as valid.

## 5. Package directory guide

```text
experiments/
├── experiment_1_gamma_B_recovery/
│   └── developmental Gamma-noise workflows
├── experiment_2_large_scale_gamma_B_recovery_HPC/
│   └── high-particle single-data-set diagnostic
├── experiment_3_gamma_B_recovery_accuracy/
│   └── earlier repeated Gamma-noise recovery study
├── experiment_4_nmif600_model_comparison/
│   ├── shared accepted data
│   ├── Gamma-noise workflow and outputs
│   └── constant-B workflow and outputs
└── experiment_5_bspline_B_recovery/
    ├── deterministic six-basis cubic B-spline workflow
    ├── validated B-spline outputs
    └── final paired three-model comparison and figures
```

Experiments 1–3 are supporting studies. Experiment 4 is the primary
two-model computational foundation. Experiment 5 is the final primary
three-model analysis.

Do not flatten these directories into a single file level. The subdirectory
structure is part of the workflow because scripts rely on relative paths.

## 6. Interpretation limits

The package supports inspection and reproduction of a controlled simulation
study. It does not demonstrate performance on real surveillance data.

The reported results are conditional on:

- the configured stochastic SIR process;
- the negative-binomial measurement model;
- the single step change in \(B(t)\);
- acceptance requiring `max(H) > 20`;
- the tested starting grids;
- the selected particle counts and `Nmif` settings;
- the 200 accepted replicates;
- the endpoint convention attached to each artifact.

Independent particle-filter likelihoods are fitting diagnostics and are not
complexity-adjusted model-selection criteria. Selected-task sampled latent
trajectories are illustrative and should not be substituted for the
observation-time Gamma filtering means used in the final primary metrics.

## 7. Version and citation

The current released version is `v1.0.0`:

- GitHub:
  https://github.com/ZeuS2U35-YX/Simulation-Based-Transmission-Rate-Recovery/releases/tag/v1.0.0
- DOI:
  https://doi.org/10.5281/zenodo.22127917

Version `v1.1.0` is not yet a public release, and no `v1.1.0` DOI is claimed.
Until it is published, a candidate full-replication ZIP should be identified
by its filename, source commit, and SHA-256 checksum rather than cited as a
released version.

See `CITATION.cff` for the current released citation metadata. See `LICENSE`
for the MIT License.

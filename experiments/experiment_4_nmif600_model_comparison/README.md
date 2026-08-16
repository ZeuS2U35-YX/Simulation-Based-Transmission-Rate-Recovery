# Experiment 4: Gamma-noise versus constant-B comparison at `Nmif = 600`

Experiment 4 is the canonical computational analysis in this repository and the primary source of quantitative evidence for the report. It replaces Experiment 3's earlier `Nmif = 100` workflow as the source of final Gamma-noise model numerical claims. Experiments 1-3 remain developmental or supporting studies.

## Research question

This experiment asks which of two fitted partially observed Markov process (POMP) models more accurately recovers the prescribed transmission-rate path `B(t)`:

1. a Gamma-noise model with latent time-varying `B(t)`;
2. a deliberately restricted constant-B comparator that cannot represent temporal variation in `B(t)`.

Both models are fitted by iterated filtering (MIF2) with `Nmif = 600` on **the same 200 accepted simulated epidemic data sets**. Acceptance requires `max(H) > 20`, so the reported performance is conditional on informative accepted outbreaks rather than unconditional over all attempted simulations. Numerical comparisons use residual sum of squares (RSS), root mean squared error (RMSE), and signed errors of the observation-time Gamma filtering mean and the repeated fitted constant-B estimate. Independent particle-filter log likelihoods are retained as descriptive fitting diagnostics, not as complexity-adjusted model-selection criteria.

## Fixed experiment settings

The shared data-generating model uses:

- `B(t) = 4` before week 5;
- `B(t) = 2` from week 5 onward;
- `mu_IR = 3`, `N = 10000`, `rho = 0.5`, `k = 10`;
- 10 weeks, 70 daily observation times, Euler step `1/30` week;
- acceptance rule `max(H) > 20`;
- simulation seeds 1001 through 1200.

Both fitted models use:

- `Nmif = 600`;
- `Np_mif = 5000`;
- five independent likelihood evaluations with `Np_eval = 50000`;
- a final particle filter with `Np_final = 50000`;
- geometric cooling with `cooling.fraction.50 = 0.5`;
- positive parameter transformations on all estimated parameters.

Gamma-noise model starts:

- `B0 in {2, 4, 6}`;
- `sigma_beta in {0.10, 0.30, 0.45}`;
- `rw.sd`: `B0 = ivp(0.20)`, `sigma_beta = 0.05`.

Constant-B starts:

- `Beta in {1, 2, 3, 4, 5, 6}`;
- `rw.sd`: `Beta = 0.05`.

## Why the data are strictly shared

`code/01_generate_shared_data_task.R` creates each accepted data set once and writes it under `shared_data/task_###/`. Both fitting scripts read the same `observed_data.csv`. MD5 checksums are verified before fitting and are copied into every model output. The final paired comparison refuses to run unless the Gamma and constant-B outputs have identical data checksums and simulation seeds for every task.

## Directory structure

```text
config/                         Single experiment configuration
code/                           Data generation, fitting, combination, analysis
hpc/                            Slurm scripts and submission wrappers
shared_data/task_###/           One immutable shared data set per task
results_raw/gamma/task_###/     Gamma task-level outputs
results_raw/constant/task_###/  Constant-B task-level outputs
results/combined/               Validated model-specific combined tables
results/comparison/             Paired model-comparison tables
results/selected_trajectory/    Selected task-1 and task-117 comparisons and provenance
figures/comparison/             Final model-comparison PDFs
figures/convergence/            Nmif = 600 convergence diagnostics
logs/                           Slurm stdout and stderr
downloads/                      Final downloadable tar.gz archive
```

Each task is written through a temporary directory and atomically renamed only after completion. A completed task contains a `COMPLETE` marker and is skipped on resubmission. Existing complete results are never overwritten automatically.

## HPC assumptions

The included Slurm files follow the environment used by the earlier experiments:

```bash
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="$HOME/packages-R4.1"
```

The Slurm scripts intentionally do **not** hard-code a partition. On the Frontenac CAC system used for this run, explicit partition names caused submission failures for this account, while submitting without `--partition` allowed the scheduler to choose an eligible partition automatically. Fitting tasks request one CPU, 12 GB memory, and a 24-hour wall-time limit.

The completed five-task pilot showed that these limits were conservative. Gamma-noise tasks took approximately 1 hour 9-11 minutes, while constant-B tasks took approximately 33-36 minutes. Peak resident memory was about 200 MB, so the requested resources should be interpreted as safety ceilings rather than expected usage.

## Recommended workflow

Run every command from the Experiment 4 root directory.

### 1. Upload and enter the folder

```bash
cd experiment_4_nmif600_model_comparison
```

### 2. Submit five full-setting pilot tasks

The pilot uses diagnostic tasks 1, 50, 100, 150, and 200. They are not a random sample. The pilot uses the full `Nmif = 600` setting, particle counts, starting grids, and likelihood evaluations. Its completed task outputs are production-quality and are reused by the full array.

```bash
bash hpc/submit_pilot.sh
```

Check status:

```bash
squeue -u "$USER"
```

After completion, inspect elapsed time and memory either directly with Slurm:

```bash
sacct -j <GAMMA_JOB_ID>,<CONSTANT_JOB_ID> \
  --format=JobID,State,Elapsed,TotalCPU,MaxRSS,AllocCPUS
```

or with the included helper:

```bash
bash hpc/report_pilot_resources.sh
cat results/pilot_resource_summary.txt
```

`submit_pilot.sh` saves the four submitted job IDs automatically in `results/pilot_job_ids.env`.

Pilot diagnostics are written to:

```text
figures/pilot/
results/pilot_comparison/
```

The pilot post-processing produces metric distributions and convergence diagnostics, but no primary `B(t)` curve: an average of five pilot tasks is not a coherent latent trajectory. The convergence PDFs contain all starts for the five selected diagnostic tasks. A dotted vertical line marks the beginning of the final 100 iterations. The script also writes tail slopes to `convergence_tail_summary.csv`. Stabilization in these selected traces provides empirical support for using `Nmif = 600`, but does not prove convergence for all 200 fitted data sets.

### 3. Submit the complete 200-task experiment

After reviewing the pilot:

```bash
bash hpc/submit_all.sh
```

The wrapper submits, in order:

1. the shared-data array;
2. the Gamma-noise model and constant-B model arrays after all shared data finish;
3. combination, validation, paired analysis, convergence plotting, and packaging after both model arrays finish.

Pilot task directories already marked `COMPLETE` are validated and skipped.

### 4. Download the completed experiment

After the post-processing job succeeds, download:

```text
downloads/experiment_4_nmif600_model_comparison_outputs.tar.gz
```

The archive includes code, configuration, shared data, task-level outputs, combined tables, figures, and logs. Its file list is stored in:

```text
downloads/experiment_4_nmif600_model_comparison_manifest.txt
```

## Main outputs

Model-specific combined outputs:

```text
results/combined/gamma/
results/combined/constant/
```

Paired comparison outputs:

- `shared_data_pair_check.csv`;
- `model_task_metrics.csv`;
- `paired_model_comparison.csv`;
- `overall_model_comparison.csv`;
- `likelihood_gaps_by_model.csv`;
- `starting_value_sensitivity.csv`.

Main figures and their estimands:

| Figure | Quantity displayed | Scientific role |
| --- | --- | --- |
| `01_selected_task_B_trajectory_comparison.pdf` | Prescribed truth, task-1 Gamma filtering mean, and the task-1 fitted constant-B estimate repeated over 70 observation times | Selected-task illustration of transmission-rate recovery |
| `02_RSS_distributions.pdf` | Across-task distributions of Gamma filtering-mean RSS and constant repeated-static-estimate RSS | Aggregate recovery comparison |
| `03_bias_distributions.pdf` | Overall, pre-switch, and post-switch signed mean errors for both estimators | Direction and timing of recovery error |
| `04_paired_RMSE_scatter.pdf` | Paired task-level RMSE values on the same simulated data sets | Direct within-task model comparison |
| `05_RMSE_distributions.pdf` | Across-task RMSE distributions | Aggregate error magnitude |
| `06_independent_loglik_difference.pdf` | Gamma-noise minus constant-B independent log likelihood | Descriptive fit diagnostic only |
| `07_task1_infectious_path_comparison.pdf` | Task-1 true infectious path and both models' filtering means | Selected-task infectious-state recovery |
| `08_task117_B_trajectory_comparison.pdf` | Prescribed truth, one task-117 ancestry-preserving Gamma trajectory, and the task-117 fitted constant-B estimate | Additional trajectory illustration |
| `09_task117_infectious_path_comparison.pdf` | Task-117 true infectious path and both models' filtering means | Additional infectious-state illustration |

Model-specific convergence diagnostic PDFs are stored in `figures/convergence/`.

The figures use a restrained publication style with a white background, no decorative grid, and color-plus-line-type encoding that remains interpretable in grayscale. Truth is shown as a solid black line, the Gamma-noise result as a blue dashed line, and the constant-B result as an orange dot-dashed line. A light-gray vertical reference marks the week-5 change point.

The Gamma curves in Figures 01, 07, and 09 are plug-in-parameter, finite-particle filtering means at the observation times. The constant-B infectious curves in Figures 07 and 09 are the corresponding filtering means from the constant-B model. Figure 08 instead preserves the original task-117 ancestry and is a finite-particle smoothing-trajectory approximation conditioned on the complete observation series. It is not a filtering mean, an exact posterior draw, or an across-task average. The constant curves in Figures 01 and 08 are fitted static estimates repeated over time, not latent trajectories.

Figures 01-05 and 07 provide the primary recovery summaries for the report. Figures 08 and 09 retain the additional task-117 illustration. Figure 06 should be presented only with its stated likelihood-interpretation boundary.

## Canonical lightweight figure regeneration

`code/07_generate_task1_comparison_figures.R` generates the task-1 artifacts, and `code/08_generate_task117_comparison_figures.R` generates the additional task-117 artifacts. Both scripts require `pomp`, `digest`, and `ggplot2`. They use `svglite` and `ragg`, when available, for editable SVG and high-resolution raster export.

The task-1 script reproduces the saved Gamma filtering mean for `B(t)` and obtains the two infectious filtering means. The task-117 script preserves the stored ancestry trajectory for Figure 08 and obtains the infectious filtering means for Figure 09. Each infectious filter uses 50,000 particles at an already fitted parameter vector. Neither script reruns MIF2 or the independent likelihood evaluations. Figure 07 uses Experiment 4 task 1 with simulation seed 1001, not the separate illustrative simulation from Experiment 1.

The canonical lightweight generators are:

- `code/07_generate_task1_comparison_figures.R` for Figures 01 and 07;
- `code/08_generate_task117_comparison_figures.R` for Figures 08 and 09;
- `code/05_compare_models.R` for Figures 02-06 and selected-artifact validation;
- `code/06_make_convergence_diagnostics.R` for the convergence PDFs and tail summary.

Run from the Experiment 4 directory:

```bash
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

`code/07_generate_selected_B_trajectory.R` is retained as a backward-compatible wrapper for the task-1 generator. The Python regeneration script is non-canonical; the R workflow above is the authoritative route for Figures 01 and 07-09. Pilot regeneration deliberately omits the selected-task figures.

## Completed-run validation and headline results

The packaged run passed the model-specific combination checks for all 200 tasks. Both models have 200 of 200 tasks present, no missing tasks, no combination problems, unique task IDs, `Nmif = 600`, and successful best-fit status. The paired data check confirms matching simulation seeds and identical observed-data MD5 checksums for the two models on every task.

The main recovery summaries from `results/comparison/overall_model_comparison.csv` compare per-task Gamma filtering means with per-task fitted constant values repeated over the 70 observation times:

- mean RSS: Gamma-noise model 24.09; constant-B 102.83;
- mean RMSE: Gamma-noise model 0.575; constant-B 1.199;
- mean absolute overall bias: Gamma-noise model 0.160; constant-B 0.606;
- mean after-switch bias: Gamma-noise model 0.139; constant-B 1.576;
- the Gamma-noise model has lower RSS and lower RMSE on all 200 paired tasks;
- the Gamma-noise model has lower absolute overall bias on 90.5% of paired tasks.

These summaries are unchanged by the selected-task figures. Figure 01 uses the task-1 Gamma filtering mean already included in the 200-task recovery metrics. Figures 07-09 are additional state or trajectory illustrations and do not enter those metrics. The independent likelihood comparison remains a descriptive fitting diagnostic rather than a complexity-adjusted model-selection test. Neither MIF2 nor the independent likelihood evaluations were rerun for the selected-task figure revision.

## Rerunning an incomplete or failed task

A complete task is skipped. An incomplete task directory causes the script to stop rather than overwrite it. Inspect the log first, then remove only the failed task directory and resubmit that array index. For example:

```bash
rm -rf results_raw/gamma/task_149
sbatch --array=149 hpc/02_gamma_array.sh
```

After repairs, rerun post-processing:

```bash
sbatch hpc/04_postprocess.sh
```

Do not delete or regenerate `shared_data/task_149` unless the shared data themselves failed validation. Both models must continue to use the same saved data set.

## Interpretation boundary

The experiment evaluates recovery only under the stated simulation model, acceptance rule, starting grids, particle counts, and `Nmif = 600`. Its performance estimand is conditional on accepted outbreaks satisfying `max(H) > 20`. Stable traces in diagnostic tasks 1, 50, 100, 150, and 200 support the computational choice for those tasks but do not prove convergence for all 200 fitted data sets. Independent particle-filter evaluations, multi-start results, and convergence diagnostics should therefore be interpreted together.

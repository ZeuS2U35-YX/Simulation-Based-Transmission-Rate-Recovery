# Experiment 4: two-model computational foundation at `Nmif = 600`

Experiment 4 is the primary computational foundation for the final study. It provides the 200 accepted data sets and the `Nmif = 600` Gamma-noise and constant-B fits. Experiment 5 fits a deterministic six-basis cubic B-spline to the same data and is the final paired three-model manuscript analysis. Experiment 4 therefore remains foundational rather than the final model comparison. Experiments 1–3 remain developmental or supporting studies.

## Research question

This experiment asks which of two fitted partially observed Markov process (POMP) models more accurately recovers the prescribed transmission-rate path `B(t)`:

1. a Gamma-noise model with latent time-varying `B(t)`;
2. a deliberately restricted constant-B comparator that cannot represent temporal variation in `B(t)`.

The B-spline model is intentionally not fitted in Experiment 4; it is added in Experiment 5 as the third fitted model.

Both models are fitted by iterated filtering (MIF2) with `Nmif = 600` on **the same 200 accepted simulated epidemic data sets**. Acceptance requires `max(H) > 20`, so the reported performance is conditional on informative accepted outbreaks rather than unconditional over all attempted simulations. For each data set, the best Gamma fit is followed by one ancestry-preserving sampled latent `B(t)` trajectory from the final particle filter. Numerical comparisons use residual sum of squares (RSS), root mean squared error (RMSE), signed mean error, and absolute overall bias (AOB) calculated from that sampled trajectory and the repeated fitted constant-B estimate. Independent particle-filter log likelihoods are retained as descriptive fitting diagnostics, not as complexity-adjusted model-selection criteria.

These sampled-trajectory results document the historical two-model Experiment 4 analysis. The final manuscript-facing three-model metrics are calculated in Experiment 5 from Gamma particle-filtering means, selected deterministic B-spline trajectories, and repeated constant-B estimates.

## Fixed experiment settings

The shared data-generating model uses:

- `B(t) = 4` before week 5;
- `B(t) = 2` from week 5 onward;
- `mu_IR = 3`, `N = 10000`, `rho = 0.5`, `k = 10`;
- 10 weeks, 70 daily observation times, Euler step `1/30` week;
- acceptance rule `max(H) > 20`;
- simulation seeds 1001 through 1200.

The Experiment 4 code and retained artifacts use the historical endpoint convention `B(t) = 4` for `t < 5` and `B(t) = 2` for `t >= 5`, including `B(5) = 2`. The final Experiment 5 three-model reporting normalizes the shared comparison truth to `B(5) = 4` and defines the low period by `t > 5`. Experiment 4 outputs remain unchanged for provenance and must not be mixed with the normalized Experiment 5 primary metrics.

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
shared_data/task_###/           One immutable shared data set per task [full bundle or generated locally; not in Git]
results_raw/gamma/task_###/     Gamma task-level outputs [full bundle or generated locally; not in Git]
results_raw/constant/task_###/  Constant-B task-level outputs [full bundle or generated locally; not in Git]
results/combined/               Validated model-specific combined tables
results/comparison/             Paired model-comparison tables
results/selected_trajectory/    Selected task-1 and task-117 comparisons and provenance
figures/comparison/             Final model-comparison PDFs
figures/convergence/            Nmif = 600 convergence diagnostics
logs/                           Slurm stdout and stderr [generated locally; not in Git]
downloads/                      Optional local HPC package output [generated locally; not in Git]
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

## Distribution boundary

A fresh Git clone contains code, configuration, validated combined tables, selected artifacts, and tracked figures. It does not contain `shared_data/`, `results_raw/`, Slurm logs, or `downloads/`. The repository-level full-replication ZIP adds the accepted shared data and required task-level raw results at their original relative paths. Use the full bundle to inspect or recombine the completed historical tasks; use the workflow below to generate a new complete run.

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

### 4. Optional local packaging

After a completed HPC run, `bash hpc/package_outputs.sh` creates the local archive
`downloads/experiment_4_nmif600_model_comparison_outputs.tar.gz` and its
manifest. These generated files are not tracked by Git and are not public
download links. The published one-download artifact is the repository-level
`v1.1.0` full-replication ZIP described in the root README.

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
| `01_selected_task_B_trajectory_comparison.pdf` | Prescribed truth, task-1 Gamma sampled latent trajectory, and the task-1 fitted constant-B estimate | Selected-task illustration of transmission-rate recovery |
| `02_RSS_distributions.pdf` | Across-task distributions of sampled-trajectory RSS and constant repeated-static-estimate RSS | Aggregate recovery comparison |
| `03_mean_error_distributions.pdf` | Overall, pre-switch, and post-switch signed mean errors for both estimators | Direction and timing of recovery error |
| `04_paired_RMSE_scatter.pdf` | Paired task-level RMSE values on the same simulated data sets | Direct within-task model comparison |
| `05_RMSE_distributions.pdf` | Across-task RMSE distributions | Aggregate error magnitude |
| `06_independent_loglik_difference.pdf` | Gamma-noise minus constant-B independent log likelihood | Descriptive fit diagnostic only |
| `08_task117_B_trajectory_comparison.pdf` | Prescribed truth, task-117 Gamma sampled latent trajectory, and the task-117 fitted constant-B estimate | Second selected-task illustration |

Model-specific convergence diagnostic PDFs are stored in `figures/convergence/`.

The figures use a restrained publication style with a white background, no decorative grid, and color-plus-line-type encoding that remains interpretable in grayscale. Truth is shown as a solid black line, the Gamma-noise result as a blue dashed line, and the constant-B result as an orange dot-dashed line. A light-gray vertical reference marks the week-5 change point.

The Gamma curves in Figures 01 and 08 are the exact sampled paths used for the corresponding task-level recovery metrics. Each is one ancestry-preserving finite-particle plug-in approximation to a smoothing trajectory conditional on the complete observation series. It is not a filtering mean, an exact posterior draw, an uncertainty interval, or an across-task average. The constant curves are fitted static estimates repeated over time, not latent trajectories. The figures include `t0`; RSS, RMSE, mean error, and AOB use only the 70 observation times.

Figures 01–05 and 08 document the historical Experiment 4 two-model recovery analysis. They are not the final three-model manuscript figures; those are produced in Experiment 5. Figure 06 should be presented only with its stated likelihood-interpretation boundary.

## Post-processing and trajectory reconstruction

`code/09_regenerate_sampled_B_trajectories.R` reconstructs the 200 canonical Gamma paths from the saved best-fit parameter records. It does not rerun MIF2 or the five independent likelihood evaluations. For each task, it runs the final 50,000-particle filter with `filter.traj = TRUE`, extracts one path with `filter_traj()`, verifies the saved final-filter log likelihood, and writes the 70 observation-time values.

`code/07_generate_task1_comparison_figures.R` and `code/08_generate_task117_comparison_figures.R` then read the canonical sampled paths used by the recovery metrics and generate the selected-task figures. They do not run another particle filter, so the displayed curve and the metric input cannot diverge.

The Experiment 4 post-processing generators are:

- `code/09_regenerate_sampled_B_trajectories.R` for the 200 sampled Gamma paths;
- `code/07_generate_task1_comparison_figures.R` for Figure 01;
- `code/08_generate_task117_comparison_figures.R` for Figure 08;
- `code/05_compare_models.R` for Figures 02-06 and selected-artifact validation;
- `code/06_make_convergence_diagnostics.R` for the convergence PDFs and tail summary.

Run from the Experiment 4 directory:

```bash
Rscript code/09_regenerate_sampled_B_trajectories.R \
  results/combined/gamma/combined_B_paths.csv \
  results/combined/gamma/sampled_B_trajectory_provenance.csv \
  4 \
  1:200

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

In a fresh clone, `code/05_compare_models.R` and `code/06_make_convergence_diagnostics.R` can use the tracked combined tables. The full 200-path reconstruction in `code/09_regenerate_sampled_B_trajectories.R` requires the accepted shared data and completed fit records supplied by the full bundle or a new run; the selected-task scripts should be treated as downstream of those validated reconstructed paths.

`code/07_generate_selected_B_trajectory.R` is retained as a backward-compatible wrapper for the task-1 generator. The Python regeneration script is non-canonical; the R workflow above is the authoritative route for Figures 01 and 08. Pilot regeneration deliberately omits the selected-task figures.

## Historical completed-run validation and two-model results

The packaged run passed the model-specific combination checks for all 200 tasks. Both models have 200 of 200 tasks present, no missing tasks, no combination problems, unique task IDs, `Nmif = 600`, and successful best-fit status. The paired data check confirms matching simulation seeds and identical observed-data MD5 checksums for the two models on every task.

The main recovery summaries from `results/comparison/overall_model_comparison.csv` compare one sampled Gamma trajectory per task with the fitted constant value repeated over the same 70 observation times:

- mean RSS: Gamma-noise model 33.054; constant-B 102.833;
- mean RMSE: Gamma-noise model 0.653; constant-B 1.199;
- mean AOB: Gamma-noise model 0.229; constant-B 0.606;
- mean post-switch signed mean error: Gamma-noise model 0.008; constant-B 1.576;
- the Gamma-noise model has lower RSS and lower RMSE on 97.5% of paired tasks;
- the Gamma-noise model has lower AOB on 88.5% of paired tasks.

Figures 01 and 08 use the same task-1 and task-117 sampled trajectories already included in the 200-task recovery metrics. The independent likelihood comparison remains a descriptive fitting diagnostic rather than a complexity-adjusted model-selection test. Neither MIF2 nor the independent likelihood evaluations were rerun for this correction.

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

# Experiment 4: Gamma-noise versus constant-B comparison at `Nmif = 600`

Experiment 4 is the canonical computational analysis in this repository and the primary source of quantitative evidence for the report. It replaces Experiment 3's earlier `Nmif = 100` workflow as the source of final Gamma-noise model numerical claims. Experiments 1-3 remain developmental or supporting studies.

## Research question

This experiment asks which of two fitted partially observed Markov process (POMP) models more accurately recovers the prescribed transmission-rate path `B(t)`:

1. a Gamma-noise model with latent time-varying `B(t)`;
2. a deliberately restricted constant-B comparator that cannot represent temporal variation in `B(t)`.

Both models are fitted by iterated filtering (MIF2) with `Nmif = 600` on **the same 200 accepted simulated epidemic data sets**. These 200 data sets came from 207 simulation attempts (96.6% acceptance). Acceptance requires `max(H) > 20`, so the reported performance is conditional on informative accepted outbreaks rather than unconditional over all attempted simulations. For each selected Gamma fit, the primary recovery estimate is the particle filtering mean of the latent `B(t)` state at each observation time, conditional on observations through that time and evaluated at the selected plug-in parameter estimate. Numerical comparisons use residual sum of squares (RSS), root mean squared error (RMSE), signed mean error, and absolute overall bias (AOB) calculated from this filtering-mean trajectory and the repeated fitted constant-B estimate. The active selected-task figures show `B(t)` and `I(t)` for the previously designated Tasks 1 and 117. One figure contains filtering means; the other contains the first seeded forward simulation from each fitted model, with the selected parameters fixed and no observation conditioning after fitting. Only the filtering-mean `B(t)` values enter the metrics. Independent particle-filter log likelihoods are retained as descriptive fitting diagnostics, not as complexity-adjusted model-selection criteria.

## Fixed experiment settings

The shared data-generating model uses:

- `B(t) = 4` before week 5;
- `B(t) = 2` from week 5 onward;
- `mu_IR = 3`, `N = 10000`, `rho = 0.5`, `k = 10`;
- 10 weeks and 70 daily observation times; maximum Euler step `1/30` week, giving five equal substeps per day and an actual step of `1/35` week;
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
| `01_selected_tasks_particle_filtering_mean_trajectories.pdf` | Task-1 and task-117 truth and fitted-model particle filtering means for `B(t)` and `I(t)` | Active report Figure 6; data-conditioned point summaries for an easier and a more challenging case |
| `09_selected_tasks_forward_simulations_at_fitted_parameters.pdf` | One seeded joint forward simulation of `B(t)`/`I(t)` from each fitted model for tasks 1 and 117 | Active report Figure 9; unconditional stochastic illustration at fixed fitted parameters only |
| `01_selected_task_B_trajectory_comparison.pdf` | Prescribed truth, task-1 Gamma ancestry path, and the task-1 fitted constant-B estimate | Legacy repository provenance; not an active report figure or metric input |
| `02_RSS_distributions.pdf` | Across-task distributions of filtering-mean RSS and constant repeated-static-estimate RSS | Aggregate recovery comparison |
| `03_mean_error_distributions.pdf` | Overall, through-week-5, and after-week-5 signed mean errors for both estimators | Direction and timing of recovery error |
| `04_paired_RMSE_scatter.pdf` | Paired task-level RMSE values on the same simulated data sets | Direct within-task model comparison |
| `05_RMSE_distributions.pdf` | Across-task RMSE distributions | Aggregate error magnitude |
| `06_independent_loglik_difference.pdf` | Gamma-noise minus constant-B independent log likelihood | Descriptive fit diagnostic only |
| `08_task117_B_trajectory_comparison.pdf` | Prescribed truth, task-117 Gamma ancestry path, and the task-117 fitted constant-B estimate | Legacy repository provenance; not used in the report or metrics |

Model-specific convergence diagnostic PDFs are stored in `figures/convergence/`.

The figures use a restrained publication style with a white background, no decorative grid, and color-plus-line-type encoding that remains interpretable in grayscale. Truth is shown as a solid black line, the Gamma-noise result as a blue dashed line, and the constant-B result as an orange dot-dashed line. A light-gray vertical reference marks the week-5 change point.

In the active filtering-mean figure, each observation-time value is a weighted particle average conditional on `Y_1:n`; the week-0 markers are fitted or fixed initial values rather than filtering means. In the forward-simulation figure, the Gamma `B(t)` and `I(t)` curves within a task come from one joint simulation initialized at the fitted `B0` and evolved with the fitted `sigma_beta`; the constant-model `I(t)` curve is simulated with the fitted static `B`. No observations are assimilated after the fitted parameters are fixed. Neither figure integrates parameter uncertainty. RSS, RMSE, mean error, and AOB use only the 70 observation-time `B(t)` filtering means.

Figures 02-05 provide the primary recovery summaries. The selected-task filtering-mean and forward-simulation figures provide complementary data-conditioned point-summary and model-generative views. Task 1 is a comparatively well-recovered case; Task 117 is a more challenging case and was already designated as the second illustration before this figure revision. They are not claimed to be representative of all 200 tasks. The historical B-only ancestry-path Figures 01 and 08 remain repository artifacts, and Figure 06 should be presented only with its stated likelihood-interpretation boundary.

## Canonical lightweight regeneration and input requirements

`code/10_regenerate_filtering_mean_B_paths.R` reconstructs the 200 primary Gamma filtering-mean trajectories from the saved best-fit parameter records. It does not rerun MIF2 or the five independent likelihood evaluations. For each task, it runs the final 50,000-particle filter with `filter.mean = TRUE`, verifies the saved final-filter log likelihood, and writes the 70 observation-time filtering means. The output records that the estimates condition on `Y_1:n`, use plug-in parameter estimates, and do not integrate parameter uncertainty.

`code/11_generate_selected_task_filtering_forward_figures.R` builds the two active Task 1/117 manuscript figures and their source-data tables. It uses the retained canonical Gamma `B(t)` filtering means and reruns the seeded final fitted Gamma and constant-model filters with filtering-mean storage enabled to recover `I(t)` means. It then fixes each model at its selected parameters and uses the first forward realization from a prespecified seed offset, without assimilating observations or screening trajectories. The provenance table records the shared-data checksums, fitted parameters, filter and forward-simulation seeds, filter likelihoods, runtime versions, retained-versus-current `B(t)` differences, empirical recovery-error percentiles, conditioning semantics, and the no-screening selection rule.
The figure exporter requires `ggplot2`, `patchwork`, `svglite`, `ragg`, and Ghostscript in addition to the experiment's POMP runtime.

`code/09_regenerate_sampled_B_trajectories.R` separately reconstructs the sampled Gamma paths retained for illustration. It runs the final filter with trajectory storage enabled and extracts one path with `filter_traj()`. These sampled paths are not used by `code/05_compare_models.R` for the primary recovery metrics.

This reconstruction is lightweight relative to the full IF2 workflow, but it
is not self-contained from the compact combined results alone. For every task,
the script reads the exact fitted case series from
`shared_data/task_###/observed_data.csv` and checks its MD5 value against the
saved best-fit record. The full workflow creates `shared_data/` in stage 1,
and the completed-experiment download archive includes it. The directory is not
committed in this compact repository milestone. Before running the command,
either restore `shared_data/` from the archive or run the shared-data generation
stage. If `shared_data/` is unavailable, retain the existing canonical
`results/combined/gamma/combined_B_filtering_means.csv` and run only the
downstream comparison and convergence commands. The retained sampled paths may
still be used to rebuild the illustration figures.

`code/07_generate_task1_comparison_figures.R` and `code/08_generate_task117_comparison_figures.R` read the retained sampled paths and generate the illustration figures. They do not run another particle filter.

The canonical lightweight generators are:

- `code/10_regenerate_filtering_mean_B_paths.R` for the 200 primary Gamma filtering-mean trajectories;
- `code/11_generate_selected_task_filtering_forward_figures.R` for the active Task 1/117 filtering-mean and forward-simulation figures;
- `code/09_regenerate_sampled_B_trajectories.R` for optional sampled Gamma illustrations;
- `code/07_generate_task1_comparison_figures.R` for Figure 01;
- `code/08_generate_task117_comparison_figures.R` for Figure 08;
- `code/05_compare_models.R` for Figures 02-06 and selected-artifact validation;
- `code/06_make_convergence_diagnostics.R` for the convergence PDFs and tail summary.

Run from the Experiment 4 directory:

```bash
# Requires shared_data/task_###/observed_data.csv for tasks 1-200.
Rscript code/10_regenerate_filtering_mean_B_paths.R \
  results/combined/gamma/combined_B_filtering_means.csv \
  results/combined/gamma/filtering_mean_provenance.csv \
  4 \
  1:200

# Active Task 1/117 manuscript figures and source-data tables.
Rscript code/11_generate_selected_task_filtering_forward_figures.R shared_data

# Optional sampled trajectories used only for illustration.
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

`code/07_generate_selected_B_trajectory.R` is retained as a backward-compatible wrapper for the task-1 generator. The Python regeneration script is non-canonical; the R workflow above is the authoritative route for Figures 01 and 08. Pilot regeneration deliberately omits the selected-task figures.

## Completed-run validation and headline results

The packaged run passed the model-specific combination checks for all 200 tasks. Both models have 200 of 200 tasks present, no missing tasks, no combination problems, unique task IDs, `Nmif = 600`, and successful best-fit status. The paired data check confirms matching simulation seeds and identical observed-data MD5 checksums for the two models on every task.

The main recovery summaries from `results/comparison/overall_model_comparison.csv` compare one Gamma filtering-mean trajectory per task with the fitted constant value repeated over the same 70 observation times. Endpoint truth is aligned to the transmission rate that drove the final Euler substep, so the week-5 target is `B = 4`:

- mean RSS: Gamma-noise model 20.055; constant-B 100.528;
- mean RMSE: Gamma-noise model 0.522; constant-B 1.186;
- mean AOB: Gamma-noise model 0.160; constant-B 0.579;
- mean signed error through week 5: Gamma-noise model -0.087; constant-B -0.424;
- mean signed error after week 5: Gamma-noise model 0.086; constant-B 1.576;
- the Gamma-noise model has lower RSS and lower RMSE on all 200 paired tasks;
- the Gamma-noise model has lower AOB on 180 of 200 paired tasks (90.0%).

The active selected-task figures distinguish observation-time particle filtering means from unconditional stochastic forward simulations at fixed fitted parameters. Figures 01 and 08 retain the older B-only task-1 and task-117 ancestry paths as historical repository artifacts. The independent likelihood comparison remains a descriptive fitting diagnostic rather than a complexity-adjusted model-selection test. Neither MIF2 nor the independent likelihood evaluations were rerun for this update; the seeded final plug-in filters were rerun only to recover `I(t)` means, and the displayed forward paths were generated separately from prespecified seeds.

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

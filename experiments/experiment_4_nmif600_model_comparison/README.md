# Experiment 4: two-model computational foundation at `Nmif = 600`

Experiment 4 is the primary computational foundation for the final study. It provides the 200 accepted data sets and the `Nmif = 600` Gamma-noise and constant-B fits. Experiment 5 fits a deterministic six-basis cubic B-spline to the same data and is the final paired three-model manuscript analysis. Experiment 4 therefore remains foundational rather than the final model comparison. Experiments 1–3 remain developmental or supporting studies.

> **Analysis correction in the post-v1.0.0 candidate.** Experiment 4 recovery
> metrics were recomputed from per-task particle filtering means rather than
> ancestry-sampled latent trajectories, and reporting truth at week 5 was
> normalized to `B(5)=4`. The simulated observations, fitted parameters, and
> independent likelihood evaluations were not changed or rerun. The v1.0.0
> release remains archived for provenance.

## Research question

This experiment asks which of two fitted partially observed Markov process (POMP) models more accurately recovers the prescribed transmission-rate path `B(t)`:

1. a Gamma-noise model with latent time-varying `B(t)`;
2. a deliberately restricted constant-B comparator that cannot represent temporal variation in `B(t)`.

The B-spline model is intentionally not fitted in Experiment 4; it is added in Experiment 5 as the third fitted model.

Both models are fitted by iterated filtering (MIF2) with `Nmif = 600` on **the same 200 accepted simulated epidemic data sets**. Acceptance requires `max(H) > 20`, so the reported performance is conditional on informative accepted outbreaks rather than unconditional over all attempted simulations. For each data set, the current workflow reconstructs the best Gamma fit with the observation-time particle filtering mean. Numerical comparisons use residual sum of squares (RSS), root mean squared error (RMSE), signed mean error, and absolute overall bias (AOB) calculated from that filtering mean and the repeated fitted constant-B estimate. Independent particle-filter log likelihoods are retained as descriptive fitting diagnostics, not as complexity-adjusted model-selection criteria.

One ancestry-preserving sampled Gamma trajectory is retained only for each selected-task illustration. It is not used in the current two-model metrics or in the final Experiment 5 three-model analysis. The older `v1.0.0` Experiment 4 tables and figures that used sampled trajectories remain available in the `v1.0.0` tag and release as historical provenance.

## Fixed experiment settings

The shared stochastic SIR process retains the step implementation:

- `B(t) = 4` for process times `t < 5`;
- `B(t) = 2` for process times `t >= 5`;
- `mu_IR = 3`, `N = 10000`, `rho = 0.5`, `k = 10`;
- 10 weeks, 70 daily observation times, Euler step `1/30` week;
- acceptance rule `max(H) > 20`;
- simulation seeds 1001 through 1200.

For the discrete observation-time recovery tables, current post-processing uses the endpoint-aligned reporting truth `B(5) = 4`: `B(t) = 4` for `t <= 5` and `B(t) = 2` for `t > 5`. When historical raw outputs are recombined, `code/04_combine_results.R` accepts the documented legacy `B(5) = 2` label, normalizes only the combined truth label, and recomputes the recovery metrics. It does not alter the simulated observations, fitted parameters, or fitted model objects.

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

Large array ranges are supplied explicitly by the guarded submission wrappers,
not embedded in the worker headers. A directly submitted worker therefore
fails unless an explicit `--array` range is provided.

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
bash hpc/submit_pilot.sh --dry-run
EXP4_CONFIRM_PILOT_SUBMIT=YES bash hpc/submit_pilot.sh
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
bash hpc/submit_all.sh --dry-run
EXP4_CONFIRM_FULL_SUBMIT=YES bash hpc/submit_all.sh
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
download links. The planned public one-download artifact is the repository-level
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
| `01_selected_task_B_trajectory_comparison.pdf` | Prescribed truth, task-1 example-only Gamma sampled trajectory, and fitted constant-B | Selected-task illustration; sampled path is not a metric input |
| `02_RSS_distributions.pdf` | Across-task distributions of Gamma filtering-mean RSS and constant repeated-static-estimate RSS | Aggregate recovery comparison |
| `03_mean_error_distributions.pdf` | Overall, through-week-5, and after-week-5 signed mean errors | Direction and timing of recovery error |
| `04_paired_RMSE_scatter.pdf` | Paired filtering-mean and constant-B RMSE values | Direct within-task model comparison |
| `05_RMSE_distributions.pdf` | Across-task filtering-mean and constant-B RMSE distributions | Aggregate error magnitude |
| `06_independent_loglik_difference.pdf` | Gamma-noise minus constant-B independent log likelihood | Descriptive fit diagnostic only |
| `08_task117_B_trajectory_comparison.pdf` | Prescribed truth, task-117 example-only Gamma sampled trajectory, and fitted constant-B | Second selected-task illustration; sampled path is not a metric input |

Model-specific convergence diagnostic PDFs are stored in `figures/convergence/`.

The figures use a restrained publication style with a white background, no decorative grid, and color-plus-line-type encoding that remains interpretable in grayscale. Truth is shown as a solid black line, the Gamma-noise result as a blue dashed line, and the constant-B result as an orange dot-dashed line. A light-gray vertical reference marks the week-5 change point.

The Gamma curves in Figures 01 and 08 are example-only ancestry-preserving sampled paths from the final plug-in particle-filter approximation. They are not filtering means, exact posterior draws, uncertainty intervals, across-task averages, or inputs to the current recovery metrics. The aggregate Figures 02–05 use the Gamma observation-time filtering means. Constant curves are fitted static estimates repeated over time.

Experiment 4 remains a two-model foundation. The final three-model manuscript figures are produced in Experiment 5. Figure 06 should be presented only with its stated likelihood-interpretation boundary.

## Post-processing and trajectory reconstruction

The current full and pilot post-processing wrappers perform these steps in
order:

1. combine and validate Gamma-noise and constant-B task outputs;
2. run `code/10_regenerate_filtering_means.R` to reconstruct the canonical
   Gamma particle filtering means without rerunning MIF2 or the independent
   likelihood evaluations;
3. use `code/09_regenerate_sampled_B_trajectories.R` only for the task-1 and
   task-117 example trajectories;
4. rebuild the two-model tables, aggregate figures, and convergence
   diagnostics.

Run the equivalent commands from the Experiment 4 directory:

```bash
Rscript code/04_combine_results.R \
  gamma results_raw/gamma results/combined/gamma 1:200

Rscript code/04_combine_results.R \
  constant results_raw/constant results/combined/constant 1:200

Rscript code/10_regenerate_filtering_means.R \
  results/combined/gamma/combined_B_paths.csv \
  results/combined/gamma/filtering_mean_provenance.csv \
  1 1:200

Rscript code/09_regenerate_sampled_B_trajectories.R \
  results/selected_trajectory/sampled_gamma_B_paths.csv \
  results/selected_trajectory/sampled_gamma_B_path_provenance.csv \
  1 1,117

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

`code/09_regenerate_sampled_B_trajectories.R` refuses to overwrite the
canonical filtering-mean file. `code/05_compare_models.R` likewise refuses to
run unless the Gamma paths are labelled `particle_filtering_mean` and the
combined truth encodes `B(5) = 4`.

The accepted shared data and saved fit records required for this post-processing
are supplied by the full-replication ZIP or a new complete run. A fresh clone
supports the final Experiment 5 figure-only workflow documented at the
repository root; it does not contain the raw inputs needed to reconstruct all
200 filtering means.

## Completed-run validation and headline results

The packaged run passed the model-specific combination checks for all 200 tasks. Both models have 200 of 200 tasks present, no missing tasks, no combination problems, unique task IDs, `Nmif = 600`, and successful best-fit status. The paired data check confirms matching simulation seeds and identical observed-data MD5 checksums for the two models on every task.

The main recovery summaries in `results/comparison/overall_model_comparison.csv` compare one Gamma particle filtering mean per task with the fitted constant value repeated over the same 70 observation times:

- mean RSS: Gamma-noise 20.055; constant-B 100.528;
- mean RMSE: Gamma-noise 0.522; constant-B 1.186;
- mean AOB: Gamma-noise 0.160; constant-B 0.579;
- mean error through week 5: Gamma-noise -0.087; constant-B -0.424;
- mean error after week 5: Gamma-noise 0.086; constant-B 1.576;
- Gamma-noise has lower RSS and RMSE on all 200 paired tasks and lower AOB on 90%.

Week 5 belongs to the `through week 5` period and has true `B = 4`; the `after week 5` period uses only `week > 5`. Figures 01 and 08 instead use sampled trajectories solely as selected-task illustrations. The independent likelihood comparison remains a descriptive fitting diagnostic rather than a complexity-adjusted model-selection test. Regenerating filtering means reruns only the final particle filters; it does not rerun MIF2 or the independent likelihood evaluations.

## Rerunning an incomplete or failed task

A complete task is skipped. An incomplete task directory causes the script to stop rather than overwrite it. Inspect the log first and move the failed task directory to a timestamped quarantine location before resubmitting that array index. Do not delete results during diagnosis. For example:

```bash
mv results_raw/gamma/task_149 results_raw/gamma/task_149.incomplete.$(date +%Y%m%dT%H%M%S)
sbatch --array=149 hpc/02_gamma_array.sh
```

After repairs, rerun post-processing:

```bash
sbatch hpc/04_postprocess.sh
```

Do not delete or regenerate `shared_data/task_149` unless the shared data themselves failed validation. Both models must continue to use the same saved data set.

## Interpretation boundary

The experiment evaluates recovery only under the stated simulation model, acceptance rule, starting grids, particle counts, and `Nmif = 600`. Its performance estimand is conditional on accepted outbreaks satisfying `max(H) > 20`. Stable traces in diagnostic tasks 1, 50, 100, 150, and 200 support the computational choice for those tasks but do not prove convergence for all 200 fitted data sets. Independent particle-filter evaluations, multi-start results, and convergence diagnostics should therefore be interpreted together.

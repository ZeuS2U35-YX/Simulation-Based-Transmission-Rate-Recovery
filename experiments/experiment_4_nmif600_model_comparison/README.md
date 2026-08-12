# Experiment 4: Gamma-noise versus constant-B comparison at Nmif = 600

Experiment 4 is the canonical final computational analysis in this repository and the primary source of quantitative evidence for the report, which is currently in progress. Experiments 1-3 are developmental or supporting studies. In particular, Experiment 3 is an earlier `Nmif = 100` recovery-accuracy study and is superseded by Experiment 4 for final Gamma-noise model numerical claims.

## Research question

This experiment compares the ability of two fitted POMP models to recover the prescribed transmission-rate trajectory `B(t)`:

1. a Gamma-noise model with latent time-varying `B(t)`;
2. a deliberately restricted constant-B comparator that cannot represent temporal variation in `B(t)`.

Both models use `Nmif = 600` and are fitted to **exactly the same 200 accepted simulated epidemic data sets**. Because acceptance requires `max(H) > 20`, all reported performance is conditional on informative accepted outbreaks satisfying this rule; it is not unconditional performance over all simulated trajectories. The numerical comparison uses observation-time errors of the Gamma filtering mean and the repeated fitted static constant-B estimate. Independent particle-filter log likelihoods are also retained as descriptive fitting diagnostics; they are not, by themselves, a complexity-adjusted model-selection criterion.

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
results/selected_trajectory/    Task-117 trajectory and provenance
figures/comparison/             Final model-comparison PDFs
figures/convergence/            Nmif=600 convergence diagnostics
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

The completed five-task pilot showed that these limits are very conservative: Gamma-noise model tasks took about 1 hour 9-11 minutes and constant-B tasks about 33-36 minutes, with peak resident memory around 0.2 GB. The resource limits are therefore safety ceilings, not expected runtimes.

## Recommended workflow

Run every command from the Experiment 4 root directory.

### 1. Upload and enter the folder

```bash
cd experiment_4_nmif600_model_comparison
```

### 2. Submit five full-setting pilot tasks

The pilot uses selected diagnostic tasks 1, 50, 100, 150, and 200. They are not described as a random sample. The pilot runs the real `Nmif=600`, particle numbers, starts, and likelihood evaluations. The outputs are production-quality and will be reused by the full array.

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

Main figures:

- `01_selected_task_B_trajectory_comparison.pdf`: prescribed truth, one prespecified task-117 ancestry-preserving Gamma trajectory, and task 117's single fitted constant-B estimate repeated over time; no series is averaged across tasks.
- `02_RSS_distributions.pdf`: comparison of observation-time point-estimator RSS distributions using common histogram bins.
- `03_bias_distributions.pdf`: a 2 × 3 panel summarizing overall, pre-switch, and post-switch mean estimation error for both models; zero is unbiased and the dotted line marks the sample mean within each panel.
- `04_paired_RMSE_scatter.pdf`: each point is one shared simulated data set, comparing Gamma filtering-mean RMSE with repeated-static-estimate RMSE; points below the diagonal favor the Gamma-noise estimator.
- `05_RMSE_distributions.pdf`: observation-time point-estimator RMSE distributions for the two models.
- `06_independent_loglik_difference.pdf`: Gamma-noise minus constant-B independent log likelihood; positive values favor the Gamma-noise model in raw likelihood, but this is not a complexity-adjusted model-selection criterion.
- model-specific convergence diagnostic PDFs in `figures/convergence/`.

The figure style follows a deliberately restrained applied-mathematics convention: black/gray comparison plots, minimal ornament, no decorative grid, and simple line-type or panel distinctions. Figure 01 uses modest color to distinguish the three scientific objects. Its Gamma curve is a finite-particle, plug-in-parameter approximation to a smoothing trajectory conditioned on task 117's complete observation series; it is not an exact posterior draw, a filtering mean, or an across-task average. The constant curve is a fitted static parameter, not a latent trajectory.

For the report, Figures 01-05 are the most directly interpretable recovery summaries. Figure 06 is a descriptive likelihood diagnostic and should be accompanied by interpretation rather than shown without explanation.

## Canonical lightweight figure regeneration

The selected trajectory is generated separately by `code/07_generate_selected_B_trajectory.R`, which runs one final 50,000-particle filter at the already fitted task-117 parameters. The canonical lightweight figure generators are `code/05_compare_models.R` for the six comparison PDFs (reading the saved selected trajectory for Figure 01) and `code/06_make_convergence_diagnostics.R` for the convergence PDFs and tail summary. Run from the Experiment 4 directory:

```bash
Rscript code/07_generate_selected_B_trajectory.R

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

`regenerate_exp4_figures.py` is retained as an optional alternative generator. It reads the same selected task-117 trajectory for the primary figure and cannot generate the retired across-task curve. Pilot regeneration deliberately omits Figure 01.


## Completed-run validation and headline results

The packaged completed run passed the model-specific combination checks for all 200 tasks: both models have 200/200 tasks present, no missing tasks, no combination problems, unique task IDs, `Nmif = 600`, and successful best-fit status. The paired data check verifies the same simulation seed and identical observed-data MD5 checksum for the two models on every task.

The main recovery summaries from `results/comparison/overall_model_comparison.csv` compare per-task Gamma filtering means with per-task fitted constant values repeated over the 70 observation times:

- mean RSS: Gamma-noise model 24.09; constant-B 102.83;
- mean RMSE: Gamma-noise model 0.575; constant-B 1.199;
- mean absolute overall bias: Gamma-noise model 0.160; constant-B 0.606;
- mean after-switch bias: Gamma-noise model 0.139; constant-B 1.576;
- the Gamma-noise model has lower RSS and lower RMSE on all 200 paired tasks;
- the Gamma-noise model has lower absolute overall bias on 90.5% of paired tasks.

These summaries are unchanged by the illustrative task-117 trajectory, which is not used to calculate any metric. The independent likelihood comparison is retained as a descriptive fitting diagnostic and should not be interpreted as a complexity-adjusted model-selection test. Neither MIF2 nor the independent likelihood evaluations were rerun for the trajectory-figure repair.

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

The experiment tests recovery under the stated simulation model, acceptance rule, starting grids, particle counts, and `Nmif=600`. Its performance estimand is conditional on informative accepted outbreaks satisfying `max(H) > 20`. Stable traces in selected diagnostic tasks 1, 50, 100, 150, and 200 provide empirical support for the computational choice, but do not prove convergence for all 200 fitted data sets. The independent pfilter evaluations and multi-start results should be reviewed together with the convergence diagnostics.

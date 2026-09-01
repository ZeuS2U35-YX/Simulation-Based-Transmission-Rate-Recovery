# Experiment 5: B-spline recovery and final paired three-model comparison

## Design invariant

Experiment 5 is one B-spline simulation-analysis batch with exactly 200 paired
simulation replicates. It does **not** simulate a new set of observations.
Instead, it reads the already accepted Experiment 4 data directly from:

```text
../experiment_4_nmif600_model_comparison/shared_data/task_001
...
../experiment_4_nmif600_model_comparison/shared_data/task_200
```

The resulting design is:

```text
200 accepted Experiment 4 observed data sets
  -> 200 Gamma-noise particle-filtering means
  -> 200 selected deterministic cubic B-spline trajectories
  -> 200 repeated constant-B estimates
  -> exactly 200 paired three-model comparison rows
```

The comparison contains exactly three fitted models:

| Fitted model | Representation | Path used in final comparison |
| --- | --- | --- |
| Gamma-noise | Positive latent stochastic transmission process | Observation-time particle filtering mean |
| B-spline | Deterministic non-periodic cubic B-spline for `log B(t)` with six basis coefficients | Selected deterministic trajectory |
| Constant-B | One fitted scalar for the complete epidemic | Selected scalar repeated over all observation times |

The six spline coefficients form one B-spline model; they are not six separate
models. The stochastic SIR simulator is the shared data-generating process and
is not a fourth fitted model.

There is no second or third group of 200 simulations. If `data/`, `results/debug`, `results/formal`, or historical figure directories are present in an older or local workspace, they are legacy single-data-set artifacts and are not inputs to the paired batch.
`code/01_generate_observed_data.R` is deliberately disabled to prevent an
accidental second data-generating workflow.

## Unit of replication and starting values

`task_id` is the simulation-replicate identifier. For every task, the B-spline
model is fitted to the same `observed_data.csv` used by the existing
Gamma-noise result with that `task_id`.

The configured ten MIF2 starting values are internal optimization attempts for
one task. They are not simulation replicates and do not create additional data
sets. Each start has independent MIF2 and particle-filter evaluation seeds.
Every start must finish all configured independent likelihood evaluations to be
eligible. The final selected task fit is the candidate with the largest
independently evaluated log likelihood.

Each completed `results/bspline/task_###/` directory therefore contains:

- `start_selection_audit.csv`: internal seeds, statuses, independent
  likelihoods, and the selected flag; alternate fitted coefficient vectors
  are not retained;
- `likelihood_evaluations.csv`: repeated independent PF evaluations;
- `best_fit_summary.csv`: exactly one final selected fit;
- `B_path.csv`: the selected deterministic B-spline trajectory at 70
  observation times;
- `best_mif2.rds`: the only retained fitted MIF2 object;
- `run_config.csv` and `COMPLETE`.

Internal audit rows are never treated as simulation replicates or as final
model-comparison rows. Only the selected coefficient vector and selected MIF2
object are retained as fitted results.

## Fixed computation

The paired batch uses:

- `n_tasks = 200`;
- six non-periodic cubic B-spline coefficients;
- `Nmif = 600`;
- `Np_mif = 5,000`;
- ten internal starting values per task;
- five independent likelihood evaluations per start;
- `Np_eval = 50,000`.

The fitted model is

\[
\log B(t)=\sum_{j=1}^{6} b_j\xi_j(t), \qquad
B(t)=\exp\left\{\sum_{j=1}^{6}b_j\xi_j(t)\right\}.
\]

## Hard pairing checks

Before fitting, `code/01_validate_paired_inputs.R` verifies that:

- the source has exactly `task_001` through `task_200`;
- every source task is accepted and complete;
- each observed-data and simulated-data MD5 matches Experiment 4's checksum
  file;
- each `simulated_data.csv` has a finite `H` column;
- `acceptance_threshold` is exactly 20 and `accepted` is exactly `TRUE`;
- the independently recomputed `max(H)` is greater than 20 and exactly equals
  `simulation_metadata.csv`'s `max_H`;
- the existing Gamma table has exactly one successful row per task;
- `task_id`, simulation seed, and observed-data MD5 agree between the
  accepted data and Gamma result.

The resulting `results/paired_input_manifest.csv` contains exactly 200 rows
and records both `acceptance_threshold` and the independently calculated
`recomputed_max_H` for every paired input.

After fitting, `code/03_combine_results.R` refuses to combine unless every task
has exactly one correctly selected final fit. `code/04_compare_with_gamma.R`
then repeats the pairing checks and writes:

```text
results/comparison/paired_gamma_bspline_comparison.csv
```

That generated legacy two-model table must contain exactly 200 rows. It is not the tracked final reporting table; the final three-model outputs are under `results/comparison_three_models/`.

The paired recovery metrics use the Gamma model's particle filtering mean,
labelled `path_semantics=particle_filtering_mean`; ancestry-preserving sampled
Gamma trajectories are example-figure artifacts only. The truth convention is
`B(t)=4` through week 5, including `week == 5`, and `B(t)=2` after week 5.
Accordingly, period summaries use `week <= 5` for `through week 5` and
`week > 5` for `after week 5`. Historical raw B-spline task files that contain
the former single-point `B(5)=2` label are accepted only if every other truth
value is consistent; `code/03_combine_results.R` normalizes the combined truth
column and recomputes `B_rmse` without changing any fitted coefficients or MIF2
objects.

## Fresh-clone and full-bundle boundary

A fresh Git clone contains the tracked paired-input manifest, validated combined
B-spline tables, final three-model summaries, figure source data, and plotting
environment. It does not contain the Experiment 4 `shared_data/` tree,
`results/bspline/`, `results_raw/`, or `results_pilot/`. Consequently,
the final figures can be reproduced from a fresh clone, but
`code/01_validate_paired_inputs.R`, `code/02_fit_bspline_B.R`, and
`code/03_combine_results.R` require the full replication bundle or locally
generated raw inputs.

## Run order

Validation only (no simulation and no model fitting):

```bash
Rscript code/01_validate_paired_inputs.R
```

### Independent task 1 production pilot

`hpc/03_task1_production_pilot.sh` is a single-task pilot using the production
settings `Nmif=600`, `Np_mif=5000`, `n_start=10`, `Np_eval=50000`, and
`n_pf_evals=5`. It writes only to `results_pilot/bspline/task_001`; it does not
write to `results/bspline` and is not called by `submit_all.sh`.

After synchronizing this experiment to the HPC, submit the pilot manually:

```bash
mkdir -p logs/pilot results_pilot/bspline
sbatch hpc/03_task1_production_pilot.sh
```

Creating or validating the pilot files does not submit the job automatically.

### Validate and promote the completed task-1 pilot

The formal 199-task array does not include task 1. Before submitting or running
post-processing, promote the completed production pilot into the canonical
result tree. The operation must preserve the pilot and must never overwrite an
existing formal task.

Run these commands from the Experiment 5 directory on the machine that holds
the completed results:

```bash
pilot=results_pilot/bspline/task_001
formal=results/bspline/task_001
stage=results/bspline/.task_001_promote_$$

Rscript code/00_validate_completed_bspline_task.R \
  1 "$pilot" results/paired_input_manifest.csv
test ! -e "$formal"
test ! -e "$stage"
mkdir -p results/bspline
cp -a "$pilot" "$stage"
diff -qr "$pilot" "$stage"
Rscript code/00_validate_completed_bspline_task.R \
  1 "$stage" results/paired_input_manifest.csv
mv "$stage" "$formal"
Rscript code/00_validate_completed_bspline_task.R \
  1 "$formal" results/paired_input_manifest.csv
```

`cp -a` leaves `results_pilot/bspline/task_001` unchanged. The two `test`
guards prevent overwriting either an existing formal result or an abandoned
staging directory. If any check fails, stop and inspect; do not remove or
replace the existing result. `submit_all.sh` independently validates the
promoted formal task before it submits anything.

Fit one task:

```bash
Rscript code/02_fit_bspline_B.R 1
```

Combine all 200 tasks and compare with the existing Gamma fits:

```bash
Rscript code/03_combine_results.R
Rscript code/04_compare_with_gamma.R
```

Submit the single HPC batch:

```bash
bash hpc/submit_all.sh
```

`hpc/submit_all.sh` first requires a valid promoted
`results/bspline/task_001`, then submits only:

1. one B-spline array with task IDs `2-200` (199 tasks);
2. one dependent combination/comparison job.

It never submits a data-generation job or another Gamma-noise fitting job.
The post-processing script reads the canonical `results/bspline` tree and
requires exactly one valid result for every task ID 1 through 200, so the
intended total is one promoted pilot plus 199 array tasks.

The guarded local wrapper is available only for intentional full local runs:

```bash
EXP5_CONFIRM_200_TASKS=YES ./run_all.sh
```

Environment overrides (`EXP5_NMIF`, `EXP5_NP_MIF`, `EXP5_N_START`,
`EXP5_NP_EVAL`, and `EXP5_N_PF_EVALS`) are intended only for isolated smoke
tests. Results produced with overrides cannot pass the canonical 200-task
combination checks unless they equal the configured production values.

## Completed local three-model analysis

For a clean, dependency-locked reproduction of the final comparison figures,
follow [`REPRODUCE_FIGURES.md`](REPRODUCE_FIGURES.md). The accompanying
`renv.lock` records the exact R package versions used for the archived figure
exports; reproducing the figures does not require the raw HPC task directories.

In the author's archival workspace and the full replication bundle, the completed production B-spline task tree is stored under `results_raw/bspline/task_001` through `task_200`. The downloaded raw task files remain unchanged and are not tracked in Git. Full-bundle post-processing uses them explicitly:

```bash
Rscript code/03_combine_results.R \
  results_raw/bspline \
  results/combined/bspline \
  results/paired_input_manifest.csv
Rscript code/05_compare_three_models.R
Rscript code/06_plot_three_model_comparison.R
```

`code/03_combine_results.R` accepts the historical raw `B(5)=2` label only
when all other truth values and provenance checks pass. It writes `B(5)=4` in
the combined analysis table and recomputes every recovery metric without
changing the fitted B-spline coefficients or MIF2 objects.

The three-model analysis pairs the B-spline results with the existing
Experiment 4 Gamma particle-filtering means and repeated constant-B estimates.
It requires exactly 200 matching task IDs, simulation seeds, observed-data
checksums, and 70 observation times per model and task. Primary reporting uses
RMSE, signed mean error, and paired RMSE win counts; RSS and AOB are retained as
secondary diagnostics.

The tracked outputs are:

- `results/combined/bspline/`: validated combined B-spline tables;
- `results/comparison_three_models/`: 200-task paired metrics and summaries;
- `figures/comparison_three_models/`: PDF, SVG, PNG and 600-dpi TIFF figures,
  source-data CSV files, and visual QA notes.

The completed comparison gives mean RMSE values of 0.5218 for Gamma-noise,
0.7128 for B-spline and 1.1858 for Constant-B. B-spline has lower RMSE than
Constant-B in 188 of 200 paired replicates; Gamma-noise has lower RMSE than
B-spline in 162 of 200 paired replicates. These are path-recovery results for
the configured step-change simulation and are not predictive-performance
comparisons.

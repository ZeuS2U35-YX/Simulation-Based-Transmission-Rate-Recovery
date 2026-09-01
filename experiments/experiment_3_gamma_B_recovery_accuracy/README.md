# Experiment 3: Transmission-Rate Recovery Accuracy

## Overview

This earlier `Nmif = 100` supporting study evaluates observation-time recovery of a prescribed time-varying epidemic transmission rate from simulated case-report data. Experiment 4 supersedes it for the final Gamma-noise fits, and Experiment 5 is the final paired Gamma-noise/B-spline/constant-B manuscript analysis.

The experiment contains 200 independently generated, accepted simulation replicates. For each replicate, the same fitting workflow is applied: the model is fitted from nine starting-value combinations, each fitted parameter vector is evaluated using repeated particle filters, the run with the largest independently evaluated log likelihood is selected, and a final particle filter estimates the filtering mean of `B(t)` at each observation time.

Recovery accuracy is summarized using errors of each task's observation-time filtering mean: residual sum of squares (RSS), root mean squared error (RMSE), and mean transmission-rate error. These are point-estimator metrics, not errors of sampled latent trajectories. The gap between the largest and second-largest independently evaluated log likelihoods is retained as a diagnostic of the multi-start search.

## Research question

**Can the Gamma-noise POMP model recover the true transmission-rate trajectory from simulated epidemic data?**

The analysis focuses on three aspects:

- how closely the observation-time filtering mean follows the true `B(t)` values;
- whether recovery shows systematic over- or under-estimation before and after the change in transmission rate;
- how variable recovery accuracy is across independently simulated data sets.

## Study design

The data-generating transmission rate is piecewise constant:

- `B(t) = 4` for `t < 5` weeks;
- `B(t) = 2` for `t >= 5` weeks.

Each data set covers 10 weeks and contains 70 observation times separated by `1/7` week. A simulated trajectory is accepted only if the accumulated number of new infections satisfies `max(H) > 20` for at least one observation interval. The reported recovery results are therefore conditional on this acceptance rule.

The Experiment 3 code and retained tables use the historical point convention `B(t) = 4` for `t < 5` and `B(t) = 2` for `t >= 5`, including `B(5) = 2`. The final Experiment 5 three-model analysis uses an endpoint-aligned reporting convention with `B(5) = 4` and the low period defined by `t > 5`. Experiment 3 outputs remain unchanged for provenance and should not be mixed with the normalized Experiment 5 primary metrics.

The data-generating process is a stochastic SIR model with Euler step `1/30` week, latent states `S`, `I`, `R`, and the incidence accumulator `H`, and initial state `S = 9990`, `I = 10`, `R = 0`, `H = 0`. The fixed values are `mu_IR = 3`, population size `N = 10000`, reporting fraction `rho = 0.5`, and negative-binomial size `k = 10`. Reports follow

\[
Y_n\mid H_n \sim \operatorname{NegBin}(\text{mean}=\rho H_n,\text{size}=k).
\]

The fitted Gamma-noise transmission-rate model uses the same SIR and measurement components but replaces the prescribed piecewise path with a positive latent `B(t)` process. Under this model, the one-step transition distribution for `B(t)` is Gamma with conditional mean equal to its previous value. The fitting model estimates only `B0` and `sigma_beta`; `mu_IR`, `N`, `rho`, `k`, and the initial epidemic state remain fixed at their data-generating values. Both estimated parameters use log transformations.

For each accepted data set:

- MIF2 is run from the nine combinations of `B0 in {2, 4, 6}` and `sigma_beta in {0.10, 0.30, 0.45}`;
- each MIF2 run uses 100 iterations and 5,000 particles;
- MIF2 uses geometric cooling with `cooling.fraction.50 = 0.5` and random-walk standard deviations `B0 = ivp(0.20)` and `sigma_beta = 0.05`;
- each fitted parameter vector is evaluated using five particle filters with 50,000 particles each;
- the repeated likelihood estimates are combined using `pomp::logmeanexp()`;
- the run with the largest independently evaluated log likelihood is selected;
- a final 50,000-particle filter is used to obtain the filtered mean of `B(t)`.

Within a simulation replicate, the same task-specific MIF2 seed and the same sequence of evaluation seeds are reused across the nine starting-value combinations. This keeps the random-number settings fixed across the multi-start comparison but does not remove Monte Carlo error.

The committed `results/paramlist.csv` records all task-level seeds. For task `i`, the simulation seed is `1000 + i`, the shared MIF2 seed is `20260728 + 1000000*i`, the evaluation seed base is `20260800 + 1000000*i` (with repetitions using base plus 1 through 5), and the final-filter seed is `20260900 + 1000000*i`.

## Workflow

Each of the 200 simulation replicates follows the same sequence:

1. generate one accepted simulated epidemic data set;
2. run MIF2 from nine starting-value combinations;
3. evaluate every fitted parameter vector using repeated particle filters;
4. select the run with the largest independently evaluated log likelihood;
5. run a final particle filter at the selected parameter vector;
6. calculate observation-time filtering-mean errors and summary metrics.

## Directory structure

```text
code/                  Simulation, fitting, combination, analysis, and reconstruction scripts
figures/               Recovery figures and the reconstructed MIF2 diagnostic
results/
    paramlist.csv      Task-specific random seeds
    combined/          Combined filtering-mean and fit outputs from 200 replicates
    selected_trajectory/ Prespecified task-145 trajectory and provenance
    recreated_mif2/    Metadata and regenerated data for one reconstructed run
README.md              Experiment documentation
```

### `code/`

- `01_create_paramlist.R` creates the task-specific seed table.
- `02_run_hpc_task.R` runs one complete simulation replicate and writes task-level outputs.
- `03_combine_results.R` validates and combines the 200 task-level result directories.
- `04_analyze_results.R` regenerates the filtering-mean error summaries, the selected-task primary figure, and aggregate metric figures from stored CSV files. It does **not** rerun MIF2.
- `05_recreate_global_best_mif2.R` reconstructs one selected MIF2 run for diagnostic purposes.
- `06_generate_selected_B_trajectory.R` runs one final task-145 particle filter with ancestry tracking and writes the selected trajectory, provenance, and primary figure. It does not rerun MIF2.
- `run_simulation_array.sh` submits the 200-replicate Slurm array.

The canonical summary-figure generator is `code/04_analyze_results.R`. It reads the selected trajectory artifact for Figure 01 and the combined filtering-mean results for the metric distributions. It never averages combined paths to construct the primary trajectory figure. The optional Matplotlib generator follows the same object distinction. The separate `code/05_recreate_global_best_mif2.R` script is the source of the two-page reconstructed MIF2 diagnostic only.

## Main outputs

### Combined numerical results

`results/combined/` contains:

- `combined_mif2_results.csv`;
- `combined_best_fit_summary.csv`;
- `combined_filtered_B_paths.csv`;
- `combined_simulated_data.csv`;
- `task_level_error_summary.csv`;
- `overall_summary.csv`;
- `likelihood_gaps_by_task.csv`.

### Figures

`figures/` contains:

- `01_selected_task_B_trajectory.pdf`: prescribed truth and one prespecified task-145 ancestry-preserving particle-filter trajectory, including `t0`;
- `rss_distribution.pdf`: distribution of observation-time filtering-mean RSS across 200 tasks;
- `rmse_distribution.pdf`: distribution of observation-time filtering-mean RMSE;
- `bias_before_after.pdf`: aligned distributions of per-task filtering-mean signed error overall, before the switch, and from week 5;
- `likelihood_gap_distribution.pdf`: distribution of the gap between the best and second-best independently evaluated log likelihoods;
- `mif2_diagnostic_task_149_run_09.pdf`: reconstructed POMP/MIF2 diagnostic for one selected run.

The reconstructed MIF2 diagnostic is a repository diagnostic, not a summary figure for the 200-replicate recovery study.

## Headline recovery results

The stored 200-replicate results have 200 successful selected fits and 200 successful final particle filters. The following values from `results/combined/overall_summary.csv` summarize errors of per-task filtering means, not the selected task-145 trajectory:

- mean RSS: **24.14** (median 21.77);
- mean RMSE: **0.576** (median 0.558);
- mean overall `B(t)` error: **0.024**;
- mean pre-switch error: **-0.090**;
- mean post-switch error: **0.132**.

The best-versus-second-best evaluated log-likelihood gap has a median of approximately **0.011 log-likelihood units**. This quantity is retained as a multi-start diagnostic; a small gap should not be interpreted as proof that the best start is scientifically distinct from the alternatives, because particle-filter likelihood estimates contain Monte Carlo error.

## Figure design and visual review

The aggregate figures use a deliberately conservative applied-mathematics style:

- no decorative grid or background;
- no large in-figure titles when the caption can carry that information;
- consistent axis labeling and units;
- grayscale distributions with thin black/gray outlines;
- restrained color only where it materially improves line identification;
- reference lines used only when they have a clear inferential meaning.

Figure 01 uses black for the truth and blue for one task-145 ancestry-preserving particle trajectory. This curve is a finite-particle, plug-in-parameter approximation to a smoothing trajectory conditioned on the complete selected observation series; it is not an exact posterior draw, a filtering mean, or an across-task average. It was prespecified by a SHA-256 rule before its values were inspected and is not used in the 200-task metric tables.

Before using a figure in a report or manuscript, perform the following check:

1. render the PDF at approximately its intended printed size;
2. verify that labels and line types remain legible in grayscale;
3. verify that panels use comparable scales when the comparison requires them;
4. remove redundant titles, legends, colors, or annotations;
5. compare the visual density and clarity with published POMP / applied-mathematics figures, especially King, Nguyen, and Ionides (2016), and simplify further if the figure is noticeably more decorative or harder to read.

Reference for the visual standard:

> King, A. A., Nguyen, D., & Ionides, E. L. (2016). *Statistical Inference for Partially Observed Markov Processes via the R Package pomp*. Journal of Statistical Software, 69(12), 1-43. https://doi.org/10.18637/jss.v069.i12

## Fresh-clone and full-bundle boundary

A fresh clone includes `results/combined/` and `results/selected_trajectory/`, but not `Results/task_*` or logs. `code/04_analyze_results.R` works from tracked files; `code/03_combine_results.R` requires the full replication bundle or a new Slurm run.

## Reproducibility

Run all commands from the Experiment 3 root directory.

### Regenerate the prespecified trajectory (only when intentionally required)

The selected task and fitted parameters are already stored, so no MIF2 rerun is required. This command runs one new 50,000-particle final filter:

```bash
Rscript code/06_generate_selected_B_trajectory.R
```

It writes the 71-row task-145 trajectory and provenance under `results/selected_trajectory/` and the primary PDF under `figures/`.

### Regenerate summaries and figures

The stored combined outputs are sufficient; no MIF2 rerun is required.

```bash
Rscript code/04_analyze_results.R
```

This command regenerates `task_level_error_summary.csv`, `likelihood_gaps_by_task.csv`, `overall_summary.csv`, the selected-task primary PDF from its saved trajectory CSV, and the metric-distribution PDFs. It does not modify the original combined CSVs and does not rerun a particle filter.

### Reconstruct the selected MIF2 diagnostic (optional)

Create the seed table if necessary:

```bash
Rscript code/01_create_paramlist.R results/paramlist.csv
```

Then run:

```bash
Rscript code/05_recreate_global_best_mif2.R \
  results/combined/combined_mif2_results.csv \
  results/paramlist.csv \
  results/recreated_mif2 \
  figures
```

The reconstruction uses the stored simulation seed, MIF2 seed, starting values, and original fitting settings for task 149, run 9. Exact bit-for-bit agreement can depend on the R, `pomp`, compiler, and computing environment.

### Rerun the full simulation study on Slurm

A full rerun requires R, `pomp`, and Slurm.

```bash
Rscript code/01_create_paramlist.R results/paramlist.csv
sbatch code/run_simulation_array.sh
```

The array tasks write to `Results/task_001/` through `Results/task_200/`. After all tasks finish:

```bash
Rscript code/03_combine_results.R \
  Results \
  results/combined

Rscript code/04_analyze_results.R
```

## Notes on the reconstructed diagnostic

The original 200-replicate study stored numerical summaries but did not save the selected MIF2 objects. The reconstruction script therefore recreates one run after the fact.

It selects the stored run with the largest evaluated log likelihood among all 1,800 MIF2 runs: task 149, run 9. This choice identifies a run for diagnostic reconstruction only. It is not intended to represent a typical replicate, the smallest-error replicate, or evidence of convergence across all 200 tasks.

The reconstructed object and the selected task-145 illustrative trajectory do not replace or modify the original filtering means, RSS, RMSE, signed errors, or likelihood summaries used in the 200-replicate analysis.

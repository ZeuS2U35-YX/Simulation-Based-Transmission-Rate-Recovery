# Experiment 3/4 B-trajectory audit

> **Historical audit record.** This audit predates the current Experiment 4 sampled-trajectory analysis. Its filtering-mean findings and values, including the earlier RMSE and win-count summaries, are preserved as evidence of the issue that prompted correction; they are not active Experiment 4 results. Use the Experiment 4 README, `results/comparison/`, and `results/selected_trajectory/` for the current estimands, values, and provenance.

## 1. Executive verdict: FAIL

The saved Gamma-noise curves in both experiments are **not** one coherent latent
particle-ancestry trajectory.  The formal Experiment 4 Figure 01, its optional
Python duplicate, and its pilot analogue plot an across-replicate mean of
per-task filtering means.  This conflicts with Felicia's binding requirement
for the illustrative Gamma-noise curve.  Experiment 3's analogous figure has
the same mathematical object.  No defect was found in the stored MIF2 fits,
best-start selection, independent likelihood evaluations, or arithmetic of the
saved filtering-mean metrics.

## 2. Severity counts

| Critical | Major | Minor | Optional |
|---:|---:|---:|
| 1 | 3 | 2 | 1 |

### Audit safety record

- Applicable `AGENTS.md`: none found in the repository or its checked ancestor
  directories.
- `git rev-parse --show-toplevel`:
  `<local-repository-checkout>`
- Initial `git status --short`: `?? overleaf_model_code/` (preserved; not
  inspected as an audit artifact and not modified).
- `git branch --show-current`: `main`.
- `git log -1 --oneline`: `7067d76 Add MIT license`.
- `git remote -v`: `origin https://github.com/ZeuS2U35-YX/2026_summer_epidemic_project.git`
  for both fetch and push.

No destructive Git command, full experiment, Slurm submission, commit, push,
stash, reset, or overwrite was performed.  This report is the sole audit
artifact created in the repository.

## 3. Confirmed mathematical distinction

The installed package is `pomp` 6.4 (`Rscript` query on 2026-08-11; Experiment
4 run configurations also record 6.4).  Its official 6.4 documentation defines
`filter_mean` as the expectation of
\(X(t_k)\mid Y(t_1),\ldots,Y(t_k)\), whereas `filter_traj` requires
`pfilter(..., filter.traj=TRUE)` and extracts one full particle trajectory.
The latter documentation calls the target a smoothing distribution,
\(X(t_k)\mid Y(t_1),\ldots,Y(t_N)\), and says trajectories from independent
PF runs are **likelihood-weighted** samples.  Thus an extracted trajectory is
not an exact posterior draw: it is a finite-particle, plug-in-parameter,
particle-filter approximation to a smoothing trajectory.  It preserves one
particle ancestry; independently selecting a marginal particle at each time
does not.

Primary sources: [pomp 6.4 `filter_traj` documentation](https://kingaa.github.io/manuals/pomp/html/filter_traj.html)
and [pomp 6.4 `filter_mean` documentation](https://kingaa.github.io/manuals/pomp/html/filter_mean.html).

The taxonomy used below is:

1. **Truth (T):** \(B_{true}(t)=4\) for \(t<5\), and 2 for \(t\ge5\).
2. **Static fitted parameter (S):** \(\widehat B_0\),
   \(\widehat\sigma_B\), or \(\widehat B_{constant}\).
3. **Per-task filtering mean (FM):**
   \(\widehat m_r(t_n)\approx E_{\widehat\theta_r}[B(t_n)\mid Y^{(r)}_{1:n}]\).
4. **Across-replicate mean of filtering means (AFM):**
   \(\bar m(t_n)=200^{-1}\sum_r\widehat m_r(t_n)\).
5. **Required coherent trajectory (CT):** one ancestry-preserving `filter_traj`
   result from the final PF at one selected \(\widehat\theta_r\), interpreted
   as the qualified, plug-in smoothing approximation above.
6. **Unconditional simulation (US):** `simulate()` from the fitted model,
   without conditioning on the selected observations.  None of the plotted
   Gamma curves is this object.
7. **Repeated static estimate (RS):** \(\widehat B_{constant,r}\) copied to
   each observation time.

Neither FM nor AFM is CT.  FM uses observations only through time \(t_n\);
CT returned by `filter_traj` is conditioned on the full observation series.
Both condition on a selected fitted parameter vector, not parameter
uncertainty.

## 4. Complete artifact data-lineage table

| Experiment | Artifact / variable | Source and transformation | Object | Current description / status |
|---|---|---|---|---|
| 3 | `simulated_data`, `B_true` | `02_run_hpc_task.R:160-169, 175-273, 306-334, 853-857` | T (data are a stochastic SIR realization under prescribed T) | Correct: truth is prescribed, not a Gamma trajectory. |
| 3 | `B0_hat`, `sigma_beta_hat` | `02_run_hpc_task.R:775-790, 859-875` | S | Correct fitted static parameters. |
| 3 | task `filtered_B_path.csv`; combined `B_filtered_mean` | final `pfilter(filter.mean=TRUE)` then `filter_mean(pf_best)["B",]`, `02_run_hpc_task.R:792-850, 957-999`; row bind/sort, `03_combine_results.R:419-488` | FM | Filename says “filtered”; it is not a CT. |
| 3 | `mean_filtered_B_path.pdf` Gamma line | `aggregate(... B_filtered_mean ... mean)`, `04_analyze_results.R:174-208` | AFM | Correctly says “mean” but fails CT illustration requirement. |
| 3 | `task_level_error_summary.csv`, `overall_summary.csv` | errors of `B_filtered_mean-B_true`, `04_analyze_results.R:88-105, 125-155` | errors of FM; across-task summaries | Numerically valid for FM recovery, not CT error. |
| 3 | likelihood-gap CSV/PDF | rank per-task evaluated MIF2 likelihoods, `04_analyze_results.R:108-123, 271-284` | another: start-selection diagnostic | Unaffected by FM/CT distinction. |
| 3 | reconstructed MIF2 diagnostic | regenerated selected MIF2 object plotted, `05_recreate_global_best_mif2.R:583-640` | another: MIF2 trace | Not a B-state curve; unaffected. |
| 4 | shared data and `B_true` | `01_generate_shared_data_task.R:53-107`; `model_components.R:29-81, 206-212` | T | Correct.  Both models use matching observed-data MD5s. |
| 4 | Gamma `B0_hat`, `sigma_beta_hat`; constant `Beta_hat` | selected finite likelihood row, `02_run_gamma_task.R:252-255`; `03_run_constant_task.R:234-236` | S | Correct static fits. |
| 4 | Gamma `combined_B_paths.csv:B_estimate` | final PF with `filter.mean=TRUE`, extraction, save/combine: `02_run_gamma_task.R:257-281, 330-370`; `04_combine_results.R:68-139` | FM | Misleading generic name `B_estimate`; not CT. |
| 4 | constant `combined_B_paths.csv:B_estimate` | selected `Beta` repeated at 70 times, `03_run_constant_task.R:238-256` | RS | Correct implementation; same column name overloads its meaning. |
| 4 | formal/pilot Figure 01 Gamma line | `aggregate(B_estimate~week,mean)`, `05_compare_models.R:232-268`; pilot routing `hpc/05_postprocess_pilot.sh:16-33` | AFM | Violates CT illustration requirement. |
| 4 | formal/pilot Figure 01 constant line | aggregate of RS, same source | across-replicate mean of RS | Correctly horizontal, but not a same-dataset constant estimate. |
| 4 | metrics/comparison tables | `B_estimate-B_true`, task RSS/RMSE/mean errors, `05_compare_models.R:62-89, 169-205` | errors of FM (Gamma), errors of RS (constant) | Valid under those estimands; not CT error. |
| 4 | likelihood, best-start, convergence | independent PF evaluations/start selection `02_run_gamma_task.R:167-230`; convergence traces/plots `06_make_convergence_diagnostics.R:47-170` | another: likelihood/optimization diagnostics | Unaffected. |

## 5. Findings ordered by severity

### C-01 — Formal Experiment 4 illustration is the wrong object

- **Experiment / evidence:** 4; `code/02_run_gamma_task.R:257-281` creates FM;
  `code/05_compare_models.R:232-268` averages it; optional
  `regenerate_exp4_figures.py:40-67, 135-145` independently repeats the same
  aggregation.  The pilot uses the same R generator through
  `hpc/05_postprocess_pilot.sh:16-33`.
- **Affected output:** `figures/comparison/01_mean_recovered_B_paths.pdf`, its
  copied report figure, and `figures/pilot/01_mean_recovered_B_paths.pdf`.
- **Actual object:** Gamma AFM; constant line is the across-task mean of RS.
- **Currently claimed object:** filename/README call it “mean recovered B
  paths” (`README.md:179-191`); the report accurately calls it an average of
  final-filter means (`Reports/.../sections/05_results.tex:12-25`).  Neither
  is a single latent CT.
- **Why it matters:** an average masks within-data-set stochastic behavior and
  cannot preserve particle ancestry.  It directly fails the supervisor's
  specified illustrative object even though the existing average is honestly
  described in the report.
- **Needed action:** code modification, regenerate this figure and its report
  copy/caption; **no metric recomputation and no MIF2 rerun**.

### M-01 — Experiment 3's corresponding curve is also AFM, not CT

- **Experiment / evidence:** 3; final FM construction is
  `code/02_run_hpc_task.R:792-850`; combination is
  `code/03_combine_results.R:419-488`; aggregation/plotting is
  `code/04_analyze_results.R:174-208`.  The optional Python generator repeats
  it at `regenerate_exp3_figures.py:40-55, 81-100`.
- **Affected output:** `figures/mean_filtered_B_path.pdf`.
- **Actual object:** AFM of 200 FMs.
- **Currently claimed object:** “mean filtered transmission-rate recovery”
  (`README.md:102-109`), which is accurate as far as it goes; it is not a CT.
- **Why it matters:** it cannot be used as the required individual
  Gamma-noise latent trajectory and “trajectory RSS/RMSE” language nearby
  invites conflation with CT.
- **Needed action:** figure/code/caption change only if this supporting figure
  remains an illustrative Gamma-noise curve; no metrics, likelihood, or MIF2
  recomputation.

### M-02 — Recovery metrics are valid filtering-mean/static-estimate metrics, but generic “trajectory recovery” overstates their object

- **Experiment / evidence:** 3 `code/04_analyze_results.R:88-105, 211-268`;
  4 `code/05_compare_models.R:62-89, 169-205, 271-418`.
- **Affected output:** E3 `task_level_error_summary.csv`, `overall_summary.csv`,
  `rss_distribution.pdf`, `rmse_distribution.pdf`, `bias_before_after.pdf`;
  E4 `model_task_metrics.csv`, `paired_model_comparison.csv`,
  `overall_model_comparison.csv`, Figures 02–05 and their pilot duplicates.
- **Actual object:** E3 and E4 Gamma RSS/RMSE/signed mean errors are errors of
  FM against T at 70 observation times.  E4 constant metrics are errors of RS
  against T.  Across-task quantities are means/distributions of those
  task-level values.
- **Currently claimed object:** generic “B(t)”, “trajectory”, “path recovery”
  in E3 `README.md:5-19, 104-107` and E4 `README.md:7-12, 181-186, 217-226`,
  and in report `sections/05_results.tex:29-83`,
  `sections/06_discussion.tex:6-19`.
- **Why it matters:** these results do **not** measure errors of sampled
  individual latent trajectories.  Their numerical values remain appropriate
  under squared-error loss for the stated point estimator; replacing them by
  errors of one sampled CT would change the estimand and inject
  trajectory-sampling variance.
- **Needed action:** wording correction throughout; no recomputation unless
  Felicia changes the desired metric estimand.  Preferred replacement:
  “observation-time filtering-mean recovery metric” for Gamma and
  “observation-time repeated-static-estimate metric” for constant-B.

### M-03 — Report gives an unsupported unweighted particle-mean formula

- **Experiment / evidence:** 4 report source
  `Reports/Summer_Report_revision_01/sections/03_method.tex:104-123` calls the
  result an arithmetic mean and displays \(N_p^{-1}\sum_jB^{(j)}\); actual
  code is `code/02_run_gamma_task.R:257-281`.
- **Affected output:** report method section, Figure 01 caption interpretation,
  and any prose derived from this definition.
- **Actual object:** FM, an approximation to the `pomp` filtering-distribution
  expectation at plug-in \(\widehat\theta\), conditioned through \(Y_{1:n}\).
- **Currently claimed object:** an equally weighted arithmetic mean of 50,000
  filtering particles.
- **Why it matters:** the official `pomp` 6.4 contract identifies the
  filtering distribution and its expectation, not the report's unqualified
  equal-weight formula.  It also obscures the distinction from CT.
- **Needed action:** wording/equation correction only; no figure, metric, or
  MIF2 recomputation.  Replace with the conditional-expectation notation in
  Section 3 above and state that it is estimated by the final PF.

### m-01 — Reuse of `final_pf_seed` is documented for the old FM run, not a new trajectory draw

- **Experiment / evidence:** 3 `code/01_create_paramlist.R:46-61` and
  `02_run_hpc_task.R:792-802`; 4 `code/00_create_paramlist.R:17-30` and
  `02_run_gamma_task.R:257-264`.  E4 persists that seed in
  `combined_run_config.csv` (`02_run_gamma_task.R:349-370`).
- **Affected output:** any replacement CT.
- **Actual object:** existing seed initializes the FM final PF.  No CT is saved.
- **Currently claimed object:** no trajectory seed/selection protocol is
  documented.
- **Why it matters:** using a new, predeclared `trajectory_seed` avoids implying
  that the old FM run selected the displayed trajectory and supports an
  auditable non-appearance-based draw.
- **Needed action:** future code/config documentation and figure generation;
  no MIF2 rerun.

### m-02 — Exact historical Monte Carlo reproduction is not fully pinned

- **Experiment / evidence:** 3 records `pomp_6.4` only in the reconstruction
  session file `results/recreated_mif2/sessionInfo.txt:1-21`; 4 records R and
  `pomp` in `code/io_helpers.R:169-188` and saved run configurations.  The
  report itself notes unrecorded `RNGkind` and no locked environment at
  `sections/D_reproducibility.tex:56-64`.
- **Affected output:** any claim of bit-for-bit reproduction of a CT.
- **Actual object:** a valid new PF can still be run at the saved selected
  parameters and data; exact original random stream is not guaranteed.
- **Currently claimed object:** no incorrect numerical figure claim, but an
  exact-reproduction claim would be unsupported.
- **Why it matters:** it limits reproducibility precision, not reconstruction
  feasibility or the existing results.
- **Needed action:** record versions/RNG kind/trajectory seed with repair;
  no MIF2 rerun.

### O-01 — Pilot figures and Python alternatives must not reintroduce the issue

- **Experiment / evidence:** 4 pilot outputs contain five selected tasks;
  `hpc/05_postprocess_pilot.sh:16-33` routes them through the canonical R
  generator.  Python produces formal figures at
  `regenerate_exp4_figures.py:176-189` when its pilot opt-in is set.
- **Affected output:** eight pilot PDFs and optional Python duplicates.
- **Actual object:** pilot Figure 01 is AFM over five FMs; pilot Figures 02–05
  are FM/RS metrics; pilot convergence is traces.
- **Currently claimed object:** pilot is a selected diagnostic, not a random
  sample (`README.md:92-95, 122-129`).
- **Why it matters:** noncanonical figures can silently recreate the non-CT
  line after the formal repair.
- **Needed action:** update/retire optional generator only if pilot figures are
  published; no MIF2 rerun.

## 6. Experiment 3 impact assessment

The suspected lineage is confirmed exactly.  The final PF uses
`filter.mean=TRUE` (E3 `02_run_hpc_task.R:796-802`), extracts
`filter_mean` (840-850), writes `B_filtered_mean` (957-999), combines it into
`combined_filtered_B_paths.csv` (E3 `03_combine_results.R:419-488`), and
averages it by week (E3 `04_analyze_results.R:174-208`).  RSS, RMSE, and mean
errors use each task's 70 FMs (88-105), not a latent CT.

Figure-by-figure:

| Figure | Canonical / alternate generator; inputs | Actual object and impact |
|---|---|---|
| `mean_filtered_B_path.pdf` | canonical `04_analyze_results.R:174-208`; optional `regenerate_exp3_figures.py:81-100`; combined path CSV columns `week,B_filtered_mean,B_true`; `aggregate(...,mean)` | T plus AFM. Affected (M-01). |
| `rss_distribution.pdf` | canonical lines 211-224; Python 102-112; `task_level_error_summary:RSS` derived from FM | distribution of FM RSS. Valid after relabeling. |
| `rmse_distribution.pdf` | canonical 227-240; Python 114-124; `RMSE` | distribution of FM RMSE. Valid after relabeling. |
| `bias_before_after.pdf` | canonical 243-268; Python 126-149; `bias_*` | distributions of within-task FM signed mean errors. Valid after relabeling. |
| `likelihood_gap_distribution.pdf` | canonical 271-284; Python 151-161; MIF2 evaluated `logLik` | best-minus-second-best likelihood diagnostic. Unaffected. |
| `mif2_diagnostic_task_149_run_09.pdf` | `05_recreate_global_best_mif2.R:583-640`; recreated MIF2 object | two-page MIF2 optimization diagnostic, not B state path. Unaffected. |

The figure filename and README correctly say “mean filtered”, but terms such as
“trajectory RSS/RMSE” should be replaced by “filtering-mean RSS/RMSE at the 70
observation times.”  The `simulate()` calls at `02_run_hpc_task.R:306-312`
generate the accepted *data-generating* epidemic only; they do not generate the
plotted Gamma curve.

## 7. Experiment 4 impact assessment

The suspected lineage is confirmed.  Gamma final PF FM extraction is at
`02_run_gamma_task.R:257-281`, save at 330-370, combine at
`04_combine_results.R:68-139`.  Figure 01 aggregates Gamma and constant
`B_estimate` by week at `05_compare_models.R:238-257`.  The constant code
repeats its fitted `Beta` at `03_run_constant_task.R:252-256`.  Metrics use
per-task `B_estimate-B_true` at `05_compare_models.R:65-89`.

| Figure family | Canonical / alternate generator; inputs | Actual object and impact |
|---|---|---|
| formal `01_mean_recovered_B_paths.pdf` | R `05_compare_models.R:232-268`; Python `regenerate_exp4_figures.py:40-67`; `combined_B_paths.csv:week,B_estimate,B_true`; mean aggregation | T + Gamma AFM + across-task RS. Critical, requires CT replacement. |
| formal `02_RSS_distributions.pdf` | R 271-307; Python 135-145; `paired_model_comparison:RSS_*` | distributions of task FM RSS (Gamma) / RS RSS (constant); relabel only. |
| formal `03_bias_distributions.pdf` | R 309-357; Python 86-108; `bias_*` | distributions of signed within-task FM/RS mean errors; relabel only. |
| formal `04_paired_RMSE_scatter.pdf` | R 359-382; Python 110-120; `RMSE_*` | paired FM versus RS RMSE on shared data; relabel only. |
| formal `05_RMSE_distributions.pdf` | R 384-418; Python 135-145; `RMSE_*` | FM/RS task RMSE distributions; relabel only. |
| formal `06_independent_loglik_difference.pdf` | R 420-441; Python 122-145; `delta_logLik_*` | paired post-fit likelihood diagnostic; unaffected. |
| formal convergence PDFs | R `06_make_convergence_diagnostics.R:100-170`; optional Python 147-180; `combined_mif2_traces.csv` | selected-start MIF2 traces (tasks 1,50,100,150,200); unaffected. |
| eight pilot PDFs | R via `hpc/05_postprocess_pilot.sh:16-33`; Python opt-in 182-189; five-task pilot tables | same objects as formal counterparts but only selected tasks; Figure 01 affected if displayed, convergence unaffected. |

The formal data tables verify the arithmetic: Gamma mean RSS 24.08688857,
mean RMSE 0.57507459; constant 102.83305685 and 1.19892370.  All 200 formal
tasks have 70 rows, unique `(task_id,week)`, and valid T; the code/data checks
and this audit found no row-order or duplication issue.

## 8. Report/documentation impact assessment

Report source **was available and audited** at
`<local-report-root>/`.
It copied Figure 01 unchanged (`sections/D_reproducibility.tex:6-18`) and calls
it an average (`sections/05_results.tex:12-25`), so the report does not falsely
call the old line a single trajectory.  It nevertheless needs a replacement
illustration and caption under the binding requirement.

Locations requiring correction (do not edit in this audit):

| Location | Replacement terminology |
|---|---|
| E3 README `5-19, 104-107, 136` | “per-task observation-time filtering mean” and “across-replicate mean of filtering means”; reserve “trajectory” for CT only. |
| E4 README `7-12, 179-191, 217-226` | Figure 01: “one selected-data-set, ancestry-preserving particle-filter smoothing trajectory at plug-in fitted parameters”; metrics: “filtering-mean/static-estimate recovery metrics.” |
| root README `3, 15-24` | “filtering-mean recovery study/comparison” unless explicitly discussing the new CT display. |
| report method `03_method.tex:104-123` | replace equal-weight formula with \(E_{\hat\theta_d}[B(t_n)\mid Y_{1:n}^{(d)}]\), identified as a PF estimate; separately define CT/full-data smoothing. |
| report results `05_results.tex:10-25, 29-83` and discussion `06_discussion.tex:6-19` | distinguish the new one-dataset CT display from numerical FM/RS metrics; use “observation-time filtering-mean RMSE” rather than “trajectory RMSE.” |
| report reproducibility `D_reproducibility.tex:32-36, 56-64` | map `B_estimate` separately for Gamma (FM) and constant (RS), and add task-selection/trajectory-seed/PF version rules. |

## 9. Unaffected results

- data generation and acceptance rule;
- MIF2 iterations and fitted static parameters;
- five-replicate independent likelihood evaluation and best-start selection;
- all likelihood-gap figures/tables;
- MIF2 convergence traces, pages, and tail summaries;
- the numerical values in existing FM/RS RSS, RMSE, and signed mean-error
  tables, provided their estimands are named correctly.

## 10. Results that remain valid only after relabeling

E3 RSS/RMSE/bias figures/tables and E4 Figures 02–05, their pilot duplicates,
and comparison tables remain numerically and scientifically correct **for**
their point-estimator estimands.  Relabel Gamma quantities as
“per-task filtering-mean errors at observation times” and constant quantities
as “per-task repeated fitted-static-estimate errors.”  Across-task means are
means of those task quantities.  They are numerically correct but
mislabelled/overinterpreted wherever the prose makes them errors of individual
stochastic latent trajectories.

## 11. Results requiring recomputation

No existing RSS, RMSE, bias, likelihood, fitted-parameter, or convergence
number requires recomputation for the stated FM/RS estimands.  The only required
new calculation is the new CT for the illustration, followed by a replacement
Figure 01 (and E3's mean figure only if it is retained as an illustrative
Gamma-noise curve).  Recomputing RSS/RMSE/bias from one CT would be a different
scientific estimand and would add particle-trajectory sampling variation; do so
only if Felicia explicitly requests that new estimand.

## 12. Minimum rerun required

**No full Experiment 3 or Experiment 4 MIF2 rerun is required.**  The minimum
is one new final 50,000-particle PF for one preselected task/model at its saved
selected parameters, with `filter.traj=TRUE`, followed by `filter_traj()`.

Feasibility:

| Requirement | E3 | E4 |
|---|---|---|
| selected task ID | not yet selected; 1–200 IDs are stored | not yet selected; 1–200 IDs are stored |
| observed data | `combined_simulated_data.csv` includes `week,reports` | full 200 `shared_data/task_###/observed_data.csv` are available under the requested Reports directory; MD5 matches saved model outputs |
| selected parameters | `combined_best_fit_summary.csv` (`B0_hat`,`sigma_beta_hat`) | Gamma/constant `combined_best_fit_summary.csv` |
| model/time grid | E3 `02_run_hpc_task.R:160-169, 349-454` | `model_components.R:84-136`, config `15-36` |
| particle count/seeds | `Np_final=50000`; `paramlist.csv` has `final_pf_seed` | `Np_final=50000`; model-specific final seeds and 6.4 are in combined run configs |
| version evidence | `pomp_6.4` in recreated session; current installed version 6.4 | saved run configs record 6.4; current installed version 6.4 |

E3's accepted-data regeneration is also available for task 149, but a future
CT need not regenerate data: it can use the saved observed columns directly.
Because `RNGkind`/full historical environment were not pinned, the replacement
is reproducible under its newly documented environment, not necessarily
bit-identical to the old FM PF.

## 13. Exact proposed repair plan, in dependency order

1. Pre-register a task-selection rule **before viewing any CT** (for example,
   a named seed and `sample(seq_len(200),1)`) and record its result.  Do not
   use task 149 merely because it is already a diagnostic, and do not select
   by curve appearance.
2. For E4 load that one shared observed-data file; verify its MD5 against both
   saved selected-fit rows.  Load the Gamma \(\widehat B_0,
   \widehat\sigma_B\) and the constant \(\widehat B_{constant}\) from the
   same task, plus common T.
3. Reconstruct only the Gamma model at its saved selected parameters and run a
   new PF using `Np=50000, filter.traj=TRUE` and a distinct documented
   `trajectory_seed` (not `final_pf_seed`).  Extract exactly one `B` series
   using `filter_traj(..., vars="B", format="data.frame")`; keep its inherited
   ancestry, do not row-average or resample timewise.
4. Plot the CT at the 70 observation times with T and the same-task RS constant
   line.  State that the CT is a finite-particle, plug-in-parameter smoothing
   approximation conditioned on all observations; do not call it exact posterior
   or one Gamma transition per displayed point.
5. Keep all existing FM/RS metric CSVs unchanged; relabel their axes, captions,
   README, and report text.  Update both R and optional Python generators (or
   clearly retire the alternative) before regenerating official outputs.
6. If E3's curve is shown as an illustration, perform the analogous single-task
   CT operation; otherwise retain it only as a clearly named AFM summary.
7. Copy only verified regenerated figures into the report and preserve the
   selection rule, seed, `pomp` version, R version, `RNGkind`, task ID, selected
   parameters, Np, extraction call, and checksum in a provenance record.

## 14. Validation criteria for the later repair

- `pomp` version is recorded and compatible with its `filter_traj` documentation.
- final PF call visibly contains `filter.traj=TRUE`; extraction calls
  `filter_traj`, not `filter_mean`, `saved_states`, or `simulate`.
- output contains exactly one CT, 70 unique sorted observation times, and no
  averaging across particles/replicates/times.
- selection is reproducible and documented before trajectory inspection; the
  trajectory seed differs from `final_pf_seed`.
- E4 Gamma CT, constant RS, data MD5, task ID, and T all reference the same
  selected data set.
- captions explicitly distinguish full-data smoothing CT from FM
  \(Y_{1:n}\)-conditioned numerical metrics and from US.
- recheck that formal and optional/Python/pilot generators cannot recreate an
  AFM line under a CT label.
- regenerate and visually inspect the changed PDFs; recompute no existing
  metrics unless a new CT-error estimand is formally approved.

## 15. Remaining scientific questions for Felicia

1. Must Experiment 3's supporting AFM figure be replaced, or may it remain as
   a clearly labelled across-replicate filtering-mean summary?
2. Should the formal E4 comparison retain an AFM supplementary panel in
   addition to the required single CT, or should it be retired?
3. Are RSS/RMSE/signed mean error to remain filtering-mean/static-estimate
   recovery metrics (recommended), or is a separate CT-error analysis desired?
4. What pre-registered task-selection seed/rule should govern the illustrative
   task?  The rule must be fixed before curve inspection.
5. Does “trajectory” in narrative conclusions refer only to the display, or is
   it intended to redefine the numerical performance estimand?

## Compact impact matrix

| Component | Impact | Required action |
|---|---|---|
| Data generation | Unaffected | None |
| MIF2 fitting | Unaffected | No rerun |
| Independent likelihood evaluation | Unaffected | None |
| Best-start selection | Unaffected | Reuse saved selected row |
| Final PF output | FM exists; CT absent | One new PF for selected task |
| Trajectory CSVs | E3 FM / E4 Gamma FM, constant RS | Add a separately named CT artifact; do not overwrite |
| RSS/RMSE/bias | Valid FM/RS metrics | Relabel; do not recompute |
| Comparison tables | Valid FM/RS tables | Relabel surrounding prose only |
| E3 figure family | Figure 1 AFM; others FM metrics/diagnostics | Replace/retire Figure 1 if illustrative; relabel others |
| E4 formal/pilot figure family | Figure 01 AFM; Figures 02–05 FM/RS metrics; 06/convergence unaffected | Replace Figure 01; relabel metric figures |
| README | Terminology can conflate objects | Wording correction |
| Report text | Figure copied as AFM and method formula is unsupported | Replace figure/caption; correct formula/terminology |
| Headline scientific conclusion | Supported only for FM/RS recovery metrics, not CT errors | Relabel conclusion; no numeric rerun |

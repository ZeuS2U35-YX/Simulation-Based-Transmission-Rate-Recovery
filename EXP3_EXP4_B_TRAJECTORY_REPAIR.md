# Experiment 3 and 4 B-Trajectory Repair

## Verdict

**PASS.** The Experiment 3 and Experiment 4 primary transmission-rate figures now use one prespecified, ancestry-preserving particle-filter trajectory from one selected fitted dataset. The Experiment 4 figure uses the fitted constant-B estimate from that same dataset. The protected numerical results are byte-for-byte unchanged, the revised report compiles, and the corrected outputs passed programmatic and visual checks.

## Original problem and corrected mathematical objects

The retired primary curves were not individual latent trajectories. Their Gamma lines were across-replicate means of per-task filtering means, and the retired Experiment 4 constant line was an average of fitted static estimates. The repair uses the following terms consistently:

| Abbreviation | Mathematical object | Conditioning and interpretation |
|---|---|---|
| FM | Per-task filtering mean, `m_hat_d(t_n) ~= E_{theta_hat_d}[B(t_n) | Y_1:n]` | A point estimate at an observation time, extracted with `filter_mean()`; not a sampled path. |
| AFM | Across-replicate mean of the task-specific FMs | A task-level aggregate that cannot preserve particle ancestry. |
| CT | One coherent, ancestry-preserving trajectory extracted with `filter_traj()` from `pfilter(..., filter.traj=TRUE)` | A finite-particle, plug-in-parameter particle-filter approximation to a smoothing trajectory conditioned on the complete selected observation series. It is not an exact posterior draw and does not integrate parameter uncertainty. |
| RS | One selected fitted static constant-B estimate repeated over plotted times | A repeated scalar point estimate; not a stochastic latent trajectory. |

The existing RSS, RMSE, and signed-error summaries remain errors of FM for Gamma and RS for constant-B at the 70 observation times. They were not recalculated from CT.

## Scope of computation

- No MIF2 fit, 200-task experiment, likelihood evaluation, start selection, accepted-data simulation, fitted-parameter table, convergence analysis, or task-level metric was rerun.
- One retained 50,000-particle trajectory extraction was completed for each experiment using already selected fitted parameters and saved observations.
- Three particle-filter executions occurred in total. The first Experiment 3 execution completed filtering but the sandbox blocked creation of its output directory, so it retained no artifact. A deterministic replacement Experiment 3 execution beyond the original two-call cap was run only after explicit user authorization. Experiment 4 required one execution. No additional particle-filter calls were made.

## Prespecified selections and input validation

The selection rule was applied to each exact UTF-8 label without a trailing newline: interpret the complete SHA-256 digest as an unsigned big-endian integer `H`, then select `1 + (H mod 200)`. Selection preceded trajectory inspection.

| Check | Experiment 3 | Experiment 4 |
|---|---|---|
| Label | `experiment3-trajectory-2026-08-12` | `experiment4-trajectory-2026-08-12` |
| SHA-256 | `a7d32e608ff0f1c1f20a36f04d68503a949ccda4e39753819f89db11c4f080f8` | `da75fae89d9e42061b02d114495261d348cc0e1527553ad0b9402f7e942eb88c` |
| Selected task | 145 | 117 |
| New trajectory seed | 900000145 | 900000117 |
| Historical FM final-PF seed | 165260900 | 137260900 |
| Simulation seed | 1145 | 1117 |
| Selected run(s) | Gamma run 9 | Gamma run 6; constant-B run 3 |
| Selected Gamma parameters | `B0=4.08137604772541`; `sigma_beta=0.318461288709797` | `B0=3.61026637785311`; `sigma_beta=0.345283889826574` |
| Selected constant parameter | Not applicable | `Beta=3.38780993855149` |
| Selected-data checksum | SHA-256 `d12699591399f1eef7701aba44524f2758690ee1f154c904e0a1f98cacba0df8` for the exact serialized 70-row `week,reports` input | MD5 `64dffb15867fda5ef262e2caf0e46bbf`; SHA-256 `f7767d357c3b871c9b737d5b619cbf1bbe368a4aab3e750c130070b4c632ab75` |
| Same-data validation | Saved task-145 observations used directly | Shared-data file, Gamma record, and constant-B record all match task 117 and the expected MD5 |

The Experiment 3 source combined-data SHA-256 is `ee159ff6d3da6dab04f5acc3f24b16c4cc97226cc260f029adbe4072b35f9976`.

## Runtime and extraction provenance

Both retained extractions used:

- R `4.5.2 (2025-10-31)`;
- `pomp` `6.4`;
- `RNGkind("Mersenne-Twister", "Inversion", "Rejection")`;
- `Np=50000`;
- `set.seed()` immediately before `pfilter()`;
- `pfilter(..., filter.traj=TRUE)` followed by `filter_traj(..., vars="B", format="data.frame")`.

The installed `pomp` 6.4 help/source describes `filter_traj()` under drawing from the smoothing distribution and reconstructs the ancestry of a terminal filter particle. The saved CT is therefore described conservatively as a finite-particle, plug-in smoothing-trajectory approximation conditioned on the full selected series, not an exact posterior draw. The original FM metrics instead condition on observations through each time `t_n`.

## New artifacts

### Experiment 3

- Generator: `experiments/experiment_3_gamma_B_recovery_accuracy/code/06_generate_selected_B_trajectory.R`
- Trajectory CSV: `experiments/experiment_3_gamma_B_recovery_accuracy/results/selected_trajectory/experiment3_task145_B_trajectory.csv`
- CSV SHA-256: `ffb94657433f58f8421462dc45403b977afb653e7a4c019ceb6506a05f67252a`
- Provenance: `experiments/experiment_3_gamma_B_recovery_accuracy/results/selected_trajectory/experiment3_task145_B_trajectory_provenance.txt`
- Corrected primary PDF: `experiments/experiment_3_gamma_B_recovery_accuracy/figures/01_selected_task_B_trajectory.pdf`
- PDF SHA-256: `9ad61182224b2ff061652c32cdd4b1be634ace56c2c3cf979d5f57fb9bb05e67`

The CSV contains one task, one `B_trajectory` series, and 71 unique, sorted, finite state times: `t0` plus 70 observation times. Truth is 4 for `t<5` and 2 for `t>=5`, including the switch at exactly week 5.

### Experiment 4

- Generator: `experiments/experiment_4_nmif600_model_comparison/code/07_generate_selected_B_trajectory.R`
- Comparison CSV: `experiments/experiment_4_nmif600_model_comparison/results/selected_trajectory/experiment4_task117_B_trajectory_comparison.csv`
- CSV SHA-256: `927431361081b49766e94bdb4c53c2a3d3d855f28eb47fc61ff9f2f034a72ad3`
- Provenance: `experiments/experiment_4_nmif600_model_comparison/results/selected_trajectory/experiment4_task117_B_trajectory_provenance.txt`
- Corrected primary PDF: `experiments/experiment_4_nmif600_model_comparison/figures/comparison/01_selected_task_B_trajectory_comparison.pdf`
- PDF SHA-256: `f571dcd5878d95b5c707b1c14499ec79a20e59067fc4446aeb62f2276dd360c6`

The CSV contains one task, one Gamma `B_trajectory`, the prescribed truth, and the same scalar task-117 constant estimate repeated at all 71 unique, sorted, finite state times. It includes `t0` plus 70 observation times. No 200-task or five-task average enters the primary figure.

## Retired aggregate figures and generator safeguards

- Retired Experiment 3 primary AFM figure: `experiments/experiment_3_gamma_B_recovery_accuracy/figures/mean_filtered_B_path.pdf`.
- Retired formal Experiment 4 AFM/mean-static figure: `experiments/experiment_4_nmif600_model_comparison/figures/comparison/01_mean_recovered_B_paths.pdf`.
- Retired pilot five-task aggregate: `experiments/experiment_4_nmif600_model_comparison/figures/pilot/01_mean_recovered_B_paths.pdf`.
- `code/04_analyze_results.R` and `code/05_compare_models.R` now read the selected-trajectory artifacts for their official primary plots.
- The pilot postprocessing script disables primary-trajectory generation because a prespecified formal selected artifact is not a pilot aggregate.
- Both optional Python generators were updated so they cannot silently recreate an AFM curve under the corrected primary label. Their syntax was checked; runtime rendering was not used because the available Python environments lacked `numpy` or `matplotlib` in the same interpreter. Canonical PDFs were generated with R.

## Protected numerical results

A pre-repair SHA-256 manifest covered 27 canonical Experiment 3 and Experiment 4 result CSVs, including combined FM/static paths, fitted-parameter tables, likelihood tables, task metrics, paired comparisons, and overall summaries. After all figure generation, this command returned `OK` for all 27 files:

```text
shasum -a 256 -c /Users/zeus/Documents/Codex/2026-08-11/you-are-conducting-an-independent-evidence/work/repair/immutable_csv_sha256_before.txt
```

Consequently, the saved RSS, RMSE, bias/signed-error, likelihood, fitted-parameter, comparison-table, and convergence inputs remain byte-for-byte unchanged. In particular, the Experiment 4 mean FM RMSE remains `0.575074590343979`, the mean RS RMSE remains `1.19892369846335`, and Gamma has lower paired FM-versus-RS RMSE in 200/200 tasks.

## Repository documentation and plotting changes

The following active documentation/generator families were corrected:

- root `README.md`;
- Experiment 3 `README.md`, `code/04_analyze_results.R`, and `regenerate_exp3_figures.py`;
- Experiment 4 `README.md`, `code/05_compare_models.R`, `hpc/05_postprocess_pilot.sh`, and `regenerate_exp4_figures.py`.

Metric plots now identify their estimands as observation-time filtering-mean errors or repeated-static-estimate errors. Searches of active source and documentation found no current primary reference to the retired filenames, no primary generator averaging combined task paths, and no primary caption presenting an AFM as a trajectory. Historical references in the independent audit are intentionally retained as audit evidence.

## LaTeX report repair

Report root: `/Users/zeus/Desktop/Math-Resources/Research_2026summer/Reports/Summer_Report_revision_01` (not Git-backed).

Updated report material:

- `sections/03_method.tex`
- `sections/05_results.tex`
- `sections/06_discussion.tex`
- `sections/A_developmental_experiments.tex`
- `sections/C_extended_code.tex`
- `sections/D_reproducibility.tex`
- `REVISION_EVIDENCE.md`
- `REVISION_NOTES.md`
- `code/experiment4_source/05_compare_models.R`
- `code/experiment4_source/07_generate_selected_B_trajectory.R`
- refreshed selected comparison/error/likelihood figure copies under `figures/experiment4/`

The old report copy `figures/experiment4/01_mean_recovered_B_paths.pdf` was retired and replaced by `figures/experiment4/01_selected_task_B_trajectory_comparison.pdf`.

Canonical build command:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error report.tex
```

The command exited successfully. The resulting `/Users/zeus/Desktop/Math-Resources/Research_2026summer/Reports/Summer_Report_revision_01/report.pdf` has 32 letter-size pages. The final log contains no LaTeX/package warnings, undefined references/citations, overfull boxes, or underfull boxes.

## Commands used for retained trajectory filters

From each experiment root:

```text
Rscript code/06_generate_selected_B_trajectory.R
Rscript code/07_generate_selected_B_trajectory.R
```

Each generator asserts its complete selection digest, task, seed, saved observations, fitted values, time grid, and particle count before filtering. Output directories are created before the PF, and the resulting artifact is validated before it is saved.

## Programmatic and visual validation

- Experiment 3: task 145; seed 900000145 distinct from 165260900; 70 observations; 71 extracted state times; exactly one finite, sorted, unique Gamma trajectory; correct truth at all times; no `filter_mean()` or aggregation in trajectory construction.
- Experiment 4: task 117; seed 900000117 distinct from 137260900; expected data MD5; matching Gamma/constant task records; 70 observations; 71 state times; exactly one finite Gamma trajectory; constant line exactly `3.38780993855149` throughout; correct truth at all times; no task average.
- Both new R generators contain `filter.traj=TRUE` and call `filter_traj()`.
- R scripts and Python generators passed syntax checks.
- The Experiment 3 primary PDF, Experiment 4 primary PDF, and the affected report pages were rendered and visually inspected. Truth steps at week 5, `t0`, line mappings, axes, legends, captions, and nearby results text were correct, readable, and unclipped.
- Protected result CSVs passed 27/27 checksum comparisons.

## Remaining limitations

- Each CT is illustrative for one prespecified dataset and must not be used as a 200-task performance summary.
- `filter_traj()` provides a finite-particle approximation at a selected plug-in parameter vector. It is not an exact posterior draw and does not propagate fitted-parameter uncertainty.
- Values are recorded at `t0` and observation times; underlying process evolution uses Euler substeps of 1/30 week. Consecutive displayed points are not individual Gamma-process substeps.
- Conditioning and resampling select an ancestry-preserving particle history relative to the complete observed series; the CT is not an unconditional simulation or the trajectory that generated the data.
- The original HPC environment was not locked completely, so bit-for-bit historical reproduction of old particle filters is not claimed. The corrective runtime, RNG kinds, seeds, inputs, parameters, and package version are recorded.
- The initial sandbox-blocked Experiment 3 execution and explicitly authorized replacement are disclosed above; only one retained trajectory artifact exists for that experiment.

## Final conclusion

The repaired primary figures now match the binding scientific requirement. The 200-task performance conclusion remains supported by unchanged FM-versus-RS metrics and is explicitly separated from the two single-dataset CT illustrations.

# Experiment 4 selected-task trajectory-figure QA

## Figure contract

- **Figure 6 conclusion:** observation-time particle filtering means recover the prescribed transmission change and corresponding infectious state in a comparatively well-recovered case (Task 1) and a more challenging case (Task 117).
- **Figure 9 conclusion:** fixing the selected parameters still produces stochastic model realizations; one joint forward simulation shows how `B(t)` and the corresponding `I(t)` evolve without observation conditioning after fitting.
- **Panel map:** panels a/b show `B(t)` for Tasks 1/117; panels c/d show the corresponding `I(t)` curves.
- **Comparison:** black solid = data-generating truth; blue dashed = Gamma-noise model; orange dot-dashed = constant-B model.

## Selection and interpretation safeguards

- Tasks 1 and 117 were already designated as the primary and second illustrations before this figure revision; they were not selected by screening the new plots.
- Task 1 has Gamma filtering-mean RMSE 0.446 (empirical percentile 0.275 across the 200 tasks), while Task 117 has RMSE 0.560 (percentile 0.700).
- The constant-B RMSE percentiles are 0.930 for Task 1 and 0.265 for Task 117, so Task 117 is a deliberately less favorable comparison rather than a showcase selected to maximize the Gamma advantage.
- The illustrations are not claimed to represent the full task distribution; aggregate recovery evidence remains in Figures 4 and 5.

## Data and trajectory integrity

- Each task uses its exact shared observation series, verified against the fitted-record checksum.
- Figure 6 Gamma `B(t)` uses the retained primary-metric filtering means. `I(t)` means come from the rerun 50,000-particle final filters.
- Figure 9 fixes each model at its selected plug-in parameters and uses the first realization from a prespecified seed equal to the saved final-filter seed plus 100,000. Within each task, the Gamma `B(t)` and `I(t)` curves come from the same joint forward simulation.
- No path was screened or selected for agreement with the truth, visual appearance, or recovery performance.
- The trajectories are unconditional simulations after the fitted parameters are fixed. They are not filtering means, fitted paths, posterior draws, or uncertainty bands, and parameter uncertainty is not integrated.
- Week-0 open markers are fitted or fixed initial states rather than filtering summaries.

## Export and visual QA targets

- Backend: R (`ggplot2` + `patchwork`), 183 mm x 150 mm.
- Panels use common `B(t)` and common `I(t)` scales within each figure.
- Vector exports: PDF and SVG; raster exports: 300-dpi PNG and 600-dpi TIFF.
- Required visual checks: shared legend, panel order, task headings, week-5 markers, open week-0 markers, line styles, common row scales, unclipped labels, and readable final-size text.
- Static source preflight: 18 pass, 2 reviewed warnings, 0 fail. The syntax warning requests an R parse check, which passed. The uncertainty warning is not applicable: each filtering mean is a conditional state summary within one data set, and each forward path is explicitly one seeded realization rather than a seed aggregate.
- R source parses successfully. Each source-data table has 718 data rows and the provenance table has two task rows.
- PDF glyph audit: PASS for both figures; minimum detected text size 6.117 pt (required minimum 5 pt).
- TIFF metadata: 600 x 600 dpi, 4322 x 3543 pixels for both figures.

## Source data

- `experiment4_selected_tasks_filtering_mean_trajectories.csv`
- `experiment4_selected_tasks_forward_simulation_trajectories.csv`
- `experiment4_selected_tasks_filtering_forward_provenance.csv`

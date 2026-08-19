# Experiment 4 Task 1 trajectory-figure QA

## Figure contract

- **Figure 6 conclusion:** the Gamma-noise particle filtering mean follows the prescribed change in transmission rate and gives an infectious-state estimate conditioned on reports through each observation time.
- **Figure 9 conclusion:** one forward realization at the fitted parameters can depart substantially from both the truth and the filtering mean; it is an illustration of process variability, not a state estimate or uncertainty interval.
- **Panel map:** panel a is transmission rate \(B(t)\); panel b is infectious population \(I(t)\).
- **Comparison:** black solid = data-generating truth; blue dashed = Gamma-noise model; orange dot-dashed = constant-\(B\) model.

## Data and provenance checks

- Dataset: Experiment 4, Task 1 (`observed_data_md5 = 64a1b5c02bdecc100b37ee56029391d3`).
- Observation grid: 70 times from week 1/7 through week 10; open markers at week 0 denote fitted or fixed initial values and are not filtering means.
- Figure 6 Gamma \(B(t)\): the retained primary-metric filtering-mean output.
- Figure 6 \(I(t)\): final 50,000-particle filters rerun at the fitted parameters with recorded seeds. Because the current `pomp` runtime does not exactly reproduce the retained Monte Carlo output, the retained-vs-rerun Gamma \(B\) maximum absolute difference (0.0119827) and both likelihood versions are recorded in the provenance CSV.
- Figure 9 seeds: 21,360,900 (Gamma-noise) and 22,260,628 (constant-\(B\)), defined before plotting as the final-filter seed plus 100,000. The first realization was retained; there was no visual or performance-based selection.
- Parameter uncertainty is not integrated into either figure.

## Quantitative range checks

| Figure | Quantity | Truth | Gamma-noise | Constant-\(B\) |
|---|---:|---:|---:|---:|
| 6 | \(B(t)\) | 2--4 | 2.072--4.927 | 4.045 |
| 6 | \(I(t)\) | 0--464 | 0.003--566.458 | 0.002--509.954 |
| 9 | \(B(t)\) | 2--4 | 3.365--8.683 | 4.045 |
| 9 | \(I(t)\) | 0--464 | 7--1,072 | 0--10 |

The large Figure 9 discrepancy is retained because it is the prespecified first forward realization and is central to the interpretation boundary.

## Export and visual QA

- Backend: R (`ggplot2` + `patchwork`), 183 mm x 150 mm.
- Vector exports: PDF 1.5 and SVG; raster exports: PNG and 600-dpi TIFF.
- Static preflight: 18 pass, 2 non-blocking warnings, 0 fail. The uncertainty warning is a false positive caused by the terms `seed` and `mean`; no across-seed aggregate is plotted.
- PDF text audit: minimum detected text size 6.97 pt in both figures; no text below 5 pt.
- Visual checks: panel order, legends, week-5 marker, open week-0 markers, line styles, axis labels, and unclipped content verified.

## Source files

- `experiment4_task1_filtering_mean_trajectories.csv`
- `experiment4_task1_forward_simulation_trajectories.csv`
- `experiment4_task1_filtering_forward_provenance.csv`
- `../../code/11_generate_task1_filtering_mean_forward_figures.R`

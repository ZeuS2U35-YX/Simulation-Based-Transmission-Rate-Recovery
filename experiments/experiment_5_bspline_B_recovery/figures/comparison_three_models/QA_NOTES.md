# Three-model figure QA

## Figure contract

- Core conclusion: allowing time variation improves recovery relative to a
  constant transmission rate, while the Gamma-noise model recovers the
  configured step-change path more accurately than the deterministic B-spline
  model in this simulation scenario.
- Archetype: quantitative comparison grid.
- Replicate unit: one accepted simulation task; `n = 200` matched tasks.
- Backend: R (`ggplot2` and `patchwork`) for plotting and every export.
- Final size: 183 mm wide; 190 mm high for the RMSE figure and 155 mm high for
  the mean-error figure.
- Exclusions: none.

## Data and statistics checks

- All three models contain task IDs 1 through 200 exactly once.
- Every model-task path contains 70 ordered observation times.
- Simulation seeds and observed-data MD5 values match across models.
- Gamma paths are particle-filtering means; sampled latent trajectories are not
  used in these figures.
- B-spline paths are deterministic selected B-spline trajectories.
- Truth is `B = 4` through week 5 (including week 5) and `B = 2` afterward.
- Figure 1 uses common RMSE histogram bins and paired RMSE scatterplots with
  equal axes and a `y = x` reference line.
- Figure 2 uses common mean-error bins and a common x-axis across all nine
  panels. Dashed lines mark zero and dotted lines mark panel means.

## Rendered visual checks

- No panel, title, axis label, strip label, caption or data mark is clipped.
- Model names and capitalization are consistent across figures.
- Gamma-noise, B-spline and Constant-B retain the same blue, orange and grey
  identities across panels.
- The plotted source data include all 200 replicates per model and panel.
- PDF and SVG retain editable text; TIFF is exported at 600 dpi and PNG at
  300 dpi.
- The R plotting source passes syntax parsing and the static figure preflight
  with no substantive failure. The PDF scanner reports Cairo's transformed
  one-point text operators; visual inspection and the SVG text sizes confirm
  that displayed labels are readable at the declared final dimensions.

## Interpretation boundary

The figures compare recovery of the latent transmission-rate path in one
controlled step-change simulation design. They do not compare out-of-sample
prediction and do not establish universal superiority of one representation.

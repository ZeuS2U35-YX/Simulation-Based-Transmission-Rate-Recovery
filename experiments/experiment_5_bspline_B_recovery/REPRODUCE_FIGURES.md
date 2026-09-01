# Reproducing the Experiment 5 figures

This guide reproduces the two three-model comparison figures from the tracked,
post-processed source data. It does **not** rerun MIF2, particle filtering, or
any HPC job.

## Requirements

- R 4.5.2 (the version recorded in `renv.lock`)
- Internet access for the first package restore
- Git, or a ZIP download of the repository

The locked plotting environment contains `ggplot2`, `patchwork`, `scales`,
`svglite`, and `ragg`, together with their recursive dependencies.
It is intentionally scoped to `code/06_plot_three_model_comparison.R` through
`.renvignore`; it is not the dependency environment for the MIF2, particle
filtering, or other analysis scripts in this directory.

## 1. Obtain the repository or replication bundle

```bash
git clone https://github.com/ZeuS2U35-YX/Simulation-Based-Transmission-Rate-Recovery.git
cd Simulation-Based-Transmission-Rate-Recovery/experiments/experiment_5_bspline_B_recovery
```

If Git is unavailable, download the repository source ZIP from GitHub,
extract it, and open a terminal in
`experiments/experiment_5_bspline_B_recovery`.

After the planned `v1.1.0` release is published, the custom
`Simulation-Based-Transmission-Rate-Recovery-v1.1.0-full-replication.zip`
asset provides the same tracked plotting inputs plus the accepted shared data
and raw task-level outputs. The custom full-replication ZIP is distinct from
GitHub's automatically generated source-code ZIP. Raw task directories are not
needed for this figure-only workflow.

## 2. Restore the locked R environment

The repository includes the standard `renv` activation files. Run the following
command once from the Experiment 5 directory:

```bash
Rscript -e 'renv::restore(prompt = FALSE)'
```

The lockfile records R 4.5.2 and the exact package versions used to generate
the archived figures. A newer R installation may work, but R 4.5.2 is
recommended when matching the archived rendering as closely as possible.
On first use, the included activation script bootstraps `renv` 1.2.4 and
creates a project-specific package library. It does not replace packages in the
user's global R library.

To confirm that the plotting environment is synchronized with the lockfile,
run:

```bash
Rscript -e 'renv::status()'
```

## 3. Generate the figures

After the restore completes, run:

```bash
Rscript code/06_plot_three_model_comparison.R
```

The script validates the three-model input tables before drawing. It requires
exactly 200 matched simulation replicates for each of the Gamma-noise,
B-spline, and Constant-B models. No plotted value is manually edited or
excluded.

## 4. Check the outputs

The command produces both figures in:

```text
figures/comparison_three_models/
```

Expected files include:

```text
01_three_model_RMSE_comparison.pdf
01_three_model_RMSE_comparison.svg
01_three_model_RMSE_comparison.png
01_three_model_RMSE_comparison.tiff
02_three_model_mean_error_comparison.pdf
02_three_model_mean_error_comparison.svg
02_three_model_mean_error_comparison.png
02_three_model_mean_error_comparison.tiff
02_three_model_mean_error_legend.txt
```

The script also regenerates the figure source-data CSV files under
`figures/comparison_three_models/source_data/`. A successful run ends with:

```text
Generated two R figures in PDF, PNG, SVG, and 600-dpi TIFF formats, with source-data CSV files.
```

For an additional check, run:

```bash
Rscript -e 'x <- read.csv("figures/comparison_three_models/source_data/figure_2_mean_error_histograms_source_data.csv"); stopifnot(nrow(x) == 1800L, length(unique(x$task_id)) == 200L); cat("Reproduction checks passed.\n")'
```

## Reproducibility notes

- The quantitative inputs used by the plotting script are tracked in
  `results/comparison_three_models/`; the raw HPC task directories are not
  required for figure reproduction.
- The figures are generated entirely in R. No Canva, Photoshop, Illustrator,
  or manual graphical adjustment is used.
- PDF and SVG are the vector deliverables; SVG text remains editable. PNG is a
  300-dpi preview, and TIFF is exported at 600 dpi.
- The intended final widths are 183 mm. Figure 1 is 183 x 190 mm and Figure 2
  is 183 x 155 mm.
- Helvetica is requested by the script. If it is unavailable, the operating
  system may substitute a metrically similar sans-serif font, causing small
  text-spacing differences without changing any plotted data.
- On systems for which CRAN does not provide a binary package, installation
  from source may require the usual R compilation tools (Xcode Command Line
  Tools on macOS or Rtools on Windows).

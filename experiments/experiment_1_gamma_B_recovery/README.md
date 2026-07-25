# Experiment 1: Transmission-Rate Recovery Using Particle Filtering

## Goal

This experiment investigates whether a partially observed Markov process (POMP) model can recover a known epidemic transmission-rate path from simulated reported-case data.

The data-generating model uses a specified transmission-rate path, while the fitted model treats the transmission rate as a latent Gamma-noise process. The recovered transmission-rate and infectious-state trajectories are compared with the corresponding simulated trajectories.

## Experimental workflow

1. Specify a known transmission-rate path.
2. Simulate the latent epidemic process and reported cases.
3. Fit a Gamma-transition POMP model to the simulated observations.
4. Use MIF2 from multiple starting values to estimate the unknown parameter or parameters.
5. Evaluate each fitted parameter vector with repeated particle filters.
6. Select the fit with the largest Monte Carlo estimate of the log likelihood.
7. Run a final particle filter using the selected parameter estimates.
8. Compare the filtered transmission-rate and infectious-state paths with the simulated truth.

## Code

### `code/Fit_Gamma_B_to_Constant_B4_Data.R`

Simulates epidemic data with a constant transmission rate, \(B(t)=4\), and fits the Gamma-transition model. Both \(B_0\) and `sigma_beta` are included in the MIF2 search. The script saves the MIF2 convergence diagnostic and the filtered-path comparisons as PDF files.

### `code/MIF2_Search_and_Likelihood_Surface.R`

Runs the original multiple-start MIF2 search for a changing transmission-rate scenario and compares the true and Gamma-filtered transmission-rate and infectious-state paths. Despite the historical filename, the current script does not construct the Part D likelihood-surface figure; that PDF is retained as an existing experiment output.

### `code/MIF2_Search_Sigma_Only_Fixed_B0_Save_Plots.R`

Treats \(B_0=4\) as known and fixed, and estimates only `sigma_beta`. Multiple initial values of `sigma_beta` are used to reduce dependence on a single MIF2 starting point. The final transmission-rate and infectious-state plots are saved as PDF files.

## Output paths

All scripts use paths relative to this experiment folder rather than machine-specific absolute paths. Therefore, after the repository is cloned, figures are written to:

```text
experiment_1_gamma_B_recovery/figures/
```

and numerical or fitted-object outputs can be written to:

```text
experiment_1_gamma_B_recovery/results/
```

The path helper works when a script is run from the repository, from the experiment folder, or directly from the `code/` folder. This avoids committing paths tied to one computer or user account.

## Figure format

All tracked figures in `figures/` are stored as PDF files. PDF is used because these plots are intended for LaTeX reports and GitHub version control, and because plots regenerated directly from R remain vector graphics.

Existing figures that were originally available only as PNG files were converted to PDF. Those converted files retain the original raster content inside a PDF container; rerunning the corresponding R script is preferable when a true vector PDF is required.

## Folder structure

```text
experiment_1_gamma_B_recovery/
├── code/       # R scripts
├── figures/    # PDF figures generated or retained for the experiment
├── notes/      # Experiment notes and unresolved questions
├── results/    # Numerical summaries and fitted objects
├── .gitignore
└── README.md
```

## Reproducibility notes

- Random seeds are specified in the scripts.
- MIF2 fits are evaluated using repeated particle filters.
- The repository does not track `.Rhistory`, `.RData`, `.DS_Store`, or RStudio session-state files.
- Empty `notes/` and `results/` folders are retained with `.gitkeep` files so that the intended project structure is preserved on GitHub.

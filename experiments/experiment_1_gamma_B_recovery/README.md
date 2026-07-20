# Experiment 1: Transmission-Rate Recovery Using Particle Filtering

## Goal
The goal of this expriment is to investigate whether a particle-filtering algorithm can recover can recover a known transmission rate from simulated epidemic reported-case data.

The recovered transimission rate is compared with the true transmission path.


## Method

1. Define a known transmission-rate path.
2. Simulate the hidden epidemic process and generate the corresponding reported cases.
3. Fit the simulated data using a particle-filtering model that assumes a Gamma-noise process for the transmission rate.
4. Use MIF2 to estimate the unknown parameters.
5. Select the parameter estimates with the largest estimated log likelihood.
6. Evaluate the selected parameter estimates using repeated particle filtering.
7. Obtain the estimated transmission-rate trajectory.
8. Compare the estimated transmission-rate trajectory with the true trajectory.


## Experiment Design

The experiment includes both constant and changing transmission-rate scenarios.

For each scenario, the epidemic data are simulated using known parameter values. The Gamma-noise filtering model is then fitted to the simulated
reported cases.

Details to be reconstructed from the original notes and code.

## Files

### `Fit_Gamma_B_to_Constant_B4_Data.R`

This script simulates data using a constant transmission rate of \(B=4\)
and fits the Gamma-noise model to the simulated data.

It runs MIF2 from multiple starting points, selects the run with the
largest estimated log likelihood, and plots the recovered trajectories.

### `MIF2_Search_and_Likelihood_Surface.R`

This script performs the MIF2 parameter search and evaluates the
likelihood surface over selected parameter values.



## Folder Structure

- `code/`: R scripts used in the experiment.
- `figures/`: Figures generated from the analysis. (MIF2 convergence diagnostics were examined for Part C only.)
- `results/`: Numerical summaries and fitted results.
- `notes/`: Detailed experiment notes, decisions, and unresolved questions.



## Status

- The two main R scripts have been checked and run successfully.
- The experiment description is being reconstructed from the original code
  and research notes.
- Figures and numerical results still need to be organized.



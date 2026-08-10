# Experiment 1: Gamma-noise model recovery on simulated epidemic data

## Overview

Experiment 1 is an early developmental simulation study of a Gamma-noise POMP model for recovering a time-varying epidemic transmission rate. It contains three controlled, single-dataset scenarios:

1. constant true transmission rate, with `B0` and `sigma_beta` estimated;
2. piecewise true transmission rate, with `B0 = 4` fixed and `sigma_beta` estimated;
3. the same piecewise dataset, with both `B0` and `sigma_beta` estimated.

The purpose is to illustrate the fitting workflow, inspect starting-value sensitivity, and compare filtered latent trajectories with known simulated truth. Because each scenario uses only one simulated epidemic, Experiment 1 is **not** a repeated-simulation assessment of bias, RMSE, coverage, or robustness. Those questions belong to later experiments.

Experiment 1 is developmental evidence, not a co-equal final analysis. Experiment 4 is the canonical computational analysis and the primary source of quantitative evidence for the report.

## Model

The latent epidemic process contains `S`, `I`, `R`, the incidence accumulator `H`, and - in the fitted model - the latent transmission-rate state `B`.

The initial epidemic state is

\[
S(0)=9990,\qquad I(0)=10,\qquad R(0)=0,\qquad H(0)=0.
\]

Reported cases follow

\[
Y_n\mid H_n \sim \operatorname{NegBin}(\text{mean}=\rho H_n,\text{size}=k).
\]

Fixed parameters are:

- `mu_IR = 3`
- `N = 10000`
- `rho = 0.5`
- `k = 10`

There are 70 observation times over 10 weeks, spaced `1/7` week apart (daily observations with time expressed in weeks). The latent process uses an Euler step of `1/30` week.

All simulations use seed `20260527`.

## Scenarios

### 1A. Constant true transmission rate

Script: `code/01_fit_constant_B4.R`

\[
B(t)=4,\qquad 0\le t\le10.
\]

Estimated parameters:

- `B0`
- `sigma_beta`

Starting grid:

\[
B_0\in\{2,4,6\},\qquad
\sigma_\beta\in\{0.10,0.30,0.45\}.
\]

Selected fit from the stored run:

- `B0_hat = 3.9893`
- `sigma_beta_hat = 0.00253`
- evaluated log likelihood `= -250.6249`
- Monte Carlo SE `= 0.0574`

The near-zero fitted `sigma_beta` is consistent with the constant data-generating transmission rate in this illustrative dataset.

### 1B. Piecewise transmission rate with known `B0`

Script: `code/02_fit_piecewise_B_fixed_B0.R`

\[
B(t)=
\begin{cases}
4,&t<5,\\
2,&t\ge5.
\end{cases}
\]

`B0` is fixed at 4. MIF2 estimates only `sigma_beta`, starting from `0.10`, `0.30`, and `0.45`.

Selected fit:

- `B0 = 4` fixed
- `sigma_beta_hat = 0.1969`
- evaluated log likelihood `= -197.7006`
- Monte Carlo SE `= 0.0499`

### 1C. Piecewise transmission rate with `B0` and `sigma_beta` estimated

Script: `code/03_fit_piecewise_B_estimate_B0_sigma.R`

This scenario uses the **same simulated piecewise dataset as Scenario 1B** because the data-generating model and simulation seed are identical.

The starting grid is the same 3 x 3 grid used in Scenario 1A.

Selected fit:

- `B0_hat = 4.5365`
- `sigma_beta_hat = 0.2206`
- evaluated log likelihood `= -197.3708`
- Monte Carlo SE `= 0.0462`

The selected run is simply the largest Monte Carlo-evaluated likelihood among the tested starts; it is not proof of a unique global maximum.

## MIF2 and particle-filter settings

Unless a scenario fixes `B0`, the scripts use the same computational structure:

- `Nmif = 100`
- `Np_mif = 1000`
- five independent particle-filter likelihood evaluations per fitted parameter vector
- `Np_eval = 5000`
- `Np_final = 50000`
- geometric cooling with `cooling.fraction.50 = 0.5`
- MIF2 seed reset to `20260628` for each starting point
- evaluation seeds `20260801` through `20260805`
- final particle-filter seed `999`

When both parameters are estimated:

- `B0 = ivp(0.20)`
- `sigma_beta = 0.05`
- both are log-transformed during estimation

When `B0` is fixed, only `sigma_beta` is perturbed and log-transformed.

The five likelihood estimates for each candidate are combined using `pomp::logmeanexp(..., se = TRUE)`.

## Reproducibility

Requirements:

- R
- `pomp`
- `tidyverse` / `ggplot2`

From the experiment root, run all three scenarios with:

```bash
chmod +x run_all.sh
./run_all.sh
```

or run them separately:

```bash
Rscript code/01_fit_constant_B4.R
Rscript code/02_fit_piecewise_B_fixed_B0.R
Rscript code/03_fit_piecewise_B_estimate_B0_sigma.R
```

Each script independently regenerates its simulated dataset, runs the multi-start fit, evaluates the candidates, selects the best stored candidate, runs a final particle filter, writes numerical outputs to `results/`, and writes PDF figures to `figures/`.

The three scenario scripts are the authoritative full-analysis workflows and are computationally nontrivial because they rerun MIF2 and particle filters. To reproduce only the nine report-facing PDFs from the committed compact CSV results, run:

```bash
Rscript code/04_regenerate_figures.R
```

Run this command from the Experiment 1 directory. It does not simulate data or rerun inference.

## Outputs

For each scenario, `results/` contains:

- `*_mif2_results.csv` - all starting-point fits and evaluated likelihoods;
- `*_best_fit.csv` - selected candidate;
- `*_filtered_B_path.csv` - true and filtered mean transmission-rate paths;
- `*_filtered_infectious_path.csv` - true and filtered mean infectious paths;
- `*_best_mif2.rds` - regenerable fitted object, ignored by Git.

The nine canonical, lightweight-regenerable PDF figures are:

- `constant_B4_filtered_B_path.pdf`
- `constant_B4_filtered_infectious_path.pdf`
- `constant_B4_starting_value_loglik.pdf`
- `piecewise_fixed_B0_filtered_B_path.pdf`
- `piecewise_fixed_B0_filtered_infectious_path.pdf`
- `piecewise_fixed_B0_starting_value_loglik.pdf`
- `piecewise_estimated_B0_filtered_B_path.pdf`
- `piecewise_estimated_B0_filtered_infectious_path.pdf`
- `piecewise_estimated_B0_starting_value_loglik.pdf`

The repository diagnostics are:

- `constant_B4_mif2_diagnostic.pdf`
- `piecewise_fixed_B0_mif2_diagnostic.pdf`
- `piecewise_estimated_B0_mif2_diagnostic.pdf`

All nine canonical PDFs are generated by `code/04_regenerate_figures.R`. Their exact source files are the matching `*_filtered_B_path.csv`, `*_filtered_infectious_path.csv`, `*_mif2_results.csv`, and `*_best_fit.csv` files under `results/`.

The MIF2 diagnostic PDFs are generated by the corresponding full scenario scripts (`code/01_fit_constant_B4.R`, `code/02_fit_piecewise_B_fixed_B0.R`, and `code/03_fit_piecewise_B_estimate_B0_sigma.R`) from the selected fitted objects. They use the standard `pomp` diagnostic plot: the first page shows final filtering diagnostics and the second page shows IF2 convergence diagnostics. They are retained as supporting computational diagnostics, not as canonical report figures; regenerating them requires the full fit or the ignored local `*_best_mif2.rds` objects.

## Figure conventions

The main figures use the same restrained visual language as the later experiments:

- PDF output only;
- no decorative grid or large in-figure title;
- true trajectories in black;
- filtered means in light blue;
- piecewise true `B(t)` drawn as two separate horizontal segments with no solid vertical connector at week 5;
- week 5 shown, when relevant, with a light-gray dashed reference line;
- `B(t)` panels use a common `0-6` vertical range;
- starting-value plots use open black points, a filled point for the selected candidate, and approximately `+/- 2` Monte Carlo standard-error bars.

These conventions are intended to keep the repository figures readable at report/manuscript scale and visually consistent with Experiments 2-4.

## Repository structure

```text
experiment_1_gamma_B_recovery/
├── code/
│   ├── 01_fit_constant_B4.R
│   ├── 02_fit_piecewise_B_fixed_B0.R
│   ├── 03_fit_piecewise_B_estimate_B0_sigma.R
│   └── 04_regenerate_figures.R
├── figures/
├── notes/
├── results/
├── .gitignore
├── run_all.sh
└── README.md
```

## Interpretation and limitations

Experiment 1 supports qualitative statements about fitting behavior on controlled examples and about sensitivity to the tested starting values. It does not support population-level statements about systematic recovery performance. For repeated-simulation recovery accuracy and model comparison, use the later experiments.

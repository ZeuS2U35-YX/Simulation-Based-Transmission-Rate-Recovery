# Simulation-Based Transmission-Rate Recovery

This repository contains the simulation, fitting, recovery, and figure-generation
workflows for a controlled comparison of three representations of the
time-varying epidemic transmission rate $B(t)$:

- a Gamma-noise process, reconstructed using observation-time particle
  filtering means;
- a deterministic non-periodic cubic B-spline for $\log B(t)$;
- a constant-$B$ model using one fitted value over the complete epidemic.

All three models are compared on the same 200 accepted simulated case series.
The simulation study uses a shared stochastic SIR process, negative-binomial
measurement model, multi-start iterated filtering framework, and 70 observation
times per replicate.

## Version

The manuscript-supporting GitHub release is version **1.0.0**. The release
records the code, retained aggregate outputs, figure source data, plotting
scripts, software documentation, seeds, and checksum-based pairing information
available in this repository. A DOI-backed archival record will be added
separately after the GitHub release has been deposited in a permanent
repository.

The repository is distributed under the [MIT License](LICENSE). Citation
metadata are provided in [CITATION.cff](CITATION.cff).

## Primary analysis

Experiment 5 extends the accepted Experiment 4 data with paired B-spline fits.
It does not generate another set of observations. The final design is:

~~~text
200 accepted Experiment 4 data sets
  -> 200 Gamma-noise filtering-mean reconstructions
  -> 200 selected B-spline trajectories
  -> 200 repeated constant-B estimates
  -> exactly 200 paired three-model comparison rows
~~~

The final three-model comparison gives mean replicate-level RMSE values of
0.5218 for Gamma-noise, 0.7128 for B-spline, and 1.1858 for constant-$B$.
Gamma-noise has lower RMSE than B-spline in 162 of 200 paired replicates, and
B-spline has lower RMSE than constant-$B$ in 188 of 200 paired replicates.
These results concern latent transmission-path recovery under the configured
single step-change simulation. They are not predictive-performance comparisons
and do not establish a general model ranking.

Primary reporting uses RMSE and signed mean error. The independent analysis
unit is the accepted simulation replicate. Multi-start fits, repeated
particle-filter likelihood evaluations, particles, Euler substeps, and
observation times are computational or within-replicate components rather than
independent replicates.

## Experiment hierarchy

| Experiment | Role | Status |
| --- | --- | --- |
| [Experiment 1](experiments/experiment_1_gamma_B_recovery/) | Developmental Gamma-noise fitting scenarios and likelihood-surface diagnostics | Supporting |
| [Experiment 2](experiments/experiment_2_large_scale_gamma_B_recovery_HPC/) | High-particle single-data-set Gamma-noise diagnostic | Supporting |
| [Experiment 3](experiments/experiment_3_gamma_B_recovery_accuracy/) | Earlier repeated Gamma-noise recovery study | Supporting; superseded for final numerical claims |
| [Experiment 4](experiments/experiment_4_nmif600_model_comparison/) | Accepted data, Gamma-noise fits, and constant-$B$ fits for 200 paired replicates | Primary computational foundation |
| [Experiment 5](experiments/experiment_5_bspline_B_recovery/) | B-spline extension and final paired three-model comparison | Primary manuscript analysis |

Experiments 1--3 document development history and supporting diagnostics.
Experiments 4 and 5 provide the paired computational evidence used in the
current manuscript.

## Final tracked outputs

The final three-model analysis is documented under
[Experiment 5](experiments/experiment_5_bspline_B_recovery/). Important tracked
outputs include:

- experiments/experiment_5_bspline_B_recovery/results/combined/bspline/
- experiments/experiment_5_bspline_B_recovery/results/comparison_three_models/
- experiments/experiment_5_bspline_B_recovery/figures/comparison_three_models/
- experiments/experiment_5_bspline_B_recovery/results/paired_input_manifest.csv

The figure directory includes PDF, SVG, PNG, and 600-dpi TIFF files, figure
source-data CSV files, and visual quality-control notes. Raw HPC task
directories are not tracked in Git; the retained combined results and source
data are sufficient for the documented lightweight figure reproduction.

## Reproducing the final figures

The Experiment 5 plotting environment is locked with renv. Follow
[REPRODUCE_FIGURES.md](experiments/experiment_5_bspline_B_recovery/REPRODUCE_FIGURES.md)
for the complete instructions. The abbreviated workflow is:

~~~bash
git clone https://github.com/ZeuS2U35-YX/Simulation-Based-Transmission-Rate-Recovery.git
cd Simulation-Based-Transmission-Rate-Recovery/experiments/experiment_5_bspline_B_recovery
Rscript -e 'renv::restore(prompt = FALSE)'
Rscript code/06_plot_three_model_comparison.R
~~~

This workflow regenerates the final comparison figures and their source-data
files. It does not rerun MIF2, the complete particle-filtering analyses, or the
HPC task arrays.

## Full computational workflows

The complete fitting workflows are computationally expensive and use Slurm.
Each experiment README records its scientific purpose, configuration,
execution order, retained outputs, and limitations. Before rerunning a full
workflow, consult the corresponding README and validate all input manifests,
simulation seeds, and checksums.

Software and computing-environment information is summarized in
[SOFTWARE.md](SOFTWARE.md).

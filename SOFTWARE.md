# Software and computing environment

This file distinguishes the dependency-locked environment used to reproduce
the final Experiment 5 comparison figures from the partially recorded
historical environments used for the full HPC fitting runs.

## Archived v1.0.0 scope

Repository version 1.0.0 contains a dependency lockfile managed by renv for
the final three-model figure-reproduction workflow under:

experiments/experiment_5_bspline_B_recovery/

This lockfile records the R and package versions used to regenerate the
archived Experiment 5 figures and figure source-data files. It is intentionally
scoped to the plotting workflow and does not claim exact restoration of the
historical environments used for the complete MIF2 and particle-filtering HPC
runs.

## Distribution modes

The project has three reproducibility modes with different software claims:

| Distribution or workflow | Included evidence | Software boundary |
| --- | --- | --- |
| Fresh Git clone or GitHub source archive | Tracked final three-model summaries, figure source data, plotting code, and `renv.lock` | Restores the final Experiment 5 plotting environment |
| Planned `v1.1.0` full-replication ZIP | The tracked source plus Experiment 4 accepted shared data, Experiment 4 Gamma/constant raw outputs, and Experiment 5 B-spline raw outputs | Supports inspection, validation, and recombination of retained task results; it does not reconstruct missing historical package versions |
| New full HPC rerun | Newly generated simulations and three fitted-model workflows | Requires Slurm, a compatible compiler, `pomp`, and a deliberately provisioned analysis environment |

Adding raw data and task outputs to the full-replication ZIP does not turn the
partially recorded historical HPC environment into a dependency-locked or
containerized environment.

## Locked Experiment 5 figure environment

The final comparison figures were generated with:

- R 4.5.2;
- renv 1.2.4;
- the package versions recorded in
  experiments/experiment_5_bspline_B_recovery/renv.lock.

The locked plotting environment includes ggplot2, patchwork, scales, svglite,
and ragg with their recursive dependencies. Restore it from the Experiment 5
directory with:

~~~bash
Rscript -e 'renv::restore(prompt = FALSE)'
Rscript -e 'renv::status()'
Rscript code/06_plot_three_model_comparison.R
~~~

The lockfile is scoped through the Experiment 5 configuration. It is not the
dependency environment for MIF2, particle filtering, or the other experiment
directories. Experiment 5 fitting, validation, and post-processing entry
points therefore invoke `Rscript --vanilla`; only the figure workflow relies
on the directory's automatic renv activation.

## Recorded historical HPC configuration

The Slurm scripts for Experiments 2--4 explicitly load:

- StdEnv/2020;
- R 4.1.2;
- the user library path recorded in the batch scripts.

This verifies the intended HPC R module, but it does not establish the exact
historical version of every installed R package. In particular, the complete
historical package environment for all 200-task fitting runs was not preserved
as a lockfile or container image.

Experiment 3 also contains a session record for a later single-run diagnostic
reconstruction. That reconstruction used:

- R 4.5.2 on Apple silicon macOS;
- pomp 6.4;
- deSolve 1.42;
- coda 0.19-4.1;
- mvtnorm 1.3-7;
- data.table 1.18.2.1;
- digest 0.6.39;
- lattice 0.22-7.

These versions document the later diagnostic reconstruction rather than the
original 200-task HPC environment.

## Currently tested analysis software

The repository has been validated with:

- R 4.5.2;
- pomp 6.4;
- tidyverse 2.0.0.

These are tested versions, not verified minimum requirements and not proof of
the complete historical HPC package set.

The principal analysis dependencies are:

- pomp for POMP construction, particle filtering, MIF2, filtering summaries,
  and Monte Carlo likelihood aggregation;
- tidyverse for developmental Experiment 1 scripts;
- ggplot2, patchwork, scales, svglite, and ragg for the final Experiment 5
  figures;
- base R packages including stats, graphics, grDevices, utils, and tools.

Compiling the embedded C snippets used by pomp requires a working compiler
toolchain compatible with the installed R distribution.

## Optional Python utility

experiments/experiment_4_nmif600_model_comparison/regenerate_exp4_figures.py
is an optional alternative figure generator. It is not the canonical final
three-model workflow.

The utility was tested with:

- Python 3.11.7;
- NumPy 1.26.4;
- pandas 2.1.4;
- Matplotlib 3.8.0.

These are tested versions rather than historical or minimum requirements.

## Command-line and HPC tools

- Rscript for local and batch execution;
- Bash for workflow wrappers;
- Slurm commands including sbatch and squeue for the full HPC workflows;
- an environment-modules implementation for the recorded module commands;
- Git for version control.

The lightweight Experiment 5 figure workflow does not require Slurm. The full
simulation, MIF2, and particle-filtering workflows are computationally
expensive and should be run only after consulting the experiment-specific
documentation.

## Reproducibility boundary

A fresh clone supports restoration of the final figure-generation environment
and regeneration of the tracked three-model comparison figures from retained
source data. The planned `v1.1.0` full-replication ZIP additionally preserves
the accepted shared data and required raw task outputs at their original
relative paths, enabling validation and recombination without claiming exact
restoration of the original HPC dependencies. A new full fit remains a
computationally expensive Slurm workflow whose complete historical package
environment was not preserved as a lockfile or container image. These
distinctions should be retained in manuscript reproducibility statements.

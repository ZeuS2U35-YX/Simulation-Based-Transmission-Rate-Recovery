# Experiment 4 shared simulated data

This directory contains the exact 200 accepted simulated epidemic data sets used by both Experiment 4 fitting models. The Gamma-noise and constant-`B` POMP models receive the same `observed_data.csv` within each task. These are synthetic simulation data, not observations from human participants or surveillance systems.

## Organization

Each `task_###/` directory contains:

- `observed_data.csv`: the 70 reported-case observations supplied to both fitting models;
- `simulated_data.csv`: the corresponding simulated latent `S`, `I`, `R`, and infection-accumulator `H` states together with reported cases;
- `simulation_metadata.csv`: simulation seed, accepted attempt, acceptance check, prescribed transmission-rate parameters, and fixed model parameters;
- `data_checksums.csv`: MD5 checksums for the three CSV files;
- `COMPLETE`: the atomic-workflow completion record used by the experiment scripts.

[`MANIFEST.csv`](MANIFEST.csv) provides one row per task so that readers can inspect seeds, acceptance information, and file checksums without opening 200 directories.

## Scientific role

The data-generating transmission-rate path is prescribed deterministically as `B(t) = 4` before week 5 and `B(t) = 2` from week 5 onward. The stochastic SIR process and negative-binomial observation model generate each accepted data set. Acceptance requires `max(H) > 20`; consequently, the reported comparison concerns accepted informative outbreaks rather than all attempted simulations.

The Gamma-noise and constant-`B` fitting models are not given the prescribed `B(t)` path or its change point. They read only the reported cases in `observed_data.csv`, together with the fixed model quantities documented in the Experiment 4 README and configuration.

## Reproduction and validation

Run commands from the Experiment 4 directory. The existing scripts verify the shared-data checksums before fitting. The retained data also allow the post-fit sampled trajectories to be reconstructed at the saved best-fit parameters without rerunning IF2 or the five independent likelihood evaluations.

```bash
Rscript code/07_generate_task1_comparison_figures.R
Rscript code/08_generate_task117_comparison_figures.R
```

The complete 200-task reconstruction command is documented in the parent [README](../README.md).

## Licence status

The root MIT License covers the repository software. A separate reuse licence for these simulated data must be confirmed before the final archival release; no additional data licence is asserted here.

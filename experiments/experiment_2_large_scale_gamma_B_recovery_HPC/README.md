# Experiment 2: Large-scale Gamma-B Recovery on HPC

## Experiment

One fixed simulated epidemic dataset is generated with

\[
B(t)=
\begin{cases}
4, & t<5,\\
2, & t\ge 5.
\end{cases}
\]

The Gamma-transition model is fitted from nine starting points:

- \(B_0 \in \{2,4,6\}\)
- \(\sigma_\beta \in \{0.10,0.30,0.45\}\)

Each starting point is assigned to one SLURM array task.

Formal settings:

- `Nmif = 100`
- `Np_mif = 50000`
- `Np_eval = 50000`
- `n_pf_evals = 5`
- `Np_final = 50000`

## Seed design

- Dataset simulation seed: `20260527`
- MIF2 seed shared by all nine starts: `20260628`
- Likelihood-evaluation seeds shared by all nine starts:
  `20260801` to `20260805`
- Final particle-filter seed: `999`

## HPC environment

- `StdEnv/2020`
- `r/4.1.2`
- `R_LIBS_USER=$HOME/packages-R4.1`

## Upload location

Upload this entire folder directly to:

```text
/global/home/hpc6245/
```

The resulting full path must be:

```text
/global/home/hpc6245/experiment2_large_scale_gamma_B_recovery_HPC
```

## Run

```bash
cd ~/experiment2_large_scale_gamma_B_recovery_HPC
chmod +x run_experiment.sh
chmod +x hpc/*.sh
./run_experiment.sh
```

This single command:

1. generates the fixed dataset;
2. submits nine MIF2 array tasks;
3. automatically submits a dependent finishing job;
4. combines results after all nine tasks succeed;
5. selects the best fit;
6. runs the final particle filter;
7. saves final data and figures.

## Outputs

### Fixed data

```text
data/fixed_piecewise_B_dataset.csv
```

### Per-task outputs

```text
results/array_output/task_001/
...
results/array_output/task_009/
```

Each task contains:

```text
mif2_result.csv
pfilter_evaluations.csv
mif2_object.rds
```

### Combined results

```text
results/combined_mif2_results.csv
results/best_fit.csv
results/best_mif2_object.rds
results/final_filtered_B_path.csv
results/final_filtered_infectious_path.csv
```

### Figures

```text
figures/best_B_path.png
figures/best_infectious_path.png
figures/best_mif2_trace.png
```

### Logs

```text
logs/slurm-mif2-<JOBID>_<TASKID>.out
logs/slurm-mif2-<JOBID>_<TASKID>.err
logs/slurm-finish-<JOBID>.out
logs/slurm-finish-<JOBID>.err
```

## Check progress

```bash
squeue -u "$USER"
```

Count completed task results:

```bash
find results/array_output -name "mif2_result.csv" | wc -l
```

View task 1 log:

```bash
tail -50 logs/slurm-mif2-<JOBID>_1.out
```

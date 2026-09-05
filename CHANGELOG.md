# Changelog

This file records user-visible analysis, reproducibility, and documentation
changes. The archived v1.0.0 release remains unchanged.

## Unreleased

### Analysis correction

- Reconstructed the Gamma-noise aggregate recovery path as the
  observation-time particle filtering mean. The archived v1.0.0 Experiment 4
  aggregate used an ancestry-sampled latent trajectory, so this correction
  changes the estimand and numerical summary; it is not documentation-only.
- Kept ancestry-preserving sampled Gamma paths as separately labelled
  selected-task illustration artifacts so they cannot be mistaken for primary
  metric inputs.
- Normalized final observation-grid reporting to $B(5)=4$ while retaining the
  simulator's process-time step at $t=5$. The correction does not change the
  simulated observations or fitted model objects.

### Added

- Provenance fields and path-semantics checks for combined recovery paths.
- A base-R release gate that checks syntax, the 200-task manifest, pairing,
  path semantics, full-window summaries, and retained headline values.
- A reproducible week-8 sensitivity analysis and task-level output.
- GitHub Actions validation for lightweight checks on pushes and pull requests.
- Dry-run modes and explicit confirmation guards for Experiment 4 and 5 Slurm
  submission wrappers.
- Explicit array ranges in guarded wrappers, so directly submitting a worker
  cannot silently launch a 199- or 200-task array.
- A guarded preview-and-submit wrapper for the Experiment 5 task-1 pilot.

### Documentation and environment

- Clarified the five-experiment hierarchy, interpretation limits, and
  fresh-clone versus full-HPC reproducibility boundary.
- Isolated the Experiment 5 plotting lockfile from fitting and post-processing
  entry points, which now use `Rscript --vanilla`.
- Added a concise project overview, quick verification path, supervision
  acknowledgment, and explicit unreleased-correction notice.

## 1.0.0 — 2026-08-27

- First archived public release.
- DOI: [10.5281/zenodo.22127917](https://doi.org/10.5281/zenodo.22127917).

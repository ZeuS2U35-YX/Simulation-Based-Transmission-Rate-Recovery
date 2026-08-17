# LaTeX report

This directory contains the complete report source, bibliography, scientific
figures, and the Experiment 4 code snapshot printed in the appendix.

The latest compiled copy is [report.pdf](report.pdf).

Build from this directory with:

```sh
latexmk -pdf report.tex
```

The Introduction is preserved from the full revision-04 report. The Model,
Inference and evaluation, Results, and limited Discussion statements are
synchronized with the current Experiment 4 definition: likelihood estimates
are aggregated with `pomp::logmeanexp`, the best fit is retained, and recovery
metrics use one ancestry-preserving sampled Gamma trajectory rather than a
filtering mean.

The canonical executable analysis and retained 200-task data/results live in
`../experiments/experiment_4_nmif600_model_comparison/`. The files under
`code/experiment4_source/` are a report-local snapshot so the code appendix is
self-contained.

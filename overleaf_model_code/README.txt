OVERLEAF REPORT CODE FILES

Purpose
-------
This directory contains report-facing code excerpts mapped to the Model and
Theory sections. The canonical executable implementation remains in

  experiments/experiment_4_nmif600_model_comparison/code/
  experiments/experiment_4_nmif600_model_comparison/config/

The excerpts retain the model definitions used by the experiment while
removing workflow, error-handling, and file-I/O code that would obscure the
scientific correspondence between the equations and their implementation.

Main-text mapping
-----------------

Section 2.5: Shared Measurement Model
  05_measurement_model.R

Section 2.6.1: Shared Stochastic SIR Component
  01_piecewise_latent_process.R, using the transition probabilities, binomial
  event draws, and S-I-R-H updates as one representative implementation of the
  shared stochastic transition mechanism.

Section 2.6.2: Prescribed Transmission-Rate Path
  01_piecewise_latent_process.R, using the t < t_switch branch.

Section 2.6.3: Gamma-Noise Transmission-Rate Process
  03_gamma_latent_process.R, using shape_B, scale_B, and rgamma.

Section 2.6.4: Constant-B Transmission-Rate Model
  09_constant_latent_process.R, where Beta enters the infection probability
  directly and is not updated as a latent state.

Section 2.7: Initial States and Fixed Parameters for Data Generation
  Present the values as a table in the main text. The implementation excerpts
  are 02_piecewise_initial_values.R, 04_gamma_initial_values.R,
  10_constant_initial_values.R, and 08_model_parameters.R.

Section 3.4: Multi-start Initialization for IF2
  Present the starting values as a table in the main text. The corresponding
  configuration, representative loop, independent particle-filter likelihood
  aggregation, and best-fit selection are in 12_multistart_initialization.R.

Section 3.5: Absolute Task-Level Mean Error and Residual Sum of Squares
  Keep the mathematical definitions in the main text. The implementation is
  shown in 13_recovery_metrics.R. Each Gamma-noise task uses the saved sampled
  trajectory; no filtering mean enters these recovery metrics.

Section 3.6: Sampled Latent Transmission-Rate Trajectory
  14_sampled_B_trajectory.R shows best-fit selection, the final particle filter
  with ancestry storage, extraction with filter_traj(), and removal of time zero
  before task-level metric calculation.

Supporting files
----------------

06_data_generating_pomp_object.R
  POMP construction for the prescribed-path data-generating model.

07_gamma_filtering_pomp_object.R
  POMP construction for the Gamma-noise fitted model.

11_constant_filtering_pomp_object.R
  POMP construction for the constant-B fitted model.

listings_setup.tex
  Add with \input{overleaf_model_code/listings_setup} in the preamble.

main_text_code_insertions.tex
  Copy the relevant blocks into Sections 2.5 and 2.6.1--2.6.4. This file is a
  placement guide and is not intended to be input wholesale.

appendix_code_listings.tex
  After \appendix, include with
  \input{overleaf_model_code/appendix_code_listings}.

Terminology
-----------

Use "multi-start initialization" for the overall IF2 strategy,
"starting-value combination" for one Gamma-noise pair, and "starting value"
for one constant-B initialization. Do not use "starting-value grid" in the
report prose.

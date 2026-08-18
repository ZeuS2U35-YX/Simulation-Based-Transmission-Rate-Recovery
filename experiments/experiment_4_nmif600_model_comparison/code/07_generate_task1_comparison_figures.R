# Generate the prespecified Experiment 4 task-1 sampled B-trajectory figure.
# Run from the Experiment 4 root:
#   Rscript code/07_generate_task1_comparison_figures.R

options(stringsAsFactors = FALSE, digits = 17)

source(file.path("code", "selected_B_trajectory_figure_helpers.R"))

generate_selected_B_trajectory_figure(
  task_id = 1L,
  expected_data_md5 = "64a1b5c02bdecc100b37ee56029391d3",
  expected_B0_hat = 4.35032977730012,
  expected_sigma_beta_hat = 0.206674783764282,
  expected_constant_B = 4.04534940404491,
  figure_stem = "01_selected_task_B_trajectory_comparison",
  selection_basis = "task 1 requested for the manuscript comparison figure"
)

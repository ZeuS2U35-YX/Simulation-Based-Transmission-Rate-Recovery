# Generate the prespecified Experiment 4 task-117 sampled B-trajectory figure.
# Run from the Experiment 4 root:
#   Rscript code/08_generate_task117_comparison_figures.R

options(stringsAsFactors = FALSE, digits = 17)

source(file.path("code", "selected_B_trajectory_figure_helpers.R"))

generate_selected_B_trajectory_figure(
  task_id = 117L,
  expected_data_md5 = "64dffb15867fda5ef262e2caf0e46bbf",
  expected_B0_hat = 3.6102663778531099,
  expected_sigma_beta_hat = 0.34528388982657399,
  expected_constant_B = 3.3878099385514902,
  figure_stem = "08_task117_B_trajectory_comparison",
  selection_basis = "task 117 requested as the second manuscript illustration"
)

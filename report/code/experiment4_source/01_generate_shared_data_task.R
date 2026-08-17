# ============================================================
# Generate one accepted shared data set for Experiment 4
#
# Both fitted models read these exact CSV files. They do not
# regenerate the epidemic independently.
#
# Usage:
# Rscript code/01_generate_shared_data_task.R \
#   <task_id> <shared_data_root> <paramlist_file>
# ============================================================

library(pomp)
options(stringsAsFactors = FALSE)

source("config/experiment_config.R")
source("code/model_components.R")
source("code/io_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript code/01_generate_shared_data_task.R <task_id> <shared_data_root> <paramlist_file>")
}

task_id <- suppressWarnings(as.integer(args[[1]]))
if (is.na(task_id) || task_id < 1L || task_id > experiment_config$n_tasks) {
  stop("task_id must be between 1 and ", experiment_config$n_tasks, ".")
}

shared_data_root <- if (length(args) >= 2L) args[[2]] else "shared_data"
paramlist_file <- if (length(args) >= 3L) args[[3]] else file.path("results", "paramlist.csv")

seed_row <- read_paramlist_task(
  paramlist_file,
  task_id,
  c("task_id", "simulation_seed")
)
simulation_seed <- as.integer(seed_row$simulation_seed[[1]])

atomic <- start_atomic_task(shared_data_root, task_id)
if (atomic$skip) quit(save = "no", status = 0L)
on.exit({
  if (dir.exists(atomic$temp_dir)) unlink(atomic$temp_dir, recursive = TRUE, force = TRUE)
}, add = TRUE)

cat(
  "Experiment = ", experiment_config$experiment_id,
  "\nTask ID = ", task_id,
  "\nSimulation seed = ", simulation_seed,
  "\nAcceptance rule = max(H) > ", experiment_config$acceptance_threshold,
  "\n", sep = ""
)

data_model <- make_data_generating_model(experiment_config)
set.seed(simulation_seed)

accepted <- FALSE
simulation_attempt <- 0L
simulated_data <- NULL

while (!accepted && simulation_attempt < experiment_config$max_simulation_attempts) {
  simulation_attempt <- simulation_attempt + 1L
  simulated_data <- simulate(
    data_model,
    params = experiment_config$true_parameters,
    nsim = 1,
    format = "data.frame",
    include.data = FALSE
  )

  accepted <- is.finite(max(simulated_data$H)) &&
    max(simulated_data$H) > experiment_config$acceptance_threshold
}

if (!accepted) {
  stop(
    "Task ", task_id, " did not generate an accepted simulation after ",
    experiment_config$max_simulation_attempts, " attempts."
  )
}

observed_data <- simulated_data[, c("week", "reports"), drop = FALSE]

observed_file <- file.path(atomic$temp_dir, "observed_data.csv")
simulated_file <- file.path(atomic$temp_dir, "simulated_data.csv")
metadata_file <- file.path(atomic$temp_dir, "simulation_metadata.csv")

write.csv(observed_data, observed_file, row.names = FALSE)
write.csv(simulated_data, simulated_file, row.names = FALSE)

metadata <- data.frame(
  experiment_id = experiment_config$experiment_id,
  task_id = task_id,
  simulation_seed = simulation_seed,
  simulation_attempt = simulation_attempt,
  accepted = TRUE,
  acceptance_threshold = experiment_config$acceptance_threshold,
  max_H = max(simulated_data$H),
  Beta_high = experiment_config$true_parameters[["Beta_high"]],
  Beta_low = experiment_config$true_parameters[["Beta_low"]],
  t_switch = experiment_config$true_parameters[["t_switch"]],
  mu_IR = experiment_config$true_parameters[["mu_IR"]],
  N = experiment_config$true_parameters[["N"]],
  rho = experiment_config$true_parameters[["rho"]],
  k = experiment_config$true_parameters[["k"]],
  stringsAsFactors = FALSE
)
write.csv(metadata, metadata_file, row.names = FALSE)

checksums <- data.frame(
  file = c("observed_data.csv", "simulated_data.csv", "simulation_metadata.csv"),
  md5 = c(file_md5(observed_file), file_md5(simulated_file), file_md5(metadata_file)),
  stringsAsFactors = FALSE
)
write.csv(
  checksums,
  file.path(atomic$temp_dir, "data_checksums.csv"),
  row.names = FALSE
)

commit_atomic_task(atomic$temp_dir, atomic$final_dir)
cat(
  "Accepted after ", simulation_attempt, " attempt(s); max(H) = ",
  max(simulated_data$H), ".\nSaved shared data to ", atomic$final_dir, "\n",
  sep = ""
)

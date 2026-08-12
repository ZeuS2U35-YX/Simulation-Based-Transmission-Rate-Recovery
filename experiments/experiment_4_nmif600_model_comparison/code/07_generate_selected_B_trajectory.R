# ============================================================
# Generate the prespecified Experiment 4 illustrative B(t) trajectory
#
# Run from the Experiment 4 root directory:
#   Rscript code/07_generate_selected_B_trajectory.R
#
# This script performs exactly one final 50,000-particle Gamma-model
# filter at the already selected task-117 parameter vector.  It compares
# the ancestry-preserving trajectory with the prescribed truth and the
# same task's fitted static constant-B estimate.  It does not rerun MIF2,
# evaluate likelihoods, or alter the existing 200-task metric inputs.
# ============================================================

options(stringsAsFactors = FALSE, digits = 17)

suppressPackageStartupMessages({
  library(pomp)
})

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The digest package is required for provenance checks.")
}

source(file.path("config", "experiment_config.R"))
source(file.path("code", "model_components.R"))

assert_close <- function(actual, expected, label, tolerance = 1e-10) {
  if (length(actual) != 1L || !is.finite(actual) ||
      abs(actual - expected) > tolerance) {
    stop(label, " mismatch: expected ", expected, ", found ", actual, ".")
  }
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

digest_task <- function(hex_digest, n_tasks = 200L) {
  bytes <- strtoi(substring(hex_digest, seq(1L, 64L, by = 2L),
                           seq(2L, 64L, by = 2L)), base = 16L)
  remainder <- 0L
  for (byte in bytes) remainder <- (remainder * 256L + byte) %% n_tasks
  1L + remainder
}

selection_label <- "experiment4-trajectory-2026-08-12"
expected_digest <- "da75fae89d9e42061b02d114495261d348cc0e1527553ad0b9402f7e942eb88c"
expected_data_md5 <- "64dffb15867fda5ef262e2caf0e46bbf"
task_id <- 117L
trajectory_seed <- 900000117L
Np <- 50000L
expected_old_final_pf_seed <- 137260900L

selection_digest <- digest::digest(
  selection_label,
  algo = "sha256",
  serialize = FALSE
)
if (!identical(selection_digest, expected_digest)) {
  stop("Selection-label SHA-256 mismatch.")
}
if (digest_task(selection_digest) != task_id) {
  stop("Complete-digest big-endian modulo rule did not select task 117.")
}
if (trajectory_seed == expected_old_final_pf_seed) {
  stop("The new trajectory seed must differ from the historical final-PF seed.")
}
if (experiment_config$Np_final != Np) {
  stop("Canonical Experiment 4 Np_final is not 50000.")
}

gamma_best_path <- file.path(
  "results", "combined", "gamma", "combined_best_fit_summary.csv"
)
constant_best_path <- file.path(
  "results", "combined", "constant", "combined_best_fit_summary.csv"
)
paramlist_path <- file.path("results", "paramlist.csv")
data_path <- file.path(
  "shared_data", sprintf("task_%03d", task_id), "observed_data.csv"
)
required_paths <- c(gamma_best_path, constant_best_path, paramlist_path, data_path)
if (any(!file.exists(required_paths))) {
  stop("Missing required input(s): ",
       paste(required_paths[!file.exists(required_paths)], collapse = ", "))
}

gamma_best <- read.csv(gamma_best_path, check.names = FALSE)
constant_best <- read.csv(constant_best_path, check.names = FALSE)
paramlist <- read.csv(paramlist_path, check.names = FALSE)
observed_data <- read.csv(data_path, check.names = FALSE)

gamma_row <- gamma_best[gamma_best$task_id == task_id, , drop = FALSE]
constant_row <- constant_best[constant_best$task_id == task_id, , drop = FALSE]
seed_row <- paramlist[paramlist$task_id == task_id, , drop = FALSE]
if (nrow(gamma_row) != 1L || nrow(constant_row) != 1L || nrow(seed_row) != 1L) {
  stop("Task 117 must have one Gamma fit, one constant-B fit, and one seed row.")
}
if (!all(c("week", "reports") %in% names(observed_data))) {
  stop("Task-117 saved data lack week or reports.")
}
observed_data <- observed_data[, c("week", "reports"), drop = FALSE]
if (nrow(observed_data) != 70L) stop("Task 117 must have exactly 70 observations.")
if (any(!is.finite(observed_data$week)) || any(!is.finite(observed_data$reports))) {
  stop("Task-117 observations contain non-finite values.")
}
if (is.unsorted(observed_data$week, strictly = TRUE) || anyDuplicated(observed_data$week)) {
  stop("Task-117 observation times must be strictly increasing and unique.")
}
expected_observation_times <- seq(
  from = experiment_config$observation_interval,
  to = experiment_config$n_weeks,
  by = experiment_config$observation_interval
)
if (length(expected_observation_times) != 70L ||
    max(abs(observed_data$week - expected_observation_times)) > 1e-12) {
  stop("Task-117 observations do not use the configured 70-time grid.")
}

selected_data_md5 <- unname(tools::md5sum(data_path))
if (!identical(selected_data_md5, expected_data_md5)) {
  stop("Task-117 observed-data MD5 mismatch.")
}
if (!identical(as.character(gamma_row$observed_data_md5[[1]]), expected_data_md5) ||
    !identical(as.character(constant_row$observed_data_md5[[1]]), expected_data_md5)) {
  stop("Gamma and constant selected-fit records do not reference task-117 data.")
}
if (gamma_row$simulation_seed[[1]] != constant_row$simulation_seed[[1]] ||
    gamma_row$simulation_seed[[1]] != seed_row$simulation_seed[[1]]) {
  stop("Gamma, constant-B, and seed records do not identify the same dataset.")
}
assert_close(gamma_row$B0_hat[[1]], 3.61026637785311, "Gamma B0_hat")
assert_close(gamma_row$sigma_beta_hat[[1]], 0.345283889826574,
             "Gamma sigma_beta_hat")
assert_close(constant_row$Beta_hat[[1]], 3.38780993855149,
             "constant Beta_hat")
if (seed_row$gamma_final_pf_seed[[1]] != expected_old_final_pf_seed) {
  stop("Task-117 historical Gamma final-PF seed mismatch.")
}

output_dir <- file.path("results", "selected_trajectory")
figure_dir <- file.path("figures", "comparison")
if (!dir.exists(output_dir) &&
    !dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)) {
  stop("Cannot create selected-trajectory output directory before filtering.")
}
if (!dir.exists(figure_dir) &&
    !dir.create(figure_dir, recursive = TRUE, showWarnings = TRUE)) {
  stop("Cannot create comparison-figure directory before filtering.")
}
trajectory_path <- file.path(output_dir, "experiment4_task117_B_trajectory_comparison.csv")
provenance_path <- file.path(output_dir, "experiment4_task117_B_trajectory_provenance.txt")
figure_path <- file.path(figure_dir, "01_selected_task_B_trajectory_comparison.pdf")

gamma_model <- make_gamma_model(observed_data, experiment_config)
theta_fitted <- gamma_baseline_parameters(experiment_config)
theta_fitted[["B0"]] <- gamma_row$B0_hat[[1]]
theta_fitted[["sigma_beta"]] <- gamma_row$sigma_beta_hat[[1]]

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
rng_kind <- RNGkind()
set.seed(trajectory_seed)
pf_trajectory <- pfilter(
  gamma_model,
  params = theta_fitted,
  Np = Np,
  filter.traj = TRUE
)
extracted <- filter_traj(
  pf_trajectory,
  vars = "B",
  format = "data.frame"
)

if (!all(c("time", "value") %in% names(extracted))) {
  stop("Unexpected filter_traj() data-frame schema: ",
       paste(names(extracted), collapse = ", "))
}
if ("name" %in% names(extracted) &&
    !identical(unique(as.character(extracted$name)), "B")) {
  stop("filter_traj() returned a state other than B.")
}
trajectory_id_columns <- setdiff(names(extracted), c("name", "time", "value"))
if (length(trajectory_id_columns) > 0L &&
    any(vapply(extracted[trajectory_id_columns], function(x) length(unique(x)),
               integer(1)) != 1L)) {
  stop("filter_traj() returned more than one trajectory.")
}
if (nrow(extracted) != 71L) stop("Expected t0 plus 70 trajectory times.")
if (any(!is.finite(extracted$time)) || any(!is.finite(extracted$value))) {
  stop("Extracted trajectory contains non-finite values.")
}
if (is.unsorted(extracted$time, strictly = TRUE) || anyDuplicated(extracted$time)) {
  stop("Extracted trajectory times must be strictly increasing and unique.")
}
if (abs(extracted$time[[1]]) > 1e-12 ||
    max(abs(extracted$time[-1L] - observed_data$week)) > 1e-12) {
  stop("Extracted times do not equal t0 followed by the 70 observation times.")
}

trajectory <- data.frame(
  experiment = 4L,
  task_id = task_id,
  week = as.numeric(extracted$time),
  B_trajectory = as.numeric(extracted$value),
  B_true = true_B_at_times(extracted$time, experiment_config),
  B_constant = rep(constant_row$Beta_hat[[1]], nrow(extracted)),
  is_time_zero = abs(extracted$time) < 1e-12,
  stringsAsFactors = FALSE
)
if (length(unique(trajectory$B_constant)) != 1L) {
  stop("The constant-B display line is not one repeated fitted value.")
}
if (any(trajectory$B_true[trajectory$week < 5] != 4) ||
    any(trajectory$B_true[trajectory$week >= 5] != 2)) {
  stop("The plotted truth does not implement the week-5 switch.")
}

write.csv(trajectory, trajectory_path, row.names = FALSE)
trajectory_sha256 <- sha256_file(trajectory_path)

repo_root <- normalizePath(file.path(getwd(), "..", ".."), mustWork = TRUE)
repo_commit <- system2("git", c("-C", repo_root, "rev-parse", "HEAD"),
                       stdout = TRUE)
generated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
provenance <- c(
  "experiment: 4",
  paste0("selection_label: ", selection_label),
  paste0("selection_sha256: ", selection_digest),
  "hash_to_task_rule: UTF-8 label without trailing newline; full SHA-256 interpreted as one unsigned big-endian integer H; task=1+(H mod 200)",
  paste0("selected_task_id: ", task_id),
  "task_selected_before_trajectory_inspection: yes",
  paste0("trajectory_seed: ", trajectory_seed),
  paste0("historical_filtering_mean_final_pf_seed: ", expected_old_final_pf_seed),
  paste0("particle_count: ", Np),
  paste0("R_version: ", R.version.string),
  paste0("pomp_version: ", as.character(packageVersion("pomp"))),
  paste0("RNGkind: ", paste(rng_kind, collapse = ", ")),
  paste0("generated_at_utc: ", generated_at_utc),
  paste0("repository_base_commit: ", repo_commit),
  paste0("selected_data_path: ", data_path),
  paste0("selected_data_md5: ", selected_data_md5),
  paste0("selected_data_sha256: ", sha256_file(data_path)),
  "Gamma_and_constant_records_match_selected_data_md5: yes",
  paste0("simulation_seed: ", gamma_row$simulation_seed[[1]]),
  paste0("Gamma_best_run: ", gamma_row$best_run[[1]]),
  paste0("constant_best_run: ", constant_row$best_run[[1]]),
  paste0("Gamma_B0_hat: ", format(gamma_row$B0_hat[[1]], digits = 17)),
  paste0("Gamma_sigma_beta_hat: ",
         format(gamma_row$sigma_beta_hat[[1]], digits = 17)),
  paste0("constant_Beta_hat: ",
         format(constant_row$Beta_hat[[1]], digits = 17)),
  "fixed_parameters: mu_IR=3; N=10000; rho=0.5; k=10",
  "truth: B_true(t)=4 for t<5; B_true(t)=2 for t>=5",
  "process_grid: Euler delta.t=1/30 week",
  "observation_grid: 70 times, 1/7 through 10 weeks",
  "PF_method: pomp::pfilter(..., Np=50000, filter.traj=TRUE) at the selected Gamma plug-in parameter vector",
  "extraction_method: pomp::filter_traj(pf_trajectory, vars='B', format='data.frame')",
  "trajectory_semantics: one terminal particle's ancestry-preserving finite-particle plug-in approximation to a smoothing trajectory conditioned on the complete selected observation series",
  "constant_line_semantics: task-117 selected static constant-B estimate repeated over the 71 plotted state times",
  "exact_posterior_draw: no",
  "parameter_uncertainty_integrated: no",
  "filtering_mean_or_across_task_average: no",
  "number_of_observations: 70",
  "number_of_extracted_trajectory_times: 71",
  "t0_included: yes",
  paste0("trajectory_csv_sha256: ", trajectory_sha256)
)
writeLines(provenance, provenance_path, useBytes = TRUE)

truth_x <- c(0, 5, 5, 10)
truth_y <- c(4, 4, 2, 2)
y_range <- range(c(0, 6, trajectory$B_trajectory, trajectory$B_constant),
                 finite = TRUE)
y_pad <- 0.04 * diff(y_range)
pdf(figure_path, width = 6.2, height = 4.2, useDingbats = FALSE)
par(bty = "l", las = 1, lend = "butt", tcl = -0.25,
    mgp = c(2.2, 0.65, 0), mar = c(3.8, 4.1, 0.8, 0.7),
    cex.axis = 0.92, cex.lab = 1.03)
plot(trajectory$week, trajectory$B_trajectory,
     type = "l", lwd = 1.45, lty = 1, col = "#0072B2",
     xlim = c(0, 10), ylim = c(y_range[[1]] - y_pad, y_range[[2]] + y_pad),
     xaxs = "i", xlab = "Week", ylab = "Transmission rate B(t)")
lines(truth_x, truth_y, lwd = 1.65, lty = 2, col = "black")
lines(trajectory$week, trajectory$B_constant,
      lwd = 1.35, lty = 4, col = "#A05A4A")
abline(v = 5, lty = 3, lwd = 0.75, col = "grey65")
lines(trajectory$week, trajectory$B_trajectory,
      lwd = 1.45, lty = 1, col = "#0072B2")
legend("topright",
       legend = c("True B(t)", "Gamma-noise trajectory (task 117)",
                  "Fitted constant B (task 117)"),
       col = c("black", "#0072B2", "#A05A4A"),
       lty = c(2, 1, 4), lwd = c(1.65, 1.45, 1.35),
       bty = "n", cex = 0.82, seg.len = 3.0)
dev.off()

cat("Experiment 4 selected trajectory generated successfully.\n",
    "Trajectory: ", trajectory_path, "\n",
    "Provenance: ", provenance_path, "\n",
    "Figure: ", figure_path, "\n", sep = "")

# ============================================================
# Generate the prespecified Experiment 3 illustrative B(t) trajectory
#
# Run from the Experiment 3 root directory:
#   Rscript code/06_generate_selected_B_trajectory.R
#
# This script performs exactly one final 50,000-particle filter at the
# already selected task-145 parameter vector.  It preserves particle
# ancestry via filter.traj=TRUE and filter_traj().  It does not rerun
# MIF2, evaluate likelihoods, or alter the existing filtering-mean files.
# ============================================================

options(stringsAsFactors = FALSE, digits = 17)

suppressPackageStartupMessages({
  library(pomp)
})

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The digest package is required for provenance checks.")
}

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

selection_label <- "experiment3-trajectory-2026-08-12"
expected_digest <- "a7d32e608ff0f1c1f20a36f04d68503a949ccda4e39753819f89db11c4f080f8"
task_id <- 145L
trajectory_seed <- 900000145L
Np <- 50000L
expected_old_final_pf_seed <- 165260900L

selection_digest <- digest::digest(
  selection_label,
  algo = "sha256",
  serialize = FALSE
)
if (!identical(selection_digest, expected_digest)) {
  stop("Selection-label SHA-256 mismatch.")
}
if (digest_task(selection_digest) != task_id) {
  stop("Complete-digest big-endian modulo rule did not select task 145.")
}
if (trajectory_seed == expected_old_final_pf_seed) {
  stop("The new trajectory seed must differ from the historical final-PF seed.")
}

best_path <- file.path("results", "combined", "combined_best_fit_summary.csv")
data_path <- file.path("results", "combined", "combined_simulated_data.csv")
paramlist_path <- file.path("results", "paramlist.csv")
required_paths <- c(best_path, data_path, paramlist_path)
if (any(!file.exists(required_paths))) {
  stop("Missing required input(s): ",
       paste(required_paths[!file.exists(required_paths)], collapse = ", "))
}

best <- read.csv(best_path, check.names = FALSE)
paramlist <- read.csv(paramlist_path, check.names = FALSE)
simulated <- read.csv(data_path, check.names = FALSE)

best_row <- best[best$task_id == task_id, , drop = FALSE]
seed_row <- paramlist[paramlist$task_id == task_id, , drop = FALSE]
task_data <- simulated[simulated$task_id == task_id, , drop = FALSE]
if (nrow(best_row) != 1L || nrow(seed_row) != 1L) {
  stop("Task 145 must have exactly one selected-fit row and one seed row.")
}
if (nrow(task_data) != 70L) stop("Task 145 must have exactly 70 observations.")
if (!all(c("week", "reports") %in% names(task_data))) {
  stop("Task-145 saved data lack week or reports.")
}
observed_data <- task_data[, c("week", "reports"), drop = FALSE]
if (any(!is.finite(observed_data$week)) || any(!is.finite(observed_data$reports))) {
  stop("Task-145 observations contain non-finite values.")
}
if (is.unsorted(observed_data$week, strictly = TRUE) || anyDuplicated(observed_data$week)) {
  stop("Task-145 observation times must be strictly increasing and unique.")
}
expected_observation_times <- seq(from = 1 / 7, to = 10, by = 1 / 7)
if (length(expected_observation_times) != 70L ||
    max(abs(observed_data$week - expected_observation_times)) > 1e-12) {
  stop("Task-145 observations do not use the expected 70-time grid.")
}

if (best_row$simulation_seed[[1]] != 1145L || seed_row$simulation_seed[[1]] != 1145L) {
  stop("Task-145 simulation seed mismatch.")
}
if (best_row$best_run[[1]] != 9L) stop("Task-145 selected run mismatch.")
assert_close(best_row$B0_hat[[1]], 4.08137604772541, "B0_hat")
assert_close(best_row$sigma_beta_hat[[1]], 0.318461288709797,
             "sigma_beta_hat")
if (seed_row$final_pf_seed[[1]] != expected_old_final_pf_seed) {
  stop("Task-145 historical final-PF seed mismatch.")
}

observation_snapshot <- tempfile(fileext = ".csv")
on.exit(unlink(observation_snapshot), add = TRUE)
write.csv(observed_data, observation_snapshot, row.names = FALSE)
selected_observation_sha256 <- sha256_file(observation_snapshot)
source_data_sha256 <- sha256_file(data_path)

output_dir <- file.path("results", "selected_trajectory")
figure_dir <- "figures"
if (!dir.exists(output_dir) &&
    !dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)) {
  stop("Cannot create selected-trajectory output directory before filtering.")
}
if (!dir.exists(figure_dir) &&
    !dir.create(figure_dir, recursive = TRUE, showWarnings = TRUE)) {
  stop("Cannot create figure directory before filtering.")
}
trajectory_path <- file.path(output_dir, "experiment3_task145_B_trajectory.csv")
provenance_path <- file.path(output_dir, "experiment3_task145_B_trajectory_provenance.txt")
figure_path <- file.path(figure_dir, "01_selected_task_B_trajectory.pdf")

sir_step_gamma <- Csnippet("
  double shape_B;
  double scale_B;
  double dN_SI;
  double dN_IR;
  double p_SI;
  double p_IR;

  if (sigma_beta > 0.0) {
    shape_B = 1.0 / (sigma_beta * sigma_beta * dt);
    scale_B = B * sigma_beta * sigma_beta * dt;
    B = rgamma(shape_B, scale_B);
  }

  p_SI = 1.0 - exp(-B * I / N * dt);
  p_IR = 1.0 - exp(-mu_IR * dt);
  dN_SI = rbinom(S, p_SI);
  dN_IR = rbinom(I, p_IR);
  S = S - dN_SI;
  I = I + dN_SI - dN_IR;
  R = R + dN_IR;
  H = H + dN_SI;
")

sir_rinit_gamma <- Csnippet("
  S = N - 10;
  I = 10;
  R = 0;
  H = 0;
  B = B0;
")

sir_rmeas <- Csnippet("
  reports = rnbinom_mu(k, rho * H);
")

sir_dmeas <- Csnippet("
  lik = dnbinom_mu(reports, k, rho * H, give_log);
")

gamma_model <- pomp(
  data = observed_data,
  times = "week",
  t0 = 0,
  rinit = sir_rinit_gamma,
  rprocess = euler(sir_step_gamma, delta.t = 1 / 30),
  rmeasure = sir_rmeas,
  dmeasure = sir_dmeas,
  accumvars = "H",
  statenames = c("S", "I", "R", "H", "B"),
  paramnames = c("B0", "sigma_beta", "mu_IR", "N", "rho", "k"),
  partrans = parameter_trans(log = c("B0", "sigma_beta"))
)

theta_fitted <- c(
  B0 = best_row$B0_hat[[1]],
  sigma_beta = best_row$sigma_beta_hat[[1]],
  mu_IR = 3,
  N = 10000,
  rho = 0.5,
  k = 10
)

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
  experiment = 3L,
  task_id = task_id,
  week = as.numeric(extracted$time),
  B_trajectory = as.numeric(extracted$value),
  B_true = ifelse(extracted$time < 5, 4, 2),
  is_time_zero = abs(extracted$time) < 1e-12,
  stringsAsFactors = FALSE
)

write.csv(trajectory, trajectory_path, row.names = FALSE)
trajectory_sha256 <- sha256_file(trajectory_path)

repo_root <- normalizePath(file.path(getwd(), "..", ".."), mustWork = TRUE)
repo_commit <- system2("git", c("-C", repo_root, "rev-parse", "HEAD"),
                       stdout = TRUE)
generated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
provenance <- c(
  "experiment: 3",
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
  paste0("selected_data_path: ", data_path, " [task_id=145; columns=week,reports]"),
  paste0("selected_observation_sha256: ", selected_observation_sha256),
  "selected_observation_checksum_method: SHA-256 of write.csv(row.names=FALSE) serialization of the exact 70-row week,reports input",
  paste0("source_combined_data_sha256: ", source_data_sha256),
  paste0("simulation_seed: ", best_row$simulation_seed[[1]]),
  paste0("best_run: ", best_row$best_run[[1]]),
  paste0("B0_hat: ", format(best_row$B0_hat[[1]], digits = 17)),
  paste0("sigma_beta_hat: ", format(best_row$sigma_beta_hat[[1]], digits = 17)),
  "fixed_parameters: mu_IR=3; N=10000; rho=0.5; k=10",
  "truth: B_true(t)=4 for t<5; B_true(t)=2 for t>=5",
  "process_grid: Euler delta.t=1/30 week",
  "observation_grid: 70 times, 1/7 through 10 weeks",
  "PF_method: pomp::pfilter(..., Np=50000, filter.traj=TRUE) at the selected plug-in parameter vector",
  "extraction_method: pomp::filter_traj(pf_trajectory, vars='B', format='data.frame')",
  "trajectory_semantics: one terminal particle's ancestry-preserving finite-particle plug-in approximation to a smoothing trajectory conditioned on the complete selected observation series",
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
y_range <- range(c(0, 6, trajectory$B_trajectory), finite = TRUE)
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
abline(v = 5, lty = 3, lwd = 0.75, col = "grey65")
lines(trajectory$week, trajectory$B_trajectory,
      lwd = 1.45, lty = 1, col = "#0072B2")
legend("topright",
       legend = c("True B(t)", "Gamma-noise trajectory (task 145)"),
       col = c("black", "#0072B2"), lty = c(2, 1),
       lwd = c(1.65, 1.45), bty = "n", cex = 0.86, seg.len = 3.0)
dev.off()

cat("Experiment 3 selected trajectory generated successfully.\n",
    "Trajectory: ", trajectory_path, "\n",
    "Provenance: ", provenance_path, "\n",
    "Figure: ", figure_path, "\n", sep = "")

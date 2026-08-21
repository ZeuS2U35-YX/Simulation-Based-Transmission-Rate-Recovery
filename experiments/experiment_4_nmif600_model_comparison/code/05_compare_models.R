# ============================================================
# Gamma-noise versus constant-B comparison
# on exactly the same shared data sets
#
# Produces validated paired summary tables and a set of figures
# designed for direct interpretation in the report.
# ============================================================

options(stringsAsFactors = FALSE)
source("config/experiment_config.R")

true_B_driver_at_endpoints <- function(times, config) {
  ifelse(
    times <= config$true_parameters[["t_switch"]],
    config$true_parameters[["Beta_high"]],
    config$true_parameters[["Beta_low"]]
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop("Usage: Rscript code/05_compare_models.R <gamma_combined_dir> <constant_combined_dir> <output_dir> <figures_dir> [include_selected_trajectory_figure]")
}

gamma_dir <- args[[1]]
constant_dir <- args[[2]]
output_dir <- args[[3]]
figures_dir <- args[[4]]
include_selected_trajectory_figure <- if (length(args) >= 5L) {
  identical(tolower(args[[5]]), "true")
} else {
  TRUE
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

read_combined <- function(dir, path_file) {
  list(
    best = read.csv(file.path(dir, "combined_best_fit_summary.csv"), check.names = FALSE),
    paths = read.csv(path_file, check.names = FALSE),
    runs = read.csv(file.path(dir, "combined_mif2_results.csv"), check.names = FALSE)
  )
}

gamma_metric_path <- Sys.getenv(
  "EXP4_GAMMA_METRIC_PATH",
  unset = file.path(gamma_dir, "combined_B_filtering_means.csv")
)
constant_metric_path <- file.path(constant_dir, "combined_B_paths.csv")
if (!file.exists(gamma_metric_path)) {
  stop(
    "Missing primary Gamma recovery input: ", gamma_metric_path, ". ",
    "Run code/10_regenerate_filtering_mean_B_paths.R first."
  )
}

gamma <- read_combined(gamma_dir, gamma_metric_path)
constant <- read_combined(constant_dir, constant_metric_path)

if (!("estimate_semantics" %in% names(gamma$paths)) ||
    !all(gamma$paths$estimate_semantics == "particle_filtering_mean")) {
  stop("Primary Gamma recovery input is not labelled as particle filtering means.")
}

gamma_tasks <- sort(unique(gamma$best$task_id))
constant_tasks <- sort(unique(constant$best$task_id))
if (!identical(gamma_tasks, constant_tasks)) {
  stop("Gamma-noise and constant-B combined outputs do not contain the same task IDs.")
}

# ------------------------------------------------------------
# Verify that both models used exactly the same data and Nmif.
# ------------------------------------------------------------

pair_check <- merge(
  gamma$best[, c("task_id", "simulation_seed", "observed_data_md5", "Nmif")],
  constant$best[, c("task_id", "simulation_seed", "observed_data_md5", "Nmif")],
  by = "task_id",
  suffixes = c("_gamma", "_constant")
)

pair_check$same_simulation_seed <- pair_check$simulation_seed_gamma == pair_check$simulation_seed_constant
pair_check$same_observed_data_md5 <- pair_check$observed_data_md5_gamma == pair_check$observed_data_md5_constant
pair_check$both_nmif_600 <- pair_check$Nmif_gamma == experiment_config$Nmif & pair_check$Nmif_constant == experiment_config$Nmif
write.csv(pair_check, file.path(output_dir, "shared_data_pair_check.csv"), row.names = FALSE)

if (!all(pair_check$same_simulation_seed) || !all(pair_check$same_observed_data_md5) || !all(pair_check$both_nmif_600)) {
  stop("The two models did not use identical shared data or identical Nmif=600 settings for every task.")
}

# ------------------------------------------------------------
# Recovery metrics.
# ------------------------------------------------------------

# B_true in older retained path tables used the right-continuous path value.
# Recovery estimates at t_n are endpoint states that drove the final Euler
# substep ending at t_n, so derive the aligned interval-driving target from the
# time column here. This changes only the evaluation target and does not rerun
# a particle filter or alter any stored B estimate.
align_endpoint_driver_truth <- function(paths) {
  paths$B_true <- true_B_driver_at_endpoints(
    paths$week,
    experiment_config
  )
  paths
}

gamma$paths <- align_endpoint_driver_truth(gamma$paths)
constant$paths <- align_endpoint_driver_truth(constant$paths)

calculate_metrics <- function(paths, model_name, estimate_semantics) {
  required_columns <- c("task_id", "week", "B_estimate", "B_true")
  if (!all(required_columns %in% names(paths))) {
    stop(model_name, " paths lack required recovery-metric columns.")
  }
  expected_driver_truth <- true_B_driver_at_endpoints(
    paths$week,
    experiment_config
  )
  stopifnot(all(abs(paths$B_true - expected_driver_truth) <= 1e-12))
  paths$error <- paths$B_estimate - paths$B_true
  paths$squared_error <- paths$error^2
  split_paths <- split(paths, paths$task_id)

  out <- do.call(rbind, lapply(split_paths, function(x) {
    expected_times <- seq(
      from = experiment_config$observation_interval,
      to = experiment_config$n_weeks,
      by = experiment_config$observation_interval
    )
    if (nrow(x) != 70L || length(expected_times) != 70L ||
        max(abs(x$week - expected_times)) > 1e-12 ||
        is.unsorted(x$week, strictly = TRUE) || anyDuplicated(x$week)) {
      stop(
        model_name, " task ", x$task_id[[1]],
        " does not contain exactly the 70 configured observation times."
      )
    }
    mean_error <- mean(x$error)
    data.frame(
      task_id = x$task_id[[1]],
      model = model_name,
      estimate_semantics = estimate_semantics,
      RSS = sum(x$squared_error),
      RMSE = sqrt(mean(x$squared_error)),
      mean_error = mean_error,
      AOB = abs(mean_error),
      mean_error_through_5 = mean(
        x$error[x$week <= experiment_config$true_parameters[["t_switch"]]]
      ),
      mean_error_after_5 = mean(
        x$error[x$week > experiment_config$true_parameters[["t_switch"]]]
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

gamma_metrics <- calculate_metrics(
  gamma$paths,
  "gamma_noise",
  "particle_filtering_mean"
)
constant_metrics <- calculate_metrics(
  constant$paths,
  "constant_B",
  "repeated_static_estimate"
)
model_task_metrics <- rbind(gamma_metrics, constant_metrics)
write.csv(model_task_metrics, file.path(output_dir, "model_task_metrics.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Likelihood gaps and starting-value sensitivity.
# ------------------------------------------------------------

likelihood_gap <- function(runs, model_name) {
  finite <- runs[is.finite(runs$logLik), , drop = FALSE]
  split_runs <- split(finite, finite$task_id)
  out <- do.call(rbind, lapply(split_runs, function(x) {
    x <- x[order(x$logLik, decreasing = TRUE), , drop = FALSE]
    data.frame(
      task_id = x$task_id[[1]], model = model_name,
      best_run = x$run[[1]], best_logLik = x$logLik[[1]],
      second_best_logLik = if (nrow(x) >= 2L) x$logLik[[2]] else NA_real_,
      logLik_gap = if (nrow(x) >= 2L) x$logLik[[1]] - x$logLik[[2]] else NA_real_,
      n_successful_starts = nrow(x),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

gamma_gaps <- likelihood_gap(gamma$runs, "gamma_noise")
constant_gaps <- likelihood_gap(constant$runs, "constant_B")
write.csv(rbind(gamma_gaps, constant_gaps), file.path(output_dir, "likelihood_gaps_by_model.csv"), row.names = FALSE)

start_sensitivity_gamma <- do.call(rbind, lapply(split(gamma$runs[is.finite(gamma$runs$logLik), ], gamma$runs$task_id[is.finite(gamma$runs$logLik)]), function(x) {
  data.frame(
    task_id = x$task_id[[1]], model = "gamma_noise",
    estimate_range_1 = diff(range(x$B0_hat)),
    estimate_range_2 = diff(range(x$sigma_beta_hat)),
    stringsAsFactors = FALSE
  )
}))
start_sensitivity_constant <- do.call(rbind, lapply(split(constant$runs[is.finite(constant$runs$logLik), ], constant$runs$task_id[is.finite(constant$runs$logLik)]), function(x) {
  data.frame(
    task_id = x$task_id[[1]], model = "constant_B",
    estimate_range_1 = diff(range(x$Beta_hat)),
    estimate_range_2 = NA_real_,
    stringsAsFactors = FALSE
  )
}))
write.csv(
  rbind(start_sensitivity_gamma, start_sensitivity_constant),
  file.path(output_dir, "starting_value_sensitivity.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Paired task-level comparison table.
# ------------------------------------------------------------

paired <- merge(
  gamma_metrics,
  constant_metrics,
  by = "task_id",
  suffixes = c("_gamma", "_constant")
)
paired <- merge(
  paired,
  gamma$best[, c("task_id", "logLik", "logLik_se", "B0_hat", "sigma_beta_hat", "status")],
  by = "task_id"
)
names(paired)[names(paired) %in% c("logLik", "logLik_se", "B0_hat", "sigma_beta_hat", "status")] <- paste0(
  names(paired)[names(paired) %in% c("logLik", "logLik_se", "B0_hat", "sigma_beta_hat", "status")],
  "_gamma"
)
paired <- merge(
  paired,
  constant$best[, c("task_id", "logLik", "logLik_se", "Beta_hat", "status")],
  by = "task_id",
  suffixes = c("", "_constant")
)
names(paired)[names(paired) == "logLik"] <- "logLik_constant"
names(paired)[names(paired) == "logLik_se"] <- "logLik_se_constant"
names(paired)[names(paired) == "Beta_hat"] <- "Beta_hat_constant"
names(paired)[names(paired) == "status"] <- "status_constant"

paired$delta_RSS_gamma_minus_constant <- paired$RSS_gamma - paired$RSS_constant
paired$delta_RMSE_gamma_minus_constant <- paired$RMSE_gamma - paired$RMSE_constant
paired$delta_AOB_gamma_minus_constant <- paired$AOB_gamma - paired$AOB_constant
paired$delta_logLik_gamma_minus_constant <- paired$logLik_gamma - paired$logLik_constant
paired$gamma_lower_RSS <- paired$RSS_gamma < paired$RSS_constant
paired$gamma_lower_RMSE <- paired$RMSE_gamma < paired$RMSE_constant
paired$gamma_lower_AOB <- paired$AOB_gamma < paired$AOB_constant
write.csv(paired, file.path(output_dir, "paired_model_comparison.csv"), row.names = FALSE)

overall <- data.frame(
  quantity = c(
    "RSS", "RMSE", "AOB (absolute task-level mean error)",
    "mean error through week 5", "mean error after week 5",
    "independent log-likelihood", "Gamma-noise model win proportion by RSS",
    "Gamma-noise model win proportion by RMSE",
    "Gamma-noise model win proportion by AOB"
  ),
  gamma_mean_or_proportion = c(
    mean(paired$RSS_gamma), mean(paired$RMSE_gamma), mean(paired$AOB_gamma),
    mean(paired$mean_error_through_5_gamma),
    mean(paired$mean_error_after_5_gamma), mean(paired$logLik_gamma),
    mean(paired$gamma_lower_RSS), mean(paired$gamma_lower_RMSE),
    mean(paired$gamma_lower_AOB)
  ),
  constant_mean = c(
    mean(paired$RSS_constant), mean(paired$RMSE_constant), mean(paired$AOB_constant),
    mean(paired$mean_error_through_5_constant),
    mean(paired$mean_error_after_5_constant), mean(paired$logLik_constant),
    NA_real_, NA_real_, NA_real_
  ),
  median_paired_difference_gamma_minus_constant = c(
    median(paired$delta_RSS_gamma_minus_constant),
    median(paired$delta_RMSE_gamma_minus_constant),
    median(paired$delta_AOB_gamma_minus_constant),
    median(paired$mean_error_through_5_gamma -
           paired$mean_error_through_5_constant),
    median(paired$mean_error_after_5_gamma -
           paired$mean_error_after_5_constant),
    median(paired$delta_logLik_gamma_minus_constant),
    NA_real_, NA_real_, NA_real_
  ),
  stringsAsFactors = FALSE
)
write.csv(overall, file.path(output_dir, "overall_model_comparison.csv"), row.names = FALSE)

# ============================================================
# Figures
# ============================================================
# ============================================================
# Figures
# ============================================================
#
# Figure style follows conventional epidemiological/applied-math
# graphics: monochrome/grayscale, thin axes, no decorative grid,
# shared scales for direct comparisons, and minimal in-panel text.
# Captions in the report should carry most interpretation.

set_figure_par <- function(mar = c(4.2, 4.5, 0.5, 0.5)) {
  par(
    mar = mar,
    family = "Helvetica",
    las = 1,
    mgp = c(2.5, 0.7, 0),
    tcl = -0.25,
    cex.axis = 0.90,
    cex.lab = 1.00,
    bty = "o"
  )
}

# ------------------------------------------------------------
# Figures 1 and 8 show illustrative sampled Gamma-noise latent trajectories
# for tasks 1 and 117. They are not inputs to the primary recovery metrics.
# This block validates their source data and confirms that the corresponding
# PDFs exist.
# Pilot post-processing disables this validation because selected-task
# figures are not part of the five-task pilot summary.
# ------------------------------------------------------------
if (include_selected_trajectory_figure) {
  selected_B_files <- file.path(
    "results", "selected_trajectory",
    c(
      "experiment4_task1_B_trajectory_comparison.csv",
      "experiment4_task117_B_trajectory_comparison.csv"
    )
  )
  selected_figure_files <- c(
    file.path(figures_dir, "01_selected_task_B_trajectory_comparison.pdf"),
    file.path(figures_dir, "08_task117_B_trajectory_comparison.pdf")
  )
  required_selected_files <- c(selected_B_files, selected_figure_files)
  if (any(!file.exists(required_selected_files))) {
    stop(
      "Missing selected-trajectory artifact(s). Run the task-1 and task-117 ",
      "comparison-figure scripts first."
    )
  }
  required_B_columns <- c(
    "experiment", "task_id", "week", "B_true",
    "B_gamma_sampled_trajectory", "B_constant", "is_time_zero"
  )
  for (j in seq_along(selected_B_files)) {
    selected_B <- read.csv(selected_B_files[[j]], check.names = FALSE)
    expected_task <- c(1L, 117L)[[j]]
    valid_selected <-
      all(required_B_columns %in% names(selected_B)) &&
      nrow(selected_B) == 71L &&
      identical(unique(selected_B$task_id), expected_task) &&
      !is.unsorted(selected_B$week, strictly = TRUE) &&
      !anyDuplicated(selected_B$week) &&
      isTRUE(selected_B$is_time_zero[[1]]) &&
      !any(selected_B$is_time_zero[-1L]) &&
      all(is.finite(unlist(selected_B[c(
        "week", "B_true", "B_gamma_sampled_trajectory", "B_constant"
      )]))) &&
      length(unique(selected_B$B_constant)) == 1L
    if (!valid_selected) {
      stop("The selected task-", expected_task, " trajectory failed validation.")
    }
  }
}

# ------------------------------------------------------------
# Figure 2: RSS distributions in vertically stacked panels.
# Filled grayscale bars avoid the visually awkward overlapping
# hollow histograms and both panels share exactly the same x scale.
# ------------------------------------------------------------

rss_max <- max(260, ceiling(max(c(paired$RSS_gamma, paired$RSS_constant), na.rm = TRUE) / 20) * 20)
rss_breaks <- seq(0, rss_max, by = 20)

cairo_pdf(
  file.path(figures_dir, "02_RSS_distributions.pdf"),
  width = 6.15, height = 4.60, pointsize = 10, family = "Helvetica"
)
layout(matrix(1:2, ncol = 1))
set_figure_par(c(2.0, 4.5, 0.5, 0.5))
hist(
  paired$RSS_gamma,
  breaks = rss_breaks,
  col = "grey82",
  border = "black",
  xlim = c(0, rss_max),
  xaxt = "n",
  xlab = "",
  ylab = "Frequency",
  main = ""
)
text(rss_max * 0.98, par("usr")[[4]] * 0.86, "Gamma-noise model", adj = c(1, 0.5), cex = 0.88)

set_figure_par(c(4.2, 4.5, 0.5, 0.5))
hist(
  paired$RSS_constant,
  breaks = rss_breaks,
  col = "grey82",
  border = "black",
  xlim = c(0, rss_max),
  xlab = "Observation-time point-estimator RSS",
  ylab = "Frequency",
  main = ""
)
text(rss_max * 0.98, par("usr")[[4]] * 0.86, "Constant-B", adj = c(1, 0.5), cex = 0.88)
dev.off()

# ------------------------------------------------------------
# Figure 3: signed mean-error distributions in a 2 x 3 panel figure.
# Columns are overall, through week 5, and after week 5.
# Rows are Gamma-noise model and constant-B. The dashed vertical
# line is zero error; the dotted line is the across-task mean error.
# ------------------------------------------------------------

mean_error_sets <- list(
  list(
    title = "Overall",
    gamma = paired$mean_error_gamma,
    constant = paired$mean_error_constant
  ),
  list(
    title = "Through week 5",
    gamma = paired$mean_error_through_5_gamma,
    constant = paired$mean_error_through_5_constant
  ),
  list(
    title = "After week 5",
    gamma = paired$mean_error_after_5_gamma,
    constant = paired$mean_error_after_5_constant
  )
)
all_mean_errors <- unlist(lapply(
  mean_error_sets, function(z) c(z$gamma, z$constant)
))
mean_error_min <- min(
  -1.25, floor(min(all_mean_errors, na.rm = TRUE) * 4) / 4
)
mean_error_max <- max(
  2.75, ceiling(max(all_mean_errors, na.rm = TRUE) * 4) / 4
)
mean_error_breaks <- seq(mean_error_min, mean_error_max, by = 0.25)

cairo_pdf(
  file.path(figures_dir, "03_mean_error_distributions.pdf"),
  width = 7.15, height = 4.55, pointsize = 10, family = "Helvetica"
)
par(mfrow = c(2, 3), oma = c(0.5, 0.5, 0.2, 2.0))
for (row_name in c("gamma", "constant")) {
  for (j in seq_along(mean_error_sets)) {
    z <- mean_error_sets[[j]]
    values <- z[[row_name]]
    set_figure_par(c(if (row_name == "constant") 4.0 else 1.7, if (j == 1) 3.8 else 2.0, 1.8, 0.5))
    hist(
      values,
      breaks = mean_error_breaks,
      col = "grey82",
      border = "black",
      xlim = c(mean_error_min, mean_error_max),
      xaxt = if (row_name == "gamma") "n" else "s",
      xlab = if (row_name == "constant") "Mean estimation error" else "",
      ylab = if (j == 1) "Frequency" else "",
      main = if (row_name == "gamma") z$title else ""
    )
    abline(v = 0, lty = 2, lwd = 0.8)
    abline(v = mean(values, na.rm = TRUE), lty = 3, lwd = 0.9, col = "grey35")
    if (j == 3) {
      mtext(
        if (row_name == "gamma") "Gamma-noise model" else "Constant-B",
        side = 4,
        line = 0.7,
        las = 3,
        cex = 0.85
      )
    }
  }
}
dev.off()

# ------------------------------------------------------------
# Figure 4: paired RMSE scatter.
# Open circles and the 1:1 line follow a conventional diagnostic
# style. Points below the diagonal have lower RMSE under Gamma.
# ------------------------------------------------------------

all_rmse <- c(paired$RMSE_gamma, paired$RMSE_constant)
rmse_limits <- c(
  floor(min(all_rmse, na.rm = TRUE) * 10) / 10,
  ceiling(max(all_rmse, na.rm = TRUE) * 10) / 10
)
cairo_pdf(
  file.path(figures_dir, "04_paired_RMSE_scatter.pdf"),
  width = 4.85, height = 4.75, pointsize = 10, family = "Helvetica"
)
set_figure_par(c(4.2, 4.5, 0.5, 0.5))
plot(
  paired$RMSE_constant,
  paired$RMSE_gamma,
  xlim = rmse_limits,
  ylim = rmse_limits,
  xaxs = "i",
  yaxs = "i",
  asp = 1,
  pch = 1,
  cex = 0.75,
  xlab = "Repeated-static-estimate RMSE",
  ylab = "Particle-filtering-mean RMSE"
)
abline(0, 1, lty = 2, lwd = 0.9)
dev.off()

# ------------------------------------------------------------
# Figure 5: RMSE distributions, again using vertically stacked
# filled histograms with a common x scale.
# ------------------------------------------------------------

rmse_breaks <- seq(rmse_limits[[1]], rmse_limits[[2]], by = 0.10)
cairo_pdf(
  file.path(figures_dir, "05_RMSE_distributions.pdf"),
  width = 6.15, height = 4.60, pointsize = 10, family = "Helvetica"
)
layout(matrix(1:2, ncol = 1))
set_figure_par(c(2.0, 4.5, 0.5, 0.5))
hist(
  paired$RMSE_gamma,
  breaks = rmse_breaks,
  col = "grey82",
  border = "black",
  xlim = rmse_limits,
  xaxt = "n",
  xlab = "",
  ylab = "Frequency",
  main = ""
)
text(
  rmse_limits[[2]] - 0.01 * diff(rmse_limits),
  par("usr")[[4]] * 0.86,
  "Gamma-noise model", adj = c(1, 0.5), cex = 0.88
)

set_figure_par(c(4.2, 4.5, 0.5, 0.5))
hist(
  paired$RMSE_constant,
  breaks = rmse_breaks,
  col = "grey82",
  border = "black",
  xlim = rmse_limits,
  xlab = "Observation-time point-estimator RMSE",
  ylab = "Frequency",
  main = ""
)
text(
  rmse_limits[[2]] - 0.01 * diff(rmse_limits),
  par("usr")[[4]] * 0.86,
  "Constant-B", adj = c(1, 0.5), cex = 0.88
)
dev.off()

# ------------------------------------------------------------
# Figure 6: descriptive independent-likelihood difference.
# ------------------------------------------------------------

lik_breaks <- seq(
  floor(min(paired$delta_logLik_gamma_minus_constant, na.rm = TRUE) / 3) * 3,
  ceiling(max(paired$delta_logLik_gamma_minus_constant, na.rm = TRUE) / 3) * 3,
  by = 3
)
cairo_pdf(
  file.path(figures_dir, "06_independent_loglik_difference.pdf"),
  width = 6.15, height = 4.00, pointsize = 10, family = "Helvetica"
)
set_figure_par(c(4.4, 4.5, 0.5, 0.5))
hist(
  paired$delta_logLik_gamma_minus_constant,
  breaks = lik_breaks,
  col = "grey82",
  border = "black",
  xlab = "Independent log likelihood difference (Gamma-noise - constant-B)",
  ylab = "Frequency",
  main = ""
)
abline(v = 0, lty = 2, lwd = 0.9)
dev.off()

cat("Paired model comparison complete. Outputs: ", output_dir, " and ", figures_dir, "\n", sep = "")

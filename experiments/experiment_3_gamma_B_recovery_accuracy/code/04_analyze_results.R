# ============================================================
# Post-process Experiment 3 combined results
#
# Reads the existing combined CSV files only. No MIF2 reruns.
# Produces publication-ready summary figures and updated CSV
# summaries from the stored combined results.
#
# Usage from the Experiment 3 directory:
#   Rscript code/04_analyze_results.R
# ============================================================

options(stringsAsFactors = FALSE)

results_dir <- file.path("results", "combined")
figures_dir <- "figures"
write_summary_csvs <- !identical(Sys.getenv("EXP3_FIGURES_ONLY"), "1")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
selected_trajectory_file <- file.path(
  "results", "selected_trajectory", "experiment3_task145_B_trajectory.csv"
)

best <- read.csv(file.path(results_dir, "combined_best_fit_summary.csv"))
paths <- read.csv(file.path(results_dir, "combined_filtered_B_paths.csv"))
runs <- read.csv(file.path(results_dir, "combined_mif2_results.csv"))
selected_trajectory <- read.csv(selected_trajectory_file, check.names = FALSE)

required_best_columns <- c(
  "task_id", "best_run", "logLik", "fit_success", "final_pf_success"
)
required_path_columns <- c("task_id", "week", "B_filtered_mean", "B_true")
required_run_columns <- c("task_id", "run", "logLik")
required_columns <- list(
  combined_best_fit_summary = required_best_columns,
  combined_filtered_B_paths = required_path_columns,
  combined_mif2_results = required_run_columns
)
values <- list(
  combined_best_fit_summary = best,
  combined_filtered_B_paths = paths,
  combined_mif2_results = runs
)
for (name in names(required_columns)) {
  missing_columns <- setdiff(required_columns[[name]], names(values[[name]]))
  if (length(missing_columns) > 0) {
    stop(name, ".csv is missing column(s): ", paste(missing_columns, collapse = ", "))
  }
}

expected_tasks <- seq_len(200L)
if (!identical(sort(unique(best$task_id)), expected_tasks)) {
  stop("combined_best_fit_summary.csv must contain tasks 1 through 200 exactly once.")
}
if (anyDuplicated(best$task_id) || nrow(best) != length(expected_tasks)) {
  stop("combined_best_fit_summary.csv must contain one row per task.")
}
if (!all(best$fit_success & best$final_pf_success)) {
  stop("The stored recovery summary contains an unsuccessful selected fit or final particle filter.")
}
if (!identical(sort(unique(paths$task_id)), expected_tasks) ||
    !all(table(paths$task_id) == 70L)) {
  stop("combined_filtered_B_paths.csv must contain 70 rows for each task 1 through 200.")
}
if (!identical(sort(unique(runs$task_id)), expected_tasks) ||
    !all(table(runs$task_id) == 9L)) {
  stop("combined_mif2_results.csv must contain nine rows for each task 1 through 200.")
}
best_from_runs <- do.call(
  rbind,
  lapply(split(runs[is.finite(runs$logLik), , drop = FALSE], runs$task_id), function(x) {
    x[which.max(x$logLik), c("task_id", "run", "logLik"), drop = FALSE]
  })
)
best_check <- merge(
  best[, c("task_id", "best_run", "logLik")],
  best_from_runs,
  by = "task_id",
  suffixes = c("_summary", "_runs"),
  sort = TRUE
)
if (nrow(best_check) != 200L ||
    any(best_check$best_run != best_check$run) ||
    any(abs(best_check$logLik_summary - best_check$logLik_runs) > 1e-12)) {
  stop("Stored selected fits do not match the largest evaluated likelihood within each task.")
}
if (any(!is.finite(paths$B_filtered_mean)) || any(!is.finite(paths$B_true))) {
  stop("The stored filtered paths contain non-finite transmission-rate values.")
}
if (any(paths$B_true[paths$week < 5] != 4) ||
    any(paths$B_true[paths$week >= 5] != 2)) {
  stop("The stored true transmission-rate path does not match the documented week-5 switch.")
}
required_selected_columns <- c(
  "experiment", "task_id", "week", "B_trajectory", "B_true", "is_time_zero"
)
if (!all(required_selected_columns %in% names(selected_trajectory)) ||
    nrow(selected_trajectory) != 71L ||
    !identical(unique(selected_trajectory$task_id), 145L) ||
    any(!is.finite(selected_trajectory$week)) ||
    any(!is.finite(selected_trajectory$B_trajectory)) ||
    is.unsorted(selected_trajectory$week, strictly = TRUE) ||
    anyDuplicated(selected_trajectory$week)) {
  stop("The selected task-145 trajectory artifact failed validation.")
}

paths$error <- paths$B_filtered_mean - paths$B_true
paths$squared_error <- paths$error^2
paths$period <- ifelse(paths$week < 5, "before", "after")

split_paths <- split(paths, paths$task_id)
error_summary <- do.call(
  rbind,
  lapply(split_paths, function(x) {
    data.frame(
      task_id = x$task_id[[1]],
      RSS = sum(x$squared_error),
      RMSE = sqrt(mean(x$squared_error)),
      bias_all = mean(x$error),
      bias_before = mean(x$error[x$week < 5]),
      bias_after = mean(x$error[x$week >= 5])
    )
  })
)
rownames(error_summary) <- NULL

successful_runs <- runs[is.finite(runs$logLik), , drop = FALSE]
split_runs <- split(successful_runs, successful_runs$task_id)
likelihood_gaps <- do.call(
  rbind,
  lapply(split_runs, function(x) {
    x <- x[order(x$logLik, decreasing = TRUE), , drop = FALSE]
    data.frame(
      task_id = x$task_id[[1]],
      best_run = x$run[[1]],
      best_logLik = x$logLik[[1]],
      second_best_logLik = if (nrow(x) >= 2) x$logLik[[2]] else NA_real_,
      logLik_gap = if (nrow(x) >= 2) x$logLik[[1]] - x$logLik[[2]] else NA_real_
    )
  })
)
rownames(likelihood_gaps) <- NULL

if (write_summary_csvs) {
  write.csv(error_summary, file.path(results_dir, "task_level_error_summary.csv"), row.names = FALSE)
  write.csv(likelihood_gaps, file.path(results_dir, "likelihood_gaps_by_task.csv"), row.names = FALSE)
}

summary_table <- data.frame(
  quantity = c("RSS", "RMSE", "bias_all", "bias_before", "bias_after", "logLik_gap"),
  mean = c(
    mean(error_summary$RSS),
    mean(error_summary$RMSE),
    mean(error_summary$bias_all),
    mean(error_summary$bias_before),
    mean(error_summary$bias_after),
    mean(likelihood_gaps$logLik_gap, na.rm = TRUE)
  ),
  median = c(
    median(error_summary$RSS),
    median(error_summary$RMSE),
    median(error_summary$bias_all),
    median(error_summary$bias_before),
    median(error_summary$bias_after),
    median(likelihood_gaps$logLik_gap, na.rm = TRUE)
  ),
  sd = c(
    sd(error_summary$RSS),
    sd(error_summary$RMSE),
    sd(error_summary$bias_all),
    sd(error_summary$bias_before),
    sd(error_summary$bias_after),
    sd(likelihood_gaps$logLik_gap, na.rm = TRUE)
  )
)
if (write_summary_csvs) {
  write.csv(summary_table, file.path(results_dir, "overall_summary.csv"), row.names = FALSE)
}

# ----------------------------
# Plotting helpers
# ----------------------------
set_clean_plot_style <- function() {
  par(
    bty = "l",
    las = 1,
    lend = "butt",
    tcl = -0.25,
    mgp = c(2.0, 0.6, 0),
    mar = c(3.6, 3.9, 1.0, 0.8),
    cex.axis = 0.95,
    cex.lab = 1.1
  )
}

# ----------------------------
# Figure 1: one prespecified ancestry-preserving particle trajectory.
# This is not a filtering mean and is not aggregated across tasks.
# ----------------------------
truth_x <- c(0, 5, 5, 10)
truth_y <- c(4, 4, 2, 2)
selected_ylim <- range(c(0, 6, selected_trajectory$B_trajectory), finite = TRUE)
selected_ylim <- selected_ylim + c(-1, 1) * 0.04 * diff(selected_ylim)
pdf(file.path(figures_dir, "01_selected_task_B_trajectory.pdf"),
    width = 6.2, height = 4.2, useDingbats = FALSE)
set_clean_plot_style()
plot(
  selected_trajectory$week,
  selected_trajectory$B_trajectory,
  type = "l",
  lwd = 1.45,
  lty = 1,
  col = "#0072B2",
  xlim = c(0, 10),
  ylim = selected_ylim,
  xlab = "Week",
  ylab = "Transmission rate B(t)",
  xaxt = "n",
  yaxt = "n"
)
axis(1, at = seq(0, 10, by = 2))
axis(2)
abline(v = 5, col = "gray70", lty = 2, lwd = 1)
lines(truth_x, truth_y, lwd = 1.65, lty = 2, col = "black")
lines(selected_trajectory$week, selected_trajectory$B_trajectory,
      lwd = 1.45, lty = 1, col = "#0072B2")
legend(
  "topright",
  legend = c("True B(t)", "Gamma-noise trajectory (task 145)"),
  lwd = c(1.65, 1.45),
  lty = c(2, 1),
  col = c("black", "#0072B2"),
  bty = "n",
  x.intersp = 0.8,
  seg.len = 3
)
dev.off()

# ----------------------------
# Figure 2: RSS distribution
# ----------------------------
pdf(file.path(figures_dir, "rss_distribution.pdf"), width = 6.2, height = 4.2)
set_clean_plot_style()
hist(
  error_summary$RSS,
  breaks = seq(0, max(ceiling(error_summary$RSS / 5) * 5, 60), by = 5),
  col = "gray85",
  border = "gray20",
  main = "",
  xlab = "Observation-time filtering-mean RSS",
  ylab = "Frequency"
)
dev.off()

# ----------------------------
# Figure 3: RMSE distribution
# ----------------------------
pdf(file.path(figures_dir, "rmse_distribution.pdf"), width = 6.2, height = 4.2)
set_clean_plot_style()
hist(
  error_summary$RMSE,
  breaks = seq(0.3, max(1.1, ceiling(max(error_summary$RMSE) * 20) / 20), by = 0.05),
  col = "gray85",
  border = "gray20",
  main = "",
  xlab = "Observation-time filtering-mean RMSE",
  ylab = "Frequency"
)
dev.off()

# ----------------------------
# Figure 4: bias distributions
# ----------------------------
pdf(file.path(figures_dir, "bias_before_after.pdf"), width = 7.4, height = 3.2)
par(mfrow = c(1, 3))
set_clean_plot_style()
par(mar = c(3.6, 3.9, 1.5, 0.8))
bias_limits <- range(c(error_summary$bias_all, error_summary$bias_before, error_summary$bias_after, -1.2, 1.2))
bias_breaks <- seq(floor(bias_limits[1] / 0.2) * 0.2, ceiling(bias_limits[2] / 0.2) * 0.2, by = 0.2)
for (j in seq_along(list(error_summary$bias_all, error_summary$bias_before, error_summary$bias_after))) {
  current_data <- list(error_summary$bias_all, error_summary$bias_before, error_summary$bias_after)[[j]]
  current_title <- c("Overall", "Before week 5", "Week 5 onward")[j]
  hist(
    current_data,
    breaks = bias_breaks,
    col = "gray85",
    border = "gray20",
    main = "",
    xlab = "Mean estimation error",
    ylab = if (j == 1) "Frequency" else "",
    xlim = range(bias_breaks)
  )
  abline(v = 0, lty = 2, lwd = 1)
  abline(v = mean(current_data), lty = 3, lwd = 1, col = "gray35")
  title(current_title, cex.main = 0.9, font.main = 1, line = 0.35)
}
dev.off()

# ----------------------------
# Figure 5: likelihood-gap distribution
# ----------------------------
pdf(file.path(figures_dir, "likelihood_gap_distribution.pdf"), width = 6.2, height = 4.2)
set_clean_plot_style()
hist(
  likelihood_gaps$logLik_gap,
  breaks = seq(0, max(0.15, ceiling(max(likelihood_gaps$logLik_gap, na.rm = TRUE) * 100) / 100), by = 0.005),
  col = "gray85",
  border = "gray20",
  main = "",
  xlab = "Best - second-best evaluated log likelihood",
  ylab = "Frequency"
)
dev.off()

cat("Post-processing complete.\n")

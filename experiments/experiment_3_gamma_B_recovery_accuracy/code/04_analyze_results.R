# ============================================================
# Post-process Experiment 3 combined results
#
# Reads the existing combined CSV files only. No MIF2 reruns.
# Produces RSS, RMSE, bias summaries, likelihood gaps, and PDFs.
# Time is stored in weeks; observations are spaced 1/7 week apart.
#
# Usage from the Experiment 3 directory:
#   Rscript code/04_analyze_results.R
# ============================================================

options(stringsAsFactors = FALSE)

results_dir <- file.path("results", "combined")
figures_dir <- "figures"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

best <- read.csv(file.path(results_dir, "combined_best_fit_summary.csv"))
paths <- read.csv(file.path(results_dir, "combined_filtered_B_paths.csv"))
runs <- read.csv(file.path(results_dir, "combined_mif2_results.csv"))

required_path_columns <- c("task_id", "week", "B_filtered_mean", "B_true")
if (!all(required_path_columns %in% names(paths))) {
  stop("combined_filtered_B_paths.csv is missing required columns.")
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

write.csv(
  error_summary,
  file.path(results_dir, "task_level_error_summary.csv"),
  row.names = FALSE
)
write.csv(
  likelihood_gaps,
  file.path(results_dir, "likelihood_gaps_by_task.csv"),
  row.names = FALSE
)

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
write.csv(
  summary_table,
  file.path(results_dir, "overall_summary.csv"),
  row.names = FALSE
)

pdf(file.path(figures_dir, "rss_distribution.pdf"), width = 7, height = 5)
hist(
  error_summary$RSS,
  breaks = 20,
  main = "Distribution of transmission-path RSS",
  xlab = "RSS across 70 observation times"
)
dev.off()

pdf(file.path(figures_dir, "bias_before_after.pdf"), width = 7, height = 5)
boxplot(
  error_summary$bias_before,
  error_summary$bias_after,
  names = c("Before change", "After change"),
  ylab = "Mean filtered-path error",
  main = "Bias before and after the transmission-rate change"
)
abline(h = 0, lty = 2)
dev.off()

mean_path <- aggregate(
  cbind(B_filtered_mean, B_true) ~ week,
  data = paths,
  FUN = mean
)

pdf(file.path(figures_dir, "mean_filtered_B_path.pdf"), width = 7, height = 5)
plot(
  mean_path$week,
  mean_path$B_filtered_mean,
  type = "l",
  lwd = 2,
  xlab = "Week",
  ylab = "Transmission rate",
  main = "Mean filtered transmission-rate path across 200 datasets"
)
lines(mean_path$week, mean_path$B_true, lty = 2, lwd = 2)
legend(
  "topright",
  legend = c("Mean filtered B", "True B"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)
dev.off()

cat("Post-processing complete.\n")

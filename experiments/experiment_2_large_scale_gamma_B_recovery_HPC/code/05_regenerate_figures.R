#!/usr/bin/env Rscript

# Rebuild Experiment 2 figures from committed compact results.
# This script does not simulate data, run MIF2, or run a particle filter.

library(pomp)

script_args <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_flag) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_flag[[1]]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
experiment_dir <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir

input_paths <- c(
  combined = file.path(experiment_dir, "results", "combined_mif2_results.csv"),
  best = file.path(experiment_dir, "results", "best_fit.csv"),
  B_path = file.path(experiment_dir, "results", "final_filtered_B_path.csv"),
  I_path = file.path(experiment_dir, "results", "final_filtered_infectious_path.csv"),
  best_object = file.path(experiment_dir, "results", "best_mif2_object.rds")
)
missing_paths <- input_paths[!file.exists(input_paths)]
if (length(missing_paths) > 0) stop("Missing required inputs: ", paste(missing_paths, collapse = ", "))

combined <- read.csv(input_paths[["combined"]], check.names = FALSE)
best <- read.csv(input_paths[["best"]], check.names = FALSE)
B_path <- read.csv(input_paths[["B_path"]], check.names = FALSE)
I_path <- read.csv(input_paths[["I_path"]], check.names = FALSE)

required_columns <- list(
  combined = c("task_id", "start_B0", "start_sigma_beta", "logLik", "logLik_se"),
  best = "task_id",
  B_path = c("week", "B_true", "B_filtered_mean"),
  I_path = c("week", "gamma_infectious", "true_infectious")
)
values <- list(combined = combined, best = best, B_path = B_path, I_path = I_path)
for (name in names(required_columns)) {
  missing_columns <- setdiff(required_columns[[name]], names(values[[name]]))
  if (length(missing_columns) > 0) {
    stop("Missing columns in ", name, ": ", paste(missing_columns, collapse = ", "))
  }
}
if (nrow(best) != 1) stop("Expected exactly one selected best fit.")
if (!best$task_id[[1]] %in% combined$task_id) stop("Selected task is absent from combined results.")

figures_dir <- file.path(experiment_dir, "figures")
results_dir <- file.path(experiment_dir, "results")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

combined$is_selected_best <- combined$task_id == best$task_id[[1]]
starting_summary <- combined[order(-combined$logLik), , drop = FALSE]
write.csv(starting_summary, file.path(results_dir, "starting_value_summary.csv"), row.names = FALSE)

open_pdf <- function(filename, width, height, bottom = 4.3, left = 4.5) {
  grDevices::pdf(file.path(figures_dir, filename), width = width, height = height)
  graphics::par(
    mar = c(bottom, left, 0.8, 0.8),
    mgp = c(2.7, 0.8, 0),
    tcl = -0.25,
    las = 1,
    bty = "l",
    family = "sans"
  )
}

open_pdf("best_B_path.pdf", 6.2, 4.2)
graphics::plot(
  B_path$week,
  B_path$B_filtered_mean,
  type = "n",
  xlim = c(0, 10),
  ylim = c(0, 6),
  xlab = "Week",
  ylab = "Transmission rate, B(t)",
  axes = FALSE
)
graphics::axis(1, at = seq(0, 10, by = 2))
graphics::axis(2, at = seq(0, 6, by = 1))
graphics::segments(c(0, 5), c(4, 2), c(5, 10), c(4, 2), lwd = 1.5)
graphics::abline(v = 5, lty = 3, col = "gray75")
graphics::lines(B_path$week, B_path$B_filtered_mean, col = "#56B4E9", lwd = 1.2)
graphics::legend(
  "topright",
  legend = c("True B(t)", "Gamma-noise fit"),
  col = c("black", "#56B4E9"),
  lwd = c(1.5, 1.2),
  bty = "n",
  cex = 0.9
)
grDevices::dev.off()

y_limit <- range(c(I_path$gamma_infectious, I_path$true_infectious), finite = TRUE)
open_pdf("best_infectious_path.pdf", 6.2, 4.2)
graphics::plot(
  I_path$week,
  I_path$true_infectious,
  type = "l",
  xlim = c(0, 10),
  ylim = c(0, y_limit[[2]] * 1.05),
  xlab = "Week",
  ylab = "Infectious individuals",
  axes = FALSE,
  lwd = 1.3
)
graphics::axis(1, at = seq(0, 10, by = 2))
graphics::axis(2)
graphics::abline(v = 5, lty = 3, col = "gray75")
graphics::lines(I_path$week, I_path$gamma_infectious, col = "#56B4E9", lwd = 1.1)
graphics::legend(
  "topright",
  legend = c("True I(t)", "Gamma-noise fit"),
  col = c("black", "#56B4E9"),
  lwd = c(1.3, 1.1),
  bty = "n",
  cex = 0.9
)
grDevices::dev.off()

plot_values <- combined[order(combined$start_sigma_beta, combined$start_B0), , drop = FALSE]
x <- seq_len(nrow(plot_values))
labels <- paste0(
  "B0=", plot_values$start_B0,
  ", sigma=", sprintf("%.2f", plot_values$start_sigma_beta)
)
lower <- plot_values$logLik - 2 * plot_values$logLik_se
upper <- plot_values$logLik + 2 * plot_values$logLik_se
y_plot <- range(c(lower, upper), finite = TRUE)
padding <- max(diff(y_plot) * 0.08, 0.01)

open_pdf("starting_value_loglik.pdf", 6.6, 4.5, bottom = 7.2, left = 5.4)
graphics::plot(
  x,
  plot_values$logLik,
  type = "n",
  ylim = y_plot + c(-padding, padding),
  xlab = "",
  ylab = "",
  axes = FALSE
)
graphics::axis(1, at = x, labels = labels, las = 2, cex.axis = 0.75)
graphics::axis(2)
graphics::mtext("Evaluated log likelihood", side = 2, line = 4.1, las = 0)
graphics::arrows(x, lower, x, upper, angle = 90, code = 3, length = 0.04)
graphics::points(
  x,
  plot_values$logLik,
  pch = ifelse(plot_values$task_id == best$task_id[[1]], 19, 1),
  cex = 0.9
)
graphics::mtext("MIF2 starting values", side = 1, line = 5.9)
grDevices::dev.off()

mif_best <- readRDS(input_paths[["best_object"]])
grDevices::pdf(
  file.path(figures_dir, "best_mif2_trace.pdf"),
  width = 8,
  height = 5.5,
  useDingbats = FALSE
)
plot(mif_best, pars = c("B0", "sigma_beta"))
grDevices::dev.off()

message("Regenerated four Experiment 2 PDFs from committed compact results.")

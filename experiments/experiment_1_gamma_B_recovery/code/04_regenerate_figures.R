#!/usr/bin/env Rscript

# Rebuild Experiment 1 report-facing figures from compact committed CSV files.
# This script does not simulate data, run MIF2, or run particle filters.

script_args <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_flag) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_flag[[1]]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

experiment_dir <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
results_dir <- file.path(experiment_dir, "results")
figures_dir <- file.path(experiment_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

read_required_csv <- function(filename, required_columns) {
  path <- file.path(results_dir, filename)
  if (!file.exists(path)) stop("Required input is missing: ", path)
  value <- read.csv(path, check.names = FALSE)
  missing_columns <- setdiff(required_columns, names(value))
  if (length(missing_columns) > 0) {
    stop("Missing columns in ", filename, ": ", paste(missing_columns, collapse = ", "))
  }
  value
}

open_pdf <- function(filename, width, height) {
  grDevices::pdf(file.path(figures_dir, filename), width = width, height = height)
  graphics::par(
    mar = c(4.3, 4.4, 0.8, 0.8),
    mgp = c(2.6, 0.8, 0),
    tcl = -0.25,
    las = 1,
    bty = "l",
    family = "sans"
  )
}

plot_B_path <- function(prefix, piecewise) {
  value <- read_required_csv(
    paste0(prefix, "_filtered_B_path.csv"),
    c("week", "B_filtered_mean", "B_true")
  )

  open_pdf(paste0(prefix, "_filtered_B_path.pdf"), 6.2, 4.2)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(
    value$week,
    value$B_filtered_mean,
    type = "n",
    xlim = c(0, 10),
    ylim = c(0, 6),
    xlab = "Week",
    ylab = "Transmission rate, B(t)",
    axes = FALSE
  )
  graphics::axis(1, at = seq(0, 10, by = 2))
  graphics::axis(2, at = seq(0, 6, by = 1))
  if (piecewise) {
    graphics::segments(c(0, 5), c(4, 2), c(5, 10), c(4, 2), lwd = 1.5)
    graphics::abline(v = 5, lty = 3, col = "gray75")
  } else {
    graphics::segments(0, 4, 10, 4, lwd = 1.5)
  }
  graphics::lines(value$week, value$B_filtered_mean, col = "#56B4E9", lwd = 1.2)
  graphics::legend(
    "topright",
    legend = c("True B(t)", "Gamma-noise fit"),
    col = c("black", "#56B4E9"),
    lwd = c(1.5, 1.2),
    bty = "n",
    cex = 0.9
  )
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
}

plot_infectious_path <- function(prefix, piecewise) {
  value <- read_required_csv(
    paste0(prefix, "_filtered_infectious_path.csv"),
    c("week", "gamma_infectious", "true_infectious")
  )
  y_limit <- range(c(value$gamma_infectious, value$true_infectious), finite = TRUE)
  y_limit <- c(0, y_limit[[2]] * 1.05)

  open_pdf(paste0(prefix, "_filtered_infectious_path.pdf"), 6.2, 4.2)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(
    value$week,
    value$true_infectious,
    type = "l",
    xlim = c(0, 10),
    ylim = y_limit,
    xlab = "Week",
    ylab = "Infectious individuals",
    axes = FALSE,
    lwd = 1.3
  )
  graphics::axis(1, at = seq(0, 10, by = 2))
  graphics::axis(2)
  if (piecewise) graphics::abline(v = 5, lty = 3, col = "gray75")
  graphics::lines(value$week, value$gamma_infectious, col = "#56B4E9", lwd = 1.1)
  graphics::legend(
    "topright",
    legend = c("True I(t)", "Gamma-noise fit"),
    col = c("black", "#56B4E9"),
    lwd = c(1.3, 1.1),
    bty = "n",
    cex = 0.9
  )
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
}

plot_starting_values <- function(prefix, fixed_B0 = FALSE) {
  results <- read_required_csv(
    paste0(prefix, "_mif2_results.csv"),
    c("run", "start_sigma_beta", "logLik", "logLik_se")
  )
  best <- read_required_csv(paste0(prefix, "_best_fit.csv"), c("run"))
  if (nrow(best) != 1) stop("Expected exactly one selected fit for ", prefix)

  labels <- if (fixed_B0) {
    paste0("sigma=", sprintf("%.2f", results$start_sigma_beta))
  } else {
    if (!"start_B0" %in% names(results)) stop("Missing start_B0 for ", prefix)
    paste0("B0=", results$start_B0, ", sigma=", sprintf("%.2f", results$start_sigma_beta))
  }
  selected <- results$run == best$run[[1]]
  x <- seq_len(nrow(results))
  lower <- results$logLik - 2 * results$logLik_se
  upper <- results$logLik + 2 * results$logLik_se
  y_limit <- range(c(lower, upper), finite = TRUE)
  padding <- max(diff(y_limit) * 0.08, 0.05)

  grDevices::pdf(file.path(figures_dir, paste0(prefix, "_starting_value_loglik.pdf")), 6.4, 4.4)
  graphics::par(
    mar = c(if (length(labels) > 3) 7.1 else 5.0, 4.5, 0.8, 0.8),
    mgp = c(2.7, 0.8, 0),
    tcl = -0.25,
    las = 1,
    bty = "l",
    family = "sans"
  )
  graphics::plot(
    x,
    results$logLik,
    type = "n",
    ylim = y_limit + c(-padding, padding),
    xlab = "",
    ylab = "Evaluated log likelihood",
    axes = FALSE
  )
  graphics::axis(1, at = x, labels = labels, las = 2, cex.axis = 0.75)
  graphics::axis(2)
  graphics::arrows(x, lower, x, upper, angle = 90, code = 3, length = 0.04)
  graphics::points(x, results$logLik, pch = ifelse(selected, 19, 1), cex = 0.9)
  graphics::mtext("MIF2 starting values", side = 1, line = if (length(labels) > 3) 5.8 else 3.8)
  grDevices::dev.off()
}

scenarios <- data.frame(
  prefix = c("constant_B4", "piecewise_fixed_B0", "piecewise_estimated_B0"),
  piecewise = c(FALSE, TRUE, TRUE),
  fixed_B0 = c(FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(scenarios))) {
  plot_B_path(scenarios$prefix[[i]], scenarios$piecewise[[i]])
  plot_infectious_path(scenarios$prefix[[i]], scenarios$piecewise[[i]])
  plot_starting_values(scenarios$prefix[[i]], scenarios$fixed_B0[[i]])
}

message("Regenerated nine Experiment 1 report-facing PDFs from committed CSV results.")

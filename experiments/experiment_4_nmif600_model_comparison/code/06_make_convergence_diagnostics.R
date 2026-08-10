# ============================================================
# Create convergence diagnostics from traces saved during the
# formal Nmif=600 runs for selected diagnostic tasks.
#
# Graphical convention: all starts are shown in light gray and
# the run with the best independent likelihood is highlighted in
# black. The vertical dashed line marks iteration 500, i.e. the
# beginning of the final 100 MIF iterations.
# ============================================================

options(stringsAsFactors = FALSE)
source("config/experiment_config.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript code/06_make_convergence_diagnostics.R <gamma_combined_dir> <constant_combined_dir> <figures_dir>")
}

gamma_dir <- args[[1]]
constant_dir <- args[[2]]
figures_dir <- args[[3]]
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

read_trace <- function(dir) {
  path <- file.path(dir, "combined_mif2_traces.csv")
  if (!file.exists(path)) return(NULL)
  read.csv(path, check.names = FALSE)
}

read_runs <- function(dir) {
  path <- file.path(dir, "combined_mif2_results.csv")
  if (!file.exists(path)) return(NULL)
  read.csv(path, check.names = FALSE)
}

gamma <- read_trace(gamma_dir)
constant <- read_trace(constant_dir)
gamma_runs <- read_runs(gamma_dir)
constant_runs <- read_runs(constant_dir)

find_loglik_col <- function(x) {
  candidates <- c("loglik", "logLik", ".loglik")
  found <- candidates[candidates %in% names(x)]
  if (length(found) == 0L) NA_character_ else found[[1]]
}

best_run_by_task <- function(runs) {
  if (is.null(runs) || nrow(runs) == 0L) return(integer(0))
  pieces <- split(runs[is.finite(runs$logLik), , drop = FALSE], runs$task_id[is.finite(runs$logLik)])
  out <- vapply(pieces, function(x) x$run[[which.max(x$logLik)]], integer(1))
  out
}

set_diag_par <- function() {
  par(
    family = "sans",
    las = 1,
    mgp = c(2.3, 0.65, 0),
    tcl = -0.22,
    cex.axis = 0.82,
    cex.lab = 0.90,
    bty = "o"
  )
}

plot_trace_panel <- function(x, value_col, ylab, best_run, show_task = FALSE) {
  finite <- is.finite(x[[value_col]]) & is.finite(x$mif_iteration)
  x <- x[finite, , drop = FALSE]
  if (nrow(x) == 0L) {
    plot.new()
    return(invisible(NULL))
  }

  yr <- range(x[[value_col]], finite = TRUE)
  plot(
    range(x$mif_iteration), yr,
    type = "n",
    xlab = "MIF iteration",
    ylab = ylab,
    xaxs = "i"
  )

  runs <- sort(unique(x$run))
  for (run_id in runs[runs != best_run]) {
    z <- x[x$run == run_id, , drop = FALSE]
    z <- z[order(z$mif_iteration), , drop = FALSE]
    lines(z$mif_iteration, z[[value_col]], col = "grey70", lwd = 0.55)
  }
  if (best_run %in% runs) {
    z <- x[x$run == best_run, , drop = FALSE]
    z <- z[order(z$mif_iteration), , drop = FALSE]
    lines(z$mif_iteration, z[[value_col]], col = "black", lwd = 1.15)
  }
  abline(v = experiment_config$Nmif - 100L, lty = 2, lwd = 0.75, col = "grey35")
  if (show_task) {
    text(par("usr")[[1]], par("usr")[[4]], paste("Task", x$task_id[[1]]), adj = c(0.05, 1.15), cex = 0.85)
  }
}

if (!is.null(gamma) && nrow(gamma) > 0L) {
  gamma_loglik <- find_loglik_col(gamma)
  gamma_best <- best_run_by_task(gamma_runs)
  pdf(file.path(figures_dir, "gamma_convergence_diagnostics.pdf"), width = 7.2, height = 2.85, useDingbats = FALSE)
  for (task_id in sort(unique(gamma$task_id))) {
    x <- gamma[gamma$task_id == task_id, , drop = FALSE]
    best <- as.integer(gamma_best[[as.character(task_id)]])
    par(mfrow = c(1, 3), mar = c(3.8, 3.8, 0.5, 0.5))
    set_diag_par()
    if (!is.na(gamma_loglik)) plot_trace_panel(x, gamma_loglik, "Internal log likelihood", best, TRUE) else plot.new()
    set_diag_par(); plot_trace_panel(x, "B0", "B0", best)
    set_diag_par(); plot_trace_panel(x, "sigma_beta", "sigma_beta", best)
    legend("topright", legend = c("Best run", "Other starts"), col = c("black", "grey70"), lwd = c(1.15, 0.55), bty = "n", cex = 0.72)
  }
  dev.off()
}

if (!is.null(constant) && nrow(constant) > 0L) {
  constant_loglik <- find_loglik_col(constant)
  constant_best <- best_run_by_task(constant_runs)
  pdf(file.path(figures_dir, "constant_B_convergence_diagnostics.pdf"), width = 7.2, height = 3.05, useDingbats = FALSE)
  for (task_id in sort(unique(constant$task_id))) {
    x <- constant[constant$task_id == task_id, , drop = FALSE]
    best <- as.integer(constant_best[[as.character(task_id)]])
    par(mfrow = c(1, 2), mar = c(3.8, 3.8, 0.5, 0.5))
    set_diag_par()
    if (!is.na(constant_loglik)) plot_trace_panel(x, constant_loglik, "Internal log likelihood", best, TRUE) else plot.new()
    set_diag_par(); plot_trace_panel(x, "Beta", "Beta", best)
    legend("topright", legend = c("Best run", "Other starts"), col = c("black", "grey70"), lwd = c(1.15, 0.55), bty = "n", cex = 0.75)
  }
  dev.off()
}

calculate_tail_slopes <- function(x, model, variables) {
  if (is.null(x) || nrow(x) == 0L) return(data.frame())
  output <- list()
  index <- 0L
  for (task_id in sort(unique(x$task_id))) {
    for (run_id in sort(unique(x$run[x$task_id == task_id]))) {
      z <- x[x$task_id == task_id & x$run == run_id, , drop = FALSE]
      cutoff <- max(z$mif_iteration, na.rm = TRUE) - 99L
      z <- z[z$mif_iteration >= cutoff, , drop = FALSE]
      for (variable in variables[variables %in% names(z)]) {
        ok <- is.finite(z$mif_iteration) & is.finite(z[[variable]])
        slope <- if (sum(ok) >= 2L) coef(lm(z[[variable]][ok] ~ z$mif_iteration[ok]))[[2]] * 100 else NA_real_
        index <- index + 1L
        output[[index]] <- data.frame(
          model = model, task_id = task_id, run = run_id, variable = variable,
          n_tail_points = sum(ok), slope_per_100_iterations = slope,
          tail_mean = if (any(ok)) mean(z[[variable]][ok]) else NA_real_,
          tail_sd = if (sum(ok) >= 2L) sd(z[[variable]][ok]) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, output)
}

gamma_vars <- if (is.null(gamma)) character(0) else c(find_loglik_col(gamma), "B0", "sigma_beta")
constant_vars <- if (is.null(constant)) character(0) else c(find_loglik_col(constant), "Beta")
tail_parts <- Filter(
  function(x) is.data.frame(x) && nrow(x) > 0L,
  list(
    calculate_tail_slopes(gamma, "gamma_noise", gamma_vars[!is.na(gamma_vars)]),
    calculate_tail_slopes(constant, "constant_B", constant_vars[!is.na(constant_vars)])
  )
)
tail_summary <- if (length(tail_parts) > 0L) do.call(rbind, tail_parts) else data.frame()
if (nrow(tail_summary) > 0L) {
  write.csv(tail_summary, file.path(figures_dir, "convergence_tail_summary.csv"), row.names = FALSE)
}

cat("Convergence diagnostics complete: ", figures_dir, "\n", sep = "")

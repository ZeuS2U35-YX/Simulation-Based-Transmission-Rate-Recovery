# Shared helper for the prespecified Experiment 4 trajectory figures.

suppressPackageStartupMessages(library(ggplot2))

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The digest package is required for provenance checks.")
}

source(file.path("config", "experiment_config.R"))

# The plotted observation-time target is the rate that drove the final Euler
# substep ending at each endpoint. It is the left limit of the prescribed path,
# so the high rate is retained at the switch endpoint itself.
true_B_at_times_for_figure <- function(times) {
  ifelse(
    times <= experiment_config$true_parameters[["t_switch"]],
    experiment_config$true_parameters[["Beta_high"]],
    experiment_config$true_parameters[["Beta_low"]]
  )
}

assert_close <- function(actual, expected, label, tolerance = 1e-10) {
  if (length(actual) != 1L || !is.finite(actual) ||
      abs(actual - expected) > tolerance) {
    stop(label, " mismatch: expected ", expected, ", found ", actual, ".")
  }
}

assert_grid <- function(actual, expected, label, tolerance = 1e-12) {
  if (length(actual) != length(expected) || any(!is.finite(actual)) ||
      max(abs(actual - expected)) > tolerance ||
      is.unsorted(actual, strictly = TRUE) || anyDuplicated(actual)) {
    stop(label, " does not match the configured time grid.")
  }
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

read_one_task <- function(path, task_id, label) {
  x <- read.csv(path, check.names = FALSE)
  x <- x[x$task_id == task_id, , drop = FALSE]
  if (nrow(x) != 1L) {
    stop(label, " must contain exactly one task-", task_id, " row.")
  }
  x
}

trajectory_theme <- function() {
  theme_classic(base_size = 8.5, base_family = "Helvetica") +
    theme(
      axis.line = element_line(colour = "#262626", linewidth = 0.38),
      axis.ticks = element_line(colour = "#262626", linewidth = 0.34),
      axis.ticks.length = grid::unit(2.5, "pt"),
      axis.title.x = element_text(size = 9.2, margin = margin(t = 6)),
      axis.title.y = element_text(size = 9.2, margin = margin(r = 7)),
      axis.text = element_text(size = 8.0, colour = "#262626"),
      legend.position = "top",
      legend.justification = "left",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 7.8, colour = "#262626"),
      legend.key.width = grid::unit(18, "pt"),
      legend.key.height = grid::unit(8, "pt"),
      legend.spacing.x = grid::unit(7, "pt"),
      legend.margin = margin(t = 0, r = 0, b = 3, l = 0),
      panel.grid = element_blank(),
      plot.margin = margin(t = 2, r = 5, b = 3, l = 3),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

export_trajectory_figure <- function(
  plot, stem, width_mm = 183, height_mm = 108, dpi = 600
) {
  width <- width_mm / 25.4
  height <- height_mm / 25.4

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width, height = height,
    family = "Helvetica", pointsize = 8.5, fallback_resolution = 600
  )
  print(plot)
  dev.off()

  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(
      paste0(stem, ".svg"), width = width, height = height,
      pointsize = 8.5, system_fonts = list(sans = "Helvetica")
    )
    print(plot)
    dev.off()
  } else {
    svg(
      paste0(stem, ".svg"), width = width, height = height,
      pointsize = 8.5, onefile = TRUE, family = "Helvetica"
    )
    print(plot)
    dev.off()
  }

  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(
      paste0(stem, ".png"), width = width, height = height,
      units = "in", res = 300, background = "white", pointsize = 8.5
    )
    print(plot)
    dev.off()
    ragg::agg_tiff(
      paste0(stem, ".tiff"), width = width, height = height,
      units = "in", res = dpi, compression = "lzw",
      background = "white", pointsize = 8.5
    )
    print(plot)
    dev.off()
  } else {
    png(
      paste0(stem, ".png"), width = width, height = height,
      units = "in", res = 300, type = "cairo", pointsize = 8.5
    )
    print(plot)
    dev.off()
    tiff(
      paste0(stem, ".tiff"), width = width, height = height,
      units = "in", res = dpi, compression = "lzw", type = "cairo",
      pointsize = 8.5
    )
    print(plot)
    dev.off()
  }
}

generate_selected_B_trajectory_figure <- function(
  task_id,
  expected_data_md5,
  expected_B0_hat,
  expected_sigma_beta_hat,
  expected_constant_B,
  figure_stem,
  selection_basis
) {
  output_root <- Sys.getenv("EXP4_OUTPUT_ROOT", unset = "")
  output_path <- function(...) {
    relative <- file.path(...)
    if (nzchar(output_root)) file.path(output_root, relative) else relative
  }

  gamma_best_path <- file.path(
    "results", "combined", "gamma", "combined_best_fit_summary.csv"
  )
  constant_best_path <- file.path(
    "results", "combined", "constant", "combined_best_fit_summary.csv"
  )
  gamma_B_path <- Sys.getenv(
    "EXP4_GAMMA_B_PATH",
    unset = file.path(
      "results", "combined", "gamma", "combined_B_paths.csv"
    )
  )
  constant_B_path <- file.path(
    "results", "combined", "constant", "combined_B_paths.csv"
  )
  observed_path <- file.path(
    "shared_data", sprintf("task_%03d", task_id), "observed_data.csv"
  )

  required_paths <- c(
    gamma_best_path, constant_best_path, gamma_B_path,
    constant_B_path, observed_path
  )
  if (any(!file.exists(required_paths))) {
    stop(
      "Missing required input(s): ",
      paste(required_paths[!file.exists(required_paths)], collapse = ", ")
    )
  }

  gamma_row <- read_one_task(gamma_best_path, task_id, "Gamma best-fit summary")
  constant_row <- read_one_task(
    constant_best_path, task_id, "Constant-B best-fit summary"
  )
  gamma_paths <- read.csv(gamma_B_path, check.names = FALSE)
  constant_paths <- read.csv(constant_B_path, check.names = FALSE)
  gamma_paths <- gamma_paths[gamma_paths$task_id == task_id, , drop = FALSE]
  constant_paths <- constant_paths[
    constant_paths$task_id == task_id, , drop = FALSE
  ]

  required_gamma_columns <- c(
    "task_id", "simulation_seed", "observed_data_md5", "week",
    "B_estimate", "B_true", "trajectory_seed", "path_semantics"
  )
  required_constant_columns <- c(
    "task_id", "simulation_seed", "observed_data_md5", "week",
    "B_estimate", "B_true"
  )
  if (!all(required_gamma_columns %in% names(gamma_paths)) ||
      !all(required_constant_columns %in% names(constant_paths))) {
    stop("Selected-task B paths lack required columns.")
  }

  expected_times <- seq(
    from = experiment_config$observation_interval,
    to = experiment_config$n_weeks,
    by = experiment_config$observation_interval
  )
  if (length(expected_times) != 70L) stop("Expected 70 observation times.")
  assert_grid(gamma_paths$week, expected_times, "Gamma B-path times")
  assert_grid(constant_paths$week, expected_times, "Constant-B path times")

  selected_data_md5 <- unname(tools::md5sum(observed_path))
  if (!identical(selected_data_md5, expected_data_md5)) {
    stop("Task-", task_id, " observed-data MD5 mismatch.")
  }
  record_md5 <- c(
    as.character(gamma_row$observed_data_md5[[1]]),
    as.character(constant_row$observed_data_md5[[1]]),
    unique(as.character(gamma_paths$observed_data_md5)),
    unique(as.character(constant_paths$observed_data_md5))
  )
  if (!all(record_md5 == expected_data_md5)) {
    stop("Selected-task fit and path records do not reference the shared data.")
  }
  if (!identical(as.character(gamma_row$status[[1]]), "success") ||
      !identical(as.character(constant_row$status[[1]]), "success")) {
    stop("Selected-task best-fit status is not successful for both models.")
  }

  assert_close(gamma_row$B0_hat[[1]], expected_B0_hat, "Gamma B0_hat")
  assert_close(
    gamma_row$sigma_beta_hat[[1]], expected_sigma_beta_hat,
    "Gamma sigma_beta_hat"
  )
  assert_close(
    constant_row$Beta_hat[[1]], expected_constant_B, "Constant-B estimate"
  )
  if (!all(gamma_paths$path_semantics ==
           "ancestry_preserving_sampled_latent_trajectory")) {
    stop("Gamma B path is not labelled as a sampled latent trajectory.")
  }
  if (length(unique(gamma_paths$trajectory_seed)) != 1L) {
    stop("Gamma B path does not have one trajectory seed.")
  }
  if (length(unique(constant_paths$B_estimate)) != 1L) {
    stop("Constant-B path is not a repeated static estimate.")
  }
  assert_close(
    unique(constant_paths$B_estimate), constant_row$Beta_hat[[1]],
    "Repeated constant-B estimate"
  )

  B_comparison <- data.frame(
    experiment = 4L,
    task_id = task_id,
    week = c(0, expected_times),
    B_true = true_B_at_times_for_figure(c(0, expected_times)),
    B_gamma_sampled_trajectory = c(
      gamma_row$B0_hat[[1]], gamma_paths$B_estimate
    ),
    B_constant = rep(constant_row$Beta_hat[[1]], 71L),
    trajectory_seed = rep(unique(gamma_paths$trajectory_seed), 71L),
    is_time_zero = c(TRUE, rep(FALSE, 70L)),
    stringsAsFactors = FALSE
  )

  result_dir <- output_path("results", "selected_trajectory")
  figure_dir <- output_path("figures", "comparison")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  B_csv <- file.path(
    result_dir,
    sprintf("experiment4_task%d_B_trajectory_comparison.csv", task_id)
  )
  provenance_path <- file.path(
    result_dir,
    sprintf("experiment4_task%d_B_trajectory_provenance.txt", task_id)
  )
  write.csv(B_comparison, B_csv, row.names = FALSE)

  series_levels <- c("Truth", "Gamma-noise", "Constant B")
  series_colours <- c(
    "Truth" = "#1A1A1A",
    "Gamma-noise" = "#0072B2",
    "Constant B" = "#D55E00"
  )
  series_linetypes <- c(
    "Truth" = "solid",
    "Gamma-noise" = "dashed",
    "Constant B" = "dotdash"
  )

  truth_plot <- data.frame(
    week = c(0, 5, 5, 10),
    value = c(4, 4, 2, 2),
    series = factor("Truth", levels = series_levels)
  )
  fitted_plot_data <- rbind(
    data.frame(
      week = B_comparison$week,
      value = B_comparison$B_gamma_sampled_trajectory,
      series = factor("Gamma-noise", levels = series_levels)
    ),
    data.frame(
      week = B_comparison$week,
      value = B_comparison$B_constant,
      series = factor("Constant B", levels = series_levels)
    )
  )
  all_values <- c(truth_plot$value, fitted_plot_data$value)
  y_upper <- ceiling(max(all_values, na.rm = TRUE) * 1.04 * 5) / 5
  y_lower <- min(0, floor(min(all_values, na.rm = TRUE) * 0.96 * 5) / 5)

  plot <- ggplot() +
    geom_vline(
      xintercept = experiment_config$true_parameters[["t_switch"]],
      colour = "#B8B8B8", linewidth = 0.42, linetype = "dotted"
    ) +
    geom_path(
      data = fitted_plot_data,
      aes(x = week, y = value, colour = series, linetype = series),
      linewidth = 0.84, lineend = "butt", linejoin = "round"
    ) +
    geom_path(
      data = truth_plot,
      aes(x = week, y = value, colour = series, linetype = series),
      linewidth = 0.92, lineend = "butt", linejoin = "mitre"
    ) +
    annotate(
      "text", x = 5.10, y = 0.985 * y_upper, label = "Week 5",
      hjust = 0, vjust = 1, size = 2.55, family = "Helvetica",
      colour = "#6F6F6F"
    ) +
    scale_colour_manual(
      values = series_colours,
      breaks = series_levels,
      labels = c(
        "True B(t)", "Gamma-noise sampled trajectory", "Constant B"
      ),
      drop = FALSE
    ) +
    scale_linetype_manual(
      values = series_linetypes, breaks = series_levels,
      drop = FALSE, guide = "none"
    ) +
    scale_x_continuous(
      breaks = seq(0, 10, by = 2), limits = c(0, 10),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = pretty(c(y_lower, y_upper), n = 6),
      limits = c(y_lower, y_upper),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Week", y = expression("Transmission rate  " * B(t))) +
    guides(
      colour = guide_legend(
        nrow = 1, byrow = TRUE,
        override.aes = list(
          linewidth = 0.9,
          linetype = unname(series_linetypes[series_levels])
        )
      )
    ) +
    trajectory_theme()

  output_stem <- file.path(figure_dir, figure_stem)
  export_trajectory_figure(plot, output_stem)

  generated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  provenance <- c(
    "experiment: 4",
    paste0("selected_task_id: ", task_id),
    paste0("selection_basis: ", selection_basis),
    paste0("simulation_seed: ", gamma_row$simulation_seed[[1]]),
    paste0("observed_data_md5: ", selected_data_md5),
    paste0("observed_data_sha256: ", sha256_file(observed_path)),
    paste0("best_run: ", gamma_row$best_run[[1]]),
    paste0("Gamma_B0_hat: ", format(gamma_row$B0_hat[[1]], digits = 17)),
    paste0(
      "Gamma_sigma_beta_hat: ",
      format(gamma_row$sigma_beta_hat[[1]], digits = 17)
    ),
    paste0(
      "constant_Beta_hat: ",
      format(constant_row$Beta_hat[[1]], digits = 17)
    ),
    paste0("trajectory_seed: ", unique(gamma_paths$trajectory_seed)),
    "trajectory_source: canonical Gamma combined_B_paths.csv used for the recovery metrics",
    "trajectory_semantics: one ancestry-preserving sampled latent trajectory from the final plug-in particle-filter approximation",
    "metric_times: 70 observation times; t0 excluded from RSS, RMSE, mean error, and AOB",
    "figure_times: t0 plus the 70 observation times",
    "filtering_mean_used: no",
    "across_task_average_used: no",
    "parameter_uncertainty_integrated: no",
    paste0("R_version: ", R.version.string),
    paste0("ggplot2_version: ", as.character(packageVersion("ggplot2"))),
    paste0("generated_at_utc: ", generated_at_utc),
    paste0("comparison_csv_sha256: ", sha256_file(B_csv))
  )
  writeLines(provenance, provenance_path, useBytes = TRUE)

  cat(
    "Experiment 4 task-", task_id,
    " sampled B-trajectory figure generated successfully.\n",
    "Comparison data: ", B_csv, "\n",
    "Provenance: ", provenance_path, "\n",
    "Figure: ", paste0(output_stem, ".pdf"), "\n",
    sep = ""
  )
}

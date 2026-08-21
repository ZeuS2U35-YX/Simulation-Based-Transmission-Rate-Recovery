# ============================================================
# Publication-ready three-model recovery figures.
#
# Figure 1: three RMSE histograms plus three paired RMSE scatters.
# Figure 2: 3 x 3 mean-error histogram grid.
#
# All graphics and previews are generated in R. No inferential
# tests or p values are added; the replicate is the simulation task.
# ============================================================

options(stringsAsFactors = FALSE)

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
)
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "path_helpers.R"))
experiment_directory <- get_experiment_directory()

args <- commandArgs(trailingOnly = TRUE)
comparison_dir <- if (length(args) >= 1L) args[[1]] else {
  file.path(experiment_directory, "results", "comparison_three_models")
}
figures_dir <- if (length(args) >= 2L) args[[2]] else {
  file.path(experiment_directory, "figures", "comparison_three_models")
}
source_data_dir <- file.path(figures_dir, "source_data")

required_packages <- c("ggplot2", "patchwork", "scales", "svglite", "ragg")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

library(ggplot2)
library(patchwork)

metrics <- read.csv(
  file.path(comparison_dir, "three_model_task_metrics_long.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
paired_differences <- read.csv(
  file.path(comparison_dir, "pairwise_RMSE_differences.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
run_check <- read.csv(
  file.path(comparison_dir, "three_model_run_check_summary.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
if (nrow(run_check) != 1L || !isTRUE(run_check$all_checks_pass[[1]])) {
  stop("Three-model run checks must pass before plotting.")
}

model_levels <- c("gamma_noise", "bspline_B", "constant_B")
display_levels <- c("Gamma-noise", "B-spline", "Constant-B")
model_palette <- c(
  "Gamma-noise" = "#4C78A8",
  "B-spline" = "#E39C37",
  "Constant-B" = "#7F7F7F"
)
metrics$model <- factor(metrics$model, levels = model_levels)
metrics$display_model <- factor(metrics$display_model, levels = display_levels)
metrics <- metrics[order(metrics$model, metrics$task_id), , drop = FALSE]

if (nrow(metrics) != 600L ||
    !all(table(metrics$model) == 200L) ||
    !all(table(metrics$task_id) == 3L) ||
    any(!is.finite(metrics$RMSE))) {
  stop("Plotting metrics must contain 200 finite rows for each of three models.")
}
if (nrow(paired_differences) != 600L ||
    !all(table(paired_differences$comparison) == 200L)) {
  stop("Paired RMSE source data must contain 200 rows for each comparison.")
}

theme_three_model <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#222222"),
      axis.text = element_text(size = 7.2, colour = "#222222"),
      axis.title = element_text(size = 8, colour = "#222222"),
      strip.background = element_rect(
        fill = "#F2F2F2", colour = "#BFBFBF", linewidth = 0.3
      ),
      strip.text = element_text(size = 8.1, face = "bold", colour = "#222222"),
      plot.title = element_text(size = 8.6, face = "bold", colour = "#222222"),
      plot.subtitle = element_text(size = 7.5, colour = "#444444"),
      plot.caption = element_text(size = 7.2, colour = "#444444", hjust = 0),
      panel.grid = element_blank(),
      legend.text = element_text(size = 7.2),
      legend.title = element_text(size = 7.5),
      plot.margin = margin(4, 5, 4, 5)
    )
}
theme_set(theme_three_model())

save_publication_figure <- function(
  plot, stem, width_mm, height_mm, dpi = 600L
) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(
    paste0(stem, ".svg"), width = width_in, height = height_in,
    pointsize = 8, system_fonts = list(sans = "Helvetica")
  )
  print(plot)
  dev.off()

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in,
    family = "Helvetica", pointsize = 8, fallback_resolution = dpi
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300, background = "white", pointsize = 8
  )
  print(plot)
  dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw",
    background = "white", pointsize = 8
  )
  print(plot)
  dev.off()
}

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Figure 1: common-scale RMSE histograms and paired scatters.
# ------------------------------------------------------------

rmse_bin_width <- 0.05
rmse_limit <- c(
  max(0, floor(min(metrics$RMSE) / 0.1) * 0.1),
  ceiling(max(metrics$RMSE) / 0.1) * 0.1
)
scatter_limit <- c(0, ceiling(max(metrics$RMSE) / 0.1) * 0.1)

make_rmse_histogram <- function(model_name, display_name) {
  plot_data <- metrics[metrics$model == model_name, , drop = FALSE]
  ggplot(plot_data, aes(x = RMSE)) +
    geom_histogram(
      binwidth = rmse_bin_width, boundary = 0,
      fill = unname(model_palette[[display_name]]),
      colour = "white", linewidth = 0.25
    ) +
    annotate(
      "text", x = Inf, y = Inf, label = "n = 200",
      hjust = 1.12, vjust = 1.35, size = 2.7, colour = "#333333"
    ) +
    coord_cartesian(xlim = rmse_limit, clip = "on") +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(title = display_name, x = "RMSE", y = "Simulation replicates") +
    theme_three_model()
}

comparison_display <- list(
  "Gamma-noise vs Constant-B" = c("Gamma-noise", "Constant-B"),
  "B-spline vs Constant-B" = c("B-spline", "Constant-B"),
  "B-spline vs Gamma-noise" = c("B-spline", "Gamma-noise")
)

make_rmse_scatter <- function(comparison_name) {
  plot_data <- paired_differences[
    paired_differences$comparison == comparison_name, , drop = FALSE
  ]
  labels <- comparison_display[[comparison_name]]
  ggplot(plot_data, aes(x = RMSE_model_a, y = RMSE_model_b)) +
    geom_abline(
      intercept = 0, slope = 1, linetype = "dashed",
      linewidth = 0.45, colour = "#555555"
    ) +
    geom_point(
      shape = 21, size = 1.65, stroke = 0.25, alpha = 0.72,
      fill = "#4C78A8", colour = "white"
    ) +
    annotate(
      "text", x = Inf, y = -Inf, label = "n = 200",
      hjust = 1.12, vjust = -0.55, size = 2.7, colour = "#333333"
    ) +
    coord_equal(xlim = scatter_limit, ylim = scatter_limit, expand = FALSE) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(
      title = comparison_name,
      x = paste0(labels[[1]], " RMSE"),
      y = paste0(labels[[2]], " RMSE")
    ) +
    theme_three_model()
}

figure_1 <- wrap_plots(
  make_rmse_histogram("gamma_noise", "Gamma-noise"),
  make_rmse_scatter("Gamma-noise vs Constant-B"),
  make_rmse_histogram("bspline_B", "B-spline"),
  make_rmse_scatter("B-spline vs Constant-B"),
  make_rmse_histogram("constant_B", "Constant-B"),
  make_rmse_scatter("B-spline vs Gamma-noise"),
  ncol = 2, widths = c(0.9, 1.1)
) +
  plot_annotation(
    title = "RMSE comparison across three recovery models",
    subtitle = paste0(
      "Matched simulation replicates; common histogram bin width = ",
      format(rmse_bin_width, nsmall = 2), "; dashed lines indicate y = x"
    ),
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(
        family = "Helvetica", size = 10, face = "bold", colour = "#222222"
      ),
      plot.subtitle = element_text(
        family = "Helvetica", size = 7.8, colour = "#444444"
      ),
      plot.tag = element_text(
        family = "Helvetica", size = 9, face = "bold", colour = "#222222"
      )
    )
  )

figure_1_hist_source <- metrics[, c(
  "task_id", "simulation_seed", "observed_data_md5",
  "model", "display_model", "RMSE"
)]
figure_1_scatter_source <- paired_differences[, c(
  "task_id", "simulation_seed", "observed_data_md5", "comparison",
  "model_a", "model_b", "RMSE_model_a", "RMSE_model_b",
  "paired_RMSE_difference_model_a_minus_model_b", "winner_lower_RMSE"
)]
write.csv(
  figure_1_hist_source,
  file.path(source_data_dir, "figure_1_RMSE_histograms_source_data.csv"),
  row.names = FALSE
)
write.csv(
  figure_1_scatter_source,
  file.path(source_data_dir, "figure_1_RMSE_paired_scatter_source_data.csv"),
  row.names = FALSE
)

save_publication_figure(
  figure_1,
  file.path(figures_dir, "01_three_model_RMSE_comparison"),
  width_mm = 183, height_mm = 190, dpi = 600L
)

# ------------------------------------------------------------
# Figure 2: 3 x 3 common-scale mean-error histograms.
# ------------------------------------------------------------

error_long <- rbind(
  data.frame(
    metrics[, c(
      "task_id", "simulation_seed", "observed_data_md5",
      "model", "display_model"
    )],
    error_scope = "Overall",
    mean_error_value = metrics$mean_error,
    stringsAsFactors = FALSE
  ),
  data.frame(
    metrics[, c(
      "task_id", "simulation_seed", "observed_data_md5",
      "model", "display_model"
    )],
    error_scope = "Through week 5",
    mean_error_value = metrics$mean_error_through_5,
    stringsAsFactors = FALSE
  ),
  data.frame(
    metrics[, c(
      "task_id", "simulation_seed", "observed_data_md5",
      "model", "display_model"
    )],
    error_scope = "After week 5",
    mean_error_value = metrics$mean_error_after_5,
    stringsAsFactors = FALSE
  )
)
error_long$display_model <- factor(
  error_long$display_model, levels = display_levels
)
error_long$error_scope <- factor(
  error_long$error_scope,
  levels = c("Overall", "Through week 5", "After week 5")
)
panel_means <- aggregate(
  mean_error_value ~ display_model + error_scope,
  data = error_long, FUN = mean
)
names(panel_means)[names(panel_means) == "mean_error_value"] <- "panel_mean"

error_bin_width <- 0.10
error_axis_step <- 0.25
error_limit <- c(
  floor(min(error_long$mean_error_value) / error_axis_step) * error_axis_step,
  ceiling(max(error_long$mean_error_value) / error_axis_step) * error_axis_step
)
error_breaks <- seq(
  ceiling(error_limit[[1]]),
  floor(error_limit[[2]]),
  by = 1
)

figure_2 <- ggplot(
  error_long,
  aes(x = mean_error_value, fill = display_model)
) +
  geom_histogram(
    binwidth = error_bin_width, boundary = 0,
    colour = "white", linewidth = 0.25
  ) +
  geom_vline(
    xintercept = 0, linetype = "dashed", linewidth = 0.24,
    colour = "#8C8C8C", alpha = 0.75
  ) +
  geom_vline(
    data = panel_means,
    aes(xintercept = panel_mean, colour = display_model),
    inherit.aes = FALSE, linetype = "dotted", linewidth = 0.65
  ) +
  facet_grid(rows = vars(display_model), cols = vars(error_scope)) +
  scale_fill_manual(values = model_palette, drop = FALSE) +
  scale_colour_manual(values = model_palette, drop = FALSE) +
  scale_x_continuous(breaks = error_breaks) +
  coord_cartesian(xlim = error_limit, clip = "on") +
  labs(
    x = "Mean error in B (week\u207B\u00B9)",
    y = "Simulation replicates"
  ) +
  theme_three_model(base_size = 6.5) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, colour = "#222222"),
    axis.title = element_text(size = 6.5, colour = "#222222"),
    strip.text = element_text(size = 7, face = "bold", colour = "#222222"),
    strip.text.y = element_text(size = 7, face = "bold", angle = 0),
    panel.spacing.x = grid::unit(3, "mm"),
    plot.margin = margin(2, 2, 2, 2)
  )

figure_2_legend <- paste0(
  "Mean-error distributions across recovery periods for three transmission-rate ",
  "representations. Histograms show mean estimation error in B (estimated minus ",
  "true B; week\u207B\u00B9) across n = 200 matched simulation replicates. Columns ",
  "show the overall period, through week 5 (week \u2264 5), and after week 5 ",
  "(week > 5); rows show the Gamma-noise, B-spline and Constant-B models. All ",
  "panels use a common bin width of 0.10 and a common x-axis. The grey dashed ",
  "line marks zero error, and the model-coloured dotted line marks the panel ",
  "mean. No simulation replicates were excluded."
)

write.csv(
  error_long,
  file.path(source_data_dir, "figure_2_mean_error_histograms_source_data.csv"),
  row.names = FALSE
)
write.csv(
  panel_means,
  file.path(source_data_dir, "figure_2_panel_means_source_data.csv"),
  row.names = FALSE
)
writeLines(
  figure_2_legend,
  file.path(figures_dir, "02_three_model_mean_error_legend.txt"),
  useBytes = TRUE
)
write.csv(
  data.frame(
    figure = c("Figure 1", "Figure 2"),
    histogram_bin_width = c(rmse_bin_width, error_bin_width),
    common_axis_min = c(rmse_limit[[1]], error_limit[[1]]),
    common_axis_max = c(rmse_limit[[2]], error_limit[[2]]),
    n_per_model_or_panel = c(200L, 200L),
    stringsAsFactors = FALSE
  ),
  file.path(source_data_dir, "figure_plotting_parameters.csv"),
  row.names = FALSE
)

save_publication_figure(
  figure_2,
  file.path(figures_dir, "02_three_model_mean_error_comparison"),
  width_mm = 183, height_mm = 155, dpi = 600L
)

cat(
  "Generated two R figures in PDF, PNG, SVG, and 600-dpi TIFF formats, ",
  "with source-data CSV files.\n",
  sep = ""
)

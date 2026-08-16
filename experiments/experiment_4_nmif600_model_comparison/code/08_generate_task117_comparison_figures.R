# ============================================================
# Generate Experiment 4 task-117 B(t) and infectious-path figures
#
# Run from the Experiment 4 root directory:
#   Rscript code/08_generate_task117_comparison_figures.R
#
# The B(t) panel preserves the prespecified task-117 ancestry trajectory
# already stored in results/selected_trajectory. The infectious panel uses
# filtering means from the saved task-117 final particle-filter seeds.
# ============================================================

options(stringsAsFactors = FALSE, digits = 17)

suppressPackageStartupMessages({
  library(pomp)
  library(ggplot2)
})

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The digest package is required for provenance checks.")
}

source(file.path("config", "experiment_config.R"))
source(file.path("code", "model_components.R"))

task_id <- 117L
expected_data_md5 <- "64dffb15867fda5ef262e2caf0e46bbf"
Np <- 50000L

output_root <- Sys.getenv("EXP4_OUTPUT_ROOT", unset = "")
output_path <- function(...) {
  relative <- file.path(...)
  if (nzchar(output_root)) file.path(output_root, relative) else relative
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

read_one_task <- function(path, id, label) {
  x <- read.csv(path, check.names = FALSE)
  x <- x[x$task_id == id, , drop = FALSE]
  if (nrow(x) != 1L) stop(label, " must contain exactly one task-", id, " row.")
  x
}

gamma_best_path <- file.path(
  "results", "combined", "gamma", "combined_best_fit_summary.csv"
)
constant_best_path <- file.path(
  "results", "combined", "constant", "combined_best_fit_summary.csv"
)
paramlist_path <- file.path("results", "paramlist.csv")
observed_path <- file.path(
  "shared_data", sprintf("task_%03d", task_id), "observed_data.csv"
)
simulated_path <- file.path(
  "shared_data", sprintf("task_%03d", task_id), "simulated_data.csv"
)
historical_B_path <- file.path(
  "results", "selected_trajectory",
  "experiment4_task117_B_trajectory_comparison.csv"
)

required_paths <- c(
  gamma_best_path, constant_best_path, paramlist_path,
  observed_path, simulated_path, historical_B_path
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
seed_row <- read_one_task(paramlist_path, task_id, "Parameter list")
observed_data <- read.csv(observed_path, check.names = FALSE)
simulated_data <- read.csv(simulated_path, check.names = FALSE)
B_comparison <- read.csv(historical_B_path, check.names = FALSE)

required_B <- c(
  "experiment", "task_id", "week", "B_trajectory", "B_true",
  "B_constant", "is_time_zero"
)
if (!all(c("week", "reports") %in% names(observed_data))) {
  stop("Task-117 observed data lack week or reports.")
}
if (!all(c("week", "I") %in% names(simulated_data))) {
  stop("Task-117 simulated data lack week or I.")
}
if (!all(required_B %in% names(B_comparison))) {
  stop("The task-117 trajectory artifact lacks required columns.")
}

observed_data <- observed_data[, c("week", "reports"), drop = FALSE]
simulated_data <- simulated_data[, c("week", "I"), drop = FALSE]
expected_times <- seq(
  from = experiment_config$observation_interval,
  to = experiment_config$n_weeks,
  by = experiment_config$observation_interval
)
trajectory_times <- c(0, expected_times)
if (length(expected_times) != 70L) stop("Expected 70 observation times.")
assert_grid(observed_data$week, expected_times, "Observed-data times")
assert_grid(simulated_data$week, expected_times, "Simulated-data times")
assert_grid(B_comparison$week, trajectory_times, "Task-117 trajectory times")

selected_data_md5 <- unname(tools::md5sum(observed_path))
if (!identical(selected_data_md5, expected_data_md5)) {
  stop("Task-117 observed-data MD5 mismatch.")
}
fit_md5 <- c(
  as.character(gamma_row$observed_data_md5[[1]]),
  as.character(constant_row$observed_data_md5[[1]])
)
if (!all(fit_md5 == expected_data_md5)) {
  stop("Gamma and constant-B task-117 records do not reference the shared data.")
}
simulation_seeds <- c(
  gamma_row$simulation_seed[[1]], constant_row$simulation_seed[[1]],
  seed_row$simulation_seed[[1]]
)
if (length(unique(simulation_seeds)) != 1L || simulation_seeds[[1]] != 1117L) {
  stop("Task-117 records do not share simulation seed 1117.")
}
if (!identical(as.character(gamma_row$status[[1]]), "success") ||
    !identical(as.character(constant_row$status[[1]]), "success")) {
  stop("Task-117 best-fit status is not successful for both models.")
}
if (experiment_config$Np_final != Np) {
  stop("Canonical Experiment 4 Np_final is not 50000.")
}

assert_close(gamma_row$B0_hat[[1]], 3.6102663778531099, "Gamma B0_hat")
assert_close(
  gamma_row$sigma_beta_hat[[1]], 0.34528388982657399,
  "Gamma sigma_beta_hat"
)
assert_close(
  constant_row$Beta_hat[[1]], 3.3878099385514902,
  "Constant-B Beta_hat"
)
if (any(B_comparison$task_id != task_id) ||
    any(B_comparison$experiment != 4L)) {
  stop("The historical trajectory artifact is not exclusively Experiment 4 task 117.")
}
if (!isTRUE(B_comparison$is_time_zero[[1]]) ||
    any(B_comparison$is_time_zero[-1L])) {
  stop("Task-117 trajectory time-zero indicator is malformed.")
}
if (any(B_comparison$B_true !=
        true_B_at_times(B_comparison$week, experiment_config))) {
  stop("Task-117 trajectory truth does not implement the week-5 switch.")
}
assert_close(B_comparison$B_trajectory[[1]], gamma_row$B0_hat[[1]],
             "Task-117 trajectory initial value")
if (length(unique(B_comparison$B_constant)) != 1L) {
  stop("Task-117 constant-B trajectory is not static.")
}
assert_close(
  unique(B_comparison$B_constant), constant_row$Beta_hat[[1]],
  "Task-117 repeated constant-B estimate"
)

gamma_model <- make_gamma_model(observed_data, experiment_config)
theta_gamma <- gamma_baseline_parameters(experiment_config)
theta_gamma[["B0"]] <- gamma_row$B0_hat[[1]]
theta_gamma[["sigma_beta"]] <- gamma_row$sigma_beta_hat[[1]]

constant_model <- make_constant_model(observed_data, experiment_config)
theta_constant <- constant_baseline_parameters(experiment_config)
theta_constant[["Beta"]] <- constant_row$Beta_hat[[1]]

gamma_final_pf_seed <- as.integer(seed_row$gamma_final_pf_seed[[1]])
constant_final_pf_seed <- as.integer(seed_row$constant_final_pf_seed[[1]])
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
rng_kind <- RNGkind()

set.seed(gamma_final_pf_seed)
gamma_pf <- pfilter(
  gamma_model, params = theta_gamma, Np = Np, filter.mean = TRUE
)
gamma_fm <- filter_mean(gamma_pf)

set.seed(constant_final_pf_seed)
constant_pf <- pfilter(
  constant_model, params = theta_constant, Np = Np, filter.mean = TRUE
)
constant_fm <- filter_mean(constant_pf)

if (!"I" %in% rownames(gamma_fm) || !"I" %in% rownames(constant_fm)) {
  stop("Final task-117 particle filters did not return infectious means.")
}
assert_grid(as.numeric(time(gamma_pf)), expected_times, "Gamma filter times")
assert_grid(as.numeric(time(constant_pf)), expected_times, "Constant filter times")
assert_close(
  as.numeric(logLik(gamma_pf)), gamma_row$final_pf_logLik[[1]],
  "Gamma final-particle-filter log likelihood", tolerance = 1e-8
)
assert_close(
  as.numeric(logLik(constant_pf)), constant_row$final_pf_logLik[[1]],
  "Constant-B final-particle-filter log likelihood", tolerance = 1e-8
)

infectious_comparison <- data.frame(
  experiment = 4L,
  task_id = task_id,
  week = expected_times,
  I_true = as.numeric(simulated_data$I),
  I_gamma_filtered_mean = as.numeric(gamma_fm["I", ]),
  I_constant_filtered_mean = as.numeric(constant_fm["I", ]),
  stringsAsFactors = FALSE
)
if (any(!is.finite(unlist(B_comparison[, c(
  "week", "B_trajectory", "B_true", "B_constant"
)]))) || any(!is.finite(unlist(infectious_comparison[-c(1L, 2L)])))) {
  stop("Task-117 comparison artifacts contain non-finite values.")
}
if (any(infectious_comparison[, c(
  "I_true", "I_gamma_filtered_mean", "I_constant_filtered_mean"
)] < 0)) {
  stop("Task-117 infectious paths contain negative values.")
}

result_dir <- output_path("results", "selected_trajectory")
figure_dir <- output_path("figures", "comparison")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

I_csv <- file.path(
  result_dir, "experiment4_task117_infectious_filter_comparison.csv"
)
provenance_path <- file.path(
  result_dir, "experiment4_task117_infectious_comparison_provenance.txt"
)
write.csv(infectious_comparison, I_csv, row.names = FALSE)

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

theme_task117 <- function() {
  theme_classic(base_size = 8.5, base_family = "Arial") +
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

add_series_scales <- function(labels) {
  list(
    scale_colour_manual(
      values = series_colours, breaks = series_levels, labels = labels,
      drop = FALSE
    ),
    scale_linetype_manual(
      values = series_linetypes, breaks = series_levels, labels = labels,
      drop = FALSE, guide = "none"
    )
  )
}

series_guide <- guides(
  colour = guide_legend(
    nrow = 1, byrow = TRUE,
    override.aes = list(
      linewidth = 0.9,
      linetype = unname(series_linetypes[series_levels])
    )
  )
)

make_B_plot <- function() {
  B_truth <- data.frame(
    week = c(0, 5, 5, 10), value = c(4, 4, 2, 2),
    series = factor("Truth", levels = series_levels)
  )
  B_gamma <- data.frame(
    week = B_comparison$week, value = B_comparison$B_trajectory,
    series = factor("Gamma-noise", levels = series_levels)
  )
  B_constant <- data.frame(
    week = c(0, 10), value = rep(B_comparison$B_constant[[1]], 2L),
    series = factor("Constant B", levels = series_levels)
  )
  y_upper <- ceiling(1.03 * max(B_gamma$value, B_truth$value) * 5) / 5

  ggplot() +
    geom_vline(
      xintercept = experiment_config$true_parameters[["t_switch"]],
      colour = "#B8B8B8", linewidth = 0.42, linetype = "dotted"
    ) +
    geom_path(
      data = B_constant,
      aes(x = week, y = value, colour = series, linetype = series),
      linewidth = 0.78, lineend = "butt"
    ) +
    geom_path(
      data = B_gamma,
      aes(x = week, y = value, colour = series, linetype = series),
      linewidth = 0.82, lineend = "butt", linejoin = "round"
    ) +
    geom_path(
      data = B_truth,
      aes(x = week, y = value, colour = series, linetype = series),
      linewidth = 0.92, lineend = "butt", linejoin = "mitre"
    ) +
    annotate(
      "text", x = 5.10, y = 0.985 * y_upper, label = "Week 5",
      hjust = 0, vjust = 1, size = 2.55, family = "Arial",
      colour = "#6F6F6F"
    ) +
    add_series_scales(c(
      "True B(t)", "Gamma-noise trajectory", "Constant B"
    )) +
    scale_x_continuous(
      breaks = seq(0, 10, by = 2), limits = c(0, 10),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = pretty(c(0, y_upper), n = 6), limits = c(0, y_upper),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Week", y = expression("Transmission rate  " * B(t))) +
    series_guide +
    theme_task117()
}

make_infectious_plot <- function() {
  y_upper <- ceiling(1.04 * max(
    infectious_comparison$I_true,
    infectious_comparison$I_gamma_filtered_mean,
    infectious_comparison$I_constant_filtered_mean
  ) / 100) * 100
  I_truth <- data.frame(
    week = infectious_comparison$week,
    value = infectious_comparison$I_true,
    series = factor("Truth", levels = series_levels)
  )
  I_gamma <- data.frame(
    week = infectious_comparison$week,
    value = infectious_comparison$I_gamma_filtered_mean,
    series = factor("Gamma-noise", levels = series_levels)
  )
  I_constant <- data.frame(
    week = infectious_comparison$week,
    value = infectious_comparison$I_constant_filtered_mean,
    series = factor("Constant B", levels = series_levels)
  )

  ggplot() +
    geom_vline(
      xintercept = experiment_config$true_parameters[["t_switch"]],
      colour = "#B8B8B8", linewidth = 0.42, linetype = "dotted"
    ) +
    geom_path(
      data = I_constant,
      aes(x = week, y = value, colour = series, linetype = series),
      linewidth = 0.78, lineend = "butt", linejoin = "round"
    ) +
    geom_path(
      data = I_gamma,
      aes(x = week, y = value, colour = series, linetype = series),
      linewidth = 0.82, lineend = "butt", linejoin = "round"
    ) +
    geom_path(
      data = I_truth,
      aes(x = week, y = value, colour = series, linetype = series),
      linewidth = 0.92, lineend = "butt", linejoin = "round"
    ) +
    annotate(
      "text", x = 5.10, y = 0.985 * y_upper, label = "Week 5",
      hjust = 0, vjust = 1, size = 2.55, family = "Arial",
      colour = "#6F6F6F"
    ) +
    add_series_scales(c("True I(t)", "Gamma-noise", "Constant B")) +
    scale_x_continuous(
      breaks = seq(0, 10, by = 2), limits = c(0, 10),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = seq(0, y_upper, by = 100), limits = c(0, y_upper),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Week", y = expression("Infectious individuals  " * I(t))) +
    series_guide +
    theme_task117()
}

export_figure <- function(
  plot, stem, width_mm = 183, height_mm = 108, dpi = 600
) {
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width, height = height,
    family = "Arial", pointsize = 8.5, fallback_resolution = 600
  )
  print(plot)
  dev.off()

  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(
      paste0(stem, ".svg"), width = width, height = height,
      pointsize = 8.5, system_fonts = list(sans = "Arial")
    )
    print(plot)
    dev.off()
  } else {
    svg(
      paste0(stem, ".svg"), width = width, height = height,
      pointsize = 8.5, onefile = TRUE, family = "Arial"
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

B_figure_stem <- file.path(
  figure_dir, "08_task117_B_trajectory_comparison"
)
I_figure_stem <- file.path(
  figure_dir, "09_task117_infectious_path_comparison"
)
export_figure(make_B_plot(), B_figure_stem)
export_figure(make_infectious_plot(), I_figure_stem)

repo_root <- normalizePath(file.path(getwd(), "..", ".."), mustWork = TRUE)
repo_commit <- system2(
  "git", c("-C", repo_root, "rev-parse", "HEAD"), stdout = TRUE
)
generated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
gamma_rmse <- sqrt(mean(
  (infectious_comparison$I_gamma_filtered_mean - infectious_comparison$I_true)^2
))
constant_rmse <- sqrt(mean(
  (infectious_comparison$I_constant_filtered_mean - infectious_comparison$I_true)^2
))
gamma_mae <- mean(abs(
  infectious_comparison$I_gamma_filtered_mean - infectious_comparison$I_true
))
constant_mae <- mean(abs(
  infectious_comparison$I_constant_filtered_mean - infectious_comparison$I_true
))
provenance <- c(
  "experiment: 4",
  paste0("selected_task_id: ", task_id),
  "selection_basis: prespecified task 117 retained as an additional illustration",
  paste0("simulation_seed: ", gamma_row$simulation_seed[[1]]),
  paste0("observed_data_md5: ", selected_data_md5),
  paste0("observed_data_sha256: ", sha256_file(observed_path)),
  paste0("particle_count_per_infectious_filter: ", Np),
  paste0("Gamma_final_pf_seed: ", gamma_final_pf_seed),
  paste0("constant_final_pf_seed: ", constant_final_pf_seed),
  paste0("Gamma_B0_hat: ", format(gamma_row$B0_hat[[1]], digits = 17)),
  paste0(
    "Gamma_sigma_beta_hat: ",
    format(gamma_row$sigma_beta_hat[[1]], digits = 17)
  ),
  paste0(
    "constant_Beta_hat: ",
    format(constant_row$Beta_hat[[1]], digits = 17)
  ),
  paste0("Gamma_final_pf_logLik: ", format(logLik(gamma_pf), digits = 17)),
  paste0(
    "constant_final_pf_logLik: ",
    format(logLik(constant_pf), digits = 17)
  ),
  "B_gamma_semantics: saved ancestry-preserving task-117 smoothing-trajectory approximation",
  "B_constant_semantics: fitted static task-117 estimate",
  "I_gamma_semantics: observation-time filtering mean from the Gamma-noise model",
  "I_constant_semantics: observation-time filtering mean from the constant-B model",
  "I_true_semantics: task-117 latent infectious state from the shared simulated data",
  paste0("I_gamma_RMSE: ", format(gamma_rmse, digits = 17)),
  paste0("I_constant_RMSE: ", format(constant_rmse, digits = 17)),
  paste0("I_gamma_MAE: ", format(gamma_mae, digits = 17)),
  paste0("I_constant_MAE: ", format(constant_mae, digits = 17)),
  "parameter_uncertainty_integrated: no",
  "across_task_average: no",
  "MIF2_rerun: no",
  "independent_likelihood_evaluations_rerun: no",
  paste0("number_of_observations: ", nrow(infectious_comparison)),
  paste0("R_version: ", R.version.string),
  paste0("pomp_version: ", as.character(packageVersion("pomp"))),
  paste0("RNGkind: ", paste(rng_kind, collapse = ", ")),
  paste0("generated_at_utc: ", generated_at_utc),
  paste0("repository_base_commit: ", repo_commit),
  paste0("historical_B_csv_sha256: ", sha256_file(historical_B_path)),
  paste0("infectious_comparison_csv_sha256: ", sha256_file(I_csv))
)
writeLines(provenance, provenance_path, useBytes = TRUE)

cat(
  "Experiment 4 task-117 comparison figures generated successfully.\n",
  "Infectious comparison data: ", I_csv, "\n",
  "Provenance: ", provenance_path, "\n",
  "B figure: ", paste0(B_figure_stem, ".pdf"), "\n",
  "Infectious figure: ", paste0(I_figure_stem, ".pdf"), "\n",
  sep = ""
)

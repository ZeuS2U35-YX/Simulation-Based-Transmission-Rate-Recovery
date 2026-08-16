# ============================================================
# Generate the Experiment 4 task-1 B(t) and infectious-path figures
#
# Run from the Experiment 4 root directory:
#   Rscript code/07_generate_task1_comparison_figures.R
#
# The script uses the saved task-1 best-fit parameters and final particle-
# filter seeds. It does not rerun MIF2 or the independent likelihood
# evaluations, and it does not alter the 200-task comparison metrics.
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

task_id <- 1L
expected_data_md5 <- "64a1b5c02bdecc100b37ee56029391d3"
Np <- 50000L

# For QA, EXP4_OUTPUT_ROOT can redirect generated artifacts to a staging
# directory while all scientific inputs continue to come from Experiment 4.
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
    stop(label, " does not match the configured observation grid.")
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
gamma_B_path <- file.path(
  "results", "combined", "gamma", "combined_B_paths.csv"
)
constant_B_path <- file.path(
  "results", "combined", "constant", "combined_B_paths.csv"
)
paramlist_path <- file.path("results", "paramlist.csv")
observed_path <- file.path(
  "shared_data", sprintf("task_%03d", task_id), "observed_data.csv"
)
simulated_path <- file.path(
  "shared_data", sprintf("task_%03d", task_id), "simulated_data.csv"
)

required_paths <- c(
  gamma_best_path, constant_best_path, gamma_B_path, constant_B_path,
  paramlist_path, observed_path, simulated_path
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

gamma_paths <- read.csv(gamma_B_path, check.names = FALSE)
constant_paths <- read.csv(constant_B_path, check.names = FALSE)
gamma_paths <- gamma_paths[gamma_paths$task_id == task_id, , drop = FALSE]
constant_paths <- constant_paths[
  constant_paths$task_id == task_id, , drop = FALSE
]
observed_data <- read.csv(observed_path, check.names = FALSE)
simulated_data <- read.csv(simulated_path, check.names = FALSE)

required_observed <- c("week", "reports")
required_simulated <- c("week", "I")
required_B <- c(
  "task_id", "simulation_seed", "observed_data_md5", "week", "B_estimate",
  "B_true"
)
if (!all(required_observed %in% names(observed_data))) {
  stop("Task-1 observed data lack week or reports.")
}
if (!all(required_simulated %in% names(simulated_data))) {
  stop("Task-1 simulated data lack week or I.")
}
if (!all(required_B %in% names(gamma_paths)) ||
    !all(required_B %in% names(constant_paths))) {
  stop("Task-1 combined B paths lack required columns.")
}

observed_data <- observed_data[, required_observed, drop = FALSE]
simulated_data <- simulated_data[, required_simulated, drop = FALSE]
expected_times <- seq(
  from = experiment_config$observation_interval,
  to = experiment_config$n_weeks,
  by = experiment_config$observation_interval
)
if (length(expected_times) != 70L) stop("Expected 70 observation times.")
assert_grid(observed_data$week, expected_times, "Observed-data times")
assert_grid(simulated_data$week, expected_times, "Simulated-data times")
assert_grid(gamma_paths$week, expected_times, "Gamma B-path times")
assert_grid(constant_paths$week, expected_times, "Constant-B path times")

selected_data_md5 <- unname(tools::md5sum(observed_path))
if (!identical(selected_data_md5, expected_data_md5)) {
  stop("Task-1 observed-data MD5 mismatch.")
}
fit_md5 <- c(
  as.character(gamma_row$observed_data_md5[[1]]),
  as.character(constant_row$observed_data_md5[[1]])
)
path_md5 <- c(
  unique(as.character(gamma_paths$observed_data_md5)),
  unique(as.character(constant_paths$observed_data_md5))
)
if (!all(c(fit_md5, path_md5) == expected_data_md5)) {
  stop("Gamma and constant-B task-1 records do not reference the shared data.")
}
simulation_seeds <- c(
  gamma_row$simulation_seed[[1]], constant_row$simulation_seed[[1]],
  seed_row$simulation_seed[[1]], unique(gamma_paths$simulation_seed),
  unique(constant_paths$simulation_seed)
)
if (length(unique(simulation_seeds)) != 1L || simulation_seeds[[1]] != 1001L) {
  stop("Task-1 records do not share simulation seed 1001.")
}
if (!identical(as.character(gamma_row$status[[1]]), "success") ||
    !identical(as.character(constant_row$status[[1]]), "success")) {
  stop("Task-1 best-fit status is not successful for both models.")
}
if (experiment_config$Np_final != Np) {
  stop("Canonical Experiment 4 Np_final is not 50000.")
}

assert_close(gamma_row$B0_hat[[1]], 4.35032977730012, "Gamma B0_hat")
assert_close(
  gamma_row$sigma_beta_hat[[1]], 0.206674783764282,
  "Gamma sigma_beta_hat"
)
assert_close(
  constant_row$Beta_hat[[1]], 4.04534940404491,
  "Constant-B Beta_hat"
)

if (any(gamma_paths$B_true != true_B_at_times(gamma_paths$week, experiment_config)) ||
    any(constant_paths$B_true !=
        true_B_at_times(constant_paths$week, experiment_config))) {
  stop("Task-1 saved B truth does not implement the week-5 switch.")
}
if (length(unique(constant_paths$B_estimate)) != 1L) {
  stop("Task-1 constant-B path is not a repeated static estimate.")
}
assert_close(
  unique(constant_paths$B_estimate), constant_row$Beta_hat[[1]],
  "Repeated constant-B estimate"
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
  gamma_model,
  params = theta_gamma,
  Np = Np,
  filter.mean = TRUE
)
gamma_fm <- filter_mean(gamma_pf)

set.seed(constant_final_pf_seed)
constant_pf <- pfilter(
  constant_model,
  params = theta_constant,
  Np = Np,
  filter.mean = TRUE
)
constant_fm <- filter_mean(constant_pf)

if (!all(c("B", "I") %in% rownames(gamma_fm)) ||
    !"I" %in% rownames(constant_fm)) {
  stop("Final particle filters did not return the required filtering means.")
}
assert_grid(as.numeric(time(gamma_pf)), expected_times, "Gamma filter times")
assert_grid(as.numeric(time(constant_pf)), expected_times, "Constant filter times")

gamma_B_reproduced <- as.numeric(gamma_fm["B", ])
if (max(abs(gamma_B_reproduced - gamma_paths$B_estimate)) > 1e-10) {
  stop("Reproduced task-1 Gamma B filtering mean differs from the saved path.")
}
assert_close(
  as.numeric(logLik(gamma_pf)), gamma_row$final_pf_logLik[[1]],
  "Gamma final-particle-filter log likelihood", tolerance = 1e-8
)
assert_close(
  as.numeric(logLik(constant_pf)), constant_row$final_pf_logLik[[1]],
  "Constant-B final-particle-filter log likelihood", tolerance = 1e-8
)

B_comparison <- data.frame(
  experiment = 4L,
  task_id = task_id,
  week = expected_times,
  B_true = as.numeric(gamma_paths$B_true),
  B_gamma_filtered_mean = gamma_B_reproduced,
  B_constant = as.numeric(constant_paths$B_estimate),
  stringsAsFactors = FALSE
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

if (any(!is.finite(unlist(B_comparison[-c(1L, 2L)]))) ||
    any(!is.finite(unlist(infectious_comparison[-c(1L, 2L)])))) {
  stop("Task-1 comparison artifacts contain non-finite values.")
}
if (any(infectious_comparison[, c(
  "I_true", "I_gamma_filtered_mean", "I_constant_filtered_mean"
)] < 0)) {
  stop("Task-1 infectious paths contain negative values.")
}

result_dir <- output_path("results", "selected_trajectory")
figure_dir <- output_path("figures", "comparison")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

B_csv <- file.path(result_dir, "experiment4_task1_B_filter_comparison.csv")
I_csv <- file.path(
  result_dir, "experiment4_task1_infectious_filter_comparison.csv"
)
provenance_path <- file.path(
  result_dir, "experiment4_task1_comparison_provenance.txt"
)
write.csv(B_comparison, B_csv, row.names = FALSE)
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

theme_task1 <- function() {
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

common_scales <- function(labels) {
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

make_B_plot <- function() {
  B_truth <- data.frame(
    week = c(0, 5, 5, 10),
    value = c(4, 4, 2, 2),
    series = factor("Truth", levels = series_levels)
  )
  B_gamma <- data.frame(
    week = B_comparison$week,
    value = B_comparison$B_gamma_filtered_mean,
    series = factor("Gamma-noise", levels = series_levels)
  )
  B_constant <- data.frame(
    week = c(0, 10),
    value = rep(B_comparison$B_constant[[1]], 2L),
    series = factor("Constant B", levels = series_levels)
  )

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
      "text", x = 5.10, y = 5.14, label = "Week 5",
      hjust = 0, vjust = 1, size = 2.55, family = "Arial",
      colour = "#6F6F6F"
    ) +
    common_scales(c("True B(t)", "Gamma-noise", "Constant B")) +
    scale_x_continuous(
      breaks = seq(0, 10, by = 2), limits = c(0, 10),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = seq(2, 5, by = 1), limits = c(1.6, 5.2),
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
    theme_task1()
}

make_infectious_plot <- function() {
  y_max <- max(
    infectious_comparison$I_true,
    infectious_comparison$I_gamma_filtered_mean,
    infectious_comparison$I_constant_filtered_mean
  )
  y_upper <- ceiling(1.04 * y_max / 100) * 100
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
    common_scales(c("True I(t)", "Gamma-noise", "Constant B")) +
    scale_x_continuous(
      breaks = seq(0, 10, by = 2), limits = c(0, 10),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = seq(0, y_upper, by = 100), limits = c(0, y_upper),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Week", y = expression("Infectious individuals  " * I(t))) +
    guides(
      colour = guide_legend(
        nrow = 1, byrow = TRUE,
        override.aes = list(
          linewidth = 0.9,
          linetype = unname(series_linetypes[series_levels])
        )
      )
    ) +
    theme_task1()
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
    svg(paste0(stem, ".svg"), width = width, height = height,
        pointsize = 8.5, onefile = TRUE, family = "Arial")
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
    png(paste0(stem, ".png"), width = width, height = height,
        units = "in", res = 300, type = "cairo", pointsize = 8.5)
    print(plot)
    dev.off()
    tiff(paste0(stem, ".tiff"), width = width, height = height,
         units = "in", res = dpi, compression = "lzw", type = "cairo",
         pointsize = 8.5)
    print(plot)
    dev.off()
  }
}

B_figure_stem <- file.path(
  figure_dir, "01_selected_task_B_trajectory_comparison"
)
I_figure_stem <- file.path(
  figure_dir, "07_task1_infectious_path_comparison"
)
B_plot <- make_B_plot()
I_plot <- make_infectious_plot()
export_figure(B_plot, B_figure_stem)
export_figure(I_plot, I_figure_stem)

repo_root <- normalizePath(file.path(getwd(), "..", ".."), mustWork = TRUE)
repo_commit <- system2(
  "git", c("-C", repo_root, "rev-parse", "HEAD"), stdout = TRUE
)
generated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
provenance <- c(
  "experiment: 4",
  paste0("selected_task_id: ", task_id),
  "selection_basis: task 1 requested for the manuscript comparison figures",
  paste0("simulation_seed: ", gamma_row$simulation_seed[[1]]),
  paste0("observed_data_md5: ", selected_data_md5),
  paste0("observed_data_sha256: ", sha256_file(observed_path)),
  "Gamma_and_constant_records_match_selected_data_md5: yes",
  paste0("particle_count_per_filter: ", Np),
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
  "B_gamma_semantics: observation-time filtering mean from the saved task-1 final particle-filter seed",
  "B_constant_semantics: fitted static task-1 estimate repeated over the 70 observation times",
  "I_gamma_semantics: observation-time filtering mean from the Gamma-noise model",
  "I_constant_semantics: observation-time filtering mean from the constant-B model",
  "I_true_semantics: task-1 latent infectious state from the shared simulated data",
  paste0(
    "I_gamma_RMSE: ",
    format(sqrt(mean(
      (infectious_comparison$I_gamma_filtered_mean -
         infectious_comparison$I_true)^2
    )), digits = 17)
  ),
  paste0(
    "I_constant_RMSE: ",
    format(sqrt(mean(
      (infectious_comparison$I_constant_filtered_mean -
         infectious_comparison$I_true)^2
    )), digits = 17)
  ),
  paste0(
    "I_gamma_constant_correlation: ",
    format(cor(
      infectious_comparison$I_gamma_filtered_mean,
      infectious_comparison$I_constant_filtered_mean
    ), digits = 17)
  ),
  "parameter_uncertainty_integrated: no",
  "across_task_average: no",
  "MIF2_rerun: no",
  "independent_likelihood_evaluations_rerun: no",
  paste0("number_of_observations: ", nrow(B_comparison)),
  paste0("R_version: ", R.version.string),
  paste0("pomp_version: ", as.character(packageVersion("pomp"))),
  paste0("RNGkind: ", paste(rng_kind, collapse = ", ")),
  paste0("generated_at_utc: ", generated_at_utc),
  paste0("repository_base_commit: ", repo_commit),
  paste0("B_comparison_csv_sha256: ", sha256_file(B_csv)),
  paste0("infectious_comparison_csv_sha256: ", sha256_file(I_csv))
)
writeLines(provenance, provenance_path, useBytes = TRUE)

cat(
  "Experiment 4 task-1 comparison figures generated successfully.\n",
  "B comparison data: ", B_csv, "\n",
  "Infectious comparison data: ", I_csv, "\n",
  "Provenance: ", provenance_path, "\n",
  "B figure: ", paste0(B_figure_stem, ".pdf"), "\n",
  "Infectious figure: ", paste0(I_figure_stem, ".pdf"), "\n",
  sep = ""
)

# ============================================================
# Generate Task 1 filtering-mean and forward-simulation figures
#
# Figure contract
# - Figure 6: data-conditioned particle filtering means for B(t) and I(t).
# - Figure 9: one joint forward simulation at fitted parameters for B(t)
#   and the corresponding I(t), without observation conditioning.
#
# Run from the Experiment 4 root:
#   Rscript code/11_generate_task1_filtering_mean_forward_figures.R \
#     [shared_data_root] [figure_root] [source_data_root]
# ============================================================

options(stringsAsFactors = FALSE, digits = 17)

suppressPackageStartupMessages({
  library(pomp)
  library(ggplot2)
  library(patchwork)
})

source(file.path("config", "experiment_config.R"))
source(file.path("code", "model_components.R"))
source(file.path("code", "io_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
shared_data_root <- if (length(args) >= 1L) args[[1]] else "shared_data"
figure_root <- if (length(args) >= 2L) {
  args[[2]]
} else {
  file.path("figures", "comparison")
}
source_data_root <- if (length(args) >= 3L) {
  args[[3]]
} else {
  file.path("results", "selected_trajectory")
}

task_id <- 1L
config <- experiment_config
expected_times <- seq(
  from = config$observation_interval,
  to = config$n_weeks,
  by = config$observation_interval
)

if (length(expected_times) != 70L) {
  stop("The configured observation grid must contain 70 times.")
}

read_task_best <- function(model_name) {
  path <- file.path(
    "results", "combined", model_name, "combined_best_fit_summary.csv"
  )
  data <- read.csv(path, check.names = FALSE)
  row <- data[data$task_id == task_id, , drop = FALSE]
  if (nrow(row) != 1L || !identical(as.character(row$status[[1]]), "success")) {
    stop("Expected one successful Task 1 best fit for ", model_name, ".")
  }
  row
}

filter_states <- function(model, theta, seed, expected_loglik, states) {
  set.seed(seed)
  pf <- pfilter(
    model,
    params = theta,
    Np = config$Np_final,
    filter.mean = TRUE
  )
  means <- filter_mean(pf)
  filter_times <- as.numeric(time(pf))
  missing_states <- setdiff(states, rownames(means))
  if (length(missing_states) > 0L) {
    stop("Filtering means lack state(s): ", paste(missing_states, collapse = ", "))
  }
  if (length(filter_times) != 70L ||
      max(abs(filter_times - expected_times)) > 1e-12 ||
      is.unsorted(filter_times, strictly = TRUE) || anyDuplicated(filter_times)) {
    stop("Filtering means do not use the configured 70-time grid.")
  }
  values <- means[states, , drop = FALSE]
  if (any(!is.finite(values))) {
    stop("Filtering means contain non-finite values.")
  }
  final_loglik <- as.numeric(logLik(pf))
  if (!is.finite(final_loglik)) {
    stop("Final particle-filter log likelihood is not finite.")
  }
  list(times = filter_times, means = values, loglik = final_loglik)
}

shared <- read_shared_data(shared_data_root, task_id)
observed <- shared$observed_data[, c("week", "reports"), drop = FALSE]
truth <- shared$simulated_data
if (nrow(observed) != 70L || nrow(truth) != 70L ||
    max(abs(observed$week - expected_times)) > 1e-12 ||
    max(abs(truth$week - expected_times)) > 1e-12) {
  stop("Task 1 shared data do not use the configured observation grid.")
}

gamma_best <- read_task_best("gamma")
constant_best <- read_task_best("constant")
if (!identical(
  as.character(gamma_best$observed_data_md5[[1]]),
  shared$observed_data_md5
) || !identical(
  as.character(constant_best$observed_data_md5[[1]]),
  shared$observed_data_md5
)) {
  stop("Task 1 shared data do not match the fitted-model records.")
}

paramlist <- read.csv(file.path("results", "paramlist.csv"), check.names = FALSE)
seed_row <- paramlist[paramlist$task_id == task_id, , drop = FALSE]
if (nrow(seed_row) != 1L) stop("Expected one Task 1 seed row.")

gamma_theta <- gamma_baseline_parameters(config)
gamma_theta[["B0"]] <- gamma_best$B0_hat[[1]]
gamma_theta[["sigma_beta"]] <- gamma_best$sigma_beta_hat[[1]]
constant_theta <- constant_baseline_parameters(config)
constant_theta[["Beta"]] <- constant_best$Beta_hat[[1]]

gamma_filter_seed <- as.integer(seed_row$gamma_final_pf_seed[[1]])
constant_filter_seed <- as.integer(seed_row$constant_final_pf_seed[[1]])

gamma_filtered <- filter_states(
  make_gamma_model(observed, config),
  gamma_theta,
  gamma_filter_seed,
  gamma_best$final_pf_logLik[[1]],
  c("B", "I")
)
constant_filtered <- filter_states(
  make_constant_model(observed, config),
  constant_theta,
  constant_filter_seed,
  constant_best$final_pf_logLik[[1]],
  "I"
)

# Reproduce the already-retained Task 1 B filtering mean exactly.
retained_B <- read.csv(
  file.path(
    "results", "combined", "gamma", "combined_B_filtering_means.csv"
  ),
  check.names = FALSE
)
retained_B <- retained_B[retained_B$task_id == task_id, , drop = FALSE]
retained_B_max_abs_diff <- max(
  abs(retained_B$B_estimate - gamma_filtered$means["B", ])
)
if (nrow(retained_B) != 70L ||
    max(abs(retained_B$week - gamma_filtered$times)) > 1e-12) {
  stop("Retained Task 1 B filtering mean uses an unexpected time grid.")
}

series_labels_filter <- c(
  truth = "Truth",
  gamma = "Gamma-noise filtering mean",
  constant = "Constant-B model"
)
series_labels_forward <- c(
  truth = "Truth",
  gamma = "Gamma-noise forward simulation",
  constant = "Constant-B forward simulation"
)
palette <- c(truth = "#252525", gamma = "#0072B2", constant = "#D55E00")
line_types <- c(truth = "solid", gamma = "longdash", constant = "dotdash")

make_rows <- function(week, value, quantity, series_key, semantics,
                      conditioning, time_zero_semantics = NA_character_) {
  data.frame(
    experiment = 4L,
    task_id = task_id,
    week = as.numeric(week),
    quantity = quantity,
    series_key = series_key,
    value = as.numeric(value),
    is_time_zero = abs(as.numeric(week)) < 1e-12,
    estimate_semantics = semantics,
    conditioning = conditioning,
    time_zero_semantics = time_zero_semantics,
    stringsAsFactors = FALSE
  )
}

# The duplicate week-5 rows draw the prescribed right-continuous step exactly.
true_B_display <- make_rows(
  c(0, 5, 5, config$n_weeks),
  c(
    config$true_parameters[["Beta_high"]],
    config$true_parameters[["Beta_high"]],
    config$true_parameters[["Beta_low"]],
    config$true_parameters[["Beta_low"]]
  ),
  "B", "truth", "data_generating_path", "none"
)
true_I_display <- make_rows(
  c(0, truth$week),
  c(10, truth$I),
  "I", "truth", "data_generating_path", "none"
)

filtering_source <- rbind(
  true_B_display,
  make_rows(
    0, gamma_theta[["B0"]], "B", "gamma", "fitted_initial_value",
    "Y_1:70_for_parameters_only", "fitted_B0"
  ),
  make_rows(
    retained_B$week, retained_B$B_estimate, "B", "gamma",
    "particle_filtering_mean", "Y_1:n"
  ),
  make_rows(
    c(0, expected_times), rep(constant_theta[["Beta"]], 71L),
    "B", "constant", "fitted_static_parameter", "Y_1:70"
  ),
  true_I_display,
  make_rows(
    0, 10, "I", "gamma", "fixed_initial_state", "none", "fixed_I0"
  ),
  make_rows(
    gamma_filtered$times, gamma_filtered$means["I", ], "I", "gamma",
    "particle_filtering_mean", "Y_1:n"
  ),
  make_rows(
    0, 10, "I", "constant", "fixed_initial_state", "none", "fixed_I0"
  ),
  make_rows(
    constant_filtered$times, constant_filtered$means["I", ],
    "I", "constant", "particle_filtering_mean", "Y_1:n"
  )
)
filtering_source$series_label <- unname(
  series_labels_filter[filtering_source$series_key]
)

# Seeds are fixed by a prespecified offset from the saved final-PF seeds.
# Only the first realization from each seed is used; no trajectory selection occurs.
gamma_forward_seed <- gamma_filter_seed + 100000L
constant_forward_seed <- constant_filter_seed + 100000L
simulation_template <- make_observation_template(config)

set.seed(gamma_forward_seed)
gamma_forward <- simulate(
  make_gamma_model(simulation_template, config),
  params = gamma_theta,
  nsim = 1,
  format = "data.frame",
  include.data = FALSE
)
set.seed(constant_forward_seed)
constant_forward <- simulate(
  make_constant_model(simulation_template, config),
  params = constant_theta,
  nsim = 1,
  format = "data.frame",
  include.data = FALSE
)

if (nrow(gamma_forward) != 70L || nrow(constant_forward) != 70L ||
    max(abs(gamma_forward$week - expected_times)) > 1e-12 ||
    max(abs(constant_forward$week - expected_times)) > 1e-12 ||
    any(!is.finite(gamma_forward$B)) ||
    any(!is.finite(gamma_forward$I)) ||
    any(!is.finite(constant_forward$I))) {
  stop("Forward simulations failed the observation-grid or finite-value checks.")
}

forward_source <- rbind(
  true_B_display,
  make_rows(
    c(0, gamma_forward$week),
    c(gamma_theta[["B0"]], gamma_forward$B),
    "B", "gamma", "forward_simulation_at_fitted_parameters",
    "none_after_t0", "fitted_B0"
  ),
  make_rows(
    c(0, expected_times), rep(constant_theta[["Beta"]], 71L),
    "B", "constant", "forward_simulation_at_fitted_parameters",
    "none_after_t0", "fitted_Beta"
  ),
  true_I_display,
  make_rows(
    c(0, gamma_forward$week), c(10, gamma_forward$I),
    "I", "gamma", "forward_simulation_at_fitted_parameters",
    "none_after_t0", "fixed_I0"
  ),
  make_rows(
    c(0, constant_forward$week), c(10, constant_forward$I),
    "I", "constant", "forward_simulation_at_fitted_parameters",
    "none_after_t0", "fixed_I0"
  )
)
forward_source$series_label <- unname(
  series_labels_forward[forward_source$series_key]
)

dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_root, recursive = TRUE, showWarnings = FALSE)

filtering_csv <- file.path(
  source_data_root, "experiment4_task1_filtering_mean_trajectories.csv"
)
forward_csv <- file.path(
  source_data_root, "experiment4_task1_forward_simulation_trajectories.csv"
)
provenance_csv <- file.path(
  source_data_root, "experiment4_task1_filtering_forward_provenance.csv"
)
write.csv(filtering_source, filtering_csv, row.names = FALSE)
write.csv(forward_source, forward_csv, row.names = FALSE)

provenance <- data.frame(
  experiment = 4L,
  task_id = task_id,
  observed_data_md5 = shared$observed_data_md5,
  simulated_data_md5 = shared$simulated_data_md5,
  B0_hat = gamma_theta[["B0"]],
  sigma_beta_hat = gamma_theta[["sigma_beta"]],
  constant_B_hat = constant_theta[["Beta"]],
  Np_filter = config$Np_final,
  gamma_filter_seed = gamma_filter_seed,
  constant_filter_seed = constant_filter_seed,
  gamma_filter_logLik = gamma_filtered$loglik,
  constant_filter_logLik = constant_filtered$loglik,
  saved_gamma_final_pf_logLik = gamma_best$final_pf_logLik[[1]],
  saved_constant_final_pf_logLik = constant_best$final_pf_logLik[[1]],
  retained_B_max_abs_diff_current_rerun = retained_B_max_abs_diff,
  gamma_B_filtering_mean_source = "retained_primary_metric_output",
  infectious_filtering_mean_source = "rerun_final_particle_filters",
  R_version = R.version.string,
  pomp_version = as.character(packageVersion("pomp")),
  gamma_forward_seed = gamma_forward_seed,
  constant_forward_seed = constant_forward_seed,
  forward_selection_rule = "first_realization_from_prespecified_seed_offset",
  figure6_conditioning = "Y_1:n_at_each_observation_time",
  figure9_conditioning = "none_after_fitted_parameters_are_fixed",
  parameter_uncertainty_integrated = FALSE,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
write.csv(provenance, provenance_csv, row.names = FALSE)

theme_trajectory <- function() {
  theme_classic(base_size = 8.2, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#252525"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#252525"),
      axis.text = element_text(size = 7.2, colour = "#252525"),
      axis.title = element_text(size = 8.2, colour = "#252525"),
      plot.title = element_text(size = 8.4, face = "bold", hjust = 0),
      plot.tag = element_text(size = 9, face = "bold"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 7.1),
      legend.key.width = grid::unit(12, "mm"),
      legend.margin = margin(0, 0, 2, 0),
      plot.margin = margin(3, 5, 2, 5)
    )
}

make_panel <- function(data, quantity, title, y_title, labels, show_x = TRUE) {
  panel <- data[data$quantity == quantity, , drop = FALSE]
  path_data <- panel[!(panel$is_time_zero & panel$series_key != "truth"), ]
  initial_data <- panel[panel$is_time_zero & panel$series_key != "truth", ]

  p <- ggplot(
    path_data,
    aes(
      x = week, y = value, colour = series_key,
      linetype = series_key, group = series_key
    )
  ) +
    geom_vline(
      xintercept = config$true_parameters[["t_switch"]],
      colour = "#B5B5B5", linewidth = 0.35, linetype = "dotted"
    ) +
    geom_path(linewidth = 0.72, lineend = "round") +
    geom_point(
      data = initial_data,
      aes(x = week, y = value, colour = series_key),
      inherit.aes = FALSE,
      shape = 21, fill = "white", stroke = 0.55, size = 2.1,
      show.legend = FALSE
    ) +
    annotate(
      "text", x = config$true_parameters[["t_switch"]] + 0.10,
      y = Inf, label = "Week 5", vjust = 1.35, hjust = 0,
      colour = "#707070", family = "Helvetica", size = 2.45
    ) +
    scale_colour_manual(values = palette, breaks = names(labels), labels = labels) +
    scale_linetype_manual(
      values = line_types, breaks = names(labels), labels = labels
    ) +
    scale_x_continuous(
      limits = c(0, config$n_weeks), breaks = seq(0, config$n_weeks, by = 2),
      expand = expansion(mult = c(0, 0.015))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
    labs(
      title = title,
      x = if (show_x) "Time (weeks)" else NULL,
      y = y_title,
      colour = NULL,
      linetype = NULL
    ) +
    guides(
      colour = guide_legend(nrow = 1, byrow = TRUE),
      linetype = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    coord_cartesian(clip = "off") +
    theme_trajectory()

  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  p
}

assemble_figure <- function(source, labels) {
  p_B <- make_panel(
    source, "B", "Transmission rate, B(t)", "Transmission rate, B(t)",
    labels, show_x = FALSE
  )
  p_I <- make_panel(
    source, "I", "Infectious population, I(t)", "Infectious individuals, I(t)",
    labels, show_x = TRUE
  )
  p_B / p_I +
    plot_layout(guides = "collect", heights = c(1, 1.08)) +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "top")
}

save_figure <- function(plot, stem, width_mm = 183, height_mm = 150, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  raw_pdf <- paste0(stem, "_cairo_raw.pdf")
  final_pdf <- paste0(stem, ".pdf")

  grDevices::cairo_pdf(
    raw_pdf, width = width_in, height = height_in,
    family = "Helvetica", onefile = TRUE
  )
  print(plot)
  grDevices::dev.off()

  if (!nzchar(Sys.which("gs"))) {
    stop("Ghostscript is required to normalize the publication PDF export.")
  }
  gs_status <- system2(
    "gs",
    c(
      "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
      "-dCompatibilityLevel=1.5", "-dPDFSETTINGS=/prepress",
      paste0("-sOutputFile=", shQuote(final_pdf)),
      shQuote(raw_pdf)
    )
  )
  if (!identical(gs_status, 0L) || !file.exists(final_pdf)) {
    stop("Ghostscript failed to normalize the publication PDF export.")
  }
  unlink(raw_pdf)

  svglite::svglite(
    paste0(stem, ".svg"), width = width_in, height = height_in,
    system_fonts = list(sans = "Helvetica")
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300, background = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw", background = "white"
  )
  print(plot)
  grDevices::dev.off()
}

figure6 <- assemble_figure(filtering_source, series_labels_filter)
figure9 <- assemble_figure(forward_source, series_labels_forward)

figure6_stem <- file.path(
  figure_root, "01_task1_particle_filtering_mean_trajectories"
)
figure9_stem <- file.path(
  figure_root, "09_task1_forward_simulations_at_fitted_parameters"
)
save_figure(figure6, figure6_stem)
save_figure(figure9, figure9_stem)

cat(
  "Generated Task 1 filtering-mean and forward-simulation figures.\n",
  "Figure 6 source: ", filtering_csv, "\n",
  "Figure 9 source: ", forward_csv, "\n",
  "Provenance: ", provenance_csv, "\n",
  "Figure 6 stem: ", figure6_stem, "\n",
  "Figure 9 stem: ", figure9_stem, "\n",
  sep = ""
)

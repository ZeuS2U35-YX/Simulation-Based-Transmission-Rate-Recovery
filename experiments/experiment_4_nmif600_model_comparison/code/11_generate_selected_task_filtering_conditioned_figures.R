# ============================================================
# Generate selected-task filtering-mean and conditioned-trajectory figures
#
# Figure 6: observation-time particle filtering means for B(t) and I(t)
#           in the previously designated Tasks 1 and 117.
# Figure 9: one ancestry-preserving, observation-conditioned particle
#           trajectory for B(t) and I(t) from each fitted model and task.
#
# Run from the Experiment 4 root:
#   Rscript code/11_generate_selected_task_filtering_conditioned_figures.R \
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

task_ids <- c(1L, 117L)
task_labels <- c(`1` = "Task 1", `117` = "Task 117")
task_roles <- c(
  `1` = "previously_designated_primary_illustration",
  `117` = "previously_designated_second_illustration"
)

expected_times <- seq(
  from = experiment_config$observation_interval,
  to = experiment_config$n_weeks,
  by = experiment_config$observation_interval
)
if (length(expected_times) != 70L) {
  stop("The configured observation grid must contain 70 times.")
}

gamma_best_all <- read.csv(
  file.path("results", "combined", "gamma", "combined_best_fit_summary.csv"),
  check.names = FALSE
)
constant_best_all <- read.csv(
  file.path("results", "combined", "constant", "combined_best_fit_summary.csv"),
  check.names = FALSE
)
paramlist <- read.csv(file.path("results", "paramlist.csv"), check.names = FALSE)
retained_B_all <- read.csv(
  file.path(
    "results", "combined", "gamma", "combined_B_filtering_means.csv"
  ),
  check.names = FALSE
)
task_metrics <- read.csv(
  file.path("results", "comparison", "model_task_metrics.csv"),
  check.names = FALSE
)

read_task_best <- function(data, model_name, task_id) {
  row <- data[data$task_id == task_id, , drop = FALSE]
  if (nrow(row) != 1L || !identical(as.character(row$status[[1]]), "success")) {
    stop("Expected one successful Task ", task_id, " best fit for ", model_name, ".")
  }
  row
}

extract_trajectory <- function(pf, states) {
  trajectory <- filter_traj(pf, vars = states, format = "data.frame")
  if (!all(c("name", "rep", "time", "value") %in% names(trajectory))) {
    stop("filter_traj returned an unexpected schema.")
  }
  if (!setequal(unique(as.character(trajectory$name)), states) ||
      !identical(unique(as.integer(trajectory$rep)), 1L)) {
    stop("filter_traj returned unexpected states or replicate identifiers.")
  }

  expected_with_t0 <- c(0, expected_times)
  for (state in states) {
    state_path <- trajectory[trajectory$name == state, , drop = FALSE]
    if (nrow(state_path) != 71L || any(!is.finite(state_path$value)) ||
        max(abs(state_path$time - expected_with_t0)) > 1e-12 ||
        is.unsorted(state_path$time, strictly = TRUE) ||
        anyDuplicated(state_path$time)) {
      stop("Conditioned trajectory failed validation for state ", state, ".")
    }
  }
  trajectory
}

run_final_filter <- function(model, theta, seed, mean_states, path_states) {
  set.seed(seed)
  pf <- pfilter(
    model,
    params = theta,
    Np = experiment_config$Np_final,
    filter.mean = TRUE,
    filter.traj = TRUE
  )

  means <- filter_mean(pf)
  filter_times <- as.numeric(time(pf))
  if (length(filter_times) != 70L ||
      max(abs(filter_times - expected_times)) > 1e-12 ||
      is.unsorted(filter_times, strictly = TRUE) || anyDuplicated(filter_times)) {
    stop("Filtering means do not use the configured 70-time grid.")
  }
  missing_states <- setdiff(mean_states, rownames(means))
  if (length(missing_states) > 0L) {
    stop("Filtering means lack state(s): ", paste(missing_states, collapse = ", "))
  }
  mean_values <- means[mean_states, , drop = FALSE]
  if (any(!is.finite(mean_values))) {
    stop("Filtering means contain non-finite values.")
  }

  final_loglik <- as.numeric(logLik(pf))
  if (!is.finite(final_loglik)) {
    stop("Final particle-filter log likelihood is not finite.")
  }

  list(
    times = filter_times,
    means = mean_values,
    trajectory = extract_trajectory(pf, path_states),
    loglik = final_loglik
  )
}

make_rows <- function(
    task_id, week, value, quantity, series_key, semantics, conditioning,
    time_zero_semantics = NA_character_) {
  data.frame(
    experiment = 4L,
    task_id = task_id,
    task_label = unname(task_labels[as.character(task_id)]),
    illustration_role = unname(task_roles[as.character(task_id)]),
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

make_true_B <- function(task_id) {
  make_rows(
    task_id,
    c(0, 5, 5, experiment_config$n_weeks),
    c(
      experiment_config$true_parameters[["Beta_high"]],
      experiment_config$true_parameters[["Beta_high"]],
      experiment_config$true_parameters[["Beta_low"]],
      experiment_config$true_parameters[["Beta_low"]]
    ),
    "B", "truth", "data_generating_path", "none"
  )
}

filtering_parts <- list()
conditioned_parts <- list()
provenance_parts <- list()

for (task_id in task_ids) {
  shared <- read_shared_data(shared_data_root, task_id)
  observed <- shared$observed_data[, c("week", "reports"), drop = FALSE]
  truth <- shared$simulated_data
  if (nrow(observed) != 70L || nrow(truth) != 70L ||
      max(abs(observed$week - expected_times)) > 1e-12 ||
      max(abs(truth$week - expected_times)) > 1e-12) {
    stop("Task ", task_id, " shared data do not use the configured grid.")
  }

  gamma_best <- read_task_best(gamma_best_all, "gamma", task_id)
  constant_best <- read_task_best(constant_best_all, "constant", task_id)
  if (!identical(
    as.character(gamma_best$observed_data_md5[[1]]),
    shared$observed_data_md5
  ) || !identical(
    as.character(constant_best$observed_data_md5[[1]]),
    shared$observed_data_md5
  )) {
    stop("Task ", task_id, " shared data do not match the fitted records.")
  }

  seed_row <- paramlist[paramlist$task_id == task_id, , drop = FALSE]
  if (nrow(seed_row) != 1L) {
    stop("Expected one seed row for Task ", task_id, ".")
  }

  gamma_theta <- gamma_baseline_parameters(experiment_config)
  gamma_theta[["B0"]] <- gamma_best$B0_hat[[1]]
  gamma_theta[["sigma_beta"]] <- gamma_best$sigma_beta_hat[[1]]
  constant_theta <- constant_baseline_parameters(experiment_config)
  constant_theta[["Beta"]] <- constant_best$Beta_hat[[1]]

  gamma_seed <- as.integer(seed_row$gamma_final_pf_seed[[1]])
  constant_seed <- as.integer(seed_row$constant_final_pf_seed[[1]])
  gamma_filtered <- run_final_filter(
    make_gamma_model(observed, experiment_config),
    gamma_theta,
    gamma_seed,
    c("B", "I"),
    c("B", "I")
  )
  constant_filtered <- run_final_filter(
    make_constant_model(observed, experiment_config),
    constant_theta,
    constant_seed,
    "I",
    "I"
  )

  retained_B <- retained_B_all[retained_B_all$task_id == task_id, , drop = FALSE]
  if (nrow(retained_B) != 70L ||
      max(abs(retained_B$week - gamma_filtered$times)) > 1e-12) {
    stop("Retained filtering mean has an unexpected grid for Task ", task_id, ".")
  }
  retained_B_max_abs_diff <- max(
    abs(retained_B$B_estimate - gamma_filtered$means["B", ])
  )

  true_B <- make_true_B(task_id)
  true_I <- make_rows(
    task_id, c(0, truth$week), c(10, truth$I), "I", "truth",
    "data_generating_path", "none"
  )

  filtering_parts[[as.character(task_id)]] <- rbind(
    true_B,
    make_rows(
      task_id, 0, gamma_theta[["B0"]], "B", "gamma",
      "fitted_initial_value", "Y_1:70_for_parameters_only", "fitted_B0"
    ),
    make_rows(
      task_id, retained_B$week, retained_B$B_estimate, "B", "gamma",
      "particle_filtering_mean", "Y_1:n"
    ),
    make_rows(
      task_id, c(0, expected_times), rep(constant_theta[["Beta"]], 71L),
      "B", "constant", "fitted_static_parameter", "Y_1:70",
      "fitted_Beta"
    ),
    true_I,
    make_rows(
      task_id, 0, 10, "I", "gamma", "fixed_initial_state", "none", "fixed_I0"
    ),
    make_rows(
      task_id, gamma_filtered$times, gamma_filtered$means["I", ],
      "I", "gamma", "particle_filtering_mean", "Y_1:n"
    ),
    make_rows(
      task_id, 0, 10, "I", "constant", "fixed_initial_state", "none", "fixed_I0"
    ),
    make_rows(
      task_id, constant_filtered$times, constant_filtered$means["I", ],
      "I", "constant", "particle_filtering_mean", "Y_1:n"
    )
  )

  gamma_B_path <- gamma_filtered$trajectory[
    gamma_filtered$trajectory$name == "B", , drop = FALSE
  ]
  gamma_I_path <- gamma_filtered$trajectory[
    gamma_filtered$trajectory$name == "I", , drop = FALSE
  ]
  constant_I_path <- constant_filtered$trajectory[
    constant_filtered$trajectory$name == "I", , drop = FALSE
  ]

  conditioned_parts[[as.character(task_id)]] <- rbind(
    true_B,
    make_rows(
      task_id, gamma_B_path$time, gamma_B_path$value, "B", "gamma",
      "observation_conditioned_ancestry_trajectory", "Y_1:70",
      "fitted_B0"
    ),
    make_rows(
      task_id, c(0, expected_times), rep(constant_theta[["Beta"]], 71L),
      "B", "constant", "fitted_static_parameter", "Y_1:70", "fitted_Beta"
    ),
    true_I,
    make_rows(
      task_id, gamma_I_path$time, gamma_I_path$value, "I", "gamma",
      "observation_conditioned_ancestry_trajectory", "Y_1:70", "fixed_I0"
    ),
    make_rows(
      task_id, constant_I_path$time, constant_I_path$value, "I", "constant",
      "observation_conditioned_ancestry_trajectory", "Y_1:70", "fixed_I0"
    )
  )

  gamma_metric <- task_metrics[
    task_metrics$task_id == task_id & task_metrics$model == "gamma_noise",
    , drop = FALSE
  ]
  constant_metric <- task_metrics[
    task_metrics$task_id == task_id & task_metrics$model == "constant_B",
    , drop = FALSE
  ]
  gamma_all <- task_metrics[task_metrics$model == "gamma_noise", , drop = FALSE]
  constant_all <- task_metrics[task_metrics$model == "constant_B", , drop = FALSE]

  provenance_parts[[as.character(task_id)]] <- data.frame(
    experiment = 4L,
    task_id = task_id,
    illustration_role = unname(task_roles[as.character(task_id)]),
    observed_data_md5 = shared$observed_data_md5,
    simulated_data_md5 = shared$simulated_data_md5,
    B0_hat = gamma_theta[["B0"]],
    sigma_beta_hat = gamma_theta[["sigma_beta"]],
    constant_B_hat = constant_theta[["Beta"]],
    Np_filter = experiment_config$Np_final,
    gamma_filter_seed = gamma_seed,
    constant_filter_seed = constant_seed,
    gamma_filter_logLik = gamma_filtered$loglik,
    constant_filter_logLik = constant_filtered$loglik,
    saved_gamma_final_pf_logLik = gamma_best$final_pf_logLik[[1]],
    saved_constant_final_pf_logLik = constant_best$final_pf_logLik[[1]],
    retained_B_max_abs_diff_current_rerun = retained_B_max_abs_diff,
    gamma_RMSE = gamma_metric$RMSE[[1]],
    gamma_RMSE_empirical_percentile = mean(gamma_all$RMSE <= gamma_metric$RMSE[[1]]),
    gamma_after_week5_error = gamma_metric$mean_error_after_5[[1]],
    gamma_abs_after_week5_error_empirical_percentile = mean(
      abs(gamma_all$mean_error_after_5) <= abs(gamma_metric$mean_error_after_5[[1]])
    ),
    constant_RMSE = constant_metric$RMSE[[1]],
    constant_RMSE_empirical_percentile = mean(
      constant_all$RMSE <= constant_metric$RMSE[[1]]
    ),
    filtering_mean_conditioning = "Y_1:n_at_each_observation_time",
    conditioned_trajectory_conditioning = "Y_1:70_via_final_particle_filter_genealogy",
    trajectory_selection_rule = paste0(
      "one_seeded_ancestry_trajectory_returned_by_filter_traj;",
      "no_visual_or_truth_based_screening"
    ),
    gamma_B_and_I_share_one_ancestry = TRUE,
    parameter_uncertainty_integrated = FALSE,
    R_version = R.version.string,
    pomp_version = as.character(packageVersion("pomp")),
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

filtering_source <- do.call(rbind, filtering_parts)
conditioned_source <- do.call(rbind, conditioned_parts)
provenance <- do.call(rbind, provenance_parts)
rownames(filtering_source) <- NULL
rownames(conditioned_source) <- NULL
rownames(provenance) <- NULL

series_labels_filter <- c(
  truth = "Truth",
  gamma = "Gamma-noise filtering mean",
  constant = "Constant-B model"
)
series_labels_conditioned <- c(
  truth = "Truth",
  gamma = "Gamma-noise conditioned trajectory",
  constant = "Constant-B conditioned trajectory"
)
palette <- c(truth = "#252525", gamma = "#0072B2", constant = "#D55E00")
line_types <- c(truth = "solid", gamma = "longdash", constant = "dotdash")
filtering_source$series_label <- unname(
  series_labels_filter[filtering_source$series_key]
)
conditioned_source$series_label <- unname(
  series_labels_conditioned[conditioned_source$series_key]
)

dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_root, recursive = TRUE, showWarnings = FALSE)

filtering_csv <- file.path(
  source_data_root, "experiment4_selected_tasks_filtering_mean_trajectories.csv"
)
conditioned_csv <- file.path(
  source_data_root,
  "experiment4_selected_tasks_observation_conditioned_trajectories.csv"
)
provenance_csv <- file.path(
  source_data_root,
  "experiment4_selected_tasks_filtering_conditioned_provenance.csv"
)
write.csv(filtering_source, filtering_csv, row.names = FALSE)
write.csv(conditioned_source, conditioned_csv, row.names = FALSE)
write.csv(provenance, provenance_csv, row.names = FALSE)

theme_trajectory <- function() {
  theme_classic(base_size = 7.4, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "#252525"),
      axis.ticks = element_line(linewidth = 0.32, colour = "#252525"),
      axis.text = element_text(size = 6.4, colour = "#252525"),
      axis.title = element_text(size = 7.3, colour = "#252525"),
      plot.title = element_text(size = 8.0, face = "bold", hjust = 0.5),
      plot.tag = element_text(size = 8.2, face = "bold"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 6.3),
      legend.key.width = grid::unit(9, "mm"),
      legend.margin = margin(0, 0, 1.5, 0),
      plot.margin = margin(2.5, 4, 2, 4)
    )
}

shared_limits <- function(source, quantity) {
  values <- source$value[source$quantity == quantity]
  limits <- range(values, finite = TRUE)
  span <- diff(limits)
  if (span <= 0) span <- max(abs(limits), 1)
  if (quantity == "I") {
    c(0, limits[[2]] + 0.09 * span)
  } else {
    c(limits[[1]] - 0.04 * span, limits[[2]] + 0.10 * span)
  }
}

make_panel <- function(
    source, task_id, quantity, labels, y_limits,
    show_x = TRUE, show_y_title = TRUE, show_title = FALSE) {
  panel <- source[
    source$task_id == task_id & source$quantity == quantity,
    , drop = FALSE
  ]
  path_data <- panel[!(panel$is_time_zero & panel$series_key != "truth"), ]
  initial_data <- panel[panel$is_time_zero & panel$series_key != "truth", ]

  plot <- ggplot(
    path_data,
    aes(
      x = week, y = value, colour = series_key,
      linetype = series_key, group = series_key
    )
  ) +
    geom_vline(
      xintercept = experiment_config$true_parameters[["t_switch"]],
      colour = "#B5B5B5", linewidth = 0.32, linetype = "dotted"
    ) +
    geom_path(linewidth = 0.66, lineend = "round") +
    geom_point(
      data = initial_data,
      aes(x = week, y = value, colour = series_key),
      inherit.aes = FALSE,
      shape = 21, fill = "white", stroke = 0.5, size = 1.7,
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = experiment_config$true_parameters[["t_switch"]] + 0.10,
      y = Inf, label = "Week 5", vjust = 1.30, hjust = 0,
      colour = "#707070", family = "Helvetica", size = 2.15
    ) +
    scale_colour_manual(values = palette, breaks = names(labels), labels = labels) +
    scale_linetype_manual(
      values = line_types, breaks = names(labels), labels = labels
    ) +
    scale_x_continuous(
      limits = c(0, experiment_config$n_weeks),
      breaks = seq(0, experiment_config$n_weeks, by = 2),
      expand = expansion(mult = c(0, 0.015))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    labs(
      title = if (show_title) unname(task_labels[as.character(task_id)]) else NULL,
      x = if (show_x) "Time (weeks)" else NULL,
      y = if (show_y_title) {
        if (quantity == "B") "Transmission rate, B(t)" else "Infectious individuals, I(t)"
      } else {
        NULL
      },
      colour = NULL,
      linetype = NULL
    ) +
    guides(
      colour = guide_legend(nrow = 1, byrow = TRUE),
      linetype = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    theme_trajectory()

  if (!show_x) {
    plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  plot
}

assemble_figure <- function(source, labels) {
  B_limits <- shared_limits(source, "B")
  I_limits <- shared_limits(source, "I")
  p_a <- make_panel(
    source, 1L, "B", labels, B_limits,
    show_x = FALSE, show_y_title = TRUE, show_title = TRUE
  )
  p_b <- make_panel(
    source, 117L, "B", labels, B_limits,
    show_x = FALSE, show_y_title = FALSE, show_title = TRUE
  )
  p_c <- make_panel(
    source, 1L, "I", labels, I_limits,
    show_x = TRUE, show_y_title = TRUE, show_title = FALSE
  )
  p_d <- make_panel(
    source, 117L, "I", labels, I_limits,
    show_x = TRUE, show_y_title = FALSE, show_title = FALSE
  )

  (p_a | p_b) / (p_c | p_d) +
    plot_layout(guides = "collect", heights = c(1, 1.04)) +
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
figure9 <- assemble_figure(conditioned_source, series_labels_conditioned)

figure6_stem <- file.path(
  figure_root, "01_selected_tasks_particle_filtering_mean_trajectories"
)
figure9_stem <- file.path(
  figure_root, "09_selected_tasks_observation_conditioned_trajectories"
)
save_figure(figure6, figure6_stem)
save_figure(figure9, figure9_stem)

cat(
  "Generated selected-task filtering-mean and conditioned-trajectory figures.\n",
  "Figure 6 source: ", filtering_csv, "\n",
  "Figure 9 source: ", conditioned_csv, "\n",
  "Provenance: ", provenance_csv, "\n",
  "Figure 6 stem: ", figure6_stem, "\n",
  "Figure 9 stem: ", figure9_stem, "\n",
  sep = ""
)

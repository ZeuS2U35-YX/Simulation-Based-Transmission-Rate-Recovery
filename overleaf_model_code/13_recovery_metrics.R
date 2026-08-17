# Task-level recovery metrics used in Experiment 4

calculate_metrics <- function(paths, model_name) {
  paths$error <- paths$B_estimate - paths$B_true
  paths$squared_error <- paths$error^2
  split_paths <- split(paths, paths$task_id)

  out <- do.call(rbind, lapply(split_paths, function(x) {
    mean_error <- mean(x$error)
    data.frame(
      task_id = x$task_id[[1]],
      model = model_name,
      RSS = sum(x$squared_error),
      RMSE = sqrt(mean(x$squared_error)),
      mean_error = mean_error,
      AOB = abs(mean_error),
      mean_error_before_5 = mean(x$error[
        x$week < experiment_config$true_parameters[["t_switch"]]
      ]),
      mean_error_after_5 = mean(x$error[
        x$week >= experiment_config$true_parameters[["t_switch"]]
      ])
    )
  }))

  rownames(out) <- NULL
  out
}

# AOB is computed within each task before paired comparison and aggregation.
paired$delta_AOB_gamma_minus_constant <-
  paired$AOB_gamma - paired$AOB_constant

gamma_mean_AOB <- mean(paired$AOB_gamma)
constant_mean_AOB <- mean(paired$AOB_constant)

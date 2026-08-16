# Task-level recovery metrics used in Experiment 4

calculate_metrics <- function(paths, model_name) {
  paths$error <- paths$B_estimate - paths$B_true
  paths$squared_error <- paths$error^2
  split_paths <- split(paths, paths$task_id)

  out <- do.call(rbind, lapply(split_paths, function(x) {
    data.frame(
      task_id = x$task_id[[1]],
      model = model_name,
      RSS = sum(x$squared_error),
      RMSE = sqrt(mean(x$squared_error)),
      bias_all = mean(x$error),
      bias_before = mean(x$error[
        x$week < experiment_config$true_parameters[["t_switch"]]
      ]),
      bias_after = mean(x$error[
        x$week >= experiment_config$true_parameters[["t_switch"]]
      ])
    )
  }))

  rownames(out) <- NULL
  out
}

# bias_all is the signed task-level mean error. With 70 observation times,
# abs(bias_all) = abs(sum(error)) / 70. The absolute quantity is applied before
# across-task comparison and aggregation.
paired$delta_abs_bias_gamma_minus_constant <-
  abs(paired$bias_all_gamma) - abs(paired$bias_all_constant)

gamma_absolute_mean_error <- mean(abs(paired$bias_all_gamma))
constant_absolute_mean_error <- mean(abs(paired$bias_all_constant))

# Sample one ancestry-preserving latent B(t) trajectory at the best fit

best <- valid[which.max(valid$logLik), , drop = FALSE]
theta_best <- theta_gamma
theta_best[["B0"]] <- best$B0_hat[[1]]
theta_best[["sigma_beta"]] <- best$sigma_beta_hat[[1]]

set.seed(final_pf_seed)
pf_best <- pfilter(
  gamma_model,
  params = theta_best,
  Np = config$Np_final,
  filter.traj = TRUE
)

sampled <- filter_traj(
  pf_best,
  vars = "B",
  format = "data.frame"
)

# Metrics use the 70 observation times; time zero is retained only for display.
sampled_observation_times <- sampled[-1L, , drop = FALSE]
B_path <- data.frame(
  week = sampled_observation_times$time,
  B_estimate = sampled_observation_times$value,
  B_true = true_B_at_times(sampled_observation_times$time, config),
  trajectory_seed = final_pf_seed,
  path_semantics = "ancestry_preserving_sampled_latent_trajectory"
)

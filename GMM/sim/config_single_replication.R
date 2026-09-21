# Initial beta[j] = true_beta[j] + Normal(0, start_sd[j]^2).
# Edit each of the 14 standard deviations independently as needed.
true_beta_input <- c(2, .5, .5, 1, 1, 1, .5, .3, .3, .3, .2, .2, .2, 6)
config <- list(
  seed = 1L,
  out_dir = "results/seed1",
  true_beta = true_beta_input,
  J = 100L,
  n_starts = 20L,
  loops = 20L, # Maximum beta-W iterations within EACH independent Stage 1 start.
  start_sd = rep(2, 14),
  start_max_attempts = 50L, # Includes the first proposal; initial checks only.
  maxit = 8000L,
  outer_tol = 1e-5, # Checked from optimization 2; Stage 1 max beta change.
  cb_max_iter = 20L,
  cb_tol = 1e-5, # Only abs(c_new - c_old); checked from alternating round 2.
  c_interval = c(-3, 3),
  cb_beta_maxit = 5000L
)

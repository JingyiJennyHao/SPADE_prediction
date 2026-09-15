# Edit these inputs before submitting. Bounds are absolute, coordinate-wise.
# +/-2 is an example wide-start experiment, not a change to the true beta/DGP.
true_beta_input <- c(2, .5, .5, 1, 1, 1, .5, .3, .3, .3, .2, .2, .2, 6)
config <- list(
  seed = 1L,
  out_csv = "single_replication/inference_seed1.csv",
  true_beta = true_beta_input,
  J = 100L,
  n_starts = 20L,
  loops = 6L,
  start_lower = true_beta_input - 2,
  start_upper = true_beta_input + 2,
  maxit = 8000L,
  outer_tol = 1e-5,
  cb_max_iter = 20L,
  cb_tol = 1e-5,
  c_interval = c(-3, 3),
  cb_beta_maxit = 5000L
)

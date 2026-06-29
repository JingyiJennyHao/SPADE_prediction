current_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
current_dir <- if (!is.null(current_file)) dirname(normalizePath(current_file)) else getwd()
source(file.path(current_dir, "..", "composite_utils.R"))

stack_score_i <- function(j, beta, sim_dat, Time) {
  scores <- lapply(seq_len(Time), function(t) {
    f_jt <- function(b) ll_row(j, t, b, sim_dat)
    numDeriv::grad(f_jt, beta)
  })
  as.numeric(unlist(scores, use.names = FALSE))
}

score_matrix <- function(beta, sim_dat, J, Time) {
  scores <- lapply(seq_len(J), stack_score_i, beta = beta, sim_dat = sim_dat, Time = Time)
  do.call(rbind, scores)
}

vcov_scores <- function(score_mat, center = TRUE) {
  if (!is.matrix(score_mat)) {
    score_mat <- as.matrix(score_mat)
  }
  if (nrow(score_mat) < 2L) {
    stop("At least two independent groups are required to estimate V.", call. = FALSE)
  }
  if (center) {
    return(stats::cov(score_mat))
  }
  crossprod(score_mat) / nrow(score_mat)
}

invert_weight_matrix <- function(V, ridge = 1e-6, use_ginv = TRUE) {
  V_ridge <- V + diag(ridge, nrow(V))
  W <- try(solve(V_ridge), silent = TRUE)
  if (!inherits(W, "try-error")) {
    return(W)
  }
  if (!use_ginv) {
    stop("Could not invert V. Try increasing ridge or set use_ginv = TRUE.", call. = FALSE)
  }
  MASS::ginv(V_ridge)
}

gmm_weight_from_beta0 <- function(beta0, sim_dat, J, Time, ridge = 1e-6, center = TRUE) {
  S0 <- score_matrix(beta0, sim_dat, J, Time)
  V <- vcov_scores(S0, center = center)
  W <- invert_weight_matrix(V, ridge = ridge)
  list(S0 = S0, V = V, W = W, ridge = ridge, center = center)
}

gmm_objective <- function(beta, W, sim_dat, J, Time) {
  S <- score_matrix(beta, sim_dat, J, Time)
  values <- rowSums((S %*% W) * S)
  sum(values)
}

draw_gmm_start <- function(beta0, start_scale = 0.1) {
  as.numeric(beta0) + stats::rnorm(length(beta0), mean = 0, sd = start_scale)
}

fit_gmm_multistart <- function(beta0,
                               sim_dat,
                               J,
                               Time,
                               n_starts = 100L,
                               seed = NULL,
                               start_scale = 0.1,
                               ridge = 1e-6,
                               maxit = 8000L,
                               center_v = TRUE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  beta0 <- as.numeric(beta0)
  weight <- gmm_weight_from_beta0(
    beta0 = beta0,
    sim_dat = sim_dat,
    J = J,
    Time = Time,
    ridge = ridge,
    center = center_v
  )

  starts <- vector("list", n_starts)
  fits <- vector("list", n_starts)
  best_idx <- NA_integer_
  best_value <- Inf

  for (k in seq_len(n_starts)) {
    starts[[k]] <- if (k == 1L) beta0 else draw_gmm_start(beta0, start_scale = start_scale)
    initial_value <- try(
      gmm_objective(starts[[k]], weight$W, sim_dat, J, Time),
      silent = TRUE
    )

    if (inherits(initial_value, "try-error") || !is.finite(initial_value)) {
      fits[[k]] <- list(
        convergence = 99L,
        value = Inf,
        par = rep(NA_real_, length(beta0)),
        initial_value = initial_value,
        message = "Initial GMM objective was invalid."
      )
      next
    }

    fit <- try(
      optim(
        par = starts[[k]],
        fn = gmm_objective,
        W = weight$W,
        sim_dat = sim_dat,
        J = J,
        Time = Time,
        control = list(maxit = maxit)
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      fits[[k]] <- list(
        convergence = 99L,
        value = Inf,
        par = rep(NA_real_, length(beta0)),
        initial_value = initial_value,
        message = as.character(fit)
      )
      next
    }

    fits[[k]] <- list(
      convergence = as.integer(fit$convergence),
      value = fit$value,
      par = fit$par,
      initial_value = initial_value,
      message = fit$message %||% NULL
    )

    if (is.finite(fit$value) && fit$value < best_value) {
      best_idx <- k
      best_value <- fit$value
    }
  }

  if (is.na(best_idx)) {
    best <- list(convergence = 99L, value = Inf, par = rep(NA_real_, length(beta0)), message = "All starts failed.")
  } else {
    best <- fits[[best_idx]]
  }

  list(
    beta0 = beta0,
    beta_hat = best$par,
    objective_value = best$value,
    convergence = best$convergence,
    message = best$message %||% NULL,
    best_start = best_idx,
    n_starts = n_starts,
    start_scale = start_scale,
    ridge = ridge,
    center_v = center_v,
    V = weight$V,
    W = weight$W,
    S0 = weight$S0,
    starts = starts,
    fits = fits
  )
}

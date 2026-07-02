current_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
current_dir <- if (!is.null(current_file)) dirname(normalizePath(current_file)) else getwd()
# source(file.path(current_dir, "..", "hpc_long_loop", "composite_utils.R"))

required_packages <- c(
  "dplyr",
  "numDeriv",
  "LaplacesDemon",
  "MASS",
  "pracma",
  "Rcpp",
  "MGLM",
  "VGAM",
  "extraDistr"
)

invisible(lapply(required_packages, library, character.only = TRUE))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

parse_args <- function(argv) {
  args <- list()
  i <- 1L
  while (i <= length(argv)) {
    token <- argv[[i]]
    if (!startsWith(token, "--")) {
      stop(sprintf("Unexpected argument: %s", token), call. = FALSE)
    }
    key <- sub("^--", "", token)
    if (identical(key, "help")) {
      args$help <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(argv)) {
      stop(sprintf("Missing value for argument: --%s", key), call. = FALSE)
    }
    args[[key]] <- argv[[i + 1L]]
    i <- i + 2L
  }
  args
}

require_args <- function(args, required) {
  missing <- required[!required %in% names(args)]
  if (length(missing) > 0L) {
    stop(sprintf("Missing required arguments: %s", paste(sprintf("--%s", missing), collapse = ", ")), call. = FALSE)
  }
}

as_int_arg <- function(args, name) {
  value <- suppressWarnings(as.integer(args[[name]]))
  if (is.na(value)) {
    stop(sprintf("Argument --%s must be an integer.", name), call. = FALSE)
  }
  value
}

normalize_weights <- function(weights) {
  weights / sum(weights)
}

prepare_data <- function(data_path, group_count) {
  if (!file.exists(data_path)) {
    stop(sprintf("Input data file not found: %s", data_path), call. = FALSE)
  }

  dat <- readRDS(data_path)
  dat_clean <- dat %>%
    dplyr::group_by(DSHOSPID.dizzy) %>%
    dplyr::filter(dplyr::n_distinct(year) == 3) %>%
    dplyr::mutate(hos_size = mean(hos_size, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(DSHOSPID.dizzy, year) %>%
    dplyr::mutate(hosp_idx = as.integer(factor(DSHOSPID.dizzy))) %>%
    dplyr::group_by(hosp_idx) %>%
    dplyr::arrange(year, .by_group = TRUE) %>%
    dplyr::mutate(year_idx = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ct_rate_z = as.numeric(scale(ct_rate)),
      mri_rate_z = as.numeric(scale(mri_rate)),
      hos_size_z = as.numeric(scale(hos_size))
    )

  available_groups <- sort(unique(dat_clean$hosp_idx))
  if (group_count > length(available_groups)) {
    stop(sprintf("Requested %d groups, but only %d are available.", group_count, length(available_groups)), call. = FALSE)
  }

  keep_groups <- available_groups[seq_len(group_count)]
  dat_subset <- dat_clean %>%
    dplyr::filter(hosp_idx %in% keep_groups) %>%
    dplyr::mutate(hosp_idx = match(hosp_idx, keep_groups))

  list(data = dat_subset, J = group_count, Time = 3L)
}

ll_row <- function(j, t, beta, sim_dat) {
  idx <- which(sim_dat[, "hosp_idx"] == j)[t]
  njt <- as.numeric(sim_dat[idx, c("ni1", "ni2", "ni3", "ni4", "ni5", "ni6")])
  xjt <- c(1, sim_dat$ct_rate_z[idx], sim_dat$mri_rate_z[idx], sim_dat$hos_size_z[idx])

  a1 <- exp(drop(t(xjt) %*% beta[1:4]))
  a2 <- exp(drop(t(xjt) %*% beta[5:8]))
  a3 <- exp(beta[9])
  a4 <- exp(drop(t(xjt) %*% beta[10:13]))
  a5 <- exp(drop(t(xjt) %*% beta[14:17]))
  a6 <- exp(beta[18])
  alpha <- c(a1, a2, a3, a4, a5, a6)

  MGLM::ddirmn(njt, alpha)
}

phi_i_fun <- function(j, weights, beta, sim_dat, Time) {
  sum(weights * sapply(seq_len(Time), function(t) ll_row(j, t, beta, sim_dat)))
}

neg_comp_llh <- function(beta, weights, sim_dat, J, Time) {
  -sum(sapply(seq_len(J), function(j) phi_i_fun(j, weights, beta, sim_dat, Time)))
}

precompute_grad_hess <- function(beta, sim_dat, J, Time) {
  p <- length(beta)
  G <- vector("list", J)
  H <- vector("list", J)
  for (j in seq_len(J)) {
    Gj <- matrix(0, p, Time)
    Hj <- vector("list", Time)
    for (t in seq_len(Time)) {
      f_jt <- function(b) ll_row(j, t, b, sim_dat)
      Gj[, t] <- numDeriv::grad(f_jt, beta)
      Hj[[t]] <- numDeriv::hessian(f_jt, beta)
    }
    G[[j]] <- Gj
    H[[j]] <- Hj
  }
  list(G = G, H = H, p = p, Time = Time, J = J)
}

j_mat_cached <- function(weights, cache) {
  p <- cache$p
  J <- cache$J
  U <- matrix(0, p, J)
  for (j in seq_len(J)) {
    U[, j] <- cache$G[[j]] %*% weights
  }
  g_bar <- rowMeans(U)
  g_bar %*% t(g_bar)
}

h_mat_cached <- function(weights, cache) {
  p <- cache$p
  J <- cache$J
  Hsum <- matrix(0, p, p)
  for (j in seq_len(J)) {
    for (t in seq_len(cache$Time)) {
      Hsum <- Hsum + weights[t] * (-cache$H[[j]][[t]])
    }
  }
  Hsum / J
}

var_sandwich_cached <- function(weights, cache) {
  H <- h_mat_cached(weights, cache)
  Jm <- j_mat_cached(weights, cache)
  Hinv <- solve(H)
  Hinv %*% Jm %*% t(Hinv)
}

softmax3 <- function(a1, a2) {
  den <- 1 + exp(a1) + exp(a2)
  c(exp(a1) / den, exp(a2) / den, 1 / den)
}

default_fixed_beta_start <- function() {
  c(
    2, 0.5, 0.5, 0.5,
    1, 1, 1, 1,
    0.5,
    0.3, 0.3, 0.3, 0.3,
    0.2, 0.2, 0.2, 0.2,
    6
  )
}

draw_beta_start <- function(start_family, beta_reference = NULL) {
  if (identical(start_family, "random")) {
    return(c(stats::runif(17, min = -2, max = 2), 6) + stats::rnorm(18, 0, 0.1))
  }
  if (identical(start_family, "fixed")) {
    return(default_fixed_beta_start() + stats::rnorm(18, 0, 0.1))
  }
  if (identical(start_family, "beta_ref")) {
    if (is.null(beta_reference)) {
      stop("beta_reference is required for start_family='beta_ref'.", call. = FALSE)
    }
    return(as.numeric(beta_reference) + stats::rnorm(length(beta_reference), 0, 0.1))
  }
  stop(sprintf("Unknown start family: %s", start_family), call. = FALSE)
}

load_beta_reference <- function(path) {
  ref <- readRDS(path)
  if (!is.null(ref$beta_hat)) {
    return(as.numeric(ref$beta_hat))
  }
  if (is.numeric(ref)) {
    return(as.numeric(ref))
  }
  stop(sprintf("Could not extract beta_hat from %s", path), call. = FALSE)
}

load_weight_reference <- function(path) {
  ref <- readRDS(path)
  weights <- ref$weights_hat %||% ref$weights
  if (is.null(weights)) {
    stop(sprintf("Could not extract weights from %s", path), call. = FALSE)
  }
  normalize_weights(as.numeric(weights))
}

objective_var_fun_factory <- function(cache) {
  function(alpha) {
    weights <- softmax3(alpha[1], alpha[2])
    variance <- var_sandwich_cached(weights, cache)
    sum(diag(variance))
  }
}

alpha_from_weights <- function(weights) {
  w <- normalize_weights(weights)
  c(log(w[1] / w[3]), log(w[2] / w[3]))
}

fail_fast_dir <- function(path) {
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

save_result <- function(path, payload) {
  fail_fast_dir(path)
  saveRDS(payload, path)
}


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

fit_gmm_single_start <- function(start,
                                 beta0,
                                 sim_dat,
                                 J,
                                 Time,
                                 ridge = 1e-6,
                                 maxit = 8000L,
                                 center_v = TRUE) {
  beta0 <- as.numeric(beta0)
  start <- as.numeric(start)
  weight <- gmm_weight_from_beta0(
    beta0 = beta0,
    sim_dat = sim_dat,
    J = J,
    Time = Time,
    ridge = ridge,
    center = center_v
  )

  initial_value <- try(
    gmm_objective(start, weight$W, sim_dat, J, Time),
    silent = TRUE
  )

  if (inherits(initial_value, "try-error") || !is.finite(initial_value)) {
    return(list(
      convergence = 99L,
      value = Inf,
      par = rep(NA_real_, length(start)),
      initial_value = initial_value,
      message = "Initial GMM objective was invalid.",
      V = weight$V,
      W = weight$W,
      S0 = weight$S0
    ))
  }

  fit <- try(
    optim(
      par = start,
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
    return(list(
      convergence = 99L,
      value = Inf,
      par = rep(NA_real_, length(start)),
      initial_value = initial_value,
      message = as.character(fit),
      V = weight$V,
      W = weight$W,
      S0 = weight$S0
    ))
  }

  list(
    convergence = as.integer(fit$convergence),
    value = fit$value,
    par = fit$par,
    initial_value = initial_value,
    message = fit$message %||% NULL,
    V = weight$V,
    W = weight$W,
    S0 = weight$S0
  )
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

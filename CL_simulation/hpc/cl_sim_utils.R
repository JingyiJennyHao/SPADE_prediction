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

true_beta <- c(2, 0.5, 0.5, 1, 1, 1, 0.5, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 6)

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

as_num_arg <- function(args, name) {
  value <- suppressWarnings(as.numeric(args[[name]]))
  if (is.na(value)) {
    stop(sprintf("Argument --%s must be numeric.", name), call. = FALSE)
  }
  value
}

softmax3 <- function(a1, a2) {
  den <- 1 + exp(a1) + exp(a2)
  c(exp(a1) / den, exp(a2) / den, 1 / den)
}

normalize_weights <- function(weights) {
  weights / sum(weights)
}

alpha_from_weights <- function(weights) {
  w <- normalize_weights(weights)
  c(log(w[1] / w[3]), log(w[2] / w[3]))
}

weights_from_alpha <- function(alpha) {
  softmax3(alpha[1], alpha[2])
}

ensure_parent_dir <- function(path) {
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

with_lock <- function(lock_path, expr, timeout_seconds = 600, poll_seconds = 1) {
  ensure_parent_dir(lock_path)
  start <- Sys.time()
  repeat {
    if (dir.create(lock_path, showWarnings = FALSE)) {
      on.exit(unlink(lock_path, recursive = TRUE), add = TRUE)
      return(force(expr))
    }
    if (as.numeric(difftime(Sys.time(), start, units = "secs")) > timeout_seconds) {
      stop(sprintf("Timed out waiting for lock: %s", lock_path), call. = FALSE)
    }
    Sys.sleep(poll_seconds)
  }
}

append_csv_locked <- function(row, out_csv) {
  stopifnot(is.data.frame(row), nrow(row) == 1L)
  ensure_parent_dir(out_csv)
  lock_path <- paste0(out_csv, ".lock")
  with_lock(lock_path, {
    write_header <- !file.exists(out_csv) || file.info(out_csv)$size == 0
    utils::write.table(
      row,
      file = out_csv,
      sep = ",",
      row.names = FALSE,
      col.names = write_header,
      append = !write_header,
      quote = TRUE
    )
  })
}

save_result <- function(path, payload) {
  ensure_parent_dir(path)
  saveRDS(payload, path)
}

read_best_beta <- function(path) {
  item <- readRDS(path)
  beta <- item$beta_hat %||% item$par
  if (is.null(beta)) {
    stop(sprintf("Could not find beta_hat in %s", path), call. = FALSE)
  }
  as.numeric(beta)
}

read_best_weights <- function(path) {
  item <- readRDS(path)
  weights <- item$weights_hat %||% item$weights
  if (is.null(weights)) {
    stop(sprintf("Could not find weights_hat in %s", path), call. = FALSE)
  }
  normalize_weights(as.numeric(weights))
}

read_best_alpha <- function(path) {
  item <- readRDS(path)
  alpha <- item$alpha_hat %||% item$alpha
  if (!is.null(alpha)) {
    return(as.numeric(alpha))
  }
  alpha_from_weights(read_best_weights(path))
}

gen_simdata <- function(J, seed, beta = true_beta) {
  set.seed(seed)
  Time <- 3
  mu <- c(1, 2)
  A <- diag(c(0.7, 0.5))
  Sigma_eps <- matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)

  hos_Xi <- vector("list", J)
  for (j in seq_len(J)) {
    X <- matrix(NA_real_, nrow = Time, ncol = 2)
    X[1, ] <- mu
    for (t in 2:Time) {
      eps <- MASS::mvrnorm(1, mu = c(0, 0), Sigma = Sigma_eps)
      X[t, ] <- mu + A %*% (X[t - 1, ] - mu) + eps
    }
    hos_Xi[[j]] <- rbind(x0 = 1, x1 = X[, 1], x2 = X[, 2])
  }

  beta1 <- beta[1:3]
  beta2 <- beta[4:6]
  beta3 <- beta[7]
  beta4 <- beta[8:10]
  beta5 <- beta[11:13]
  beta6 <- beta[14]

  sim_dat <- matrix(NA_real_, nrow = J * Time, ncol = 24)
  colnames(sim_dat) <- c(
    "j", "t", "Nj", paste0("pj", 1:6), paste0("nj", 1:6),
    paste0("a", 1:6), "x0", "x1", "x2"
  )

  row_idx <- 1L
  for (j in seq_len(J)) {
    Nj1 <- sample(100:2000, 1)
    Nj2 <- sample(100:2000, 1)
    Nj3 <- sample(100:2000, 1)
    Xj1 <- hos_Xi[[j]][, 1]
    Xj2 <- hos_Xi[[j]][, 2]
    Xj3 <- hos_Xi[[j]][, 3]

    alpha1 <- c(
      exp(drop(beta1 %*% Xj1)),
      exp(drop(beta2 %*% Xj1)),
      exp(beta3),
      exp(drop(beta4 %*% Xj1)),
      exp(drop(beta5 %*% Xj1)),
      exp(beta6)
    )
    alpha2 <- c(
      exp(drop(beta1 %*% Xj2)),
      exp(drop(beta2 %*% Xj2)),
      exp(beta3),
      exp(drop(beta4 %*% Xj2)),
      exp(drop(beta5 %*% Xj2)),
      exp(beta6)
    )
    alpha3 <- c(
      exp(drop(beta1 %*% Xj3)),
      exp(drop(beta2 %*% Xj3)),
      exp(beta3),
      exp(drop(beta4 %*% Xj3)),
      exp(drop(beta5 %*% Xj3)),
      exp(beta6)
    )

    p_j <- as.numeric(LaplacesDemon::rdirichlet(1, c(alpha1, alpha2, alpha3)))
    p_j1 <- p_j[1:6]
    p_j2 <- p_j[7:12]
    p_j3 <- p_j[13:18]
    p_j1_norm <- p_j1 / sum(p_j1)
    p_j2_norm <- p_j2 / sum(p_j2)
    p_j3_norm <- p_j3 / sum(p_j3)

    n_j1 <- as.numeric(stats::rmultinom(1, size = Nj1, prob = p_j1))
    n_j2 <- as.numeric(stats::rmultinom(1, size = Nj2, prob = p_j2))
    n_j3 <- as.numeric(stats::rmultinom(1, size = Nj3, prob = p_j3))

    sim_dat[row_idx, ] <- c(j, 1, Nj1, p_j1_norm, n_j1, alpha1, Xj1)
    row_idx <- row_idx + 1L
    sim_dat[row_idx, ] <- c(j, 2, Nj2, p_j2_norm, n_j2, alpha2, Xj2)
    row_idx <- row_idx + 1L
    sim_dat[row_idx, ] <- c(j, 3, Nj3, p_j3_norm, n_j3, alpha3, Xj3)
    row_idx <- row_idx + 1L
  }

  sim_dat
}

ll_row <- function(j, t, beta, sim_dat) {
  idx <- which(sim_dat[, "j"] == j)[t]
  njt <- as.numeric(sim_dat[idx, paste0("nj", 1:6)])
  xjt <- c(1, sim_dat[idx, "x1"], sim_dat[idx, "x2"])

  alpha <- c(
    exp(drop(t(xjt) %*% beta[1:3])),
    exp(drop(t(xjt) %*% beta[4:6])),
    exp(beta[7]),
    exp(drop(t(xjt) %*% beta[8:10])),
    exp(drop(t(xjt) %*% beta[11:13])),
    exp(beta[14])
  )

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

# J_mat_cached <- function(weights, cache) {
#   p <- cache$p
#   J <- cache$J
#   U <- matrix(0, p, J)
#   for (j in 1:J) U[, j] <- cache$G[[j]] %*% weights
#   g_bar <- rowMeans(U)
#   # 210 is the number of hospitals
#   res <- g_bar %*% t(g_bar)*210
#   return(res)
# }

# revised J
j_mat_cached <- function(weights, cache) {
  p <- cache$p
  J <- cache$J
  U <- matrix(0, p, J)
  for (j in 1:J) U[, j] <- cache$G[[j]] %*% weights
  return(U %*% t(U) / J)
}

h_mat_cached <- function(weights, cache) {
  p <- cache$p
  Hsum <- matrix(0, p, p)
  for (j in seq_len(cache$J)) {
    for (t in seq_len(cache$Time)) {
      Hsum <- Hsum + weights[t] * (-cache$H[[j]][[t]])
    }
  }
  Hsum / cache$J
}

var_sandwich_cached <- function(weights, cache) {
  H <- h_mat_cached(weights, cache)
  Jm <- j_mat_cached(weights, cache)
  Hinv <- solve(H)
  Hinv %*% Jm %*% t(Hinv)
}

objective_var_fun_factory <- function(cache) {
  function(alpha) {
    weights <- weights_from_alpha(alpha)
    V <- var_sandwich_cached(weights, cache)
    sum(diag(V))
  }
}

draw_beta_start <- function(seed, beta_ref = true_beta, sd = 0.1) {
  set.seed(seed)
  as.numeric(beta_ref) + stats::rnorm(length(beta_ref), mean = 0, sd = sd)
}

draw_alpha_start <- function(seed, alpha_ref = c(0, 0), sd = 0.1) {
  set.seed(seed)
  as.numeric(alpha_ref) + stats::rnorm(length(alpha_ref), mean = 0, sd = sd)
}

make_sim_data_for_round <- function(seed, J) {
  list(data = gen_simdata(J = J, seed = seed), J = J, Time = 3L)
}

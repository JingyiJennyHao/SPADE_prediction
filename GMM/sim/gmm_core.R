suppressPackageStartupMessages({
  library(dplyr)
  library(numDeriv)
  library(LaplacesDemon)
  library(MASS)
  library(pracma)
  library(Rcpp)
  library(MGLM)
  library(VGAM)
  library(extraDistr)
})

true_beta <- c(2, 0.5, 0.5, 1, 1, 1, 0.5, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 6)

write_row <- function(df_row, file, lock_timeout = 1800, stale_after = 7200) {
  stopifnot(is.data.frame(df_row), nrow(df_row) == 1)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  lock_dir <- paste0(file, ".lock")
  started <- Sys.time()
  repeat {
    if (dir.create(lock_dir, showWarnings = FALSE)) break

    lock_info <- file.info(lock_dir)
    if (!is.na(lock_info$mtime) &&
        as.numeric(difftime(Sys.time(), lock_info$mtime, units = "secs")) > stale_after) {
      unlink(lock_dir, recursive = TRUE, force = TRUE)
      next
    }

    if (as.numeric(difftime(Sys.time(), started, units = "secs")) > lock_timeout) {
      stop("Timed out waiting for CSV lock: ", lock_dir)
    }
    Sys.sleep(runif(1, 0.25, 1.25))
  }
  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)

  if (!file.exists(file)) {
    write.table(df_row, file, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    write.table(df_row, file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

softmax3 <- function(a1, a2) {
  den <- 1 + exp(a1) + exp(a2)
  c(exp(a1) / den, exp(a2) / den, 1 / den)
}

gen_simdata <- function(J, true_beta = true_beta, seed = seed) {
  set.seed(seed)
  Time <- 3
  mu <- c(1, 2)
  A <- diag(c(0.7, 0.5))
  Sigma_eps <- matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)

  hos.Xi <- vector("list", J)
  for (j in 1:J) {
    X <- matrix(NA, nrow = Time, ncol = 2)
    X[1, ] <- mu
    for (t in 2:Time) {
      eps <- mvrnorm(1, mu = c(0, 0), Sigma = Sigma_eps)
      X[t, ] <- mu + A %*% (X[t - 1, ] - mu) + eps
    }
    hos.Xi[[j]] <- rbind(x0 = 1, x1 = X[, 1], x2 = X[, 2])
  }

  true_beta1 <- true_beta[1:3]
  true_beta2 <- true_beta[4:6]
  true_beta3 <- true_beta[7]
  true_beta4 <- true_beta[8:10]
  true_beta5 <- true_beta[11:13]
  true_beta6 <- true_beta[14]

  sim_dat1 <- matrix(NA, nrow = J * Time, ncol = 24)
  colnames(sim_dat1) <- c(
    "j", "t", "Nj", "pj1", "pj2", "pj3", "pj4", "pj5", "pj6",
    "nj1", "nj2", "nj3", "nj4", "nj5", "nj6",
    "a1", "a2", "a3", "a4", "a5", "a6", "x0", "x1", "x2"
  )
  row_idx <- 1

  for (j in 1:J) {
    Nj1 <- sample(100:2000, 1)
    Nj2 <- sample(100:2000, 1)
    Nj3 <- sample(100:2000, 1)
    Xj1 <- hos.Xi[[j]][, 1]
    Xj2 <- hos.Xi[[j]][, 2]
    Xj3 <- hos.Xi[[j]][, 3]

    alpha1 <- c(
      exp(true_beta1 %*% Xj1), exp(true_beta2 %*% Xj1), exp(true_beta3),
      exp(true_beta4 %*% Xj1), exp(true_beta5 %*% Xj1), exp(true_beta6)
    )
    alpha2 <- c(
      exp(true_beta1 %*% Xj2), exp(true_beta2 %*% Xj2), exp(true_beta3),
      exp(true_beta4 %*% Xj2), exp(true_beta5 %*% Xj2), exp(true_beta6)
    )
    alpha3 <- c(
      exp(true_beta1 %*% Xj3), exp(true_beta2 %*% Xj3), exp(true_beta3),
      exp(true_beta4 %*% Xj3), exp(true_beta5 %*% Xj3), exp(true_beta6)
    )

    p_j <- LaplacesDemon::rdirichlet(1, c(alpha1, alpha2, alpha3))
    p_j1 <- p_j[1:6]
    p_j2 <- p_j[7:12]
    p_j3 <- p_j[13:18]

    p_j1_norm <- p_j1 / sum(p_j1)
    p_j2_norm <- p_j2 / sum(p_j2)
    p_j3_norm <- p_j3 / sum(p_j3)

    n_j1 <- rmultinom(1, size = Nj1, prob = p_j1)
    n_j2 <- rmultinom(1, size = Nj2, prob = p_j2)
    n_j3 <- rmultinom(1, size = Nj3, prob = p_j3)

    sim_dat1[row_idx, ] <- c(j, 1, Nj1, p_j1_norm, n_j1, alpha1, Xj1)
    row_idx <- row_idx + 1
    sim_dat1[row_idx, ] <- c(j, 2, Nj2, p_j2_norm, n_j2, alpha2, Xj2)
    row_idx <- row_idx + 1
    sim_dat1[row_idx, ] <- c(j, 3, Nj3, p_j3_norm, n_j3, alpha3, Xj3)
    row_idx <- row_idx + 1
  }

  set.seed(seed + 5000)

  sim_dat2 <- sim_dat1
  for (r in 1:nrow(sim_dat2)) {
    Nj_new <- sample(100:2000, 1)
    p_norm <- as.numeric(sim_dat1[r, paste0("pj", 1:6)])
    n_new <- as.numeric(rmultinom(1, size = Nj_new, prob = p_norm))
    sim_dat2[r, "Nj"] <- Nj_new
    sim_dat2[r, paste0("nj", 1:6)] <- n_new
  }

  list(sim_dat1 = sim_dat1, sim_dat2 = sim_dat2)
}

ll_row <- function(j, t, beta, sim_dat) {
  idx <- which(sim_dat[, "j"] == j)[t]
  njt <- as.numeric(sim_dat[idx, c("nj1", "nj2", "nj3", "nj4", "nj5", "nj6")])
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

score_ll_row <- function(j, t, beta, sim_dat) {
  f_jt <- function(b) ll_row(j, t, b, sim_dat)
  numDeriv::grad(f_jt, beta)
}

score_stack <- function(j, beta, sim_dat, Time = 3) {
  unlist(lapply(seq_len(Time), function(t) score_ll_row(j, t, beta, sim_dat)))
}

score_average <- function(beta, sim_dat, Time = 3) {
  j_list <- sort(unique(sim_dat[, "j"]))
  G <- sapply(j_list, function(j) score_stack(j, beta, sim_dat, Time = Time))
  rowMeans(G)
}

W_hat_fun <- function(beta, sim_dat, Time = 3, ridge = 1e-6) {
  j_list <- sort(unique(sim_dat[, "j"]))
  G <- sapply(j_list, function(j) score_stack(j, beta, sim_dat, Time = Time))
  G <- t(G)
  S_hat <- cov(G)
  solve(S_hat + ridge * diag(ncol(S_hat)))
}

G_hat_fun <- function(beta, sim_dat, Time = 3) {
  j_list <- sort(unique(sim_dat[, "j"]))
  p <- length(beta)
  G_sum <- matrix(0, nrow = Time * p, ncol = p)
  for (j in j_list) {
    Gj <- numDeriv::jacobian(function(b) score_stack(j, b, sim_dat, Time), beta)
    G_sum <- G_sum + Gj
  }
  G_sum / length(j_list)
}

gmm_obj <- function(beta, sim_dat, W_hat, Time = 3) {
  gbar <- score_average(beta, sim_dat, Time = Time)
  as.numeric(t(gbar) %*% W_hat %*% gbar)
}

stats_samples1 <- function(dir_param) {
  x <- rdirichlet(1000, dir_param)
  x[, 1] + x[, 2] - 1 / 9 * (x[, 4] + x[, 5])
}

stats_samples2 <- function(dir_param) {
  x <- rdirichlet(1000, dir_param)
  x[, 1] - 1 / 9 * x[, 4]
}

stats_samples3 <- function(dir_param, size) {
  x <- rdirmnom(1000, size = size, dir_param)
  x[, 1] + x[, 2] - 1 / 9 * (x[, 4] + x[, 5])
}

stats_samples4 <- function(dir_param) {
  x <- rdirichlet(1000, dir_param)
  x[, 2] - 1 / 9 * x[, 5]
}

hpd <- function(sample_list, prob = 0.95) {
  x <- sort(sample_list[is.finite(sample_list)])
  M <- length(x)
  k <- ceiling(prob * M)
  if (k >= M) return(c(min(x), max(x)))
  w <- x[k:M] - x[1:(M - k + 1)]
  i <- which.min(w)
  c(x[i], x[i + k - 1])
}

PI_betap <- function(alpha, shape1, shape2, scale) {
  ratio_samples <- sort(rbetapr(5000, shape1, shape2, scale))
  pr_samples <- dbetapr(ratio_samples, shape1, shape2, scale)
  cdf_samples <- pbetapr(ratio_samples, shape1, shape2, scale)
  mode <- which.max(pr_samples)
  if (mode != 1) {
    intervals <- numeric(mode)
    index1s <- numeric(mode)
    index2s <- numeric(mode)
    for (index in 1:mode) {
      p2_index_range <- (mode + 1):length(pr_samples)
      differences <- abs(pr_samples[p2_index_range] - pr_samples[index])
      index2 <- p2_index_range[which.min(differences)]
      intervals[index] <- cdf_samples[index2] - cdf_samples[index]
      index1s[index] <- index
      index2s[index] <- index2
    }
    closest_index <- which.min(abs(intervals - alpha))
    return(list(l = ratio_samples[index1s[closest_index]], u = ratio_samples[index2s[closest_index]]))
  }
  list(l = 0, u = qbetapr(0.95, shape1, shape2, scale))
}

run_start_point <- function(seed, loop, start_id, base_beta, out_csv,
                            J = 100, Time = 3, start_sd = 0.1, maxit = 5000) {
  datasets <- gen_simdata(J, true_beta = true_beta, seed = seed)
  sim_dat <- datasets$sim_dat1
  W_current <- if (loop <= 1) diag(Time * length(base_beta)) else W_hat_fun(base_beta, sim_dat, Time = Time)

  start_seed <- seed + loop * 100000L + start_id
  set.seed(start_seed)
  stpt <- base_beta + rnorm(length(base_beta), 0, start_sd)

  fit <- tryCatch(
    optim(stpt, gmm_obj, sim_dat = sim_dat, W_hat = W_current, Time = Time, control = list(maxit = maxit)),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    out <- data.frame(
      seed = seed, loop = loop, start_id = start_id, start_seed = start_seed,
      convergence = NA_integer_, value = NA_real_, diff_from_base = NA_real_,
      error = conditionMessage(fit)
    )
    out[paste0("b", seq_along(base_beta))] <- as.list(rep(NA_real_, length(base_beta)))
  } else {
    out <- data.frame(
      seed = seed, loop = loop, start_id = start_id, start_seed = start_seed,
      convergence = fit$convergence, value = fit$value,
      diff_from_base = max(abs(fit$par - base_beta)), error = ""
    )
    out[paste0("b", seq_along(fit$par))] <- as.list(fit$par)
  }

  write_row(out, out_csv)
  invisible(out)
}

best_beta_from_loop_file <- function(loop_csv, require_converged = TRUE) {
  if (!file.exists(loop_csv)) stop("Missing loop CSV: ", loop_csv)
  x <- read.csv(loop_csv, stringsAsFactors = FALSE)
  if (require_converged) {
    x <- x[x$convergence == 0 & is.finite(x$value), , drop = FALSE]
  } else {
    x <- x[is.finite(x$value), , drop = FALSE]
  }
  if (!nrow(x)) stop("No usable rows in ", loop_csv)
  as.numeric(unlist(x[which.min(x$value), paste0("b", 1:14)], use.names = FALSE))
}

complete_rows <- function(loop_csv, require_converged = TRUE) {
  if (!file.exists(loop_csv)) return(0L)
  x <- tryCatch(read.csv(loop_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(0L)
  if (require_converged) return(sum(x$convergence == 0 & is.finite(x$value), na.rm = TRUE))
  nrow(x)
}

run_gmm_inference <- function(beta_hat, seed, out_csv, J = 200, Time = 3,
                              interval_prob = 0.95) {
  if (length(beta_hat) != length(true_beta) || any(!is.finite(beta_hat))) {
    stop("Inference requires 14 finite beta estimates.")
  }
  if (Time != 3) {
    stop("The existing simulation and prediction-inference code requires Time = 3.")
  }

  datasets <- gen_simdata(J, true_beta = true_beta, seed = seed)
  sim_dat <- datasets$sim_dat1
  sim_dat2 <- datasets$sim_dat2

  W_hat <- W_hat_fun(beta_hat, sim_dat, Time = Time)
  G_hat <- G_hat_fun(beta_hat, sim_dat, Time = Time)
  information <- t(G_hat) %*% W_hat %*% G_hat
  Vhat <- solve(information) / J
  se <- sqrt(diag(Vhat))
  critical_value <- qnorm(1 - (1 - interval_prob) / 2)
  lower <- beta_hat - critical_value * se
  upper <- beta_hat + critical_value * se
  cover <- as.integer(true_beta >= lower & true_beta <= upper)

  prediction_names <- c(
    paste0("PI1", 1:6, "_coverage"),
    paste0("PI2", 1:6, "_coverage"),
    paste0("PI3", 1:5, "_coverage")
  )
  prediction_hits <- setNames(numeric(length(prediction_names)), prediction_names)
  in_interval <- function(value, interval) {
    as.integer(value >= interval[1] && value <= interval[2])
  }

  # Make the Monte Carlo prediction intervals reproducible without changing
  # the simulated estimation and prediction datasets.
  set.seed(seed + 9000000L)

  for (row in seq_len(nrow(sim_dat2))) {
    xi <- as.numeric(sim_dat2[row, c("x0", "x1", "x2")])
    Ni <- as.numeric(sim_dat[row, "Nj"])
    Ni_star <- as.numeric(sim_dat2[row, "Nj"])
    pi <- as.numeric(sim_dat2[row, paste0("pj", 1:6)])
    ni <- as.numeric(sim_dat[row, paste0("nj", 1:6)])
    ni_star <- as.numeric(sim_dat2[row, paste0("nj", 1:6)])

    alpha_hat1 <- c(
      exp(drop(t(xi) %*% beta_hat[1:3])),
      exp(drop(t(xi) %*% beta_hat[4:6])),
      exp(beta_hat[7]),
      exp(drop(t(xi) %*% beta_hat[8:10])),
      exp(drop(t(xi) %*% beta_hat[11:13])),
      exp(beta_hat[14])
    )
    alpha_hat2 <- alpha_hat1 + ni

    truth_1 <- pi[1] - pi[4] / 9
    truth_2 <- pi[2] - pi[5] / 9
    truth_3 <- pi[2] / pi[1]
    truth_4 <- pi[5] / pi[4]
    truth_5 <- pi[1] + pi[2] - (pi[4] + pi[5]) / 9

    prediction_hits["PI11_coverage"] <- prediction_hits["PI11_coverage"] +
      in_interval(truth_1, hpd(stats_samples2(alpha_hat1), prob = interval_prob))
    prediction_hits["PI12_coverage"] <- prediction_hits["PI12_coverage"] +
      in_interval(truth_2, hpd(stats_samples4(alpha_hat1), prob = interval_prob))
    interval <- PI_betap(interval_prob, alpha_hat1[2], alpha_hat1[1], 1)
    prediction_hits["PI13_coverage"] <- prediction_hits["PI13_coverage"] +
      in_interval(truth_3, c(interval$l, interval$u))
    interval <- PI_betap(interval_prob, alpha_hat1[5], alpha_hat1[4], 1)
    prediction_hits["PI14_coverage"] <- prediction_hits["PI14_coverage"] +
      in_interval(truth_4, c(interval$l, interval$u))
    prediction_hits["PI15_coverage"] <- prediction_hits["PI15_coverage"] +
      in_interval(truth_5, hpd(stats_samples1(alpha_hat1), prob = interval_prob))
    count_truth_1 <- ni[1] + ni[2] - (ni[4] + ni[5]) / 9
    prediction_hits["PI16_coverage"] <- prediction_hits["PI16_coverage"] +
      in_interval(count_truth_1, hpd(stats_samples3(alpha_hat1, Ni_star), prob = interval_prob))

    prediction_hits["PI21_coverage"] <- prediction_hits["PI21_coverage"] +
      in_interval(truth_1, hpd(stats_samples2(alpha_hat2), prob = interval_prob))
    prediction_hits["PI22_coverage"] <- prediction_hits["PI22_coverage"] +
      in_interval(truth_2, hpd(stats_samples4(alpha_hat2), prob = interval_prob))
    interval <- PI_betap(interval_prob, alpha_hat2[2], alpha_hat2[1], 1)
    prediction_hits["PI23_coverage"] <- prediction_hits["PI23_coverage"] +
      in_interval(truth_3, c(interval$l, interval$u))
    interval <- PI_betap(interval_prob, alpha_hat2[5], alpha_hat2[4], 1)
    prediction_hits["PI24_coverage"] <- prediction_hits["PI24_coverage"] +
      in_interval(truth_4, c(interval$l, interval$u))
    prediction_hits["PI25_coverage"] <- prediction_hits["PI25_coverage"] +
      in_interval(truth_5, hpd(stats_samples1(alpha_hat2), prob = interval_prob))
    count_truth_2 <- ni_star[1] + ni_star[2] - (ni_star[4] + ni_star[5]) / 9
    prediction_hits["PI26_coverage"] <- prediction_hits["PI26_coverage"] +
      in_interval(count_truth_2, hpd(stats_samples3(alpha_hat2, Ni_star), prob = interval_prob))

    pi14_sample <- rdirichlet(
      1000,
      c(
        alpha_hat1[1] + ni[1],
        alpha_hat1[4] + ni[4],
        sum(alpha_hat1[c(2, 3, 5, 6)]) + (Ni - ni[1] - ni[4])
      )
    )
    remaining_sample <- rdirichlet(1000, alpha_hat1[c(2, 3, 5, 6)])
    pi1_sample <- pi14_sample[, 1]
    pi4_sample <- pi14_sample[, 2]
    residual_sample <- pi14_sample[, 3]
    pi2_sample <- remaining_sample[, 1] * residual_sample
    pi5_sample <- remaining_sample[, 3] * residual_sample

    prediction_hits["PI31_coverage"] <- prediction_hits["PI31_coverage"] +
      in_interval(truth_1, hpd(pi1_sample - pi4_sample / 9, prob = interval_prob))
    prediction_hits["PI32_coverage"] <- prediction_hits["PI32_coverage"] +
      in_interval(truth_2, hpd(pi2_sample - pi5_sample / 9, prob = interval_prob))
    prediction_hits["PI33_coverage"] <- prediction_hits["PI33_coverage"] +
      in_interval(
        truth_5,
        hpd(pi1_sample + pi2_sample - (pi4_sample + pi5_sample) / 9,
            prob = interval_prob)
      )
    prediction_hits["PI34_coverage"] <- prediction_hits["PI34_coverage"] +
      in_interval(truth_3, hpd(pi2_sample / pi1_sample, prob = interval_prob))
    prediction_hits["PI35_coverage"] <- prediction_hits["PI35_coverage"] +
      in_interval(truth_4, hpd(pi5_sample / pi4_sample, prob = interval_prob))
  }

  prediction_coverage <- prediction_hits / nrow(sim_dat2)
  out <- data.frame(
    seed = seed,
    J = J,
    Time = Time,
    interval_prob = interval_prob,
    stringsAsFactors = FALSE
  )
  out[names(prediction_coverage)] <- as.list(prediction_coverage)
  out[paste0("b", seq_along(beta_hat))] <- as.list(beta_hat)
  out[paste0("v", seq_along(beta_hat))] <- as.list(diag(Vhat))
  out[paste0("se", seq_along(beta_hat))] <- as.list(se)
  out[paste0("lower", seq_along(beta_hat))] <- as.list(lower)
  out[paste0("upper", seq_along(beta_hat))] <- as.list(upper)
  out[paste0("c", seq_along(beta_hat))] <- as.list(cover)

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(out, out_csv, row.names = FALSE)
  invisible(out)
}

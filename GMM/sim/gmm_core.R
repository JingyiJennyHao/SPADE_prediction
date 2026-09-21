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
  
  Sigma_eps <- matrix(
    c(
      0.2, 0.1,
      0.1, 0.2
    ),
    2, 2
  )
  
  hos.Xi <- vector("list", J)
  
  for (j in 1:J) {
    
    X <- matrix(
      NA,
      nrow = Time,
      ncol = 2
    )
    
    ## ---------------------------------------
    ## OLD:
    ## X[1, ] <- mu
    ##
    ## NEW:
    ## allow between-subject variation at t=1
    ## ---------------------------------------
    
    X[1, ] <- MASS::mvrnorm(
      n = 1,
      mu = mu,
      Sigma = Sigma_eps
    )
    
    for (t in 2:Time) {
      
      eps <- MASS::mvrnorm(
        n = 1,
        mu = c(0, 0),
        Sigma = Sigma_eps
      )
      
      X[t, ] <- mu +
        A %*% (X[t - 1, ] - mu) +
        eps
    }
    
    hos.Xi[[j]] <- rbind(
      x0 = 1,
      x1 = X[, 1],
      x2 = X[, 2]
    )
  }
  
  ## ---------------------------------------
  ## everything below stays exactly the same
  ## as your original gen_simdata()
  ## ---------------------------------------
  
  true_beta1 <- true_beta[1:3]
  true_beta2 <- true_beta[4:6]
  true_beta3 <- true_beta[7]
  true_beta4 <- true_beta[8:10]
  true_beta5 <- true_beta[11:13]
  true_beta6 <- true_beta[14]
  
  sim_dat1 <- matrix(
    NA,
    nrow = J * Time,
    ncol = 24
  )
  
  colnames(sim_dat1) <- c(
    "j", "t", "Nj",
    "pj1", "pj2", "pj3", "pj4", "pj5", "pj6",
    "nj1", "nj2", "nj3", "nj4", "nj5", "nj6",
    "a1", "a2", "a3", "a4", "a5", "a6",
    "x0", "x1", "x2"
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
      exp(true_beta1 %*% Xj1),
      exp(true_beta2 %*% Xj1),
      exp(true_beta3),
      exp(true_beta4 %*% Xj1),
      exp(true_beta5 %*% Xj1),
      exp(true_beta6)
    )
    
    alpha2 <- c(
      exp(true_beta1 %*% Xj2),
      exp(true_beta2 %*% Xj2),
      exp(true_beta3),
      exp(true_beta4 %*% Xj2),
      exp(true_beta5 %*% Xj2),
      exp(true_beta6)
    )
    
    alpha3 <- c(
      exp(true_beta1 %*% Xj3),
      exp(true_beta2 %*% Xj3),
      exp(true_beta3),
      exp(true_beta4 %*% Xj3),
      exp(true_beta5 %*% Xj3),
      exp(true_beta6)
    )
    
    p_j <- LaplacesDemon::rdirichlet(
      1,
      c(alpha1, alpha2, alpha3)
    )
    
    p_j1 <- p_j[1:6]
    p_j2 <- p_j[7:12]
    p_j3 <- p_j[13:18]
    
    p_j1_norm <- p_j1 / sum(p_j1)
    p_j2_norm <- p_j2 / sum(p_j2)
    p_j3_norm <- p_j3 / sum(p_j3)
    
    n_j1 <- rmultinom(
      1,
      size = Nj1,
      prob = p_j1
    )
    
    n_j2 <- rmultinom(
      1,
      size = Nj2,
      prob = p_j2
    )
    
    n_j3 <- rmultinom(
      1,
      size = Nj3,
      prob = p_j3
    )
    
    sim_dat1[row_idx, ] <- c(
      j, 1, Nj1,
      p_j1_norm,
      n_j1,
      alpha1,
      Xj1
    )
    
    row_idx <- row_idx + 1
    
    sim_dat1[row_idx, ] <- c(
      j, 2, Nj2,
      p_j2_norm,
      n_j2,
      alpha2,
      Xj2
    )
    
    row_idx <- row_idx + 1
    
    sim_dat1[row_idx, ] <- c(
      j, 3, Nj3,
      p_j3_norm,
      n_j3,
      alpha3,
      Xj3
    )
    
    row_idx <- row_idx + 1
  }
  
  set.seed(seed + 5000)
  
  sim_dat2 <- sim_dat1
  
  for (r in 1:nrow(sim_dat2)) {
    
    Nj_new <- sample(100:2000, 1)
    
    p_norm <- as.numeric(
      sim_dat1[r, paste0("pj", 1:6)]
    )
    
    n_new <- as.numeric(
      rmultinom(
        1,
        size = Nj_new,
        prob = p_norm
      )
    )
    
    sim_dat2[r, "Nj"] <- Nj_new
    sim_dat2[r, paste0("nj", 1:6)] <- n_new
  }
  
  return(
    list(
      sim_dat1 = sim_dat1,
      sim_dat2 = sim_dat2
    )
  )
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

# Original numerical-gradient implementation, retained for comparison.
score_ll_row_numeric <- function(j, t, beta, sim_dat) {
  f_jt <- function(b) ll_row(j, t, b, sim_dat)
  numDeriv::grad(f_jt, beta)
}

# Closed-form score for the Dirichlet-multinomial log-likelihood.
score_ll_row_closed <- function(j, t, beta, sim_dat) {
  idx <- which(sim_dat[, "j"] == j)[t]

  y <- as.numeric(
    sim_dat[idx, c("nj1", "nj2", "nj3", "nj4", "nj5", "nj6")]
  )
  x <- c(1, sim_dat[idx, "x1"], sim_dat[idx, "x2"])

  alpha <- c(
    exp(drop(t(x) %*% beta[1:3])),
    exp(drop(t(x) %*% beta[4:6])),
    exp(beta[7]),
    exp(drop(t(x) %*% beta[8:10])),
    exp(drop(t(x) %*% beta[11:13])),
    exp(beta[14])
  )

  A <- sum(alpha)
  N <- sum(y)
  dloglik_dalpha <- (
    digamma(y + alpha) -
      digamma(alpha) +
      digamma(A) -
      digamma(N + A)
  )
  dloglik_deta <- alpha * dloglik_dalpha

  as.numeric(c(
    x * dloglik_deta[1],
    x * dloglik_deta[2],
    dloglik_deta[3],
    x * dloglik_deta[4],
    x * dloglik_deta[5],
    dloglik_deta[6]
  ))
}

# The production score used by score_stack(), score_average(), and W_hat_fun().
score_ll_row <- function(j, t, beta, sim_dat) {
  score_ll_row_closed(j, t, beta, sim_dat)
}

score_stack <- function(j, beta, sim_dat, Time = 3) {
  unlist(lapply(seq_len(Time), function(t) score_ll_row(j, t, beta, sim_dat)))
}

score_average <- function(beta, sim_dat, Time = 3) {
  j_list <- sort(unique(sim_dat[, "j"]))
  G <- sapply(j_list, function(j) score_stack(j, beta, sim_dat, Time = Time))
  rowMeans(G)
}

W_hat_fun <- function(beta, sim_dat, Time = 3, ridge = 0) {
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



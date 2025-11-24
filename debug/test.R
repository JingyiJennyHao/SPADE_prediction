### How to start: set working directory and run
### three edit places: J value (number of hospitals); opt_rep1 (number of random starts in beta-step first part); opt_ret2. All marked with "edit X"



library(dplyr)
library(numDeriv)
library(LaplacesDemon)
library(MASS)
library(pracma)
library(Rcpp)
library(MGLM)
library(VGAM)
library(extraDistr)

##setwd("~/Desktop/harm_prediction_project/01_code_simulation")
true_beta <- c(2,0.5,0.5, 1,1,1, 0.5, 0.3,0.3,0.3, 0.2,0.2,0.2, 6)

############################### functions ###############################
write_row <- function(df_row, file) {
  stopifnot(is.data.frame(df_row), nrow(df_row) == 1)
  if (!file.exists(file)) {
    write.table(df_row, file, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    write.table(df_row, file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

softmax3 <- function(a1, a2) {
  den <- 1 + exp(a1) + exp(a2)
  c(exp(a1)/den, exp(a2)/den, 1/den)
}

## dataset
gen_simdata <- function(J, true_beta=true_beta, seed=seed){
  set.seed(seed)
  Time <- 3
  mu <- c(1, 2)
  A  <- diag(c(0.7, 0.5))
  Sigma_eps <- matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)
  
  hos.Xi <- vector("list", J)
  for (j in 1:J) {
    X <- matrix(NA, nrow = Time, ncol = 2)
    X[1, ] <- mu
    for (t in 2:Time) {
      eps <- mvrnorm(1, mu = c(0,0), Sigma = Sigma_eps)
      X[t, ] <- mu + A %*% (X[t-1, ] - mu) + eps
    }
    hos.Xi[[j]] <- rbind(x0 = 1, x1 = X[,1], x2 = X[,2])
  }
  
  true_beta1 <- true_beta[1:3]
  true_beta2 <- true_beta[4:6]
  true_beta3 <- true_beta[7]
  true_beta4 <- true_beta[8:10]
  true_beta5 <- true_beta[11:13]
  true_beta6 <- true_beta[14]
  
  sim_dat <- matrix(NA, nrow = J * Time, ncol = 24)
  colnames(sim_dat) <- c("j","t","Nj","pj1","pj2","pj3","pj4","pj5","pj6",
                         "nj1","nj2","nj3","nj4","nj5","nj6",
                         "a1","a2","a3","a4","a5","a6","x0","x1","x2")
  
  row_idx <- 1
  hosp_eff <- 3
  for (j in 1:J) {
    hj <- rgamma(1, shape = hosp_eff, rate = 10)
    hj <- rep(hj, Time)
    for (t in 1:Time) {
      Xjt <- hos.Xi[[j]][, t]
      shape_jt1 <- exp(Xjt %*% true_beta1)
      shape_jt2 <- exp(Xjt %*% true_beta2)
      shape_jt3 <- exp(true_beta3)
      shape_jt4 <- exp(Xjt %*% true_beta4)
      shape_jt5 <- exp(Xjt %*% true_beta5)
      shape_jt6 <- exp(true_beta6)
      
      rij <- c(
        rgamma(1, shape_jt1, 10),
        rgamma(1, shape_jt2, 10),
        rgamma(1, shape_jt3, 10),
        rgamma(1, shape_jt4, 10),
        rgamma(1, shape_jt5, 10),
        rgamma(1, shape_jt6, 10)
      )
      eij <- rgamma(6, 10, 10)
      p_jt <- (rij + eij) / sum(rij + eij)
      
      Nj  <- sample(100:2000, 1)
      n_jt <- rmultinom(1, size = Nj, prob = p_jt)
      
      alpha_jt <- c(shape_jt1+10, shape_jt2+10, shape_jt3+10,
                    shape_jt4+10, shape_jt5+10, shape_jt6+10)
      
      sim_dat[row_idx, ] <- c(j, t, Nj,
                              p_jt, n_jt,
                              alpha_jt, Xjt)
      row_idx <- row_idx + 1
    }}
  return(sim_dat)
}

## likelihood for each j,t row
ll_row <- function(j, t, beta, sim_dat) {
  idx <- which(sim_dat[, "j"] == j)[t]
  njt <- as.numeric(sim_dat[idx, c('nj1','nj2','nj3','nj4','nj5','nj6')])
  xjt <- c(1, sim_dat[idx,'x1'], sim_dat[idx,'x2'])
  
  a1 <- exp(drop(t(xjt) %*% beta[1:3]))
  a2 <- exp(drop(t(xjt) %*% beta[4:6]))
  a3 <- exp(beta[7])
  a4 <- exp(drop(t(xjt) %*% beta[8:10]))
  a5 <- exp(drop(t(xjt) %*% beta[11:13]))
  a6 <- exp(beta[14])
  alpha <- c(a1,a2,a3,a4,a5,a6)
  
  MGLM::ddirmn(njt, alpha)
}

Phi_i_fun <- function(j, weights, beta, sim_dat, Time) {
  sum(weights * sapply(1:Time, function(t) ll_row(j, t, beta, sim_dat)))
}

neg_comp_llh <- function(beta, weights, sim_dat, J, Time) {
  -sum(sapply(1:J, function(j) Phi_i_fun(j, weights, beta, sim_dat, Time)))
}

# calculate grad and hess for each j and t in term of beta
precompute_grad_hess <- function(beta, sim_dat, J, Time) {
  p <- length(beta)
  G <- vector("list", J)
  H <- vector("list", J)
  for (j in 1:J) {
    Gj <- matrix(0, p, Time)
    Hj <- vector("list", Time)
    for (t in 1:Time) {
      f_jt <- function(b) ll_row(j, t, b, sim_dat)
      Gj[, t] <- numDeriv::grad(f_jt, beta)
      Hj[[t]] <- numDeriv::hessian(f_jt, beta)
    }
    G[[j]] <- Gj
    H[[j]] <- Hj
  }
  result <- list(G = G, H = H, p = p, Time = Time, J=J)
  return(result)
}

# J_mat_cached <- function(weights, cache) {
#   p <- cache$p
#   J <- cache$J
#   U <- matrix(0, p, J)
#   for (j in 1:J) U[, j] <- cache$G[[j]] %*% weights
#   res <- (U %*% t(U)) / (J^2)
#   return(res)
# }

# Var = H^-1 J H^-1

J_mat_cached <- function(weights, cache) {
  p <- cache$p
  J <- cache$J
  U <- matrix(0, p, J)
  for (j in 1:J) U[, j] <- cache$G[[j]] %*% weights  # g_j in column j
  g_bar <- rowMeans(U)            # (1/J) * sum_j g_j
  res   <- g_bar%*%t(g_bar)       # g_bar %*% t(g_bar)
  return(res)
}

H_mat_cached <- function(weights, cache) {
  p <- cache$p
  J <- cache$J
  Hsum <- matrix(0, p, p)
  for (j in 1:J) for (t in 1:cache$Time)
    Hsum <- Hsum + weights[t] * (- cache$H[[j]][[t]])
  return(Hsum / J)
}

var_sandwich_cached <- function(weights, cache) {
  H  <- H_mat_cached(weights, cache)
  Jm <- J_mat_cached(weights, cache)
  Hinv <- solve(H)
  Hinv %*% Jm %*% t(Hinv)
}


## -------------------- One simulation --------------------
run_one_sim <- function(seed, out_csv) {
  set.seed(seed)
  Time <- 3
  true_beta <- c(2,0.5,0.5, 1,1,1, 0.5, 0.3,0.3,0.3, 0.2,0.2,0.2, 6)
  J <- 20                                                                                               # edit 1
  sim_dat <- gen_simdata(J, true_beta = true_beta, seed = seed)
  
  fixed_weights <- c(1,1,1)  # start equal
  beta_best_val <- NA_real_
  beta_hat <- rep(NA_real_, 14)
  beta_conv_flag <- NA_integer_
  w_conv_flag <- NA_integer_
  
  for (loop in 1:2) {
    # β-step:
    opt_rep1 <- 1                                                                                       # edit 2
    optimal_list <- vector("list", opt_rep1)                                                                  
    for (rep in 1:opt_rep1) {                                                                                 
      stpt <- c(2,0.5,0.5, 1,1,1, 0.5, 0.3,0.3,0.3, 0.2,0.2,0.2, 6) + rnorm(14,0,0.1)
      #stpt <- c(runif(13, min=0, max=1), 6)
      fit <- optim(stpt, neg_comp_llh, weights = fixed_weights, sim_dat = sim_dat, J=J, Time=Time, 
                   control = list(maxit = 8000))
      optimal_list[[rep]] <- list(value = fit$value, par = fit$par, convergence = fit$convergence)
    }
    converged <- Filter(function(x) x$convergence == 0, optimal_list)
    if (length(converged) == 0) stop("No converged beta fits.")
    pick <- which.min(sapply(converged, `[[`, "value"))
    beta_hat <- converged[[pick]]$par
    beta_best_val <- converged[[pick]]$value
    beta_conv_flag <- 1L
    
    # cache at beta_hat
    cache <- precompute_grad_hess(beta_hat, sim_dat, J, Time)
    
    # w-step: 2-dim softmax param
    objective_var_fun <- function(alpha) {
      w <- softmax3(alpha[1], alpha[2])
      V <- var_sandwich_cached(w, cache)
      sum(diag(V))
    }
    ow <- optim(par = c(0,0), fn = objective_var_fun, control = list(maxit = 8000))
    w_conv_flag <- as.integer(ow$convergence == 0)
    
    fixed_weights <- softmax3(ow$par[1], ow$par[2])
    
    if (w_conv_flag != 1){
      var_ori <- var_sandwich_cached(c(1,1,1), cache)
      if (ow$value >= sum(diag(var_ori))){
        fixed_weights <- c(1,1,1)
      }
    }
  }
  
  # with final weight, find final beta
  opt_rep2 <- 1                                                                                             # edit 3 
  optimal_list <- vector("list", opt_rep2)                                                                              
  for (rep in 1:opt_rep2) {                                                                                 
    stpt <- c(2,0.5,0.5, 1,1,1, 0.5, 0.3,0.3,0.3, 0.2,0.2,0.2, 6) + rnorm(14,0,0.1)
    fit <- optim(stpt, neg_comp_llh, weights = fixed_weights, sim_dat = sim_dat, J=J, Time=Time,
                 control = list(maxit = 8000))
    optimal_list[[rep]] <- list(value = fit$value, par = fit$par, convergence = fit$convergence)
  }
  converged <- Filter(function(x) x$convergence == 0, optimal_list)
  if (length(converged) == 0) stop("No converged beta fits.")
  pick <- which.min(sapply(converged, `[[`, "value"))
  beta_hat <- converged[[pick]]$par
  beta_best_val <- converged[[pick]]$value
  beta_conv_flag <- 1
  
  # inference
  cache_final <- precompute_grad_hess(beta_hat, sim_dat, J, Time)
  Vhat <- var_sandwich_cached(fixed_weights, cache_final)
  se   <- sqrt(diag(Vhat))
  upper <- beta_hat + 1.96 * se
  lower <- beta_hat - 1.96 * se
  cover <- as.integer((true_beta >= lower) & (true_beta <= upper))
  
  # ------------ Record a single row ------------
  out <- data.frame(
    seed = seed,
    w1 = fixed_weights[1], w2 = fixed_weights[2], w3 = fixed_weights[3],
    beta_conv = beta_conv_flag, w_conv = w_conv_flag,
    beta_val = beta_best_val
  )
  
  # add beta_hat (b1..b14), var (v1..v14), cover (c1..c14)
  names_beta <- paste0("b", 1:14)
  names_var  <- paste0("v", 1:14)
  names_cov  <- paste0("c", 1:14)
  
  out[ names_beta ] <- as.list(beta_hat)
  out[ names_var  ] <- as.list(diag(Vhat))
  out[ names_cov  ] <- as.list(cover)
  
  write_row(out, out_csv)
  invisible(out)
}

##### for cluester setting #####

# args <- commandArgs(trailingOnly = TRUE)
# # usage: Rscript run_one_sim.R <out_csv> <seed>
# out_csv <- ifelse(length(args) >= 1, args[1], "sim_results.csv")
# seed    <- ifelse(length(args) >= 2, as.integer(args[2]), {
#   # try LSF job index, else default
#   ji <- Sys.getenv("LSB_JOBINDEX", unset = NA)
#   if (!is.na(ji) && nzchar(ji)) 1000 + as.integer(ji) else 12345
# })
# cat(sprintf("Running seed=%d, writing to %s\n", seed, out_csv))


# save the resultwith current time
time <- format(Sys.time(), "%m%d_%H%M%S")
seed <- 12345
file_name <- paste0("debug_result_time_", time, ".csv")
run_one_sim(seed, file_name)

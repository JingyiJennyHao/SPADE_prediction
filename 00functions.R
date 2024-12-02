# generate simulated dataset
generate_datasets <- function(beta, sample_size, seed) {
  set.seed(seed)
  
  results_df1 <- data.frame(Ji_list = numeric(),
                            pi1_list = numeric(),
                            pi2_list = numeric(),
                            sum1_list = numeric(),
                            sum2_list = numeric(),
                            sum3_list = numeric(),
                            info_uniform_list = numeric(),
                            info_binary_list = numeric())
  
  for (i in 1:sample_size) {
    Ji <- round(500 * rweibull(1, shape = 1, scale = 1.5))
    info_uniform <- runif(1, min = 0, max = 1)
    info_binary <- sample(0:1, 1)
    info <- c(1, info_uniform, info_binary)
    a1 <- exp(info %*% beta[1:3])
    a2 <- exp(info %*% beta[4:6])
    a3 <- exp(beta[7])
    alpha <- c(a1, a2, a3)
    pi <- rdirichlet(1, alpha)
    pi1 <- pi[1, 1]
    pi2 <- pi[1, 2]
    
    pairs <- sample(c("(1,0)", "(0,1)", "(0,0)"), Ji, replace = TRUE, prob = c(pi1, pi2, 1 - pi1 - pi2))
    matrix_pairs <- matrix(as.numeric(unlist(strsplit(gsub("\\(|\\)", "", pairs), ","))), ncol = 2, byrow = TRUE)
    Dij_1 <- matrix_pairs[, 1]
    Dij_2 <- matrix_pairs[, 2]
    
    sum1 <- sum(Dij_1)
    sum2 <- sum(Dij_2)
    sum3 <- sum((1 - Dij_1) * (1 - Dij_2))
    
    results_df1 <- rbind(results_df1, data.frame(Ji_list = Ji, pi1_list = pi1, pi2_list = pi2, sum1_list = sum1, sum2_list = sum2, sum3_list = sum3, info_uniform_list = info_uniform, info_binary_list = info_binary))
  }
  
  set.seed(seed+5000)
  
  results_df2 <- data.frame(Ji_list = numeric(),
                            pi1_list = numeric(),
                            pi2_list = numeric(),
                            sum1_list = numeric(),
                            sum2_list = numeric(),
                            sum3_list = numeric(),
                            info_uniform_list = numeric(),
                            info_binary_list = numeric())
  
  for (i in 1:sample_size) {
    Ji <- results_df1$Ji_list[i]
    info_uniform <- results_df1$info_uniform_list[i]
    info_binary <- results_df1$info_binary_list[i]
    info <- c(1, info_uniform, info_binary)
    pi1 <- results_df1$pi1_list[i]
    pi2 <- results_df1$pi2_list[i]
    
    pairs <- sample(c("(1,0)", "(0,1)", "(0,0)"), Ji, replace = TRUE, prob = c(pi1, pi2, 1 - pi1 - pi2))
    matrix_pairs <- matrix(as.numeric(unlist(strsplit(gsub("\\(|\\)", "", pairs), ","))), ncol = 2, byrow = TRUE)
    Dij_1 <- matrix_pairs[, 1]
    Dij_2 <- matrix_pairs[, 2]
    
    sum1 <- sum(Dij_1)
    sum2 <- sum(Dij_2)
    sum3 <- sum((1 - Dij_1) * (1 - Dij_2))
    
    results_df2 <- rbind(results_df2, data.frame(Ji_list = Ji, pi1_list = pi1, pi2_list = pi2, sum1_list = sum1, sum2_list = sum2, sum3_list = sum3, info_uniform_list = info_uniform, info_binary_list = info_binary))
  }
  
  return(list(results_df1 = results_df1, results_df2 = results_df2))
}

# likelihood function
llh.mglm <- function(dataset,beta){
  llhi <- sapply(1:nrow(dataset), function(i) {
    ni <- c(dataset[i,'sum1_list'], dataset[i,'sum2_list'], dataset[i,'sum3_list'])
    xi <- c(1,dataset[i,'info_uniform_list'], dataset[i,'info_binary_list'])
    a1 <- exp(t(xi)%*%beta[1:3])
    a2 <- exp(t(xi)%*%beta[4:6])
    a3 <- exp(beta[7])
    alpha <- c(a1,a2,a3)
    return(ddirmn(ni,alpha))})
  NegativeSum <- sum(-llhi)
  return(NegativeSum)}

# Hessian
h <- function(dataset, beta) {
  lpolya_alpha <- function(a) {
    .llh_mglm(dataset,a)
  }
  hess <- hessian(lpolya_alpha, beta)
  return(hess)
}

##### Prediction Interval_Beta Prime #####
PI_betap <- function(alpha,shape1,shape2,scale) {
  ratio_samples <- sort(rbetapr(5000, shape1, shape2, scale))
  pr_samples <- dbetapr(ratio_samples, shape1, shape2, scale)
  cdf_samples <- pbetapr(ratio_samples, shape1, shape2, scale)
  mode <- which.max(pr_samples)
  if (mode != 1){
    intervals <- numeric(mode)
    index1s <- numeric(mode)
    index2s <- numeric(mode)
    for (index in 1:mode) {
      p1 <- pr_samples[index]
      # Find the index of the value in pr_samples[mode+1:] that is closest to p1
      p2_index_range <- (mode + 1):length(pr_samples)
      differences <- abs(pr_samples[p2_index_range] - p1)
      index2 <- p2_index_range[which.min(differences)]
      # Calculate the interval
      interval <- cdf_samples[index2] - cdf_samples[index]
      # Store the interval and indices
      intervals[index] <- interval
      index1s[index] <- index
      index2s[index] <- index2
    }
    # Find the closest interval to 0.95
    closest_index <- which.min(abs(intervals - alpha))
    # Get the corresponding indices
    index <- index1s[closest_index]
    index2 <- index2s[closest_index]
    return(list(l = ratio_samples[index], u = ratio_samples[index2]))
  }
  if (mode == 1){
    upper <- qbetapr(0.95, shape1, shape2, scale)
    return(list(l = 0, u = upper))}
}

##### Prediction Interval_Beta #####
PI_beta <- function(alpha,shape1,shape2,scale) {
  ratio_samples <- sort(rbeta(5000, shape1, shape2, scale))
  pr_samples <- dbeta(ratio_samples, shape1, shape2, scale)
  cdf_samples <- pbeta(ratio_samples, shape1, shape2, scale)
  mode <- which.max(pr_samples)
  if (mode != 1){
    intervals <- numeric(mode)
    index1s <- numeric(mode)
    index2s <- numeric(mode)
    for (index in 1:mode) {
      p1 <- pr_samples[index]
      # Find the index of the value in pr_samples[mode+1:] that is closest to p1
      p2_index_range <- (mode + 1):length(pr_samples)
      differences <- abs(pr_samples[p2_index_range] - p1)
      index2 <- p2_index_range[which.min(differences)]
      # Calculate the interval
      interval <- cdf_samples[index2] - cdf_samples[index]
      # Store the interval and indices
      intervals[index] <- interval
      index1s[index] <- index
      index2s[index] <- index2
    }
    # Find the closest interval to 0.95
    closest_index <- which.min(abs(intervals - alpha))
    # Get the corresponding indices
    index <- index1s[closest_index]
    index2 <- index2s[closest_index]
    return(list(l = ratio_samples[index], u = ratio_samples[index2]))
  }
  if (mode == 1){
    upper <- qbeta(0.95, shape1, shape2, scale)
    return(list(l = 0, u = upper))}
}

##### Prediction Interval_BetaSum #####
PI_SumBeta <- function(alpha,alpha_hat){
  n_sim <- 2000
  sample <- rdirichlet(2000,alpha_hat)
  sum_samples <- sort(sample[,2]+sample[,3])
  density_estimate <- density(sum_samples)
  pr <- (approx(density_estimate$x, density_estimate$y, xout = sum_samples)$y)
  mode_index <- which.max(pr)
  
  intervals <- numeric(mode_index)
  index1s <- numeric(mode_index)
  index2s <- numeric(mode_index)
  
  if (mode_index != 1){
    for (index in 1:mode_index) {
      p1 <- pr[index]
      index2_range <- (mode_index + 1):length(sum_samples)
      differences <- abs(pr[index2_range] - p1)
      index2 <- index2_range[which.min(differences)]
      interval <- index2-index
      intervals[index] <- interval
      index1s[index] <- index
      index2s[index] <- index2
    }
    closest_index <- which.min(abs(intervals/length(sum_samples) - alpha))
    # Get the corresponding indices
    index1 <- index1s[closest_index]
    index2 <- index2s[closest_index]
    return(list(l = sum_samples[index1], u = sum_samples[index2]))
  }
  if (mode_index == 1){
    upper <- quantile(sum_samples,0.95)
    return(list(l = 0, u = upper))}
}

##### Prediction Interval_Binomial #####
PI_binomial <- function(alpha, N, p){
  upper <- N*p + 1.96*sqrt(N*p*(1-p))
  lower <- N*p - 1.96*sqrt(N*p*(1-p))
  return(list(l = lower, u = upper))
} 

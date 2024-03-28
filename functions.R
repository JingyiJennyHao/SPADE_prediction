#sourceCpp("llh.cpp")

################################ simulated dataset ################################
dataset <- function(sample_size,alpha){
  results_df <- data.frame(Ji_list = numeric(),
                           zi_list = numeric(),
                           pi1_list = numeric(),
                           pi2_list = numeric(),
                           sum1_list = numeric(),
                           sum2_list = numeric(),
                           sum3_list = numeric())
  
  for (i in 1:sample_size) {
    Ji <- round(500 * rweibull(1, shape = 1, scale = 1.5))
    Z_i <- sample(c(0, 1), 1)
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
    
    # Append the values to the data frame
    results_df <- rbind(results_df, data.frame(Ji_list = Ji, pi1_list = pi1, pi2_list = pi2, sum1_list = sum1, sum2_list = sum2, sum3_list = sum3))
  }
  
  return(results_df)
}

################################# parameter scale #################################
# to make sure alpha is in (0,inf)
logit <- function(k){
  temp <- log((k-0.5)/(150.5-k))
  return(temp)}

expit <- function(k){
  temp <- 150*exp(k)/(1+exp(k))+0.5
  return(temp)}


#################### generate complete log-likelihood function ####################
# li: dirichlet-multinomial
li <- function(ni1,ni2,ni3,p,alpha){
  temp <- log_factorial(ni1 + ni2 + ni3) - log_factorial(ni1) - log_factorial(ni2) - log_factorial(ni3) + ni1 * log(p[1]) + ni2 * log(p[2]) + ni3 * log(p[3]) + ddirichlet(p, alpha, log = TRUE)
  return(temp)
}

# complete log-likelihood
llhc_R <- function(t_alpha, t_params, simulated_dataset) {
  alpha <- expit(t_alpha)
  params <- expit(t_params)
  sum_temp_mean <- 0
  pis <- rdirichlet(1000, params)
  sum_temp_mean <- sum(sapply(1:nrow(simulated_dataset), function(i) {
    ni1 <- simulated_dataset$sum1_list[i]
    ni2 <- simulated_dataset$sum2_list[i]
    ni3 <- simulated_dataset$sum3_list[i]
    log_likelihoods <- apply(pis, 1, function(pi_row) {
      li(ni1, ni2, ni3, pi_row, alpha)})
    return(mean(log_likelihoods))
  }))
  return(-sum_temp_mean)
}

#################################### gradient ####################################
g <- function(x,alpha){
  log_ddirichlet <- function(x, alpha) {
    return(ddirichlet(x, alpha, log=TRUE))}
  gradient <- grad(function(a) log_ddirichlet(x, a), alpha)
  return(gradient)}

# li_grad <- function(ni1,ni2,ni3,p,alpha){
#   temp <-  g(p,alpha)
#   return(temp)
# }

# gradient of complete likelihood
S <- function(t_alpha, t_params, simulated_dataset) {
  alpha <- expit(t_alpha)
  params <- expit(t_params)
  #sum_temp_mean <- 0
  pis <- rdirichlet(1000, params)
  sum_temp_mean <- sapply(1:nrow(simulated_dataset), function(i) {
    ni1 <- simulated_dataset$sum1_list[i]
    ni2 <- simulated_dataset$sum2_list[i]
    ni3 <- simulated_dataset$sum3_list[i]
    log_likelihoods <- rowMeans(apply(pis, 1, function(pi_row) {
      g(pi_row, alpha)
    }))
    mat <- matrix(as.vector(log_likelihoods), nrow = 3, byrow = TRUE)
    mat_list <- mat%*%t(mat) 
    return(mat_list) # It is a list of all 3*3 matrices
  })
  matrices_list <- lapply(1:ncol(sum_temp_mean), function(i) {
    matrix(sum_temp_mean[, i], nrow = 3, byrow = TRUE)
  })
  result <- Reduce(`+`, matrices_list)
  return(result) 
}

# gradient of individual
Si <- function(t_alpha, t_params, simulated_dataset) {
  alpha <- expit(t_alpha)
  params <- expit(t_params)
  sum_temp_mean <- list()  # Initialize as a list
  pis <- rdirichlet(1000, params)
  sum_temp_mean <- lapply(1:nrow(simulated_dataset), function(i) {
    ni1 <- simulated_dataset$sum1_list[i]
    ni2 <- simulated_dataset$sum2_list[i]
    ni3 <- simulated_dataset$sum3_list[i]
    log_likelihoods <- rowMeans(apply(pis, 1, function(pi_row) {
      li_grad(ni1, ni2, ni3, pi_row, alpha)
    }))
    mat <- matrix(as.vector(log_likelihoods), nrow = 3, byrow = TRUE)
    return(mat)
  }) 
  return(sum_temp_mean) 
}

##################################### hessian #####################################
h <- function(x, alpha) {
  log_ddirichlet <- function(a) {
    return(ddirichlet(x, a, log=TRUE))
  }
  hess <- hessian(log_ddirichlet, alpha)
  return(hess)
}

# li_hessian <- function(ni1,ni2,ni3,p,alpha){
#   temp <- h(p,alpha)
#   return(temp)}

# hessian of complete likelihood
B <- function(t_alpha, t_params, simulated_dataset) {
  alpha <- expit(t_alpha)
  params <- expit(t_params)
  pis <- rdirichlet(1000, params)
  sum_temp_mean <- sapply(1:nrow(simulated_dataset), function(i) {
    ni1 <- simulated_dataset$sum1_list[i]
    ni2 <- simulated_dataset$sum2_list[i]
    ni3 <- simulated_dataset$sum3_list[i]
    
    log_likelihoods <- apply(pis, 1, function(pi_row){
      li_hess <- h(pi_row, alpha)
      return(li_hess)})
    
    log_likelihoods_list <- lapply(1:ncol(log_likelihoods), function(i) {
      matrix(log_likelihoods[, i], nrow = 3, byrow = TRUE)})
    
    mean_matrix <- Reduce(`+`, log_likelihoods_list) / length(log_likelihoods_list)
    return(mean_matrix)})
  
  matrices_list <- lapply(1:ncol(sum_temp_mean), function(i) {
    matrix(sum_temp_mean[, i], nrow = 3, byrow = TRUE)})
  
  return(Reduce(`+`, matrices_list))}



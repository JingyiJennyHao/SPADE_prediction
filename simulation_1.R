### Attempt 1225  Merry Christmas!!!
### Update: Something wrong with the likelihood function
### Updated 0106  Revise the function of likelihood
### Update 0112  Revise the function of Likelihood
### Update 0116 Still llh

# library(MASS)
# library(mvtnorm)

likelihood_function <- function(mu_hat, sig_hat) {
  # Define the integrand function
  f <- function(pi1, pi2) {
    term1 <- pi1^sum(Dij_1) * pi2^sum(Dij_2) * (1 - pi1 - pi2)^sum((1 - Dij_1) * (1 - Dij_2))
    term2 <- dmvnorm(c(pi1, pi2), mean = mu_hat, sigma = sig_hat)
    return(term1 * term2)
  }
  
  # Perform the double integral
  result <- integrate(Vectorize(function(pi1) sapply(pi1, function(pi2) f(pi1, pi2))), -Inf, Inf)
  
  # Return the result
  return(result)
}

generate_data <- function(simulation) {
  results_list <- list()
  f <- function(params) { 0 }
  
  for (i in 1:simulation) {
    J_i <- round(500 * rweibull(1, shape = 1, scale = 1.5))
    Z_i <- sample(c(0, 1), 1)
    X_i <- t(cbind(1, J_i, Z_i))
    
    # True parameter values
    beta1 <- runif(3, -10, 10)
    beta2 <- runif(3, -10, 10)
    sig1 <- runif(1)
    sig2 <- runif(1)
    rho <- runif(1)
    sig12 <- rho * sig1 * sig2
    
    # MVN
    mean_vector <- c(1 / (1 + exp(t(X_i) %*% beta1)), 1 / (1 + exp(t(X_i) %*% beta2)))
    covariance_matrix <- matrix(c(sig1^2, sig12, sig12, sig2^2), nrow = 2)
    prob_vector <- mvrnorm(n = 1, mu = mean_vector, Sigma = covariance_matrix)
    pi1 <- prob_vector[1]
    pi2 <- prob_vector[2]
    
    # Ensure pi1 and pi2 are within the [0, 1] range and pi1 + pi2 <= 1
    while (!(0 <= pi1 && pi1 <= 1) || !(0 <= pi2 && pi2 <= 1) || !((pi1 + pi2) <= 1)) {
      prob_vector <- mvrnorm(n = 1, mu = t(mean_vector), Sigma = covariance_matrix)
      pi1 <- prob_vector[1]
      pi2 <- prob_vector[2]
    }
    
    # Generate Dij_1 and Dij_2
    Dij_1 <- sample(c(0, 1), J_i, replace = TRUE)
    Dij_2 <- sample(c(0, 1), J_i, replace = TRUE)
    indice <- which(Dij_1 == 1 & Dij_2 == 1)
    Dij_1[indice] <- 0
    
    # Calculate pi1_hat and pi2_hat
    pi1_hat <- sum(Dij_1) / J_i
    pi2_hat <- sum(Dij_2) / J_i
    
    fi <- function(params) {
      mu_hat <- params[1:2]
      sig_hat <- matrix(c(params[3], params[4], params[4], params[5]), nrow = 2)
      result <- -likelihood_function(mu_hat, sig_hat)$value
      return(result)
    }
    
    # Update the f function
    f <- function(params) { f(params) + log(fi(params)) }
    
    # Store the results in a named list
    simulation_result <- list(
      J_i = J_i,
      X_i = X_i,
      beta1 = beta1,
      beta2 = beta2,
      sig1 = sig1,
      sig2 = sig2,
      rho = rho,
      sig12 = sig12,
      mean_vector = mean_vector,
      covariance_matrix = covariance_matrix,
      pi1 = pi1,
      pi2 = pi2,
      pi1_hat = pi1_hat,
      pi2_hat = pi2_hat
    )
    
    # Append the simulation result to the list
    results_list[[i]] <- simulation_result
  }
  
  # Return both the list of results and the final f function
  return(list(results_list = results_list, final_f = f))
}

Dataset <- generate_data(10)

llh <- function(params){Dataset$final_f(params)}

initial_params <- c(0, 0, 0, 0, 0)

# Use optim to find the parameter values that minimize the likelihood
result <- optim(par = initial_params, fn = llh, method = "L-BFGS-B")

# Extract the optimized parameters
optimized_params <- result$par
mu_hat_optimized <- optimized_params[1:2]
sig_hat_params_optimized <- optimized_params[3:5]
sig_hat_optimized <- matrix(c(sig_hat_params_optimized[1], sig_hat_params_optimized[2], sig_hat_params_optimized[2], sig_hat_params_optimized[3]), nrow = 2)

# Print the optimized parameters
print(mu_hat_optimized)
print(sig_hat_optimized)

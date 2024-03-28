library(numDeriv)
library(LaplacesDemon)
library(MASS)
library(pracma)
library(Rcpp)

source("functions.R")
sourceCpp("llh.cpp")

################################### EM algorithm #######################################
sample_size <- 150
real_alpha <- c(1,1,100)
dataset1 <- dataset(sample_size,real_alpha)

initial_alpha <- logit(c(1,1,100))
initial_params <- logit(c(1,1,100))
for (i in 1:100) {
  result <- optim(par = initial_alpha, fn = llhc, t_params = initial_params, simulated_dataset=dataset1, control = list(maxit=100))
  if (all(abs(initial_alpha - result$par) < 1e-2)) {
    break
  }
  initial_alpha <- result$par
  initial_params <- result$par
  print(initial_alpha)
}

################################# Statistical Inference #################################
#calculate three parts of Observed Information matrix
# reference: http://www.jstor.org/stable/2345828
params_hat <- expit(initial_alpha)
part1 <- B(params_hat, params_hat, simulated_dataset) 
part2 <- S(params_hat, params_hat, simulated_dataset) 
Si_list <- Si(params_hat, params_hat, simulated_dataset)
size_Si_list <- length(Si_list)
sum_Si_ij <- matrix(0, nrow = 3, ncol = 3)

# Loop over all pairs of indices (i, j) where i < j
for (i in 1:(size_Si_list - 1)) {
  for (j in (i + 1):size_Si_list) {
    sum_Si_ij <- sum_Si_ij + Si_list[[i]] %*% t(Si_list[[j]])
  }
}
part3 <- 2*sum_Si_ij

fisher_matrix <- -part1 - part2 - part3

sigma <- solve(fisher_matrix)
z <- qnorm(0.975)
SE <- sqrt(diag(sigma))
lower_ci <- params_hat - z * SE
upper_ci <- params_hat + z * SE

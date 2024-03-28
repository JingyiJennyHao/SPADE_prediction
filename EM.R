# EM

library(numDeriv)
library(LaplacesDemon)
library(MASS)
library(pracma)
library(Rcpp)

source("functions.R")
sourceCpp("llh.cpp")

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
  

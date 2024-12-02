library(dplyr)
library(numDeriv)
if (!requireNamespace("LaplacesDemon", quietly = TRUE)) {
  install.packages("LaplacesDemon")
}
library(LaplacesDemon)
library(MASS)
library(pracma)
if (!requireNamespace("MGLM", quietly = TRUE)) {
  install.packages("MGLM")
}
library(MGLM)
if (!requireNamespace("extraDistr", quietly = TRUE)) {
  install.packages("extraDistr")
}
library(extraDistr)
if (!requireNamespace("VGAM", quietly = TRUE)) {
  install.packages("VGAM")
}
library(VGAM)

naive_summary <- read.csv("ni_list.csv")
hos_cov <- read.csv("hos_cov.csv")

data_scaled <- scale(hos_cov[, c("ct", "mri", "num_patients")])
sum1_list <- naive_summary$sum1_list
sum2_list <- naive_summary$sum2_list
sum3_list <- naive_summary$sum3_list
data <- as.data.frame(cbind(sum1_list, sum2_list, sum3_list, data_scaled))

.llh_mglm <- function(dataset, beta, cov_vec) {
  llhi <- sapply(1:nrow(dataset), function(i) {
    ni <- c(dataset[i, 'sum1_list'], dataset[i, 'sum2_list'], dataset[i, 'sum3_list'])
    xi <- c(1, as.numeric(dataset[i, cov_vec]))
    k <- length(xi)
    beta1 <- beta[1:k]
    beta2 <- beta[(k + 1):(2 * k)]  
    beta3 <- beta[length(beta)]
    a1 <- exp(t(xi)%*%beta1)
    a2 <- exp(t(xi)%*%beta2)
    a3 <- exp(beta3)
    alpha <- c(a1, a2, a3)
    return(ddirmn(ni, alpha))
  })
  NegativeSum <- -sum(llhi)
  return(NegativeSum)
}

h <- function(dataset, beta, cov_vec) {
  lpolya_alpha <- function(a) {
    .llh_mglm(dataset, a, cov_vec)
  }
  hess <- hessian(lpolya_alpha, beta)
  return(hess)
}

replications <- 150
summary <- matrix(NA, ncol = 12, nrow = replications)

for (i in 1:replications) {
  set.seed(i)
  print(i)
  # attempt 1
  stpt <- runif(8, min = -3, max = 3)
  stpt <- c(stpt, log(3000))

  # attempt 2
  # stpt <- runif(9, min = -2, max = 2)

  # attempt 3
  # stpt <- runif(9, min = 0, max = 3)


  cov1 <- c("ct", "mri", "num_patients")
  
  # Using tryCatch to handle errors and continue to next iteration
  result <- tryCatch({
    optim(stpt, .llh_mglm, dataset = data, cov_vec = cov1, control = list(maxit = 5000))
  }, error = function(e) {
    message("Error at seed ", i, ": ", conditionMessage(e))
    return(NULL) # Return NULL in case of an error
  })
  
  # Check if result is NULL and skip to the next iteration if so
  if (is.null(result)) {
    next
  }
  
  HessMatrix <- h(data, result$par, cov1)
  eigenvalues1 <- eig(HessMatrix)
  saddle <- if (all(eigenvalues1 > 0)) 0 else 1
  summary[i, ] <- c(result$par, result$value, result$convergence, saddle)
}

saveRDS(summary, "realdata1.RDS")


library(MGLM)
library(numDeriv)
library(MASS)
library(pracma)
library(Rcpp)
library(dplyr)
library(LaplacesDemon)
library(extraDistr)
library(VGAM)
library(ggplot2)


beta_hat <- c(1.699798, 0.10270310, -0.1722666, 0.1326662, 1.948011, 0.007584371, -0.1347554, -0.3090747, 8.494690)


hos_cov <- read.csv("~/Desktop/harm_prediction_project/02_code_data/hos_cov_0526.csv")
mean(hos_cov$num_patients)
var(hos_cov$num_patients)


# generate syn_cov
set.seed(42)  

# Number of rows
rows <- 216

# Generate the data of ct, mri, num of patients
ct <- round(runif(rows, 0, 1), 3)
mri <- round(runif(rows, 0, 1), 3)
raw_num_patients <- rnorm(rows, 0, 1)
scaled_mean <- 1563
scaled_var <- 1143664
scaled_sd <- sqrt(scaled_var)
num_patients <- round((raw_num_patients - mean(raw_num_patients)) / sd(raw_num_patients) * scaled_sd + scaled_mean)
num_patients[num_patients < 0] <- 0  # Ensure all values are >= 0

syn_cov <- data.frame(ct = ct, mri = mri, num_patients = num_patients)
write.csv(syn_cov, file = "syn_cov.csv", row.names = FALSE)


# Extract beta coefficients
beta1 <- beta_hat[1:4]
beta2 <- beta_hat[5:8]
beta3 <- beta_hat[9]


syn_cov_scale <- scale(syn_cov[, c("ct", "mri", "num_patients")])

syn_ni_list <- matrix(0, nrow = nrow(syn_cov), ncol = 3)

for (i in 1:nrow(syn_cov)) {
  xi <- syn_cov_scale[i, ]
  
  # Calculate alpha values
  a1 <- exp(sum(xi * beta1))  # Sum performs matrix-vector multiplication
  a2 <- exp(sum(xi * beta2))
  a3 <- exp(beta3)
  alpha <- c(a1, a2, a3)
  
  # Generate ni using Dirichlet distribution
  ni <- rdirichlet(1, alpha)
  
  # Assign ni values scaled to the number of patients
  syn_ni_list[i, 1] <- round(ni[1] * syn_cov$num_patients[i])
  syn_ni_list[i, 2] <- round(ni[2] * syn_cov$num_patients[i])
  syn_ni_list[i, 3] <- syn_cov$num_patients[i] - syn_ni_list[i, 1] - syn_ni_list[i, 2]
}

# Convert to a data frame for easier handling
colnames(syn_ni_list) <- c("sum1_list", "sum2_list", "sum3_list")
syn_ni_list <- as.data.frame(syn_ni_list)

write.csv(syn_cov, file = "syn_ni_list.csv", row.names = FALSE)

```

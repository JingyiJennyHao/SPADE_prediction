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
source("functions.R")

##### Simulation ##### 
# seed from 1 to 1000

replications <- 1000
sample_size <- 200
summary <- matrix(NA, ncol = 39, nrow = replications)

for (i in 1:replications){
  true_beta <- c(0.5,0.5,0.5,0.5,0.5,0.5,7)
  datasets <- generate_datasets(true_beta, sample_size, i)
  results_df1 <- datasets$results_df1
  results_df2 <- datasets$results_df2
  stpt <- c(0.5,0.5,0.5,0.5,0.5,0.5,7)
  result <- optim(stpt, llh.mglm, dataset = results_df1, control = list(maxit = 3000))
  print(result$par)
  HessMatrix <- h(results_df1, result$par)
  sigma <- solve(HessMatrix)
  z <- qnorm(0.975)
  SE <- sqrt(diag(sigma))
  lower_ci <- result$par - z * SE
  upper_ci <- result$par + z * SE
  eigenvalues1 = eig(HessMatrix)
  beta_hat <- result$par
  bias <- beta_hat - true_beta
  emp_se <- SE 
  mse <- (beta_hat - true_beta) ^ 2
  # betahat sqrt(variance) ==> ESE
  coverage <- as.integer(lower_ci <= true_beta & true_beta <= upper_ci)
  convergence <- 1-result$convergence
  saddle <- if (all(eigenvalues1 > 0)) 1 else 0
  # check in the test set
  PI_rate11 <- 0
  PI_rate21 <- 0
  PI_rate31 <- 0
  PI_rate41 <- 0
  PI_rate51 <- 0
  
  PI_rate12 <- 0
  PI_rate22 <- 0
  PI_rate32 <- 0
  PI_rate42 <- 0
  PI_rate52 <- 0
  
  for (j in (1:nrow(results_df2))){
    xi <- c(1,results_df2[j,'info_uniform_list'], results_df2[j,'info_binary_list'])
    pi1_star <- results_df2[j,"pi1_list"]
    pi2_star <- results_df2[j,"pi2_list"]
    a1 <- exp(t(xi)%*%beta_hat[1:3])
    a2 <- exp(t(xi)%*%beta_hat[4:6])
    a3 <- exp(beta_hat[7])
    alpha_i <- c(a1,a2,a3)
    
    # check pi1 (under alpha)
    shape11 <- alpha_i[1] 
    shape12 <- alpha_i[2]+alpha_i[3]
    scale <- 1
    interval11 <- PI_beta(0.95,shape11,shape12,scale)
    PI_rate11_i <- as.integer((pi1_star >= interval11$l) && (pi1_star <= interval11$u))
    PI_rate11 <- PI_rate11 + PI_rate11_i
    # check pi1 (under alpha and ni)
    shape11 <- alpha_i[1] + results_df1[j,"sum1_list"]
    shape12 <- alpha_i[2] + alpha_i[3] + results_df1[j,"sum2_list"] + results_df1[j,"sum3_list"]
    scale <- 1
    interval12 <- PI_beta(0.95,shape11,shape12,scale)
    PI_rate12_i <- as.integer((pi1_star >= interval12$l) && (pi1_star <= interval12$u))
    PI_rate12 <- PI_rate12 + PI_rate12_i
    
    # check pi2 (under alpha)
    shape21 <- alpha_i[2]
    shape22 <- alpha_i[1]+alpha_i[3]
    interval21 <- PI_beta(0.95,shape21,shape22,scale)
    PI_rate21_i <- as.integer((pi2_star >= interval21$l) && (pi2_star <= interval21$u))
    PI_rate21 <- PI_rate21 + PI_rate21_i
    # check pi2 (under alpha and ni)
    shape21 <- alpha_i[2] + results_df1[j,"sum2_list"]
    shape22 <- alpha_i[1] + alpha_i[3] + results_df1[j,"sum1_list"] + results_df1[j,"sum3_list"]
    interval22 <- PI_beta(0.95,shape21,shape22,scale)
    PI_rate22_i <- as.integer((pi2_star >= interval22$l) && (pi2_star <= interval22$u))
    PI_rate22 <- PI_rate22 + PI_rate22_i
    
    # check pi1+pi2 (under alpha)
    shape31 <- alpha_i[1]+alpha_i[2] 
    shape32 <- alpha_i[3]
    interval31 <- PI_beta(0.95,shape31,shape32,scale)
    PI_rate31_i <- as.integer(((pi1_star+pi2_star) >= interval31$l) && ((pi1_star+pi2_star) <= interval31$u))
    PI_rate31 <- PI_rate31 + PI_rate31_i
    # check pi1+pi2 (under alpha and ni)
    shape31 <- alpha_i[1]+alpha_i[2]+ results_df1[j,"sum1_list"] + results_df1[j,"sum2_list"]
    shape32 <- alpha_i[3]+results_df1[j,"sum3_list"]
    interval32 <- PI_beta(0.95,shape31,shape32,scale)
    PI_rate32_i <- as.integer(((pi1_star+pi2_star) >= interval32$l) && ((pi1_star+pi2_star) <= interval32$u))
    PI_rate32 <- PI_rate32 + PI_rate32_i
    
    # check pi2/pi1 (under alpha)
    ratio_star <- results_df2[j,"pi2_list"]/results_df2[j,"pi1_list"]
    shape41 <- alpha_i[2]
    shape42 <- alpha_i[1]
    interval41 <- PI_betap(0.95,shape41,shape42,scale)
    PI_rate41_i <- as.integer((ratio_star >= interval41$l) && (ratio_star <= interval41$u))
    PI_rate41 <- PI_rate41 + PI_rate41_i
    
    # check pi2/pi1 (under alpha and ni)
    shape41 <- results_df1[j,"sum2_list"]+alpha_i[2]
    shape42 <- results_df1[j,"sum1_list"]+alpha_i[1]
    interval42 <- PI_betap(0.95,shape41,shape42,scale)
    PI_rate42_i <- as.integer((ratio_star >= interval42$l) && (ratio_star <= interval42$u))
    PI_rate42 <- PI_rate42 + PI_rate42_i
    
    # check ni2_star: Binomial 
    ni2_star <- results_df2[j,"sum2_list"]
    size <- results_df1[j,"sum2_list"]+results_df1[j,"sum3_list"]
    interval5 <- PI_binomial(0.95, size, alpha_i[2]/(alpha_i[2]+alpha_i[3]))
    PI_rate51_i <- as.integer((ni2_star >= interval5$l) && (ni2_star <= interval5$u))
    PI_rate51 <- PI_rate51 + PI_rate51_i
  } 
  PI_percentage11 <- PI_rate11 / nrow(results_df2)
  PI_percentage12 <- PI_rate12 / nrow(results_df2)
  PI_percentage21 <- PI_rate21 / nrow(results_df2)
  PI_percentage22 <- PI_rate22 / nrow(results_df2)
  PI_percentage31 <- PI_rate31 / nrow(results_df2)
  PI_percentage32 <- PI_rate32 / nrow(results_df2)
  PI_percentage41 <- PI_rate41 / nrow(results_df2)
  PI_percentage42 <- PI_rate42 / nrow(results_df2)
  PI_percentage51 <- PI_rate51 / nrow(results_df2)
  summary[i, ] <- c(bias, emp_se, mse, coverage, convergence, saddle, PI_percentage11, PI_percentage12, PI_percentage21, PI_percentage22, PI_percentage31,PI_percentage32, PI_percentage41,PI_percentage42, PI_percentage51)
}

combined_data <- summary

column_names <- c(paste("bias", 1:7, sep = ""), 
                  paste("se", 1:7, sep = ""), 
                  paste("mse", 1:7, sep = ""), 
                  paste("coverage", 1:7, sep = ""),
                  paste("convergence","_",sep=""),
                  paste("saddle","_",sep=""),
                  paste("PI_rate11","_",sep=""),
                  paste("PI_rate12","_",sep=""),
                  paste("PI_rate21","_",sep=""),
                  paste("PI_rate22","_",sep=""),
                  paste("PI_rate31","_",sep=""),
                  paste("PI_rate32","_",sep=""),
                  paste("PI_rate41","_",sep=""),
                  paste("PI_rate42","_",sep=""),
                  paste("PI_rate51","_",sep="")
                  )

# Assign column names to the data frame
colnames(combined_data) <- column_names

# Filter rows where the 29th column equals 1
filtered_data <- combined_data[combined_data[,29] == 1, ]

#Select the first 1000 rows from the filtered data
first_1000 <- filtered_data[1:1000, ]

means <- colMeans(first_1000)
print(means)

ese_PI_rate1 <- sd(first_1000[31,]) / sqrt(1000)
ese_PI_rate2 <- sd(first_1000[32,]) / sqrt(1000)
ese_PI_rate3 <- sd(first_1000[33,]) / sqrt(1000)
ese_PI_rate4 <- sd(first_1000[34,]) / sqrt(1000)
ese_PI_rate5 <- sd(first_1000[35,]) / sqrt(1000)

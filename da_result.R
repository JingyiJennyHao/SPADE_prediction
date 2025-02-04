---
title: "data_pre_0914"
author: "Jingyi Hao"
date: "2024-09-14"
output: html_document
---

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
library(ggplot2)
library(patchwork)
library(cowplot)


file1 <- readRDS("/Users/jingyihao/Desktop/harm_prediction_project/00_result/Real_Data_Analysis/result_1111/realdata1.RDS")
file2 <- readRDS("/Users/jingyihao/Desktop/harm_prediction_project/00_result/Real_Data_Analysis/result_1111/realdata2.RDS")
file3 <- readRDS("/Users/jingyihao/Desktop/harm_prediction_project/00_result/Real_Data_Analysis/result_1111/realdata3.RDS")


ordered_file1 <- file1[order(file1[, 10]), ]
ordered_file2 <- file2[order(file2[, 10]), ]
ordered_file3 <- file3[order(file3[, 10]), ]
beta_hat <- ordered_file1[1,1:9]
print(beta_hat)


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
  ratio_samples <- sort(rbeta(10000, shape1, shape2, scale))
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


naive_summary <- read.csv("~/Desktop/harm_prediction_project/02_code_data/ni_list_0526.csv")
hos_cov <- read.csv("~/Desktop/harm_prediction_project/02_code_data/hos_cov_0526.csv")
#naive_summary <- read.csv("ni_list_0526.csv")
#hos_cov <- read.csv("hos_cov_0526.csv")

data_scaled <- scale(hos_cov[, c("ct", "mri", "num_patients")])
sum1_list <- naive_summary$sum1_list
sum2_list <- naive_summary$sum2_list
sum3_list <- naive_summary$sum3_list
data <- as.data.frame(cbind(sum1_list, sum2_list, sum3_list, data_scaled))


# Confidence interval of the parameter
cov1 <- c("ct", "mri", "num_patients")
HessMatrix <- h(data, beta_hat, cov1)
sigma <- solve(HessMatrix)
z <- qnorm(0.975)
SE <- sqrt(diag(sigma))
lower_ci <- beta_hat - z * SE
upper_ci <- beta_hat + z * SE

print(lower_ci)
print(upper_ci)


## Inference Interval
cov1 <- c("ct", "mri", "num_patients")

# Initialize an empty data frame to store the results
inference_result <- data.frame()

for (i in 1:nrow(data)) {
  # Extract the counts and covariates for the current observation
  ni <- c(data[i, 'sum1_list'], data[i, 'sum2_list'], data[i, 'sum3_list'])
  xi <- c(1, as.numeric(data[i, cov1]))
  
  # Calculate the alpha parameters
  a1 <- exp(xi %*% beta_hat[1:4])
  a2 <- exp(xi %*% beta_hat[5:8])
  a3 <- exp(beta_hat[9])
  alpha_i <- c(a1, a2, a3)
  
  # Compute the intervals
  ### pi1
  shape11 <- alpha_i[1]
  shape12 <- alpha_i[2] + alpha_i[3]
  scale <- 1
  interval_pi1 <- PI_beta(0.95, shape11, shape12, scale)
  
  ### pi2
  shape21 <- alpha_i[2]
  shape22 <- alpha_i[1] + alpha_i[3]
  interval_pi2 <- PI_beta(0.95, shape21, shape22, scale)
  
  ### pi1 + pi2
  shape31 <- alpha_i[1] + alpha_i[2]
  shape32 <- alpha_i[3]
  interval_plus <- PI_beta(0.95, shape31, shape32, scale)
  
  ### ni2
  size <- data[i, "sum2_list"] + data[i, "sum3_list"]
  prob <- alpha_i[2] / (alpha_i[2] + alpha_i[3])
  interval_ni2 <- PI_binomial(0.95, size, prob)
  
  ### pi2 / pi2|Xi
  shape51 <- data[i, "sum2_list"] + alpha_i[2]
  shape52 <- data[i, "sum1_list"] + alpha_i[1]
  interval_ratio <- PI_betap(0.95, shape51, shape52, scale)
  
  # Store the results in a data frame for the current observation
  obs_result <- data.frame(
    sum1_list = data[i, 'sum1_list'],
    sum2_list = data[i, 'sum2_list'],
    sum3_list = data[i, 'sum3_list'],
    interval_pi1_l = interval_pi1$l,
    interval_pi1_u = interval_pi1$u,
    interval_pi2_l = interval_pi2$l,
    interval_pi2_u = interval_pi2$u,
    interval_plus_l = interval_plus$l,
    interval_plus_u = interval_plus$u,
    interval_ni2_l = interval_ni2$l,
    interval_ni2_u = interval_ni2$u,
    interval_ratio_l = interval_ratio$l,
    interval_ratio_u = interval_ratio$u
  )
  
  # Append the current observation's results to the inference_result data frame
  inference_result <- rbind(inference_result, obs_result)
}


print(inference_result)


g <- function(beta,dataset,hospital_i){
  ni1 <- dataset[hospital_i,"sum1_list"]
  ni2 <- dataset[hospital_i,"sum2_list"]
  ni3 <- dataset[hospital_i,"sum3_list"]
  beta1_hat <- beta[1:4]
  beta2_hat <- beta[5:8]
  beta3_hat <- beta[9]
  cov1 <- c("ct","mri","num_patients")
  xi <- c(1, as.numeric(dataset[hospital_i, cov1]))
  gxi <- (exp(t(xi)%*%beta1_hat)+exp(t(xi)%*%beta2_hat)+ni1+ni2)/(exp(t(xi)%*%beta1_hat)+exp(t(xi)%*%beta2_hat)+exp(beta3_hat)+ni1+ni2+ni3) 
  return(as.numeric(gxi))
}

g_prime <- function(beta, dataset, hospital_i) {
  temp <- function(b) g(b, dataset, hospital_i)
  grad_value <- grad(temp, beta)
  return(grad_value)
}

cov1 <- c("ct","mri","num_patients")
I <- h(data, beta_hat, cov1)
g_val <- g(beta_hat,data,1)
g_grad <- g_prime(beta_hat,data,1)
se <- t(g_grad)%*%solve(I)%*%g_grad
CI_lower <- g_val - 1.96*se
CI_upper <- g_val + 1.96*se


g_values <- nrow(data)
for (i in 1:216) {
  g_values[i] <- g(beta_hat, data, i)
  }
g_median <- median(g_values)
print(g_median)


beta1_hat <- beta_hat[1:4]
beta2_hat <- beta_hat[5:8]
beta3_hat <- beta_hat[9]
num_hospitals <- nrow(data)
CI_results <- data.frame(hospital = 1:num_hospitals, CI_lower = NA, CI_upper = NA, Epi1 = NA, Epi12 = NA, Epi2 = NA, Epi22 = NA, EpiSum = NA, EpiSum2 = NA, Epiratio= NA,Epiratio2= NA, contains_median = NA, above_median = NA)
I_matrix <- h(data, beta_hat, cov1)
I_invse <- solve(I_matrix)
for (i in 1:num_hospitals) {
  ni1 <- data[i,"sum1_list"]
  ni2 <- data[i,"sum2_list"]
  ni3 <- data[i,"sum3_list"]
  xi <- c(1, as.numeric(data[i, cov1]))
  exp1 <- exp(t(xi)%*%beta1_hat)
  exp2 <- exp(t(xi)%*%beta2_hat)
  exp3 <- exp(beta3_hat)
  # E[pi1] based on Xi
  Epi1 <- exp1/(exp1+exp2+exp3)
  # E[pi1] based on n & Xi
  Epi12 <- (ni1+exp1)/(ni1+ni2+ni3+exp1+exp2+exp3)
  # E[pi2]
  Epi2 <- exp2/(exp1+exp2+exp3)
  Epi22 <- (ni2+exp2)/(ni1+ni2+ni3+exp1+exp2+exp3)
  # E[pi1+pi2]
  EpiSum <- (exp1+exp2)/(exp1+exp2+exp3)
  EpiSum2 <- (ni1+exp1+ni2+exp2)/(ni1+ni2+ni3+exp1+exp2+exp3)
  # E[pi2/pi1]
  Epiratio <- exp2/(exp1-1)
  Epiratio2 <- (ni2+exp2)/(ni1+exp1-1)
  
  # Calculate g value and gradient for the current hospital
  g_val <- g(beta_hat, data, i)
  g_grad <- g_prime(beta_hat, data, i)

  # Calculate standard error
  se <- sqrt(t(g_grad) %*% I_invse %*% g_grad)

  # Calculate CI lower and upper bounds of E[pi1+pi2]
  CI_lower <- g_val - 1.96 * se
  CI_upper <- g_val + 1.96 * se
  contains_median <- ifelse(g_median < CI_lower, -1,  # g_median is less than CI_lower
                   ifelse(g_median >= CI_lower & g_median <= CI_upper, 0,  # g_median is within the CI range
                   ifelse(g_median > CI_upper, 1, NA)))
  above_median <- ifelse(EpiSum > g_median, 1, 0)
  
  # Store results in the data frame
  CI_results$CI_lower[i] <- CI_lower
  CI_results$CI_upper[i] <- CI_upper
  CI_results$Epi1[i] <- Epi1
  CI_results$Epi12[i] <- Epi12
  CI_results$Epi2[i] <- Epi2
  CI_results$Epi22[i] <- Epi22
  CI_results$EpiSum[i] <- EpiSum
  CI_results$EpiSum2[i] <- EpiSum2
  CI_results$Epiratio[i] <- Epiratio
  CI_results$Epiratio2[i] <- Epiratio2
  CI_results$contains_median[i] <- contains_median
  CI_results$above_median[i] <-above_median
}


result_plot <- cbind(CI_results, hos_cov$ct, hos_cov$mri, hos_cov$num_patients)
print(result_plot)


# Define x and y variables for plotting
x_vars <- c("hos_cov$ct", "hos_cov$mri", "hos_cov$num_patients")
y_vars <- c("Epi12", "Epi22", "EpiSum2", "Epiratio2")

# Define axis labels for clarity
x_labels <- c("CT Rate", "MRI Rate", "hospital size")
y_labels <- c(
  "Same hospital return rate",
  "Different hospital return rate",
  "Total return rate",
  "Cross-over ratio"
)

plot_labels <- c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3", "D1", "D2", "D3")

# Initialize empty list to store plots
plot_list <- list()

# Loop through y variables and x variables to create the 4x3 grid of plots
plot_count <- 1
for (j in seq_along(y_vars)) {
  for (i in seq_along(x_vars)) {
    
    # Create each plot without a legend
    plot <- ggplot(result_plot, aes_string(x = x_vars[i], y = y_vars[j], color = "as.factor(contains_median)")) +
      geom_point(size = 3, alpha = 0.7) +  # Scatter points
      scale_color_manual(values = c("-1" = "red", "0" = "black", "1" = "blue")) +
      labs(
        x = x_labels[i],                  # Dynamic x-axis label
        y = "Expectation", # Dynamic y-axis label
        title = paste0(plot_labels[plot_count], ": ", y_labels[j], " vs ", x_labels[i])  # Dynamic title
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = -0.05, size = 16),  # Left align the title
        axis.title.x = element_text(size = 16),             # Increase x-axis title size
        axis.title.y = element_text(size = 16),
        legend.position = "none",                           # Remove legend from individual plots
        plot.margin = margin(20, 20, 20, 20)                # Add margins around the plot
      )
    
    # Store plot in list
    plot_list[[plot_count]] <- plot
    plot_count <- plot_count + 1
  }
}

# Create a dummy plot to extract the legend
dummy_plot <- ggplot(result_plot, aes_string(x = x_vars[1], y = y_vars[1], color = "as.factor(contains_median)")) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("-1" = "red", "0" = "black", "1" = "blue"), 
                     labels = c("-1" = "Needing Improvement", "0" = "Average", "1" = "Better Performing")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  ) +
  guides(color = guide_legend(nrow = 1))  # Arrange legend items in a single row

# Extract the legend from the dummy plot
legend <- cowplot::get_legend(dummy_plot)

# Combine the 12 plots into a 4x3 grid
combined_plots <- plot_grid(plotlist = plot_list, ncol = 3, nrow = 4)

# Combine the grid with the shared legend below
final_plot <- plot_grid(combined_plots, legend, ncol = 1, rel_heights = c(1, 0.1))

# Display the final combined plot
print(final_plot)

ggsave("final_plot_with_shared_legend_1726.png", plot = final_plot, width = 20, height = 24, dpi = 300)



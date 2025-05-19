rm(list = ls())
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

#### set up parameters ####
# modifiable parameters
trt_effect_rel = 0.9
n_perwsq = 15
prop_attribution = 0

const <- 900

# relatively fixed parameters
alpha_31_90_z = 0
alpha_90_plus_z = 0
alpha_1_30_t = 0
alpha_31_90_t = 0
n_wave = 9
n_hosp_perwave = 3
n_quarter = 10
n_quarter_burn = 0
seed = 79
n_hosp_perwave_sim = 100
alpha = 0.05

#### run power calculation ####
# load proportion of variances attributable to hospital/quarter
load("propvar_AVERT_ICC.rda")

#### set up fixed parameters for simulation ####
# set up alpha0 parameters for each interval
## 1-30
alpha_1_30_0_residual <- log(mean_prop_1_30 * prop_residual_1_30)
alpha_1_30_0_hospital <- log(mean_prop_1_30 * prop_hospital_1_30)
alpha_1_30_0_quarter <- log(mean_prop_1_30 * prop_quarter_1_30)
alpha_1_30_0 <- log(mean_prop_1_30)
## 31-90
alpha_31_90_0_residual <- log(mean_prop_31_90 * prop_residual_31_90)
alpha_31_90_0_hospital <- log(mean_prop_31_90 * prop_hospital_31_90)
alpha_31_90_0_quarter <- log(mean_prop_31_90 * prop_quarter_31_90)
alpha_31_90_0 <- log(mean_prop_31_90)
## 91-180
alpha_90_plus_0_residual <- log(mean_prop_90_plus * prop_residual_90_plus)
alpha_90_plus_0_hospital <- log(mean_prop_90_plus * prop_hospital_90_plus)
alpha_90_plus_0_quarter <- log(mean_prop_90_plus * prop_quarter_90_plus)
alpha_90_plus_0 <- log(mean_prop_90_plus)
#### set up changing parameters for testing ####

# account for attribution
n_dizzyout_perwsq <- floor(n_perwsq * (1 - prop_attribution))
alpha_1_30_z <- log(1 - trt_effect_rel)
print(alpha_1_30_z)

# generate data for n_hosp_perwave_sim hospitals
# for each wave, block out the transition time
# for wave 1, block out 2:(n_quarter_burn + 1)
# for wave 2, block out 3:(n_quarter_burn + 2) etc...

## generate quarter, wave, and hospital_id
quarter_out <- NULL
trt_out <- NULL
wave_out <- NULL
hospital_id_out <- NULL
for (wave in 1:n_wave) {
      quarter_vec <- 1:n_quarter
      if (n_quarter_burn != 0) {
            quarter_vec[1:n_quarter_burn + wave] <- NA
      }
      quarter_vec <- quarter_vec[which(!is.na(quarter_vec))]
      quarter_out <- c(quarter_out, 
                       rep(quarter_vec, each = n_hosp_perwave_sim))
      
      trt_vec <- rep(1, n_quarter - n_quarter_burn)
      trt_vec[1:wave] <- 0
      
      trt_out <- c(trt_out, rep(trt_vec, each = n_hosp_perwave_sim))
      
      wave_out <- c(wave_out, 
                    rep(wave, 
                        each = n_hosp_perwave_sim * 
                              (n_quarter - n_quarter_burn)))
      hospital_id_out <- c(hospital_id_out, 
                           rep(1:n_hosp_perwave_sim + 
                                     (wave - 1) * n_hosp_perwave_sim,
                               (n_quarter - n_quarter_burn)))
}
# wave_out <- rep(1:n_wave, 
#                 each = n_hosp_perwave_sim * 
#                       (n_quarter - n_quarter_burn))


set.seed(seed)
## generate quarter-specific p component
p_1_30_quarter_vec <- rgamma(n = n_quarter,
                             shape = exp(alpha_1_30_0_quarter) * const,
                             rate = 100)
p_31_90_quarter_vec <- rgamma(n = n_quarter,
                              shape = exp(alpha_31_90_0_quarter) * const,
                              rate = 100)
p_90_plus_quarter_vec <- rgamma(n = n_quarter,
                                shape = exp(alpha_90_plus_0_quarter) * const,
                                rate = 100)
p_quarter_dat <- data.frame(quarter = 1:n_quarter,
                            p_1_30_quarter = p_1_30_quarter_vec,
                            p_31_90_quarter = p_31_90_quarter_vec,
                            p_90_plus_quarter = p_90_plus_quarter_vec)
## generate hospital-specific p component
p_1_30_hospital_vec <- rgamma(n = n_hosp_perwave_sim * n_wave,
                              shape = exp(alpha_1_30_0_hospital) * const,
                              rate = 100)
p_31_90_hospital_vec <- rgamma(n = n_hosp_perwave_sim * n_wave,
                               shape = exp(alpha_31_90_0_hospital) * const,
                               rate = 100)
p_90_plus_hospital_vec <- rgamma(n = n_hosp_perwave_sim * n_wave,
                                 shape = exp(alpha_90_plus_0_hospital) * const,
                                 rate = 100)
p_hospital_dat <- data.frame(hospital_id = 1:(n_hosp_perwave_sim * n_wave),
                             p_1_30_hospital = p_1_30_hospital_vec,
                             p_31_90_hospital = p_31_90_hospital_vec,
                             p_90_plus_hospital = p_90_plus_hospital_vec)

## generate residual p component
p_1_30_residual_vec <- rgamma(n = n_hosp_perwave_sim * n_wave * 
                                    (n_quarter - n_quarter_burn),
                              shape = exp(alpha_1_30_0_residual) * const,
                              rate = 100)
p_31_90_residual_vec <- rgamma(n = n_hosp_perwave_sim * n_wave * 
                                     (n_quarter - n_quarter_burn),
                               shape = exp(alpha_31_90_0_residual) * const,
                               rate = 100)
p_90_plus_residual_vec <- rgamma(n = n_hosp_perwave_sim * n_wave * 
                                       (n_quarter - n_quarter_burn),
                                 shape = exp(alpha_90_plus_0_residual) * const,
                                 rate = 100)
p_residual_dat <- data.frame(hospital_id = hospital_id_out,
                             quarter = quarter_out,
                             wave = wave_out,
                             trt = trt_out,
                             p_1_30_residual = p_1_30_residual_vec,
                             p_31_90_residual = p_31_90_residual_vec,
                             p_90_plus_residual = p_90_plus_residual_vec)

## merge three components into one data.frame
p_dat <- merge(p_residual_dat, p_hospital_dat, 
               by = c("hospital_id"), 
               all.x = TRUE)
# View(p_dat[which(p_dat[, "hospital_id"] == 1), ])
# View(p_dat[which(p_dat[, "hospital_id"] == 2000), ])
p_dat <- merge(p_dat, p_quarter_dat,
               by = c("quarter"), 
               all.x = TRUE)
p_dat$p_1_30_raw <- p_dat$p_1_30_hospital +
      p_dat$p_1_30_quarter + p_dat$p_1_30_residual
p_dat$p_31_90_raw <- p_dat$p_31_90_hospital +
      p_dat$p_31_90_quarter + p_dat$p_31_90_residual
p_dat$p_90_plus_raw <- p_dat$p_90_plus_hospital +
      p_dat$p_90_plus_quarter + p_dat$p_90_plus_residual

## compute probabilities after accounting time and treatment effect
p_dat$pt_1_30 <- exp(alpha_1_30_t * p_dat$quarter)
p_dat$pt_31_90 <- exp(alpha_31_90_t * p_dat$quarter)
p_dat$pt_90_plus <- exp(0  * p_dat$quarter)
p_dat$pz_1_30 <- exp(alpha_1_30_z * p_dat$trt)
p_dat$pz_31_90 <- exp(alpha_31_90_z * p_dat$trt)
p_dat$pz_90_plus <- exp(0 * p_dat$trt)
p_dat$p_1_30_unsd <- p_dat$p_1_30_raw * 
      p_dat$pt_1_30 * p_dat$pz_1_30
p_dat$p_31_90_unsd <- p_dat$p_31_90_raw *
      p_dat$pt_31_90 * p_dat$pz_31_90
p_dat$p_90_plus_unsd <- p_dat$p_90_plus_raw *
      p_dat$pt_90_plus * p_dat$pz_90_plus
## no need to standardize to sum of 1, because these are parameters of the dirichlet distribution


## generate of sum1_list, sum2_list, sum3_list for each row
## from dirichlet distribution with parameters 
## (p_1_30_unsd, p_31_90_unsd, p_90_plus_unsd)
temp <- t(sapply(1:nrow(p_dat), function(i) {
      result <- rdirichlet(n = 1,
                           alpha = c(p_dat$p_1_30_unsd[i],
                                     p_dat$p_31_90_unsd[i],
                                     p_dat$p_90_plus_unsd[i]))
      return(result)
}))
# temp <- rdirichlet(n = nrow(p_dat),
#                    alpha = cbind(p_dat$p_1_30_unsd,
#                                  p_dat$p_31_90_unsd,
#                                  p_dat$p_90_plus_unsd))
## simulate counts and then proportion using the dirichlet prob
n_1_30 <- sapply(1:nrow(p_dat), function(i) {
      rbinom(n = 1, size = n_dizzyout_perwsq, 
             prob = temp[i, 1])
})
n_31_90 <- sapply(1:nrow(p_dat), function(i) {
      rbinom(n = 1, size = n_dizzyout_perwsq, 
             prob = temp[i, 2])
})
n_90_plus <- sapply(1:nrow(p_dat), function(i) {
      rbinom(n = 1, size = n_dizzyout_perwsq, 
             prob = temp[i, 3])
})
# n_1_30 <- rbinom(nrow(p_dat), size = rep(n_dizzyout_perwsq,
#                                          nrow(p_dat)),
#                  prob = temp[, 1])
# n_31_90 <- rbinom(nrow(p_dat), size = rep(n_dizzyout_perwsq,
#                                           nrow(p_dat)),
#                   prob = temp[, 2])

p_dat$sum1_list <- n_1_30
p_dat$sum2_list <- n_31_90
p_dat$sum3_list <- n_dizzyout_perwsq - n_1_30 - n_31_90
mean(p_dat$sum1_list) # should be equal to n_dizzyout_perwsq)

# check if each hospital have quarters 1:10
n_quarter_by_hosp <- aggregate(p_dat$quarter, 
                               by = list(p_dat$hospital_id), 
                               FUN = function(x) length(unique(x)))
table(n_quarter_by_hosp$x)
# check if each hospital have treatment, and number of trt = 0 equal to wave
n_trt_by_hosp <- aggregate(p_dat$trt, 
                           by = list(p_dat$hospital_id), 
                           FUN = function(x) sum((x ==1)))
table(n_trt_by_hosp$x)
# check if each hospital belongs to one wave
n_wave_by_hosp <- aggregate(p_dat$wave, 
                            by = list(p_dat$hospital_id), 
                            FUN = function(x) length(unique(x)))
table(n_wave_by_hosp$x)

# construct function that computes log-likelihood for each row
llhi_fun <- function(beta) {
      llhi <- sapply(1:nrow(p_dat), function(i) {
            ni <- c(p_dat[i,'sum1_list'], 
                    p_dat[i,'sum2_list'], 
                    p_dat[i,'sum3_list'])
            xi <- c(1, p_dat$quarter[i], p_dat$trt[i])
            a1 <- exp(t(xi) %*% beta[1:3])
            a2 <- exp(t(xi) %*% beta[4:6])
            a3 <- exp(beta[7])
            alpha <- c(a1,a2,a3)
            return(ddirmn(ni,alpha))
      })
      return(llhi)
}

# llh_hosp <- function(beta) {
#       llhi <- llhi_fun(beta = beta)
#       p_dat$llhi <- llhi
#       temp <- aggregate(llhi ~ hospital_id, 
#                         data = p_dat, sum)
#       return(temp$llhi)
# }
# 
# llh_quarter <- function(beta) {
#       llhi <- llhi_fun(beta = beta)
#       p_dat$llhi <- llhi
#       temp <- aggregate(llhi ~ quarter, 
#                         data = p_dat, sum)
#       return(temp$llhi)
# }

# llh <- function(beta) {
#       llhi <- sapply(1:nrow(p_dat), function(i) {
#             ni <- c(p_dat$sum1_list[i], 
#                     p_dat$sum2_list[i], 
#                     p_dat$sum3_list[i])
#             xi <- c(1, p_dat$quarter[i], p_dat$trt[i])
#             a1 <- exp(sum(xi * beta[1:3]))
#             a2 <- exp(sum(xi * beta[4:6]))
#             a3 <- exp(beta[7])
#             alpha <- c(a1,a2,a3)
#             return(ddirmn(ni,alpha))
#       })
#       return(sum(llhi))
# }


llh_dat <- function(beta, dat) {
      llhi <- sapply(1:nrow(dat), function(i) {
            ni <- c(dat$sum1_list[i], 
                    dat$sum2_list[i], 
                    dat$sum3_list[i])
            xi <- c(1, dat$quarter[i], dat$trt[i])
            a1 <- exp(sum(xi * beta[1:3]))
            a2 <- exp(sum(xi * beta[4:6]))
            a3 <- exp(beta[7])
            alpha <- c(a1,a2,a3)
            return(ddirmn(ni,alpha))
      })
      return(mean(llhi))
}



# get jacobian of llhi_hosp and hessian of llh at beta_true

# known true parameters
beta_true <- c(alpha_1_30_0_residual + log(const) - log(100),
               alpha_1_30_t,
               alpha_1_30_z,
               alpha_31_90_0 + log(const) - log(100),
               alpha_31_90_t,
               alpha_31_90_z,
               alpha_90_plus_0 + log(const) - log(100))


# fit_optim <- optim(par = beta_0, 
#                    fn = llh, 
#                    method = "BFGS",
#                    control = list(fnscale = -1, 
#                                   maxit = 1000, 
#                                   reltol = 1e-8), 
#                    hessian = T)
# llh_hess <- fit_optim$hessian # check if it's positive definite
# print(is.positive.definite(-llh_hess))
# fit_optim$convergence
# beta_true <- fit_optim$par

# llh_hosp_grad <- numDeriv::jacobian(func = llh_hosp, 
#                                     x = beta_true, 
#                                     method = "Richardson")

# # compute hospital-quarter level jacobian
# llh_i_grad <- numDeriv::jacobian(func = llhi_fun, 
#                                     x = beta_true, 
#                                     method = "Richardson")

# compute the score function for each hospital
score_ij <- sapply(1:nrow(p_dat), function(i) {
      print(i)
      dat <- p_dat[i, , drop = F]
      tempfun <- function(beta) {
            result <- llh_dat(beta = beta, dat = dat)
            return(result) 
      }
      llh_i_grad <- numDeriv::jacobian(func = tempfun,
                                       x = beta_true, 
                                       method = "Richardson")
      return(llh_i_grad)
})
p_dat <- cbind(p_dat, t(score_ij))

# aggregate s_ij into s_i, so that for hospital i, s_i is the concatenation of s_ij over batches
score_i <- NULL
# j is the batch number, not necessarily the quarter number 
# quarter = (22 - j - wave) %% 10; e.g., j = 1, wave = 1, quarter = 10 (not zero)
for (i in 1:(n_hosp_perwave_sim * n_wave)) {
      score_tempi <- NULL
      for (j in 1:n_quarter) {
            # quarter_temp <- n_quarter - p_dat$wave[p_dat$hospital_id == i]
            quarter_temp <- (22 - j - unique(p_dat$wave[p_dat$hospital_id == i])) %% 10
            if (quarter_temp == 0) {
                  quarter_temp <- 10
            }
            if (any(p_dat$hospital_id == i & 
                    p_dat$quarter == quarter_temp)) {
                  score_tempi <- c(score_tempi, 
                                   unlist(p_dat[p_dat$hospital_id == i & 
                                                      p_dat$quarter == quarter_temp, 29:35]))
            } else {
                  score_tempi <- c(score_tempi, 
                                   rep(0, 7)) # 7 is the number of parameters
            }
            # if (j %in% c(1, 10)) {
            #       score_tempi <- c(score_tempi, 
            #                        unlist(p_dat[p_dat$hospital_id == i & 
            #                                           p_dat$quarter == j, 
            #                                     c(29:30, 32:33, 35)]))
            # }
      }
      score_i <- rbind(score_i, score_tempi)
}
# score_i <- score_i[, -c(3, 66)]
dim(score_i) # n_hospital_perwave_sim x (n_quarter * p=7)
Matrix::rankMatrix(score_i) # check the rank of the matrix

# temp <- t(score_i) %*% score_i  - 
#       colMeans(score_i) %*% t(colMeans(score_i))

score_i_covmat <- t(score_i) %*% score_i / nrow(score_i) - 
      colMeans(score_i) %*% t(colMeans(score_i))
# score_i_covmat <- cov(score_i)
dim(score_i_covmat) # 70 x 70
which(diag(score_i_covmat) == 0)
det(score_i_covmat) # check if the covariance matrix is singular



# compute the hessian of llh for each quarter 12.
llh_hess_quarterstacked <- NULL
for (batch in 1:n_quarter) {
      print(batch)
      dat <- NULL
      for (wave in 1:n_wave) {
            quarter_temp <- (22 - batch - wave) %% 10
            if (quarter_temp == 0) {
                  quarter_temp <- 10
            }
            tempdat <- p_dat[which(p_dat$wave == wave & 
                                         p_dat$quarter == quarter_temp), , drop = F]
            dat <- rbind(dat, tempdat)
      }
      tempfun <- function(beta) {
            result <- llh_dat(beta = beta, 
                              dat = dat)
            return(result)
      }
      hess_j <- numDeriv::hessian(
            func = tempfun,
            x = beta_true,
            method = "Richardson"
      )
      # if (j %in% c(1, 10)) {
      #       hess_j <- hess_j[-c(3, 6), ]
      # }
      llh_hess_quarterstacked <- rbind(llh_hess_quarterstacked, 
                                       hess_j)
}
dim(llh_hess_quarterstacked)
Matrix::rankMatrix(llh_hess_quarterstacked) # check the rank of the matrix



# compute optimal variance-covariance for GMM
library(MASS)
# temp = ginv(score_i_covmat) # compute the generalized inverse of the covariance matrix)
# temp <- (temp + t(temp)) / 2 # make sure it's symmetric
temp <- solve(score_i_covmat) # even if it's invertible, it can still be ill-conditioned
Sigma <- solve(t(llh_hess_quarterstacked) %*%
                     temp %*%
                     (llh_hess_quarterstacked)) / 
      (n_hosp_perwave * n_wave) 
#/ n_hosp_perwave_sim # normalizing constant?
var_eff_z <- Sigma[3, 3]
print(var_eff_z)
# Compare optimization methods locally for one GMM loop.
#
# Usage:
#   source("try_optim_methods.R")
#
# Optional settings before source():
#   test_seed <- 1
#   test_J <- 200
#   test_maxit <- 10000
#   test_base_beta_file <- "base_beta_loop2.txt"
#
# If test_base_beta_file is NULL or does not exist, true_beta is used as the
# starting beta. The script uses the same simulated data, weight matrix, and
# starting beta for both methods.

if (!exists("test_seed", inherits = FALSE)) test_seed <- 1L
if (!exists("test_J", inherits = FALSE)) test_J <- 200L
if (!exists("test_Time", inherits = FALSE)) test_Time <- 3L
if (!exists("test_maxit", inherits = FALSE)) test_maxit <- 10000L
if (!exists("test_base_beta_file", inherits = FALSE)) {
  test_base_beta_file <- "base_beta_loop2.txt"
}

script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (!is.null(script_path) && nzchar(script_path)) {
  dirname(normalizePath(script_path, mustWork = TRUE))
} else {
  getwd()
}

source(file.path(script_dir, "gmm_core.R"))

read_test_beta <- function(file) {
  if (is.null(file) || !nzchar(file) || !file.exists(file)) {
    message("Using true_beta as the starting beta.")
    return(true_beta)
  }

  beta <- scan(file, quiet = TRUE)
  if (length(beta) != length(true_beta) || any(!is.finite(beta))) {
    stop("Expected 14 finite beta values in ", file)
  }
  message("Using starting beta from: ", normalizePath(file))
  beta
}

base_beta <- read_test_beta(test_base_beta_file)
datasets <- gen_simdata(test_J, true_beta = true_beta, seed = test_seed)
sim_dat <- datasets$sim_dat1

# Match the weight update used by loop 2 and later.
W_current <- W_hat_fun(base_beta, sim_dat, Time = test_Time)

run_method <- function(method) {
  started <- Sys.time()
  fit <- tryCatch(
    optim(
      par = base_beta,
      fn = gmm_obj,
      sim_dat = sim_dat,
      W_hat = W_current,
      Time = test_Time,
      method = method,
      control = list(maxit = test_maxit, trace = 0)
    ),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  if (inherits(fit, "error")) {
    return(data.frame(
      method = method,
      convergence = NA_integer_,
      value = NA_real_,
      diff_from_start = NA_real_,
      elapsed_seconds = elapsed,
      message = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }

  result <- data.frame(
    method = method,
    convergence = fit$convergence,
    value = fit$value,
    diff_from_start = max(abs(fit$par - base_beta)),
    elapsed_seconds = elapsed,
    message = if (!is.null(fit$message)) fit$message else "",
    stringsAsFactors = FALSE
  )
  result[paste0("b", seq_along(fit$par))] <- as.list(fit$par)
  result
}

cat("\nGMM optimization comparison\n")
cat("seed=", test_seed, " J=", test_J, " maxit=", test_maxit, "\n", sep = "")
cat("Starting beta:\n")
print(base_beta)

comparison <- rbind(
  run_method("Nelder-Mead"),
  run_method("BFGS")
)

cat("\nResults:\n")
print(comparison[, c(
  "method", "convergence", "value", "diff_from_start",
  "elapsed_seconds", "message"
)], row.names = FALSE)

cat("\nConvergence codes: 0 usually means successful convergence.\n")
cat("A lower value is better only because both methods use the same data and W.\n")

invisible(comparison)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
source(file.path(script_dir, "gmm_utils.R"))

usage <- paste(
  "Usage:",
  "Rscript gmm_runner.R --groups <int> --seed <int> --data <path> --beta0-file <path> --out <path>",
  "[--n-starts <int>] [--start-scale <numeric>] [--ridge <numeric>] [--maxit <int>] [--uncentered-v <0|1>]",
  sep = "\n"
)

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

require_args(args, c("groups", "seed", "data", "beta0-file", "out"))

group_count <- as_int_arg(args, "groups")
seed <- as_int_arg(args, "seed")
data_path <- args[["data"]]
beta0_file <- args[["beta0-file"]]
output_path <- args[["out"]]
n_starts <- if (!is.null(args[["n-starts"]])) as_int_arg(args, "n-starts") else 100L
maxit <- if (!is.null(args[["maxit"]])) as_int_arg(args, "maxit") else 8000L
start_scale <- if (!is.null(args[["start-scale"]])) as.numeric(args[["start-scale"]]) else 0.1
ridge <- if (!is.null(args[["ridge"]])) as.numeric(args[["ridge"]]) else 1e-6
center_v <- !identical(args[["uncentered-v"]] %||% "0", "1")

if (!file.exists(beta0_file)) {
  stop(sprintf("beta0 file not found: %s", beta0_file), call. = FALSE)
}
if (!is.finite(start_scale) || start_scale < 0) {
  stop("Argument --start-scale must be a non-negative numeric value.", call. = FALSE)
}
if (!is.finite(ridge) || ridge < 0) {
  stop("Argument --ridge must be a non-negative numeric value.", call. = FALSE)
}
if (n_starts < 1L) {
  stop("Argument --n-starts must be at least 1.", call. = FALSE)
}

set.seed(seed)
prepared <- prepare_data(data_path, group_count)
beta0 <- load_beta_reference(beta0_file)

fit <- try(
  fit_gmm_multistart(
    beta0 = beta0,
    sim_dat = prepared$data,
    J = prepared$J,
    Time = prepared$Time,
    n_starts = n_starts,
    seed = seed,
    start_scale = start_scale,
    ridge = ridge,
    maxit = maxit,
    center_v = center_v
  ),
  silent = TRUE
)

result <- list(
  stage = "gmm",
  seed = seed,
  groups = group_count,
  data_path = data_path,
  beta0_file = beta0_file,
  n_starts = n_starts,
  start_scale = start_scale,
  ridge = ridge,
  maxit = maxit,
  center_v = center_v,
  convergence = 99L,
  objective_value = Inf,
  beta0 = beta0,
  beta_hat = rep(NA_real_, length(beta0)),
  best_start = NA_integer_,
  message = NULL
)

if (inherits(fit, "try-error")) {
  result$message <- as.character(fit)
  save_result(output_path, result)
  quit(save = "no", status = 1)
}

result$convergence <- fit$convergence
result$objective_value <- fit$objective_value
result$beta_hat <- fit$beta_hat
result$best_start <- fit$best_start
result$message <- fit$message
result$V <- fit$V
result$W <- fit$W
result$S0 <- fit$S0
result$starts <- fit$starts
result$fits <- fit$fits

save_result(output_path, result)

if (as.integer(result$convergence) != 0L || !is.finite(result$objective_value)) {
  quit(save = "no", status = 1)
}

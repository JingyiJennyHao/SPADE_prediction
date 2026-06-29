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
  "Rscript gmm_start_runner.R --groups <int> --task-id <int> --seed <int> --data <path> --beta0-file <path> --out <path>",
  "[--start-scale <numeric>] [--ridge <numeric>] [--maxit <int>] [--uncentered-v <0|1>]",
  sep = "\n"
)

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

require_args(args, c("groups", "task-id", "seed", "data", "beta0-file", "out"))

group_count <- as_int_arg(args, "groups")
task_id <- as_int_arg(args, "task-id")
seed <- as_int_arg(args, "seed")
data_path <- args[["data"]]
beta0_file <- args[["beta0-file"]]
output_path <- args[["out"]]
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

set.seed(seed)
prepared <- prepare_data(data_path, group_count)
beta0 <- load_beta_reference(beta0_file)
start <- if (task_id == 1L) as.numeric(beta0) else draw_gmm_start(beta0, start_scale = start_scale)

fit <- fit_gmm_single_start(
  start = start,
  beta0 = beta0,
  sim_dat = prepared$data,
  J = prepared$J,
  Time = prepared$Time,
  ridge = ridge,
  maxit = maxit,
  center_v = center_v
)

result <- list(
  stage = "gmm_start",
  task_id = task_id,
  seed = seed,
  groups = group_count,
  data_path = data_path,
  beta0_file = beta0_file,
  start_scale = start_scale,
  ridge = ridge,
  maxit = maxit,
  center_v = center_v,
  beta0 = beta0,
  start = start,
  convergence = fit$convergence,
  value = fit$value,
  par = fit$par,
  initial_value = fit$initial_value,
  message = fit$message,
  V = fit$V,
  W = fit$W,
  S0 = fit$S0
)

save_result(output_path, result)

if (as.integer(result$convergence) != 0L || !is.finite(result$value)) {
  quit(save = "no", status = 1)
}

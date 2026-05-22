source("composite_utils.R")

usage <- paste(
  "Usage:",
  "Rscript stage2_weight_runner.R --groups <int> --task-id <int> --seed <int>",
  "--beta-file <path> --data <path> --out <path>",
  sep = "\n"
)

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

require_args(args, c("groups", "task-id", "seed", "beta-file", "data", "out"))

group_count <- as_int_arg(args, "groups")
task_id <- as_int_arg(args, "task-id")
seed <- as_int_arg(args, "seed")
beta_file <- args[["beta-file"]]
data_path <- args[["data"]]
output_path <- args[["out"]]

if (!file.exists(beta_file)) {
  stop(sprintf("Stage 1 summary file not found: %s", beta_file), call. = FALSE)
}

set.seed(seed)
prepared <- prepare_data(data_path, group_count)
beta_summary <- readRDS(beta_file)

if (is.null(beta_summary$beta_hat)) {
  stop("Stage 1 summary file does not contain beta_hat.", call. = FALSE)
}

beta_hat <- beta_summary$beta_hat
cache <- precompute_grad_hess(beta_hat, prepared$data, prepared$J, prepared$Time)
objective_var_fun <- objective_var_fun_factory(cache)
start_weights <- normalize_weights(c(1, 1, 1) + rnorm(3, mean = 0, sd = 0.1))
start_alpha <- alpha_from_weights(start_weights)

fit <- try(
  optim(
    par = start_alpha,
    fn = objective_var_fun,
    control = list(maxit = 8000)
  ),
  silent = TRUE
)

result <- list(
  stage = "stage2",
  task_id = task_id,
  seed = seed,
  groups = group_count,
  beta_file = beta_file,
  beta_hat = beta_hat,
  start_alpha = start_alpha,
  start_weights = start_weights,
  convergence = 99L,
  objective_value = Inf,
  weights = c(NA_real_, NA_real_, NA_real_),
  alpha = c(NA_real_, NA_real_),
  baseline_trace = NA_real_,
  optimized_trace = NA_real_,
  message = NULL
)

if (inherits(fit, "try-error")) {
  result$message <- as.character(fit)
  save_result(output_path, result)
  quit(save = "no", status = 1)
}

weights_hat <- softmax3(fit$par[1], fit$par[2])
baseline_trace <- sum(diag(var_sandwich_cached(start_weights, cache)))
optimized_trace <- sum(diag(var_sandwich_cached(weights_hat, cache)))

result$convergence <- as.integer(fit$convergence)
result$objective_value <- fit$value
result$weights <- weights_hat
result$alpha <- fit$par
result$baseline_trace <- baseline_trace
result$optimized_trace <- optimized_trace
result$message <- fit$message %||% NULL

save_result(output_path, result)

if (as.integer(fit$convergence) != 0L || !is.finite(fit$value)) {
  quit(save = "no", status = 1)
}

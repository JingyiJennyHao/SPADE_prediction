source("cl_sim_utils.R")

usage <- paste(
  "Usage:",
  "Rscript weight_job.R --round <int> --task-id <int> --seed <int> --sim-seed <int>",
  "--groups <int> --beta-file <path> --out-csv <path> --maxit <int>",
  "[--alpha-ref-file <path>]",
  sep = "\n"
)

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

require_args(args, c("round", "task-id", "seed", "sim-seed", "groups", "beta-file", "out-csv", "maxit"))

round_id <- as_int_arg(args, "round")
task_id <- as_int_arg(args, "task-id")
seed <- as_int_arg(args, "seed")
sim_seed <- as_int_arg(args, "sim-seed")
J <- as_int_arg(args, "groups")
beta_file <- args[["beta-file"]]
out_csv <- args[["out-csv"]]
maxit <- as_int_arg(args, "maxit")
alpha_ref_file <- args[["alpha-ref-file"]] %||% ""

beta_hat <- read_best_beta(beta_file)
alpha_ref <- if (nzchar(alpha_ref_file)) read_best_alpha(alpha_ref_file) else c(0, 0)
start_alpha <- draw_alpha_start(seed = seed, alpha_ref = alpha_ref)
start_weights <- weights_from_alpha(start_alpha)

prepared <- make_sim_data_for_round(seed = sim_seed, J = J)
cache <- precompute_grad_hess(beta_hat, prepared$data, prepared$J, prepared$Time)
objective_var_fun <- objective_var_fun_factory(cache)
started_at <- Sys.time()

fit <- try(
  optim(
    par = start_alpha,
    fn = objective_var_fun,
    control = list(maxit = maxit)
  ),
  silent = TRUE
)

message <- ""
convergence <- 99L
objective_value <- Inf
alpha_hat <- rep(NA_real_, 2)
weights_hat <- rep(NA_real_, 3)
baseline_trace <- try(objective_var_fun(start_alpha), silent = TRUE)
optimized_trace <- NA_real_

if (inherits(fit, "try-error")) {
  message <- as.character(fit)
} else {
  convergence <- as.integer(fit$convergence)
  objective_value <- fit$value
  alpha_hat <- fit$par
  weights_hat <- weights_from_alpha(alpha_hat)
  optimized_trace <- try(objective_var_fun(alpha_hat), silent = TRUE)
  if (inherits(optimized_trace, "try-error")) {
    optimized_trace <- NA_real_
  }
  message <- fit$message %||% ""
}

if (inherits(baseline_trace, "try-error")) {
  baseline_trace <- NA_real_
}

elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

row <- data.frame(
  round = round_id,
  task_id = task_id,
  seed = seed,
  sim_seed = sim_seed,
  groups = J,
  maxit = maxit,
  beta_file = beta_file,
  alpha_ref_file = alpha_ref_file,
  convergence = convergence,
  objective_value = objective_value,
  baseline_trace = baseline_trace,
  optimized_trace = optimized_trace,
  elapsed_seconds = elapsed_seconds,
  message = message,
  stringsAsFactors = FALSE
)

for (i in seq_along(start_alpha)) {
  row[[sprintf("start_alpha_%02d", i)]] <- start_alpha[[i]]
}
for (i in seq_along(alpha_hat)) {
  row[[sprintf("alpha_%02d", i)]] <- alpha_hat[[i]]
}
for (i in seq_along(start_weights)) {
  row[[sprintf("start_weight_%02d", i)]] <- start_weights[[i]]
}
for (i in seq_along(weights_hat)) {
  row[[sprintf("weight_%02d", i)]] <- weights_hat[[i]]
}

append_csv_locked(row, out_csv)

if (convergence != 0L || !is.finite(objective_value)) {
  quit(save = "no", status = 1)
}

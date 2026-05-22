source("composite_utils.R")

usage <- paste(
  "Usage:",
  "Rscript stage1_beta_runner.R --groups <int> --task-id <int> --start-family <random|fixed|beta_ref>",
  "--seed <int> --data <path> --out <path> [--weights-file <path>] [--beta-ref-file <path>]",
  sep = "\n"
)

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

require_args(args, c("groups", "task-id", "start-family", "seed", "data", "out"))

group_count <- as_int_arg(args, "groups")
task_id <- as_int_arg(args, "task-id")
seed <- as_int_arg(args, "seed")
start_family <- args[["start-family"]]
data_path <- args[["data"]]
output_path <- args[["out"]]
weights_file <- args[["weights-file"]] %||% ""
beta_ref_file <- args[["beta-ref-file"]] %||% ""

if (!start_family %in% c("random", "fixed", "beta_ref")) {
  stop("Argument --start-family must be 'random', 'fixed', or 'beta_ref'.", call. = FALSE)
}

if (identical(start_family, "beta_ref") && identical(beta_ref_file, "")) {
  stop("Argument --beta-ref-file is required when --start-family beta_ref is used.", call. = FALSE)
}

set.seed(seed)
prepared <- prepare_data(data_path, group_count)
beta_reference <- if (!identical(beta_ref_file, "")) load_beta_reference(beta_ref_file) else NULL
weights <- if (!identical(weights_file, "")) load_weight_reference(weights_file) else c(1, 1, 1)
start <- draw_beta_start(start_family, beta_reference)

initial_value <- try(
  neg_comp_llh(start, weights, prepared$data, prepared$J, prepared$Time),
  silent = TRUE
)

result <- list(
  stage = "stage1",
  task_id = task_id,
  seed = seed,
  start_family = start_family,
  groups = group_count,
  weights = weights,
  weights_file = weights_file,
  beta_ref_file = beta_ref_file,
  start = start,
  convergence = 99L,
  value = Inf,
  par = rep(NA_real_, length(start)),
  initial_value = initial_value,
  message = NULL
)

if (inherits(initial_value, "try-error") || is.na(initial_value) || !is.finite(initial_value)) {
  result$message <- "Initial objective was invalid."
  save_result(output_path, result)
  quit(save = "no", status = 1)
}

fit <- try(
  optim(
    par = start,
    fn = neg_comp_llh,
    weights = weights,
    sim_dat = prepared$data,
    J = prepared$J,
    Time = prepared$Time,
    control = list(maxit = 8000)
  ),
  silent = TRUE
)

if (inherits(fit, "try-error")) {
  result$message <- as.character(fit)
  save_result(output_path, result)
  quit(save = "no", status = 1)
}

result$convergence <- as.integer(fit$convergence)
result$value <- fit$value
result$par <- fit$par
result$message <- fit$message %||% NULL

save_result(output_path, result)

if (as.integer(fit$convergence) != 0L || !is.finite(fit$value)) {
  quit(save = "no", status = 1)
}

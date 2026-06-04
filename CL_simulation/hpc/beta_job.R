source("cl_sim_utils.R")

usage <- paste(
  "Usage:",
  "Rscript beta_job.R --round <int> --task-id <int> --seed <int> --sim-seed <int>",
  "--groups <int> --out-csv <path> --maxit <int>",
  "[--weights-file <path>] [--beta-ref-file <path>]",
  sep = "\n"
)

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

require_args(args, c("round", "task-id", "seed", "sim-seed", "groups", "out-csv", "maxit"))

round_id <- as_int_arg(args, "round")
task_id <- as_int_arg(args, "task-id")
seed <- as_int_arg(args, "seed")
sim_seed <- as_int_arg(args, "sim-seed")
J <- as_int_arg(args, "groups")
out_csv <- args[["out-csv"]]
maxit <- as_int_arg(args, "maxit")
weights_file <- args[["weights-file"]] %||% ""
beta_ref_file <- args[["beta-ref-file"]] %||% ""

weights <- if (nzchar(weights_file)) read_best_weights(weights_file) else c(1 / 3, 1 / 3, 1 / 3)
beta_ref <- if (nzchar(beta_ref_file)) read_best_beta(beta_ref_file) else true_beta
start_beta <- draw_beta_start(seed = seed, beta_ref = beta_ref)

prepared <- make_sim_data_for_round(seed = sim_seed, J = J)
started_at <- Sys.time()

initial_value <- try(
  neg_comp_llh(start_beta, weights, prepared$data, prepared$J, prepared$Time),
  silent = TRUE
)

fit <- NULL
message <- ""
convergence <- 99L
objective_value <- Inf
beta_hat <- rep(NA_real_, length(start_beta))

if (inherits(initial_value, "try-error") || !is.finite(initial_value)) {
  message <- "Initial beta objective was invalid."
} else {
  fit <- try(
    optim(
      par = start_beta,
      fn = neg_comp_llh,
      weights = weights,
      sim_dat = prepared$data,
      J = prepared$J,
      Time = prepared$Time,
      control = list(maxit = maxit)
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    message <- as.character(fit)
  } else {
    convergence <- as.integer(fit$convergence)
    objective_value <- fit$value
    beta_hat <- fit$par
    message <- fit$message %||% ""
  }
}

elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

row <- data.frame(
  round = round_id,
  task_id = task_id,
  seed = seed,
  sim_seed = sim_seed,
  groups = J,
  maxit = maxit,
  weights_file = weights_file,
  beta_ref_file = beta_ref_file,
  convergence = convergence,
  objective_value = objective_value,
  initial_value = if (inherits(initial_value, "try-error")) NA_real_ else initial_value,
  elapsed_seconds = elapsed_seconds,
  message = message,
  stringsAsFactors = FALSE
)

for (i in seq_along(weights)) {
  row[[sprintf("weight_%02d", i)]] <- weights[[i]]
}
for (i in seq_along(start_beta)) {
  row[[sprintf("start_beta_%02d", i)]] <- start_beta[[i]]
}
for (i in seq_along(beta_hat)) {
  row[[sprintf("beta_%02d", i)]] <- beta_hat[[i]]
}

append_csv_locked(row, out_csv)

if (convergence != 0L || !is.finite(objective_value)) {
  quit(save = "no", status = 1)
}

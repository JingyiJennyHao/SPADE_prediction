source("cl_sim_utils.R")

usage <- paste(
  "Usage:",
  "Rscript summarize_beta.R --input-csv <path> --out <path>",
  sep = "\n"
)

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

require_args(args, c("input-csv", "out"))

input_csv <- args[["input-csv"]]
out_path <- args[["out"]]

if (!file.exists(input_csv)) {
  stop(sprintf("Input CSV not found: %s", input_csv), call. = FALSE)
}

results <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
ok <- results$convergence == 0 & is.finite(results$objective_value)
if (!any(ok)) {
  stop(sprintf("No converged finite beta results found in %s", input_csv), call. = FALSE)
}

converged <- results[ok, , drop = FALSE]
best <- converged[which.min(converged$objective_value), , drop = FALSE]
beta_cols <- sprintf("beta_%02d", seq_along(true_beta))
beta_hat <- as.numeric(best[1, beta_cols])
weight_cols <- sprintf("weight_%02d", 1:3)
weights <- as.numeric(best[1, weight_cols])

payload <- list(
  stage = "beta_summary",
  round = best$round[[1]],
  input_csv = input_csv,
  best_task_id = best$task_id[[1]],
  best_seed = best$seed[[1]],
  sim_seed = best$sim_seed[[1]],
  groups = best$groups[[1]],
  maxit = best$maxit[[1]],
  convergence = best$convergence[[1]],
  objective_value = best$objective_value[[1]],
  beta_hat = beta_hat,
  weights = weights,
  best_row = best
)

save_result(out_path, payload)
cat(sprintf("Saved best beta from task %s to %s\n", payload$best_task_id, out_path))

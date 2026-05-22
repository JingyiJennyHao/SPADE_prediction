source("composite_utils.R")

usage <- paste(
  "Usage:",
  "Rscript stage2_summarize.R --input-dir <dir> --out <path>",
  sep = "\n"
)

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

require_args(args, c("input-dir", "out"))

input_dir <- args[["input-dir"]]
output_path <- args[["out"]]

if (!dir.exists(input_dir)) {
  stop(sprintf("Input directory not found: %s", input_dir), call. = FALSE)
}

files <- sort(list.files(input_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0L) {
  stop(sprintf("No RDS files found in %s", input_dir), call. = FALSE)
}

results <- lapply(files, readRDS)
converged_idx <- vapply(
  results,
  function(item) as.integer(item$convergence) == 0L && is.finite(item$objective_value),
  logical(1)
)

if (!any(converged_idx)) {
  stop("No converged stage 2 results were found.", call. = FALSE)
}

converged <- results[converged_idx]
values <- vapply(converged, function(item) item$objective_value, numeric(1))
best_idx <- which.min(values)
best <- converged[[best_idx]]

summary_payload <- list(
  stage = "stage2_summary",
  files = basename(files),
  converged_files = basename(files[converged_idx]),
  best_file = basename(files[converged_idx][best_idx]),
  weights_hat = best$weights,
  objective_value = best$objective_value,
  optimized_trace = best$optimized_trace,
  baseline_trace = best$baseline_trace,
  best_seed = best$seed,
  best_task_id = best$task_id,
  groups = best$groups,
  beta_file = best$beta_file
)

save_result(output_path, summary_payload)

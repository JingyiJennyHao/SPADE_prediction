source("composite_utils.R")

usage <- paste(
  "Usage:",
  "Rscript stage1_summarize.R --input-dir <dir> --out <path>",
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
  function(item) as.integer(item$convergence) == 0L && is.finite(item$value),
  logical(1)
)

if (!any(converged_idx)) {
  stop("No converged stage 1 results were found.", call. = FALSE)
}

converged <- results[converged_idx]
values <- vapply(converged, function(item) item$value, numeric(1))
best_idx <- which.min(values)
best <- converged[[best_idx]]

summary_payload <- list(
  stage = "stage1_summary",
  files = basename(files),
  converged_files = basename(files[converged_idx]),
  best_file = basename(files[converged_idx][best_idx]),
  beta_hat = best$par,
  beta_best_val = best$value,
  best_seed = best$seed,
  best_task_id = best$task_id,
  best_start_family = best$start_family,
  groups = best$groups,
  weights = best$weights,
  weights_file = best$weights_file,
  beta_ref_file = best$beta_ref_file
)

save_result(output_path, summary_payload)

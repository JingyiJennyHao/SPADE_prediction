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
  "Rscript gmm_summarize.R --input-dir <dir> --out <path>",
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
values <- vapply(results, function(item) item$value %||% item$objective_value %||% Inf, numeric(1))
converged_idx <- vapply(
  results,
  function(item) as.integer(item$convergence) == 0L,
  logical(1)
) & is.finite(values)

if (!any(converged_idx)) {
  stop("No converged GMM results were found.", call. = FALSE)
}

converged <- results[converged_idx]
converged_values <- values[converged_idx]
best_idx <- which.min(converged_values)
best <- converged[[best_idx]]
beta_hat <- best$par %||% best$beta_hat

summary_payload <- list(
  stage = "gmm_summary",
  files = basename(files),
  converged_files = basename(files[converged_idx]),
  best_file = basename(files[converged_idx][best_idx]),
  beta_hat = as.numeric(beta_hat),
  beta_best_val = converged_values[[best_idx]],
  objective_value = converged_values[[best_idx]],
  best_seed = best$seed,
  best_task_id = best$task_id,
  groups = best$groups,
  beta0 = best$beta0,
  beta0_file = best$beta0_file,
  start = best$start,
  start_scale = best$start_scale,
  ridge = best$ridge,
  maxit = best$maxit,
  center_v = best$center_v,
  V = best$V,
  W = best$W,
  S0 = best$S0
)

save_result(output_path, summary_payload)
cat(sprintf("Saved best GMM beta from %s to %s\n", summary_payload$best_file, output_path))

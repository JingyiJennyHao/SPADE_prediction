#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

get_arg <- function(name, default = NULL, required = FALSE) {
  key <- paste0("--", name, "=")
  hit <- grep(paste0("^", key), args, value = TRUE)
  if (length(hit)) return(sub(key, "", hit[[length(hit)]], fixed = TRUE))
  pos <- tail(which(args == paste0("--", name)), 1)
  if (length(pos) && pos < length(args)) return(args[[pos + 1]])
  if (required) stop("Missing required argument --", name)
  default
}

file_arg <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "gmm_finalize_seed.R"
)
script_dir <- dirname(normalizePath(file_arg, mustWork = TRUE))
source(file.path(script_dir, "gmm_core.R"))

seed <- as.integer(get_arg("seed", required = TRUE))
loops <- as.integer(get_arg("loops", "6"))
expected_starts <- as.integer(get_arg("starts", "20"))
J <- as.integer(get_arg("J", "200"))
Time <- as.integer(get_arg("Time", "3"))
out_dir <- normalizePath(get_arg("out-dir", required = TRUE), mustWork = FALSE)
run_inference <- as.logical(as.integer(get_arg("run-inference", "1")))
inference_prob <- as.numeric(get_arg("inference-prob", "0.95"))

if (!is.finite(seed) || seed < 1 || !is.finite(loops) || loops < 1 ||
    !is.finite(expected_starts) || expected_starts < 1 || !is.finite(J) || J < 1 ||
    !is.finite(Time) || Time < 1) {
  stop("seed, loops, starts, J, and Time must be positive integers")
}
if (!is.finite(inference_prob) || inference_prob <= 0 || inference_prob >= 1) {
  stop("--inference-prob must be between 0 and 1")
}

final_results_file <- file.path(
  out_dir, sprintf("final_start_results_seed%d.csv", seed)
)
if (!file.exists(final_results_file)) {
  stop("Missing per-start final results: ", final_results_file)
}

results <- read.csv(final_results_file, stringsAsFactors = FALSE)
beta_cols <- paste0("b", seq_along(true_beta))
if (!nrow(results) || !all(c("value", "start_id", beta_cols) %in% names(results))) {
  stop("Per-start final results are empty or missing required columns: ", final_results_file)
}
completed_start_ids <- unique(results$start_id[is.finite(results$start_id)])

usable <- is.finite(results$value) &
  apply(is.finite(as.matrix(results[beta_cols])), 1, all)
if (!any(usable)) {
  stop("No usable beta/objective rows in per-start final results: ", final_results_file)
}

candidates <- results[usable, , drop = FALSE]
best_row <- candidates[which.min(candidates$value), , drop = FALSE]
best_beta <- as.numeric(unlist(best_row[beta_cols], use.names = FALSE))
best_start_id <- best_row$start_id[[1L]]
best_value <- best_row$value[[1L]]

selection <- data.frame(
  seed = seed,
  final_loop = best_row$final_loop[[1L]],
  max_loops = loops,
  expected_starts = expected_starts,
  completed_start_points = length(completed_start_ids),
  completed_rows = nrow(results),
  usable_rows = nrow(candidates),
  best_start_id = best_start_id,
  best_value = best_value,
  selection_rule = "minimum objective value among finite beta rows",
  source_file = final_results_file,
  stringsAsFactors = FALSE
)
if ("outer_converged" %in% names(best_row)) {
  selection$outer_converged <- best_row$outer_converged[[1L]]
}
selection[beta_cols] <- as.list(best_beta)
selection_file <- file.path(out_dir, sprintf("final_selection_seed%d.csv", seed))
write.csv(selection, selection_file, row.names = FALSE)

cat(sprintf(
  "seed=%d final_loop=%d max_loops=%d completed=%d usable=%d best_start_id=%d best_value=%.12g\n",
  seed, best_row$final_loop[[1L]], loops, nrow(results), nrow(candidates),
  best_start_id, best_value
))
cat("Selected beta file: ", selection_file, "\n", sep = "")

if (run_inference) {
  inference_file <- file.path(out_dir, sprintf("inference_seed%d.csv", seed))
  run_gmm_inference(
    beta_hat = best_beta,
    seed = seed,
    out_csv = inference_file,
    J = J,
    Time = Time,
    interval_prob = inference_prob
  )
  cat("Inference file: ", inference_file, "\n", sep = "")
}

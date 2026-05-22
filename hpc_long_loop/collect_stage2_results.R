source("composite_utils.R")

usage <- paste(
  "Usage:",
  "Rscript collect_stage2_results.R --input-dir <dir> --out-csv <path>",
  sep = "\n"
)

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

require_args(args, c("input-dir", "out-csv"))

input_dir <- args[["input-dir"]]
out_csv <- args[["out-csv"]]

if (!dir.exists(input_dir)) {
  stop(sprintf("Input directory not found: %s", input_dir), call. = FALSE)
}

files <- sort(list.files(input_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0L) {
  stop(sprintf("No RDS files found in %s", input_dir), call. = FALSE)
}

flatten_stage2 <- function(item, file_name) {
  weights_vec <- item$weights %||% rep(NA_real_, 3)
  alpha_vec <- item$alpha %||% rep(NA_real_, 2)
  start_weights <- item$start_weights %||% rep(NA_real_, 3)
  start_alpha <- item$start_alpha %||% rep(NA_real_, 2)

  row <- data.frame(
    file = basename(file_name),
    task_id = item$task_id %||% NA_integer_,
    seed = item$seed %||% NA_integer_,
    groups = item$groups %||% NA_integer_,
    convergence = item$convergence %||% NA_integer_,
    objective_value = item$objective_value %||% NA_real_,
    baseline_trace = item$baseline_trace %||% NA_real_,
    optimized_trace = item$optimized_trace %||% NA_real_,
    beta_file = item$beta_file %||% "",
    message = item$message %||% "",
    stringsAsFactors = FALSE
  )

  for (i in seq_along(start_weights)) {
    row[[sprintf("start_weight_%02d", i)]] <- as.numeric(start_weights[[i]])
  }
  for (i in seq_along(start_alpha)) {
    row[[sprintf("start_alpha_%02d", i)]] <- as.numeric(start_alpha[[i]])
  }
  for (i in seq_along(weights_vec)) {
    row[[sprintf("weight_%02d", i)]] <- as.numeric(weights_vec[[i]])
  }
  for (i in seq_along(alpha_vec)) {
    row[[sprintf("alpha_%02d", i)]] <- as.numeric(alpha_vec[[i]])
  }

  row
}

rows <- Map(function(path, item) flatten_stage2(item, path), files, lapply(files, readRDS))
result_df <- do.call(rbind, rows)
result_df <- result_df[order(result_df$convergence, result_df$objective_value, result_df$task_id), ]

fail_fast_dir(out_csv)
utils::write.csv(result_df, out_csv, row.names = FALSE)
cat(sprintf("Wrote %d stage 2 results to %s\n", nrow(result_df), out_csv))

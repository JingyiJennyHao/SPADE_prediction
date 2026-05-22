source("composite_utils.R")

usage <- paste(
  "Usage:",
  "Rscript collect_stage1_results.R --input-dir <dir> --out-csv <path>",
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

flatten_stage1 <- function(item, file_name) {
  par_vec <- item$par %||% rep(NA_real_, 18)
  start_vec <- item$start %||% rep(NA_real_, 18)

  row <- data.frame(
    file = basename(file_name),
    task_id = item$task_id %||% NA_integer_,
    seed = item$seed %||% NA_integer_,
    start_family = item$start_family %||% NA_character_,
    groups = item$groups %||% NA_integer_,
    convergence = item$convergence %||% NA_integer_,
    value = item$value %||% NA_real_,
    initial_value = if (inherits(item$initial_value, "try-error")) NA_real_ else (item$initial_value %||% NA_real_),
    message = item$message %||% "",
    stringsAsFactors = FALSE
  )

  for (i in seq_along(start_vec)) {
    row[[sprintf("start_%02d", i)]] <- as.numeric(start_vec[[i]])
  }
  for (i in seq_along(par_vec)) {
    row[[sprintf("beta_%02d", i)]] <- as.numeric(par_vec[[i]])
  }

  row
}

rows <- Map(function(path, item) flatten_stage1(item, path), files, lapply(files, readRDS))
result_df <- do.call(rbind, rows)
result_df <- result_df[order(result_df$convergence, result_df$value, result_df$task_id), ]

fail_fast_dir(out_csv)
utils::write.csv(result_df, out_csv, row.names = FALSE)
cat(sprintf("Wrote %d stage 1 results to %s\n", nrow(result_df), out_csv))

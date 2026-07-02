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
names(results) <- basename(files)

flatten_gmm <- function(item, file_name) {
  par_vec <- item$par %||% item$beta_hat %||% rep(NA_real_, 18)
  start_vec <- item$start %||% rep(NA_real_, length(par_vec))
  beta0_vec <- item$beta0 %||% rep(NA_real_, length(par_vec))
  objective_value <- item$value %||% item$objective_value %||% NA_real_

  row <- data.frame(
    file = basename(file_name),
    task_id = item$task_id %||% NA_integer_,
    seed = item$seed %||% NA_integer_,
    groups = item$groups %||% NA_integer_,
    convergence = item$convergence %||% NA_integer_,
    value = objective_value,
    initial_value = if (inherits(item$initial_value, "try-error")) NA_real_ else (item$initial_value %||% NA_real_),
    message = item$message %||% "",
    beta0_file = item$beta0_file %||% "",
    stringsAsFactors = FALSE
  )

  for (i in seq_along(beta0_vec)) {
    row[[sprintf("beta0_%02d", i)]] <- as.numeric(beta0

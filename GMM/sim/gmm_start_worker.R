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

file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "gmm_start_worker.R")
script_dir <- dirname(normalizePath(file_arg, mustWork = TRUE))
source(file.path(script_dir, "gmm_core.R"))

seed <- as.integer(get_arg("seed", required = TRUE))
start_id <- as.integer(get_arg("start-id", Sys.getenv("LSB_JOBINDEX", unset = "1")))
loops_arg <- get_arg("loops", NULL)
out_dir_arg <- get_arg("out-dir", NULL)
J <- as.integer(get_arg("J", "100"))
Time <- as.integer(get_arg("Time", "3"))
start_sd <- as.numeric(get_arg("start-sd", "0.1"))
maxit <- as.integer(get_arg("maxit", "8000"))
outer_tol <- as.numeric(get_arg("outer-tol", "1e-5"))

if (!is.finite(seed) || !is.finite(start_id) || start_id < 1) {
  stop("--seed and --start-id must be positive integers")
}
if (!is.finite(J) || J < 1 || !is.finite(Time) || Time < 1 ||
    !is.finite(maxit) || maxit < 1) {
  stop("--J, --Time, and --maxit must be positive")
}
if (!is.finite(outer_tol) || outer_tol <= 0) {
  stop("--outer-tol must be positive")
}

loops <- as.integer(loops_arg %||% "6")
if (!is.finite(loops) || loops < 1 || is.null(out_dir_arg)) {
  stop("All-loop mode requires a positive --loops and --out-dir")
}

out_dir <- normalizePath(out_dir_arg, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

loop_csv_file <- function(loop) {
  file.path(out_dir, sprintf("loop_%02d_start_results_seed%d.csv", loop, seed))
}

base_beta <- true_beta
last_result <- NULL
last_valid_result <- NULL
last_valid_loop <- NA_integer_
loops_attempted <- 0L
outer_converged <- 0L
stop_reason <- "max_loops"
final_csv <- file.path(out_dir, sprintf("final_start_results_seed%d.csv", seed))
cat(sprintf(
  "seed=%d start_id=%d running loops 1-%d outer_tol=%g out_dir=%s\n",
  seed, start_id, loops, outer_tol, out_dir
))

for (loop in seq_len(loops)) {
  loops_attempted <- loop
  out_csv <- loop_csv_file(loop)
  result <- run_start_point(
    seed = seed,
    loop = loop,
    start_id = start_id,
    base_beta = base_beta,
    out_csv = out_csv,
    J = J,
    Time = Time,
    start_sd = start_sd,
    maxit = maxit
  )
  last_result <- result

  beta_cols <- paste0("b", seq_along(true_beta))
  next_beta <- as.numeric(unlist(result[beta_cols], use.names = FALSE))
  if (length(next_beta) == length(true_beta) && all(is.finite(next_beta))) {
    last_valid_result <- result
    last_valid_loop <- loop
    base_beta <- next_beta
    max_abs <- as.numeric(result$diff_from_base)
    cat(sprintf("seed=%d start_id=%d loop=%d complete value=%s max_abs=%s\n",
                seed, start_id, loop, format(result$value, digits = 12),
                format(max_abs, digits = 12)))
    if (loop >= 2L && is.finite(max_abs) && max_abs < outer_tol) {
      outer_converged <- 1L
      stop_reason <- "outer_tol"
      cat(sprintf(
        "seed=%d start_id=%d outer convergence reached at loop=%d: max_abs=%g < outer_tol=%g\n",
        seed, start_id, loop, max_abs, outer_tol
      ))
      break
    }
  } else {
    # Preserve the path and finish the remaining loops even if one optimizer
    # call fails. The next loop then retries from the last usable beta.
    cat(sprintf("seed=%d start_id=%d loop=%d returned no usable beta; retaining previous beta\n",
                seed, start_id, loop), file = stderr())
  }
}

if (is.null(last_valid_result)) {
  final_row <- last_result
  final_loop <- loops_attempted
  stop_reason <- "no_usable_beta"
} else {
  final_row <- last_valid_result
  final_loop <- last_valid_loop
}
final_row$final_loop <- final_loop
final_row$loops_attempted <- loops_attempted
final_row$outer_converged <- outer_converged
final_row$stop_reason <- stop_reason
write_row(final_row, final_csv)
cat(sprintf(
  "seed=%d start_id=%d final_loop=%d loops_attempted=%d outer_converged=%d final_results=%s\n",
  seed, start_id, final_loop, loops_attempted, outer_converged, final_csv
))

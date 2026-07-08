#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

get_arg <- function(name, default = NULL, required = FALSE) {
  key <- paste0("--", name, "=")
  hit <- grep(paste0("^", key), args, value = TRUE)
  if (length(hit)) return(sub(key, "", hit[[1]], fixed = TRUE))
  pos <- match(paste0("--", name), args)
  if (!is.na(pos) && pos < length(args)) return(args[[pos + 1]])
  if (required) stop("Missing required argument --", name)
  default
}

file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "gmm_hpc_driver.R")
script_dir <- dirname(normalizePath(file_arg, mustWork = TRUE))
source(file.path(script_dir, "gmm_core.R"))

seed <- as.integer(get_arg("seed", "1"))
n_loops <- as.integer(get_arg("loops", "5"))
n_starts <- as.integer(get_arg("starts", "120"))
target <- as.integer(get_arg("target", "100"))
poll_seconds <- as.numeric(get_arg("poll-seconds", "60"))
max_wait_hours <- as.numeric(get_arg("max-wait-hours", "72"))
J <- as.integer(get_arg("J", "100"))
Time <- as.integer(get_arg("Time", "3"))
start_sd <- as.numeric(get_arg("start-sd", "0.1"))
maxit <- as.integer(get_arg("maxit", "5000"))
out_dir <- normalizePath(get_arg("out-dir", file.path(script_dir, "hpc_runs")), mustWork = FALSE)
queue <- get_arg("queue", "")
extra_bsub <- get_arg("bsub-extra", "")
stop_remaining <- as.logical(as.integer(get_arg("stop-remaining", "1")))
driver_log <- get_arg("driver-log", "")
tol <- as.numeric(get_arg("tol", "0"))

if (target > n_starts) stop("--target cannot be larger than --starts")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
logs_dir <- file.path(out_dir, "logs")
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

if (!nzchar(driver_log)) {
  driver_log <- file.path(logs_dir, sprintf("driver_seed%d.log", seed))
}

log_msg <- function(...) {
  text <- paste0(...)
  cat(text, "\n", sep = "")
  cat(format(Sys.time(), "%F %T"), " ", text, "\n", sep = "", file = driver_log, append = TRUE)
}

write_beta <- function(beta, file) {
  writeLines(paste(beta, collapse = " "), con = file)
}

parse_bsub_jobid <- function(x) {
  hit <- regmatches(x, regexpr("<[0-9]+>", x))
  if (!length(hit) || hit == "") return(NA_character_)
  gsub("[<>]", "", hit)
}

submit_loop <- function(loop, base_beta_file, loop_csv) {
  job_name <- sprintf("gmm_s%d_l%d[1-%d]", seed, loop, n_starts)
  wrapper <- file.path(script_dir, "gmm_start_worker_lsf.sh")
  out_log <- file.path(logs_dir, sprintf("loop_%02d_seed%d.out", loop, seed))
  err_log <- file.path(logs_dir, sprintf("loop_%02d_seed%d.err", loop, seed))

  cmd <- c(
    "bsub",
    if (nzchar(queue)) c("-q", shQuote(queue)) else character(),
    "-J", shQuote(job_name),
    "-o", shQuote(out_log),
    "-e", shQuote(err_log),
    if (nzchar(extra_bsub)) extra_bsub else character(),
    shQuote(wrapper),
    seed,
    loop,
    shQuote(base_beta_file),
    shQuote(loop_csv),
    J,
    Time,
    start_sd,
    maxit
  )

  cmd <- paste(cmd, collapse = " ")
  log_msg("Submitting loop ", loop, ": ", cmd)
  submit_out <- system(cmd, intern = TRUE)
  log_msg(paste(submit_out, collapse = "\n"))
  parse_bsub_jobid(paste(submit_out, collapse = "\n"))
}

base_beta <- true_beta
base_beta_file <- file.path(out_dir, "base_beta_loop1.txt")
write_beta(base_beta, base_beta_file)

summary_file <- file.path(out_dir, sprintf("loop_summary_seed%d.csv", seed))

for (loop in seq_len(n_loops)) {
  loop_csv <- file.path(out_dir, sprintf("loop_%02d_start_results_seed%d.csv", loop, seed))
  if (file.exists(loop_csv)) {
    file.rename(loop_csv, paste0(loop_csv, ".old_", format(Sys.time(), "%Y%m%d%H%M%S")))
  }

  job_id <- submit_loop(loop, base_beta_file, loop_csv)
  started <- Sys.time()

  repeat {
    done <- complete_rows(loop_csv, require_converged = TRUE)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "hours"))
    log_msg(sprintf("loop=%d converged=%d/%d elapsed=%.2fh",
                    loop, done, target, elapsed))

    if (done >= target) break
    if (elapsed > max_wait_hours) {
      stop("Timed out waiting for loop ", loop, " to reach ", target, " converged starts.")
    }
    Sys.sleep(poll_seconds)
  }

  if (stop_remaining && !is.na(job_id)) {
    kill_out <- system(paste("bkill", shQuote(job_id)), intern = TRUE)
    log_msg(paste(kill_out, collapse = "\n"))
  }

  best_beta <- best_beta_from_loop_file(loop_csv, require_converged = TRUE)
  previous_beta <- base_beta
  base_beta <- best_beta
  beta_diff <- max(abs(base_beta - previous_beta))
  next_base_file <- file.path(out_dir, sprintf("base_beta_loop%d.txt", loop + 1))
  write_beta(base_beta, next_base_file)

  loop_rows <- read.csv(loop_csv, stringsAsFactors = FALSE)
  converged_rows <- loop_rows[loop_rows$convergence == 0 & is.finite(loop_rows$value), , drop = FALSE]
  best_value <- min(converged_rows$value)
  best_start_id <- converged_rows$start_id[which.min(converged_rows$value)]

  row <- data.frame(
    seed = seed,
    loop = loop,
    job_id = job_id,
    requested_starts = n_starts,
    target_converged = target,
    completed_converged = nrow(converged_rows),
    best_start_id = best_start_id,
    best_value = best_value,
    diff_from_previous_base = beta_diff,
    converged_by_tol = as.integer(tol > 0 && beta_diff <= tol),
    loop_csv = loop_csv
  )
  row[paste0("b", seq_along(base_beta))] <- as.list(base_beta)
  write_row(row, summary_file)

  base_beta_file <- next_base_file
  log_msg(sprintf("Loop %d complete: best_start=%s best_value=%.8g beta_diff=%.8g",
                  loop, best_start_id, best_value, beta_diff))

  if (tol > 0 && beta_diff <= tol) {
    log_msg(sprintf("Stopping early: beta_diff %.8g <= tol %.8g", beta_diff, tol))
    break
  }
}

log_msg("All requested loops complete.")
log_msg("Summary: ", summary_file)
log_msg("Final beta: ", base_beta_file)

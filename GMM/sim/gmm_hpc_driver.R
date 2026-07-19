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
n_loops <- as.integer(get_arg("loops", "20"))
n_starts <- as.integer(get_arg("starts", "120"))
target <- as.integer(get_arg("target", "70"))
start_loop <- as.integer(get_arg("start-loop", "1"))
poll_seconds <- as.numeric(get_arg("poll-seconds", "60"))
max_wait_hours <- as.numeric(get_arg("max-wait-hours", "72"))
J <- as.integer(get_arg("J", "100"))
Time <- as.integer(get_arg("Time", "3"))
start_sd <- as.numeric(get_arg("start-sd", "0.1"))
maxit <- as.integer(get_arg("maxit", "5000"))
out_dir <- normalizePath(get_arg("out-dir", file.path(script_dir, "hpc_runs")), mustWork = FALSE)
queue <- get_arg("queue", "serial")
wall_time <- get_arg("wall-time", "14:00")
memory_gb <- get_arg("memory-gb", "8")
r_module <- get_arg("r-module", "R/4.4.0")
r_libs_user <- get_arg("r-libs-user", file.path(Sys.getenv("HOME"), "R/x86_64-pc-linux-gnu-library/4.4"))
stop_remaining <- as.logical(as.integer(get_arg("stop-remaining", "1")))
driver_log <- get_arg("driver-log", "")
tol <- as.numeric(get_arg("tol", "0"))
require_converged_best <- as.logical(as.integer(get_arg("require-converged-best", "1")))

if (target > n_starts) stop("--target cannot be larger than --starts")
if (start_loop < 1 || start_loop > n_loops) stop("--start-loop must be between 1 and --loops")

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

job_active <- function(job_id) {
  if (is.na(job_id) || !nzchar(job_id)) return(FALSE)
  out <- suppressWarnings(system2("bjobs", c("-noheader", job_id), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  is.null(status) || identical(status, 0L)
}

loop_csv_file <- function(loop) {
  file.path(out_dir, sprintf("loop_%02d_start_results_seed%d.csv", loop, seed))
}

base_beta_file_for_loop <- function(loop) {
  file.path(out_dir, sprintf("base_beta_loop%d.txt", loop))
}

loop_counts <- function(loop_csv) {
  if (!file.exists(loop_csv)) return(list(total = 0L, converged = 0L))
  x <- tryCatch(read.csv(loop_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(list(total = 0L, converged = 0L))
  list(
    total = nrow(x),
    converged = sum(x$convergence == 0 & is.finite(x$value), na.rm = TRUE)
  )
}

summarize_loop <- function(loop, loop_csv, previous_beta, job_id = NA_character_,
                           reason = "completed") {
  counts <- loop_counts(loop_csv)
  best_beta <- best_beta_from_loop_file(loop_csv, require_converged = require_converged_best)
  beta_diff <- max(abs(best_beta - previous_beta))

  loop_rows <- read.csv(loop_csv, stringsAsFactors = FALSE)
  candidate_rows <- loop_rows[is.finite(loop_rows$value), , drop = FALSE]
  if (require_converged_best) {
    candidate_rows <- candidate_rows[candidate_rows$convergence == 0, , drop = FALSE]
  }
  best_value <- min(candidate_rows$value)
  best_start_id <- candidate_rows$start_id[which.min(candidate_rows$value)]

  next_base_file <- base_beta_file_for_loop(loop + 1)
  write_beta(best_beta, next_base_file)

  row <- data.frame(
    seed = seed,
    loop = loop,
    job_id = job_id,
    requested_starts = n_starts,
    target_results = target,
    completed_rows = counts$total,
    converged_rows = counts$converged,
    require_converged_best = as.integer(require_converged_best),
    best_start_id = best_start_id,
    best_value = best_value,
    diff_from_previous_base = beta_diff,
    converged_by_tol = as.integer(tol > 0 && beta_diff <= tol),
    stop_reason = reason,
    loop_csv = loop_csv
  )
  row[paste0("b", seq_along(best_beta))] <- as.list(best_beta)
  write_row(row, summary_file)

  list(beta = best_beta, beta_diff = beta_diff, next_base_file = next_base_file,
       counts = counts, best_value = best_value, best_start_id = best_start_id)
}

submit_loop <- function(loop, base_beta_file, loop_csv) {
  job_name <- sprintf("gmm_s%d_l%d[1-%d]", seed, loop, n_starts)
  out_log <- file.path(logs_dir, sprintf("loop_%02d_seed%d.out", loop, seed))
  err_log <- file.path(logs_dir, sprintf("loop_%02d_seed%d.err", loop, seed))
  queue_directive <- if (nzchar(queue)) sprintf("#BSUB -q %s", queue) else ""

  job_script <- paste(
    "#!/bin/bash",
    sprintf("#BSUB -J %s", job_name),
    sprintf("#BSUB -W %s", wall_time),
    "#BSUB -n 1",
    sprintf('#BSUB -R "rusage[mem=%sGB]"', memory_gb),
    sprintf("#BSUB -o %s", out_log),
    sprintf("#BSUB -e %s", err_log),
    queue_directive,
    "set -euo pipefail",
    sprintf("cd %s", shQuote(script_dir)),
    sprintf("module load %s", shQuote(r_module)),
    sprintf("export R_LIBS_USER=%s", shQuote(r_libs_user)),
    "START_ID=${LSB_JOBINDEX:-1}",
    "if command -v Rscript >/dev/null 2>&1; then",
    sprintf("  Rscript %s \\", shQuote(file.path(script_dir, "gmm_start_worker.R"))),
    sprintf("    --seed %s \\", seed),
    sprintf("    --loop %s \\", loop),
    "    --start-id \"${START_ID}\" \\",
    sprintf("    --base-beta-file %s \\", shQuote(base_beta_file)),
    sprintf("    --out-csv %s \\", shQuote(loop_csv)),
    sprintf("    --J %s \\", J),
    sprintf("    --Time %s \\", Time),
    sprintf("    --start-sd %s \\", start_sd),
    sprintf("    --maxit %s", maxit),
    "else",
    sprintf("  R --vanilla --slave --file=%s --args \\", shQuote(file.path(script_dir, "gmm_start_worker.R"))),
    sprintf("    --seed %s \\", seed),
    sprintf("    --loop %s \\", loop),
    "    --start-id \"${START_ID}\" \\",
    sprintf("    --base-beta-file %s \\", shQuote(base_beta_file)),
    sprintf("    --out-csv %s \\", shQuote(loop_csv)),
    sprintf("    --J %s \\", J),
    sprintf("    --Time %s \\", Time),
    sprintf("    --start-sd %s \\", start_sd),
    sprintf("    --maxit %s", maxit),
    "fi",
    sep = "\n"
  )

  log_msg("Submitting loop ", loop, " with LSF script:\n", job_script)
  submit_out <- system2("bsub", input = job_script, stdout = TRUE, stderr = TRUE)
  log_msg(paste(submit_out, collapse = "\n"))
  parse_bsub_jobid(paste(submit_out, collapse = "\n"))
}

summary_file <- file.path(out_dir, sprintf("loop_summary_seed%d.csv", seed))

base_beta <- true_beta
base_beta_file <- base_beta_file_for_loop(1)
if (!file.exists(base_beta_file)) {
  write_beta(base_beta, base_beta_file)
}

if (start_loop > 1) {
  for (prior_loop in seq_len(start_loop - 1)) {
    prior_csv <- loop_csv_file(prior_loop)
    next_base <- base_beta_file_for_loop(prior_loop + 1)
    counts <- loop_counts(prior_csv)
    log_msg(sprintf(
      "Existing loop %d check: completed=%d, converged=%d, next_base_exists=%s",
      prior_loop, counts$total, counts$converged, file.exists(next_base)
    ))

    if (!file.exists(next_base)) {
      if (counts$converged == 0) {
        stop("Cannot start from loop ", start_loop,
             ": missing ", next_base, " and no converged rows in ", prior_csv)
      }
      previous_beta <- if (prior_loop == 1) true_beta else scan(base_beta_file_for_loop(prior_loop), quiet = TRUE)
      summary <- summarize_loop(prior_loop, prior_csv, previous_beta, reason = "resume_summary")
      log_msg(sprintf(
        "Loop %d already got %d completed rows and %d converged rows; wrote %s.",
        prior_loop, summary$counts$total, summary$counts$converged, next_base
      ))
    } else {
      log_msg(sprintf("Loop %d already has %s; using it for resume.", prior_loop, next_base))
    }
  }
}

base_beta_file <- base_beta_file_for_loop(start_loop)
if (!file.exists(base_beta_file)) {
  stop("Missing base beta file for start loop: ", base_beta_file)
}
base_beta <- scan(base_beta_file, quiet = TRUE)
if (length(base_beta) != length(true_beta)) {
  stop("Expected ", length(true_beta), " beta values in ", base_beta_file)
}

log_msg(sprintf(
  "Starting driver at loop %d/%d with target=%d completed rows, starts=%d.",
  start_loop, n_loops, target, n_starts
))

for (loop in start_loop:n_loops) {
  loop_csv <- loop_csv_file(loop)
  counts <- loop_counts(loop_csv)

  if (counts$total > 0L) {
    log_msg(sprintf(
      "Loop %d existing result file found: completed=%d, converged=%d.",
      loop, counts$total, counts$converged
    ))
  }

  job_id <- NA_character_
  stop_reason <- "existing_results"
  if (counts$total == 0L) {
    job_id <- submit_loop(loop, base_beta_file, loop_csv)
    stop_reason <- "target_results"
    started <- Sys.time()

    repeat {
      counts <- loop_counts(loop_csv)
      elapsed <- as.numeric(difftime(Sys.time(), started, units = "hours"))
      active <- job_active(job_id)
      log_msg(sprintf(
        "loop=%d completed=%d/%d converged=%d active=%s elapsed=%.2fh",
        loop, counts$total, target, counts$converged, active, elapsed
      ))

      if (counts$total >= target) break
      if (!active) {
        stop_reason <- "job_finished_before_target"
        log_msg(sprintf(
          "Loop %d job %s is no longer active with completed=%d/%d; proceeding.",
          loop, job_id, counts$total, target
        ))
        break
      }
      if (elapsed > max_wait_hours) {
        stop("Timed out waiting for loop ", loop, " to reach ", target,
             " completed rows or finish all jobs.")
      }
      Sys.sleep(poll_seconds)
    }
  } else {
    if (counts$total >= target) {
      log_msg(sprintf(
        "Loop %d already got %d completed results, at least target %d; start from next loop after summarizing.",
        loop, counts$total, target
      ))
    } else {
      stop_reason <- "existing_results_below_target"
      log_msg(sprintf(
        "Loop %d has %d existing completed results below target %d and no active job is tracked; treating as no jobs left and proceeding.",
        loop, counts$total, target
      ))
    }
  }

  if (stop_remaining && !is.na(job_id)) {
    kill_out <- system(paste("bkill", shQuote(job_id)), intern = TRUE)
    log_msg(paste(kill_out, collapse = "\n"))
  }

  counts <- loop_counts(loop_csv)
  finite_rows <- if (file.exists(loop_csv)) {
    loop_rows <- read.csv(loop_csv, stringsAsFactors = FALSE)
    sum(is.finite(loop_rows$value), na.rm = TRUE)
  } else {
    0L
  }
  if (require_converged_best && counts$converged == 0L) {
    stop("Loop ", loop, " has no converged rows to choose next beta from: ", loop_csv)
  }
  if (!require_converged_best && finite_rows == 0L) {
    stop("Loop ", loop, " has no finite objective rows to choose next beta from: ", loop_csv)
  }

  previous_beta <- base_beta
  summary <- summarize_loop(loop, loop_csv, previous_beta, job_id = job_id, reason = stop_reason)
  base_beta <- summary$beta
  base_beta_file <- summary$next_base_file
  log_msg(sprintf("Loop %d complete: best_start=%s best_value=%.8g beta_diff=%.8g",
                  loop, summary$best_start_id, summary$best_value, summary$beta_diff))

  if (tol > 0 && summary$beta_diff <= tol) {
    log_msg(sprintf("Stopping early: beta_diff %.8g <= tol %.8g", summary$beta_diff, tol))
    break
  }
}

log_msg("All requested loops complete.")
log_msg("Summary: ", summary_file)
log_msg("Final beta: ", base_beta_file)

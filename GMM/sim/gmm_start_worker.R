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

file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "gmm_start_worker.R")
script_dir <- dirname(normalizePath(file_arg, mustWork = TRUE))
source(file.path(script_dir, "gmm_core.R"))

seed <- as.integer(get_arg("seed", required = TRUE))
loop <- as.integer(get_arg("loop", required = TRUE))
start_id <- as.integer(get_arg("start-id", Sys.getenv("LSB_JOBINDEX", unset = "1")))
base_beta_file <- get_arg("base-beta-file", required = TRUE)
out_csv <- get_arg("out-csv", required = TRUE)
J <- as.integer(get_arg("J", "100"))
Time <- as.integer(get_arg("Time", "3"))
start_sd <- as.numeric(get_arg("start-sd", "0.1"))
maxit <- as.integer(get_arg("maxit", "5000"))

base_beta <- scan(base_beta_file, quiet = TRUE)
if (length(base_beta) != 14) {
  stop("Expected 14 beta values in ", base_beta_file, ", found ", length(base_beta))
}

cat(sprintf(
  "seed=%d loop=%d start_id=%d out_csv=%s\n",
  seed, loop, start_id, out_csv
))

run_start_point(
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

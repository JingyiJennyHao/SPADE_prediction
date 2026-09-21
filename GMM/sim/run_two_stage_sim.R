#!/usr/bin/env Rscript
# One replication per invocation. Source only definitions from the existing Rmd.
args <- commandArgs(trailingOnly = TRUE)
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
if (length(args) != 1L) stop("Usage: Rscript run_two_stage_sim.R /path/to/config.R")
config_file <- normalizePath(args[1], mustWork = TRUE)
config_env <- new.env(parent = baseenv())
sys.source(config_file, config_env)
config <- config_env$config
if (!is.list(config) || is.null(config$seed) || (is.null(config$out_dir) && is.null(config$out_csv)))
  stop("Config must define a list named config with seed and out_dir (or legacy out_csv).")
# Resolve output paths relative to the launch directory, before changing cwd.
absolute_output <- function(path) {
  path <- path.expand(path)
  if (!startsWith(path, "/")) path <- file.path(getwd(), path)
  path
}
if (!is.null(config$out_csv)) config$out_csv <- absolute_output(config$out_csv)
if (is.null(config$out_dir)) config$out_dir <- dirname(config$out_csv)
config$out_dir <- absolute_output(config$out_dir)
if (!is.null(config$diagnostics_file))
  config$diagnostics_file <- absolute_output(config$diagnostics_file)
setwd(script_dir)
source("load_sim_definitions.R")
do.call(run_one_sim, config)

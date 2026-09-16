#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: gmm_stage2_worker.R input.rds candidate_rank output_dir")
input_file <- normalizePath(args[1], mustWork = TRUE)
out_dir <- normalizePath(args[3], mustWork = TRUE)
id <- as.integer(args[2])
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(normalizePath(file_arg)))
source("load_sim_definitions.R")
input <- readRDS(input_file)
stopifnot(length(id) == 1L, !is.na(id), id >= 1L, id <= length(input$candidates))
result <- gmm_refine_candidate(input$candidates[[id]], input$sim_dat,
  input$cb_max_iter, input$cb_tol, input$c_interval, input$cb_beta_maxit)
output <- file.path(out_dir, sprintf("candidate_%04d.rds", id))
tmp <- paste0(output, ".tmp")
saveRDS(result, tmp)
if (!file.rename(tmp, output)) stop("Failed to publish result ", output)
if (!is.null(result$error)) {
  cat(sprintf("candidate=%d start=%d ERROR: %s\n", id, result$start_id, result$error))
} else {
  cat(sprintf("candidate=%d start=%d iterations=%d converged=%s objective=%.12g\n",
              id, result$start_id, result$iterations, result$converged, result$objective))
}

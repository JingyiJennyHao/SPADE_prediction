#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: gmm_stage1_worker.R input.rds start_id output_dir")
input_file <- normalizePath(args[1], mustWork = TRUE)
out_dir <- normalizePath(args[3], mustWork = TRUE)
id <- as.integer(args[2])
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(normalizePath(file_arg)))
source("load_sim_definitions.R")
input <- readRDS(input_file)
stopifnot(length(id) == 1L, !is.na(id), id >= 1L, id <= length(input$initials))
result <- gmm_fit_start(id, input$initials[[id]], input$sim_dat,
                        input$true_beta, input$maxit, input$loops, input$outer_tol)
output <- file.path(out_dir, sprintf("start_%04d.rds", id))
tmp <- paste0(output, ".tmp")
saveRDS(result, tmp)
if (!file.rename(tmp, output)) stop("Failed to publish result ", output)
cat(sprintf("start=%d iterations=%d path_converged=%s stop_reason=%s optim_convergence=%s objective=%.12g\n",
  id, result$iterations, result$path_converged, result$stop_reason,
  result$fit$convergence, result$fit$value))
if (!is.null(result$error)) cat("ERROR:", result$error, "\n")

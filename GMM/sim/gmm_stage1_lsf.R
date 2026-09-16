# One independent beta -> W(beta) -> optim(beta; W) path per start.
gmm_fit_start <- function(id, initial, sim_dat, true_beta, maxit, loops, outer_tol) {
  beta <- initial
  history <- list()
  path_converged <- FALSE
  failure <- NULL
  fit <- list(par = initial, value = Inf, convergence = NA_integer_)
  for (iteration in seq_len(loops)) {
    step <- tryCatch({
      W <- W_hat_fun(beta, sim_dat, Time = 3)
      opt <- optim(beta, gmm_obj, sim_dat = sim_dat, W_hat = W,
                   Time = 3, control = list(maxit = maxit))
      if (any(!is.finite(opt$par)) || !is.finite(opt$value))
        stop("Nonfinite optimizer result")
      list(W = W, fit = opt)
    }, error = function(e) list(error = conditionMessage(e)))
    if (!is.null(step$error)) {
      failure <- step$error
      break
    }
    fit <- step$fit
    delta <- max(abs(fit$par - beta))
    history[[iteration]] <- list(iteration = iteration, beta_before = beta,
      W_used = step$W, beta = fit$par, objective = fit$value,
      beta_change = delta, optim_convergence = fit$convergence,
      optim_message = fit$message)
    beta <- fit$par
    path_converged <- delta < outer_tol
    cat(sprintf("start=%d iteration=%d objective=%.12g beta_change=%g optim_convergence=%d path_converged=%s\n",
      id, iteration, fit$value, delta, fit$convergence, path_converged))
    flush.console()
    if (path_converged) break
  }
  list(start_id = id, initial = initial,
    initial_distance = sqrt(sum((initial - true_beta)^2)), fit = fit,
    path_converged = path_converged, iterations = length(history),
    stop_reason = if (!is.null(failure)) "error" else if (path_converged) "outer_tol" else "max_loops",
    error = failure, history = history)
}

gmm_stage1_lsf <- function(initials, sim_dat, true_beta, maxit, loops,
                           outer_tol, diagnostics_file) {
  results <- gmm_lsf_array(
    payload = list(initials = initials, sim_dat = sim_dat,
      true_beta = true_beta, maxit = maxit, loops = loops, outer_tol = outer_tol),
    n_jobs = length(initials), worker = "gmm_stage1_worker.R",
    label = "starts", output_prefix = "start", loop = 1L,
    diagnostics_file = diagnostics_file,
    validate_result = function(result, id) {
      stopifnot(is.list(result), identical(as.integer(result$start_id), as.integer(id)),
        is.logical(result$path_converged), length(result$path_converged) == 1L,
        !is.na(result$path_converged), is.list(result$fit),
        is.numeric(result$fit$par), length(result$fit$par) == length(true_beta),
        is.numeric(result$fit$value), length(result$fit$value) == 1L,
        identical(result$initial, initials[[id]]))
    },
    failed_result = function(id, reason) {
      # Keep a placeholder at the original index so start IDs never shift.
      list(start_id = id, initial = initials[[id]],
        initial_distance = sqrt(sum((initials[[id]] - true_beta)^2)),
        fit = list(par = rep(NA_real_, length(true_beta)), value = Inf,
                   convergence = NA_integer_),
        path_converged = FALSE, iterations = NA_integer_,
        stop_reason = "unavailable_result", error = reason, history = list())
    })
  stopifnot(identical(vapply(results, function(s) as.integer(s$start_id), integer(1)),
                      seq_along(initials)))
  results
}

# Shared submission and synchronization for both independent candidate arrays.
gmm_lsf_array <- function(payload, n_jobs, worker, label, output_prefix,
                          loop, diagnostics_file, validate_result = NULL,
                          failed_result = NULL) {
  if (is.null(diagnostics_file)) stop("LSF backend requires diagnostics_file")
  root <- paste0(diagnostics_file, ".", label)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  # A unique shared-filesystem directory prevents stale results being reused.
  run_dir <- tempfile(sprintf("loop_%02d_", loop), tmpdir = normalizePath(root))
  dir.create(run_dir)
  input <- file.path(run_dir, "input.rds")
  saveRDS(payload, input)
  script_dir <- normalizePath(getwd())
  module <- Sys.getenv("GMM_TWO_STAGE_R_MODULE", "R/4.4.0")
  libraries <- Sys.getenv("GMM_TWO_STAGE_R_LIBS", Sys.getenv("R_LIBS_USER"))
  queue <- Sys.getenv("GMM_TWO_STAGE_QUEUE", "serial")
  wall <- Sys.getenv("GMM_TWO_STAGE_WALL", "72:00")
  memory <- Sys.getenv("GMM_TWO_STAGE_MEMORY", "8")
  job_script <- file.path(run_dir, "worker.sh")
  writeLines(c("#!/usr/bin/env bash", "set -euo pipefail",
    paste("module load", shQuote(module)),
    paste0("export R_LIBS_USER=", shQuote(libraries)),
    paste("Rscript", shQuote(file.path(script_dir, worker)),
          shQuote(input), '"${LSB_JOBINDEX}"', shQuote(run_dir))), job_script)
  submit <- system2("bsub", c("-J", shQuote(sprintf("gmm_%s_l%d[1-%d]", label, loop, n_jobs)),
    "-q", shQuote(queue), "-W", shQuote(wall), "-n", "1",
    "-R", shQuote(sprintf("rusage[mem=%sGB]", memory)),
    "-o", shQuote(file.path(run_dir, paste0(output_prefix, "_%I.out"))),
    "-e", shQuote(file.path(run_dir, paste0(output_prefix, "_%I.err")))),
    stdin = job_script, stdout = TRUE, stderr = TRUE)
  writeLines(submit, file.path(run_dir, "submission.txt"))
  if (!is.null(attr(submit, "status"))) stop(label, " array submission failed: ", paste(submit, collapse = "\n"))
  hit <- regmatches(submit, regexec("Job <([0-9]+)>", submit))
  ids <- vapply(hit[lengths(hit) > 1L], `[`, character(1), 2L)
  if (length(ids) != 1L) stop("Cannot identify submitted array; inspect ", run_dir)
  cat(sprintf("%s array job %s submitted (%d jobs); files: %s\n", label, ids, n_jobs, run_dir))
  flush.console()
  # Wait on a dependent barrier, rather than polling the scheduler. ended()
  # releases after all array elements terminate, including failed jobs.
  barrier <- file.path(run_dir, "barrier.sh")
  writeLines(c("#!/usr/bin/env bash", "exit 0"), barrier)
  status <- system2("bsub", c("-K", "-J", shQuote(paste0("gmm_join_", ids)),
    "-w", shQuote(sprintf("ended(%s)", ids)), "-q", shQuote(queue),
    "-W", "00:05", "-n", "1", "-R", shQuote("rusage[mem=1GB]"),
    "-o", shQuote(file.path(run_dir, "barrier.out")),
    "-e", shQuote(file.path(run_dir, "barrier.err"))), stdin = barrier,
    stdout = file.path(run_dir, "barrier_submission.txt"), stderr = "")
  if (status != 0L) stop(label, " barrier failed; inspect ", run_dir)
  outputs <- file.path(run_dir, sprintf(paste0(output_prefix, "_%04d.rds"), seq_len(n_jobs)))
  issues <- list()
  results <- lapply(seq_len(n_jobs), function(id) {
    tryCatch({
      if (!file.exists(outputs[id])) stop("Missing result file")
      result <- readRDS(outputs[id])
      if (!is.null(validate_result)) validate_result(result, id)
      result
    }, error = function(e) {
      reason <- conditionMessage(e)
      issues[[length(issues) + 1L]] <<- data.frame(
        index = id, result_file = outputs[id], reason = reason,
        stringsAsFactors = FALSE)
      if (is.null(failed_result)) return(NULL)
      cat(sprintf("%s job index=%d excluded: %s\n", label, id, reason))
      failed_result(id, reason)
    })
  })
  if (length(issues)) {
    issue_table <- do.call(rbind, issues)
    write.csv(issue_table, file.path(run_dir, "unavailable_results.csv"), row.names = FALSE)
    if (is.null(failed_result)) {
      stop(label, " jobs ended with missing/invalid results for indices: ",
           paste(issue_table$index, collapse = ", "), ". Inspect ", run_dir)
    }
    cat(sprintf("%s: %d/%d result files unavailable; continuing with remaining results\n",
                label, nrow(issue_table), n_jobs))
  }
  results
}

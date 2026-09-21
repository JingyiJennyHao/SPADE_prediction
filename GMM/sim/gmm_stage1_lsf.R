# One independent beta -> W(beta) -> optim(beta; W) path per start.
gmm_fit_start <- function(id, initial, sim_dat, true_beta, maxit, loops, outer_tol,
                          candidates = list(initial)) {
  attempts <- list()
  accepted_attempt <- NA_integer_
  initial_W <- NULL
  for (attempt in seq_along(candidates)) {
    candidate <- candidates[[attempt]]
    check <- tryCatch({
      stopifnot(is.numeric(candidate), length(candidate) == length(true_beta),
                all(is.finite(candidate)))
      for (j in sort(unique(sim_dat[, "j"]))) {
        scores <- score_stack(j, candidate, sim_dat, Time = 3)
        if (length(scores) != 3L * length(true_beta) || any(!is.finite(scores)))
          stop("Nonfinite or incomplete initial score")
      }
      W <- diag(3L * length(true_beta))
      value <- gmm_obj(candidate, sim_dat, W_hat = W, Time = 3)
      if (length(value) != 1L || !is.finite(value)) stop("Nonfinite initial objective")
      list(W = W, objective = value)
    }, error = function(e) list(error = conditionMessage(e)))
    accepted <- is.null(check$error)
    attempts[[attempt]] <- list(attempt = attempt, beta = candidate,
      accepted = accepted, error = check$error,
      objective = if (accepted) check$objective else NA_real_)
    cat(sprintf("start=%d initial_attempt=%d/%d %s%s\n", id, attempt,
      length(candidates), if (accepted) "accepted" else "rejected: ",
      if (accepted) "" else check$error))
    flush.console()
    if (accepted) {
      accepted_attempt <- attempt
      initial <- candidate
      initial_W <- check$W
      break
    }
  }
  if (is.na(accepted_attempt)) {
    return(list(start_id = id, initial = rep(NA_real_, length(true_beta)),
      initial_distance = NA_real_, accepted_attempt = NA_integer_,
      initial_attempts = attempts,
      fit = list(par = rep(NA_real_, length(true_beta)), value = Inf, convergence = NA_integer_),
      path_converged = FALSE, iterations = 0L, stop_reason = "initial_checks_exhausted",
      error = sprintf("No valid initial point in %d attempts", length(attempts)), history = list()))
  }
  beta <- initial
  history <- list()
  path_converged <- FALSE
  failure <- NULL
  fit <- list(par = initial, value = Inf, convergence = NA_integer_)
  for (iteration in seq_len(loops)) {
    step <- tryCatch({
      W <- if (iteration == 1L) initial_W else W_hat_fun(beta, sim_dat, Time = 3)
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
    path_converged <- iteration >= 2L && delta < outer_tol
    cat(sprintf("start=%d iteration=%d objective=%.12g beta_change=%g optim_convergence=%d path_converged=%s\n",
      id, iteration, fit$value, delta, fit$convergence, path_converged))
    flush.console()
    if (path_converged) break
  }
  list(start_id = id, initial = initial,
    accepted_attempt = accepted_attempt, initial_attempts = attempts,
    initial_distance = sqrt(sum((initial - true_beta)^2)), fit = fit,
    path_converged = path_converged, iterations = length(history),
    stop_reason = if (!is.null(failure)) "error" else if (path_converged) "outer_tol" else "max_loops",
    error = failure, history = history)
}

gmm_stage1_lsf <- function(initials, sim_dat, true_beta, maxit, loops,
                           outer_tol, diagnostics_file, candidate_pools, start_settings) {
  results <- gmm_lsf_array(
    payload = list(initials = initials, candidate_pools = candidate_pools,
      start_settings = start_settings, sim_dat = sim_dat,
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
        is.list(result$initial_attempts), length(result$accepted_attempt) == 1L)
      a <- result$accepted_attempt
      attempts <- result$initial_attempts
      stopifnot(length(attempts) >= 1L, length(attempts) <= length(candidate_pools[[id]]))
      for (k in seq_along(attempts)) {
        stopifnot(identical(attempts[[k]]$beta, candidate_pools[[id]][[k]]),
                  identical(attempts[[k]]$attempt, k))
      }
      if (is.na(a)) {
        stopifnot(identical(result$stop_reason, "initial_checks_exhausted"),
          !result$path_converged, all(is.na(result$initial)),
          length(attempts) == length(candidate_pools[[id]]),
          all(vapply(attempts, function(x) identical(x$accepted, FALSE), logical(1))))
      } else {
        stopifnot(is.numeric(a), is.finite(a), a == as.integer(a), a >= 1L,
          a == length(attempts), identical(result$initial, candidate_pools[[id]][[a]]),
          isTRUE(attempts[[a]]$accepted))
        if (a > 1L) stopifnot(all(vapply(attempts[seq_len(a - 1L)],
          function(x) identical(x$accepted, FALSE), logical(1))))
      }
    },
    failed_result = function(id, reason) {
      # Keep a placeholder at the original index so start IDs never shift.
      list(start_id = id, initial = rep(NA_real_, length(true_beta)),
        initial_distance = NA_real_, accepted_attempt = NA_integer_, initial_attempts = list(),
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
  root <- file.path(dirname(diagnostics_file), label)
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

source("gmm_stage1_lsf.R")

# Stage 1 iterates W separately per start; Stage 2 recomputes W at every
# effective-beta evaluation, matching real_data_GMM.Rmd.
gmm_default_beta <- c(2, .5, .5, 1, 1, 1, .5, .3, .3, .3, .2, .2, .2, 6)

gmm_refine <- function(beta_init, objective, max_iter = 20L, tol = 1e-5,
                       c_interval = c(-3, 3), beta_maxit = 5000L) {
  v_c <- numeric(14)
  v_c[c(1, 4, 7, 8, 11, 14)] <- 1
  beta <- beta_init
  cc <- 0
  history <- list()
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    opt_c <- optimize(function(c) objective(beta + c * v_c), c_interval)
    c_new <- opt_c$minimum
    opt_beta <- optim(beta, function(b) objective(b + c_new * v_c),
                      method = "Nelder-Mead",
                      control = list(maxit = beta_maxit, reltol = 1e-8))
    beta_new <- opt_beta$par
    q_new <- objective(beta_new + c_new * v_c)
    beta_change <- max(abs(beta_new - beta))
    c_change <- abs(c_new - cc)
    q_change <- if (iter == 1L) NA_real_ else
      abs(q_new - history[[iter - 1L]]$Q_after_beta)
    history[[iter]] <- list(iteration = iter, c = c_new,
      Q_after_c = opt_c$objective, Q_after_beta = q_new, beta = beta_new,
      beta_effective = beta_new + c_new * v_c,
      beta_change = beta_change, c_change = c_change, Q_change = q_change,
      optim_convergence = opt_beta$convergence, message = opt_beta$message)
    converged <- iter > 1L && beta_change < tol && c_change < tol && q_change < tol
    beta <- beta_new
    cc <- c_new
    if (converged) break
  }
  list(beta = beta, c = cc, beta_effective = beta + cc * v_c,
       objective = objective(beta + cc * v_c), iterations = length(history),
       converged = converged, optim_convergence = opt_beta$convergence,
       history = history)
}

gmm_two_stage <- function(sim_dat, true_beta, n_starts = 20L, loops = 6L,
                          start_lower = NULL, start_upper = NULL,
                          start_sd = .1, maxit = 8000L, outer_tol = 1e-5,
                          cb_max_iter = 20L, cb_tol = 1e-5,
                          c_interval = c(-3, 3), cb_beta_maxit = 5000L,
                          diagnostics_file = NULL) {
  stopifnot(length(true_beta) == 14L, all(is.finite(true_beta)))
  for (x in list(n_starts, loops, maxit, cb_max_iter, cb_beta_maxit))
    stopifnot(length(x) == 1L, is.finite(x), x >= 1, x == as.integer(x))
  stopifnot(length(outer_tol) == 1L, is.finite(outer_tol), outer_tol > 0,
            length(cb_tol) == 1L, is.finite(cb_tol), cb_tol > 0,
            length(c_interval) == 2L, all(is.finite(c_interval)),
            c_interval[1] < c_interval[2])
  uniform <- !is.null(start_lower) || !is.null(start_upper)
  if (uniform) {
    stopifnot(length(start_lower) == 14L, length(start_upper) == 14L,
              all(is.finite(start_lower)), all(is.finite(start_upper)),
              all(start_lower < start_upper))
  } else {
    stopifnot(length(start_sd) %in% c(1L, 14L),
              all(is.finite(start_sd)), all(start_sd > 0))
  }
  diagnostics <- list(settings = list(true_beta = true_beta, n_starts = n_starts,
    loops = loops, start_lower = start_lower, start_upper = start_upper,
    start_sd = start_sd, maxit = maxit, outer_tol = outer_tol,
    cb_max_iter = cb_max_iter, cb_tol = cb_tol, c_interval = c_interval,
    cb_beta_maxit = cb_beta_maxit,
    workflow = "independent_stage1_paths_then_dynamic_W_refinement"))
  checkpoint <- function() {
    if (!is.null(diagnostics_file)) saveRDS(diagnostics, diagnostics_file)
  }
  cat(sprintf("Stage 1: %d independent paths, maximum %d beta-W iterations each\n", n_starts, loops))
  t1 <- proc.time()[[3L]]
  initials <- lapply(seq_len(n_starts), function(id) {
    if (uniform) runif(14, start_lower, start_upper) else
      true_beta + rnorm(14, 0, start_sd)
  })
  diagnostics$initials <- initials
  checkpoint()
  if (identical(Sys.getenv("GMM_STAGE1_BACKEND"), "lsf")) {
    starts <- gmm_stage1_lsf(initials, sim_dat, true_beta, maxit,
                             loops, outer_tol, diagnostics_file)
  } else {
    starts <- lapply(seq_len(n_starts), function(id) {
      gmm_fit_start(id, initials[[id]], sim_dat, true_beta, maxit, loops, outer_tol)
    })
  }
  stage1_seconds <- proc.time()[[3L]] - t1
  # Path convergence is the beta-change criterion. The final inner optimizer
  # code is retained separately and is not a substitute for path convergence.
  valid <- which(vapply(starts, function(s) isTRUE(s$path_converged) &&
    is.null(s$error) && all(is.finite(s$fit$par)) && is.finite(s$fit$value), logical(1)))
  diagnostics$starts <- starts
  diagnostics$unavailable_start_ids <- which(vapply(starts,
    function(s) identical(s$stop_reason, "unavailable_result"), logical(1)))
  diagnostics$failed_start_ids <- which(vapply(starts,
    function(s) !is.null(s$error), logical(1)))
  diagnostics$n_received <- n_starts - length(diagnostics$unavailable_start_ids)
  cat(sprintf("Stage 1: submitted=%d received=%d unavailable=%d failed=%d converged=%d\n",
    n_starts, diagnostics$n_received, length(diagnostics$unavailable_start_ids),
    length(diagnostics$failed_start_ids), length(valid)))
  diagnostics$n_converged <- length(valid)
  checkpoint()
  if (!length(valid)) stop("No Stage 1 paths met outer_tol; inspect diagnostics (no unconverged fallback).")
  ordered <- valid[order(vapply(starts[valid], function(s) s$fit$value, numeric(1)))]
  top <- head(ordered, 10L)
  diagnostics$top_start_ids <- top
  diagnostics$stage1 <- list(start_id = ordered[1], beta = starts[[ordered[1]]]$fit$par,
    objective = starts[[ordered[1]]]$fit$value, seconds = stage1_seconds,
    iterations = starts[[ordered[1]]]$iterations,
    objective_definition = "last optim value using W(beta_before_last_optim)")
  checkpoint()
  cat(sprintf("Stage 1: %d paths converged; Stage 2: refining %d candidates\n", length(valid), length(top)))
  t2 <- proc.time()[[3L]]
  if (identical(Sys.getenv("GMM_STAGE1_BACKEND"), "lsf")) {
    refined <- gmm_lsf_array(
      payload = list(candidates = starts[top], sim_dat = sim_dat,
        cb_max_iter = cb_max_iter, cb_tol = cb_tol, c_interval = c_interval,
        cb_beta_maxit = cb_beta_maxit),
      n_jobs = length(top), worker = "gmm_stage2_worker.R",
      label = "refinements", output_prefix = "candidate", loop = 1L,
      diagnostics_file = diagnostics_file)
    stopifnot(identical(vapply(refined, function(s) as.integer(s$start_id), integer(1)),
                        as.integer(top)))
  } else {
    refined <- lapply(starts[top], function(candidate) {
      gmm_refine_candidate(candidate, sim_dat, cb_max_iter, cb_tol,
                           c_interval, cb_beta_maxit)
    })
  }
  diagnostics$refined <- refined
  checkpoint()
  if (any(vapply(refined, function(s) !is.null(s$error) ||
      !is.finite(s$objective) || any(!is.finite(s$beta_effective)), logical(1)))) {
    stop("Candidate refinement failed; inspect diagnostics.")
  }
  pick <- which.min(vapply(refined, function(s) s$objective, numeric(1)))
  selected <- refined[[pick]]
  beta_hat <- selected$beta_effective
  diagnostics$selected_candidate_rank <- pick
  diagnostics$stage2 <- list(start_id = selected$start_id, beta = beta_hat,
    objective = selected$objective, seconds = proc.time()[[3L]] - t2,
    iterations = selected$iterations, converged = selected$converged,
    objective_definition = "Q(beta_effective; W(beta_effective))")
  # Keep the literal Stage 1 ranking value, and also F(beta) for comparable
  # before/after reporting. The additional F does not affect Stage 1 ranking.
  s1_dynamic <- refined[[1L]]$initial_dynamic_objective
  summary_rows <- lapply(c("stage1", "stage2"), function(label) {
    s <- diagnostics[[label]]
    err <- s$beta - true_beta
    row <- data.frame(stage = label, start_id = s$start_id,
      objective = s$objective, objective_definition = s$objective_definition,
      dynamic_objective = if (label == "stage1") s1_dynamic else s$objective,
      dynamic_objective_improvement = s1_dynamic - diagnostics$stage2$objective,
      seconds = s$seconds, iterations = s$iterations,
      l2_error = sqrt(sum(err^2)), n_submitted = n_starts,
      n_received = diagnostics$n_received,
      n_unavailable = length(diagnostics$unavailable_start_ids),
      n_failed = length(diagnostics$failed_start_ids),
      n_converged = length(valid), n_refined = length(top))
    row[paste0("b", 1:14)] <- as.list(s$beta)
    row[paste0("error_b", 1:14)] <- as.list(err)
    row
  })
  diagnostics$comparison <- do.call(rbind, summary_rows)
  diagnostics$outer_converged <- starts[[selected$start_id]]$path_converged
  diagnostics$selected_stage1_iterations <- starts[[selected$start_id]]$iterations
  checkpoint()
  W_final <- W_hat_fun(beta_hat, sim_dat, Time = 3)
  diagnostics$W_final <- W_final
  checkpoint()
  cat(sprintf("Selected Stage 2 candidate=%d start=%d objective=%.12g\n", pick, selected$start_id, selected$objective))
  list(beta = beta_hat, W = W_final, objective = selected$objective,
       beta_convergence = selected$optim_convergence, diagnostics = diagnostics)
}

# Every evaluation recomputes W from the effective beta (real-data Q_beta_c).
gmm_dynamic_objective <- function(beta, sim_dat) {
  W <- W_hat_fun(beta, sim_dat, Time = 3)
  gmm_obj(beta, sim_dat, W_hat = W, Time = 3)
}

gmm_refine_candidate <- function(candidate, sim_dat, cb_max_iter, cb_tol,
                                  c_interval, cb_beta_maxit) {
  objective <- function(b) gmm_dynamic_objective(b, sim_dat)
  fit <- tryCatch({
    initial_dynamic <- objective(candidate$fit$par)
    result <- gmm_refine(candidate$fit$par, objective,
      max_iter = cb_max_iter, tol = cb_tol, c_interval = c_interval,
      beta_maxit = cb_beta_maxit)
    c(list(initial_dynamic_objective = initial_dynamic), result)
  }, error = function(e) list(error = conditionMessage(e)))
  c(list(start_id = candidate$start_id,
         original_objective = candidate$fit$value), fit)
}

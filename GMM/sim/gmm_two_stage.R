# Fixed-W optimization helpers. Moments, objective and W estimator are supplied
# by gmm_sim.Rmd; no weighting updates occur inside candidate refinement.
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
  W <- diag(42)
  beta_old <- true_beta
  diagnostics <- list(settings = list(true_beta = true_beta, n_starts = n_starts,
    loops = loops, start_lower = start_lower, start_upper = start_upper,
    start_sd = start_sd, maxit = maxit, outer_tol = outer_tol,
    cb_max_iter = cb_max_iter, cb_tol = cb_tol, c_interval = c_interval,
    cb_beta_maxit = cb_beta_maxit), stages = list(), outer_converged = FALSE)
  checkpoint <- function() {
    if (!is.null(diagnostics_file)) saveRDS(diagnostics, diagnostics_file)
  }
  summary_rows <- list()
  for (loop in seq_len(loops)) {
    cat(sprintf("GMM loop %d: %d starts\n", loop, n_starts))
    objective <- function(b) gmm_obj(b, sim_dat, W_hat = W, Time = 3)
    t1 <- proc.time()[[3L]]
    starts <- lapply(seq_len(n_starts), function(id) {
      initial <- if (uniform) runif(14, start_lower, start_upper) else
        beta_old + rnorm(14, 0, start_sd)
      fit <- tryCatch(optim(initial, objective, control = list(maxit = maxit)),
                      error = function(e) list(par = rep(NA_real_, 14),
                        value = Inf, convergence = NA_integer_, message = conditionMessage(e)))
      list(start_id = id, initial = initial,
           initial_distance = sqrt(sum((initial - true_beta)^2)), fit = fit)
    })
    stage1_seconds <- proc.time()[[3L]] - t1
    valid <- which(vapply(starts, function(s) isTRUE(s$fit$convergence == 0L) &&
      all(is.finite(s$fit$par)) && is.finite(s$fit$value), logical(1)))
    stage <- list(W_used = W, starts = starts, n_converged = length(valid))
    if (!length(valid)) {
      stage$error <- "No converged starts; replication stopped."
      diagnostics$stages[[loop]] <- stage
      checkpoint()
      stop(stage$error)
    }
    ordered <- valid[order(vapply(starts[valid], function(s) s$fit$value, numeric(1)))]
    top <- head(ordered, 10L)
    stage$stage1 <- list(start_id = ordered[1], beta = starts[[ordered[1]]]$fit$par,
                         objective = starts[[ordered[1]]]$fit$value, seconds = stage1_seconds)
    cat(sprintf("GMM loop %d: %d converged, refining %d candidates\n", loop, length(valid), length(top)))
    t2 <- proc.time()[[3L]]
    refined <- lapply(top, function(id) {
      fit <- tryCatch(gmm_refine(starts[[id]]$fit$par, objective,
        max_iter = cb_max_iter, tol = cb_tol, c_interval = c_interval,
        beta_maxit = cb_beta_maxit), error = function(e) list(error = conditionMessage(e)))
      c(list(start_id = id, original_objective = starts[[id]]$fit$value), fit)
    })
    stage$refined <- refined
    # Never silently drop a candidate whose refinement errored.
    if (any(vapply(refined, function(s) !is.null(s$error) ||
        !is.finite(s$objective) || any(!is.finite(s$beta_effective)), logical(1)))) {
      stage$error <- "Candidate refinement failed; inspect diagnostics."
      diagnostics$stages[[loop]] <- stage
      checkpoint()
      stop(stage$error)
    }
    pick <- which.min(vapply(refined, function(s) s$objective, numeric(1)))
    selected <- refined[[pick]]
    beta_new <- selected$beta_effective
    stage$selected_candidate_rank <- pick
    stage$stage2 <- list(start_id = selected$start_id, beta = beta_new,
      objective = selected$objective, seconds = proc.time()[[3L]] - t2)
    stage$diff_beta <- max(abs(beta_new - beta_old))
    for (label in c("stage1", "stage2")) {
      s <- stage[[label]]
      err <- s$beta - true_beta
      row <- data.frame(loop = loop, stage = label, start_id = s$start_id,
        objective = s$objective, seconds = s$seconds, l2_error = sqrt(sum(err^2)),
        n_converged = length(valid), n_refined = length(top),
        objective_improvement = stage$stage1$objective - stage$stage2$objective)
      row[paste0("b", 1:14)] <- as.list(s$beta)
      row[paste0("error_b", 1:14)] <- as.list(err)
      summary_rows[[length(summary_rows) + 1L]] <- row
    }
    diagnostics$comparison <- do.call(rbind, summary_rows)
    diagnostics$stages[[loop]] <- stage
    checkpoint() # Preserve selected candidates even if the W update fails.
    W <- W_hat_fun(beta_new, sim_dat, Time = 3)
    diagnostics$stages[[loop]]$W_updated <- W
    beta_old <- beta_new
    diagnostics$outer_converged <- stage$diff_beta < outer_tol
    checkpoint()
    if (diagnostics$outer_converged) break
  }
  list(beta = beta_old, W = W, objective = stage$stage2$objective,
       beta_convergence = selected$optim_convergence, diagnostics = diagnostics)
}

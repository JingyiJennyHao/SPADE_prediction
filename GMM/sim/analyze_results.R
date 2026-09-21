#!/usr/bin/env Rscript

## Usage:
##   Rscript analyze_results.R /path/to/experiment_folder
##
## If no path is supplied, the script asks for one interactively.

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args)) args[[1L]] else readline("Input result folder: ")
input_dir <- normalizePath(path.expand(input_dir), mustWork = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
}
library(ggplot2)

bind_rows_base <- function(data_list) {
  data_list <- Filter(Negate(is.null), data_list)
  if (!length(data_list)) return(data.frame())

  all_names <- unique(unlist(lapply(data_list, names), use.names = FALSE))
  data_list <- lapply(data_list, function(dat) {
    missing_names <- setdiff(all_names, names(dat))
    for (nm in missing_names) dat[[nm]] <- NA
    dat[all_names]
  })
  do.call(rbind, data_list)
}

safe_read <- function(file) {
  tryCatch(
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      warning("Could not read ", file, ": ", conditionMessage(e))
      NULL
    }
  )
}

extract_number <- function(x, pattern) {
  hit <- regmatches(x, regexpr(pattern, x))
  if (!length(hit) || !nzchar(hit)) return(NA_integer_)
  as.integer(sub("[^0-9]*", "", hit))
}

is_seed_dir <- function(x) grepl("^seed[0-9]+$", basename(x))

if (is_seed_dir(input_dir)) {
  seed_dirs <- input_dir
  result_root <- dirname(input_dir)
} else {
  seed_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[is_seed_dir(seed_dirs)]
  result_root <- input_dir
}

if (!length(seed_dirs)) {
  stop("No seed folders named seed1, seed2, ... were found in: ", input_dir)
}

analysis_dir <- file.path(result_root, "analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

loop_files <- unlist(lapply(seed_dirs, function(seed_dir) {
  list.files(
    seed_dir,
    pattern = "^loop_[0-9]+_start_results_seed[0-9]+\\.csv$",
    full.names = TRUE
  )
}), use.names = FALSE)

if (!length(loop_files)) {
  stop("No loop result CSV files were found in: ", input_dir)
}

loop_data <- lapply(loop_files, function(file) {
  dat <- safe_read(file)
  if (is.null(dat) || !nrow(dat)) return(NULL)

  dat$seed <- if ("seed" %in% names(dat)) {
    as.integer(dat$seed)
  } else {
    extract_number(basename(file), "seed[0-9]+")
  }
  dat$loop <- extract_number(basename(file), "loop_[0-9]+")
  dat$source_file <- file
  dat
})

loop_results <- bind_rows_base(loop_data)
beta_cols <- names(loop_results)[grepl("^b[0-9]+$", names(loop_results))]
beta_cols <- beta_cols[order(as.integer(sub("^b", "", beta_cols)))]

if (!length(beta_cols) || !all(c("seed", "loop", "start_id", "value") %in% names(loop_results))) {
  stop("Loop files must contain seed, loop, start_id, value, and beta columns.")
}

loop_results$seed <- as.integer(loop_results$seed)
loop_results$loop <- as.integer(loop_results$loop)
loop_results$start_id <- as.integer(loop_results$start_id)
loop_results$value <- as.numeric(loop_results$value)

for (nm in beta_cols) loop_results[[nm]] <- as.numeric(loop_results[[nm]])

finite_beta <- apply(
  is.finite(as.matrix(loop_results[beta_cols])),
  1L,
  all
)
usable_loop <- is.finite(loop_results$value) & finite_beta

## The worker's diff_from_base is the outer-loop max-absolute beta change.
## For loop 1 it measures the perturbation from the initial value, so it is
## excluded from the outer-loop convergence analysis.
if ("diff_from_base" %in% names(loop_results)) {
  loop_results$max_abs_change <- as.numeric(loop_results$diff_from_base)
  loop_results$max_abs_change[loop_results$loop < 2L] <- NA_real_
} else {
  loop_results$max_abs_change <- NA_real_
  group_key <- interaction(loop_results$seed, loop_results$start_id, drop = TRUE)
  for (rows in split(seq_len(nrow(loop_results)), group_key)) {
    rows <- rows[order(loop_results$loop[rows])]
    if (length(rows) < 2L) next
    for (i in 2:length(rows)) {
      current <- rows[[i]]
      previous <- rows[[i - 1L]]
      if (usable_loop[current] && finite_beta[previous]) {
        loop_results$max_abs_change[current] <- max(abs(
          as.numeric(unlist(loop_results[current, beta_cols], use.names = FALSE)) -
            as.numeric(unlist(loop_results[previous, beta_cols], use.names = FALSE))
        ))
      }
    }
  }
}

write.csv(
  loop_results,
  file.path(analysis_dir, "loop_results_combined.csv"),
  row.names = FALSE
)

## 1. Objective value by loop and start point
objective_data <- loop_results[usable_loop, , drop = FALSE]
if (nrow(objective_data)) {
  objective_plot <- ggplot(
    objective_data,
    aes(
      x = loop,
      y = value,
      color = factor(start_id),
      group = interaction(seed, start_id)
    )
  ) +
    geom_line() +
    geom_point(size = 1.5) +
    facet_wrap(~seed, scales = "free_y") +
    labs(
      title = "Objective Value by GMM Loop",
      x = "Outer GMM loop",
      y = "Objective value",
      color = "Start point"
    ) +
    theme_bw()

  ggsave(
    file.path(analysis_dir, "objective_by_loop.png"),
    objective_plot,
    width = 11,
    height = 7,
    dpi = 300
  )
}

## 2. Maximum absolute beta change by loop
change_data <- loop_results[
  is.finite(loop_results$max_abs_change) & loop_results$loop >= 2L,
  ,
  drop = FALSE
]
write.csv(
  change_data,
  file.path(analysis_dir, "beta_change_by_seed_start_loop.csv"),
  row.names = FALSE
)

if (nrow(change_data)) {
  change_plot <- ggplot(
    change_data,
    aes(
      x = loop,
      y = max_abs_change,
      color = factor(start_id),
      group = interaction(seed, start_id)
    )
  ) +
    geom_hline(yintercept = 1e-5, linetype = "dashed", color = "black") +
    geom_line() +
    geom_point(size = 1.5) +
    facet_wrap(~seed, scales = "free_y") +
    labs(
      title = "Maximum Absolute Beta Change",
      subtitle = "Dashed line: outer-loop convergence tolerance = 1e-5",
      x = "Outer GMM loop",
      y = "max(abs(beta_new - beta_old))",
      color = "Start point"
    ) +
    theme_bw()

  ggsave(
    file.path(analysis_dir, "beta_change_by_loop.png"),
    change_plot,
    width = 11,
    height = 7,
    dpi = 300
  )
}

## 3. Beta trajectories, one plot per seed
beta_long <- do.call(rbind, lapply(beta_cols, function(nm) {
  data.frame(
    seed = loop_results$seed,
    start_id = loop_results$start_id,
    loop = loop_results$loop,
    parameter = nm,
    beta = loop_results[[nm]],
    stringsAsFactors = FALSE
  )
}))

for (seed_value in sort(unique(beta_long$seed))) {
  dat <- beta_long[beta_long$seed == seed_value & is.finite(beta_long$beta), , drop = FALSE]
  if (!nrow(dat)) next

  beta_plot <- ggplot(
    dat,
    aes(
      x = loop,
      y = beta,
      color = factor(start_id),
      group = start_id
    )
  ) +
    geom_line() +
    geom_point(size = 1) +
    facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    labs(
      title = paste("Beta Trajectories: Seed", seed_value),
      x = "Outer GMM loop",
      y = "Beta value",
      color = "Start point"
    ) +
    theme_bw()

  ggsave(
    file.path(analysis_dir, sprintf("beta_trajectories_seed%d.png", seed_value)),
    beta_plot,
    width = 12,
    height = 10,
    dpi = 300
  )
}

## 4. Final objective value by start point
final_files <- unlist(lapply(seed_dirs, function(seed_dir) {
  list.files(
    seed_dir,
    pattern = "^final_start_results_seed[0-9]+\\.csv$",
    full.names = TRUE
  )
}), use.names = FALSE)

if (length(final_files)) {
  final_results <- bind_rows_base(lapply(final_files, safe_read))
} else {
  ## Backward-compatible fallback: use the last observed loop for each path.
  ordered <- loop_results[order(loop_results$seed, loop_results$start_id, loop_results$loop), , drop = FALSE]
  key <- interaction(ordered$seed, ordered$start_id, drop = TRUE)
  final_results <- do.call(rbind, lapply(split(ordered, key), function(dat) dat[nrow(dat), , drop = FALSE]))
}

if (nrow(final_results)) {
  final_results$value <- as.numeric(final_results$value)
  final_results$seed <- as.integer(final_results$seed)
  final_results$start_id <- as.integer(final_results$start_id)
  write.csv(
    final_results,
    file.path(analysis_dir, "final_results_combined.csv"),
    row.names = FALSE
  )

  final_usable <- is.finite(final_results$value)
  if (any(final_usable)) {
    final_plot <- ggplot(
      final_results[final_usable, , drop = FALSE],
      aes(x = factor(start_id), y = value, fill = factor(seed))
    ) +
      geom_col(position = "dodge") +
      facet_wrap(~seed, scales = "free_x") +
      labs(
        title = "Final Objective Value by Start Point",
        x = "Start point",
        y = "Final objective value",
        fill = "Seed"
      ) +
      theme_bw()

    ggsave(
      file.path(analysis_dir, "final_objective_by_start.png"),
      final_plot,
      width = 11,
      height = 7,
      dpi = 300
    )
  }
}

## 5. Final inference results and overall averages across seeds
inference_files <- unlist(lapply(seed_dirs, function(seed_dir) {
  list.files(
    seed_dir,
    pattern = "^inference_seed[0-9]+\\.csv$",
    full.names = TRUE
  )
}), use.names = FALSE)

if (length(inference_files)) {
  inference_data <- lapply(inference_files, function(file) {
    dat <- safe_read(file)
    if (is.null(dat) || !nrow(dat)) return(NULL)
    if (!"seed" %in% names(dat)) {
      dat$seed <- extract_number(basename(file), "seed[0-9]+")
    }
    dat$source_file <- file
    dat
  })
  inference_all <- bind_rows_base(inference_data)

  inference_all$seed <- as.integer(inference_all$seed)
  write.csv(
    inference_all,
    file.path(analysis_dir, "inference_all_seeds.csv"),
    row.names = FALSE
  )

  numeric_names <- names(inference_all)[
    vapply(inference_all, is.numeric, logical(1))
  ]
  average_names <- setdiff(
    numeric_names,
    c("seed", "J", "Time", "interval_prob")
  )

  inference_average <- do.call(rbind, lapply(average_names, function(nm) {
    x <- inference_all[[nm]]
    data.frame(
      variable = nm,
      n_seeds = sum(is.finite(x)),
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  write.csv(
    inference_average,
    file.path(analysis_dir, "inference_overall_average.csv"),
    row.names = FALSE
  )

  ## Average beta estimates across seeds, with between-seed standard deviation.
  inference_beta_cols <- names(inference_all)[
    grepl("^b[0-9]+$", names(inference_all))
  ]
  if (length(inference_beta_cols)) {
    beta_average <- inference_average[
      inference_average$variable %in% inference_beta_cols,
      ,
      drop = FALSE
    ]
    beta_average$parameter <- factor(
      beta_average$variable,
      levels = inference_beta_cols
    )

    beta_plot <- ggplot(
      beta_average,
      aes(x = parameter, y = mean)
    ) +
      geom_hline(yintercept = 0, color = "grey70") +
      geom_errorbar(
        aes(ymin = mean - sd, ymax = mean + sd),
        width = 0.2
      ) +
      geom_point(size = 2) +
      labs(
        title = "Average Final Beta Across Seeds",
        subtitle = "Error bars show between-seed standard deviation",
        x = "Parameter",
        y = "Average beta"
      ) +
      theme_bw()

    ggsave(
      file.path(analysis_dir, "inference_beta_overall_average.png"),
      beta_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
  }

  ## Average inference coverage across seeds.
  coverage_cols <- names(inference_all)[
    grepl("coverage", names(inference_all), ignore.case = TRUE)
  ]
  if (length(coverage_cols)) {
    coverage_average <- inference_average[
      inference_average$variable %in% coverage_cols,
      ,
      drop = FALSE
    ]
    coverage_average$variable <- factor(
      coverage_average$variable,
      levels = coverage_cols
    )

    coverage_plot <- ggplot(
      coverage_average,
      aes(x = variable, y = mean)
    ) +
      geom_col(fill = "steelblue") +
      geom_hline(yintercept = 0.95, linetype = "dashed") +
      coord_cartesian(ylim = c(0, 1)) +
      labs(
        title = "Average Inference Coverage Across Seeds",
        subtitle = "Dashed line: nominal coverage level = 0.95",
        x = "Coverage measure",
        y = "Average coverage"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(
      file.path(analysis_dir, "inference_coverage_overall_average.png"),
      coverage_plot,
      width = 12,
      height = 7,
      dpi = 300
    )
  }
} else {
  warning("No final inference_seedX.csv files were found; inference analysis was skipped.")
}

cat("Analysis complete. Outputs were written to:\n", analysis_dir, "\n", sep = "")

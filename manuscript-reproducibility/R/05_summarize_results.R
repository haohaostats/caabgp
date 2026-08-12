source(file.path("R", "00_settings.R"))
source(file.path("R", "01_surfaces.R"))

standard_error <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

summarise_mean_se <- function(df, group_cols, metrics) {
  if (nrow(df) == 0) return(data.frame())
  groups <- split(df, interaction(df[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE))
  out <- lapply(groups, function(g) {
    row <- g[1, group_cols, drop = FALSE]
    for (m in metrics) {
      row[[paste0(m, "_mean")]] <- mean(g[[m]], na.rm = TRUE)
      row[[paste0(m, "_se")]] <- standard_error(g[[m]])
    }
    row$n_replicates <- length(unique(g$replicate))
    row
  })
  ans <- do.call(rbind, out)
  row.names(ans) <- NULL
  if ("method" %in% names(ans)) ans$method <- as.character(method_order(ans$method))
  ans
}

final_rows_by <- function(df, group_cols) {
  groups <- split(df, interaction(df[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE))
  out <- lapply(groups, function(g) g[g$n == max(g$n, na.rm = TRUE), , drop = FALSE][1, , drop = FALSE])
  ans <- do.call(rbind, out)
  row.names(ans) <- NULL
  ans
}

format_mean_se <- function(mean, se, digits = 3) {
  if (!is.finite(mean)) return("--")
  if (!is.finite(se)) return(sprintf(paste0("%.", digits, "f (--)" ), mean))
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean, se)
}

scenario_result_caption <- function(scenario) {
  switch(
    as.character(scenario),
    "1" = "Fixed-budget operating characteristics in Scenario 1: shared main surface with moderate subgroup deviations. Values are Monte Carlo means, with standard errors in parentheses.",
    "2" = "Fixed-budget operating characteristics in Scenario 2: rare, noisy, and clinically prioritized stratum. Values are Monte Carlo means, with standard errors in parentheses.",
    "3" = "Fixed-budget operating characteristics in Scenario 3: mixed heterogeneity, rare-stratum prioritization, and high manufacturing cost. Values are Monte Carlo means, with standard errors in parentheses.",
    "4" = "Fixed-budget operating characteristics in Scenario 4: high outcome noise and value of replication. Values are Monte Carlo means, with standard errors in parentheses.",
    "Fixed-budget operating characteristics. Values are Monte Carlo means, with standard errors in parentheses."
  )
}

allocation_by_stratum_summary <- function(allocations, final_traj, scenario, n_mc, prefix) {
  alloc_s <- allocations[allocations$scenario == scenario, ]
  K <- length(scenario_prevalence(scenario))
  methods_s <- unique(final_traj$method[final_traj$scenario == scenario])
  reps_s <- unique(final_traj$replicate[final_traj$scenario == scenario])

  if (nrow(alloc_s) == 0) return(data.frame())
  alloc_count_long <- do.call(rbind, lapply(seq_len(K), function(k) {
    data.frame(
      method = alloc_s$method,
      replicate = alloc_s$replicate,
      stratum = k,
      n_patients = alloc_s[[paste0("n_z", k)]]
    )
  }))
  alloc_count_rep <- aggregate(n_patients ~ method + stratum + replicate, alloc_count_long, sum, na.rm = TRUE)
  base <- expand.grid(method = methods_s, replicate = reps_s, stratum = seq_len(K))
  alloc_count_rep <- merge(base, alloc_count_rep, by = c("method", "replicate", "stratum"), all.x = TRUE)
  alloc_count_rep$n_patients[!is.finite(alloc_count_rep$n_patients)] <- 0

  init <- initial_counts_by_stratum(scenario)
  alloc_count_rep$n_patients <- alloc_count_rep$n_patients + init[alloc_count_rep$stratum]
  out <- summarise_mean_se(
    alloc_count_rep,
    group_cols = c("method", "stratum"),
    metrics = c("n_patients")
  )
  write.csv(out, file.path(DIR_RESULTS, sprintf("%s_allocation_by_stratum_n%d.csv", prefix, n_mc)), row.names = FALSE)
  out
}

final_operational_summary <- function(final_traj, scenario, n_mc, prefix) {
  out <- summarise_mean_se(
    final_traj[final_traj$scenario == scenario, ],
    group_cols = c("method"),
    metrics = c("n", "unique_doses", "total_cost", "budget_fraction", "regret", "rpsel")
  )
  write.csv(out, file.path(DIR_RESULTS, sprintf("%s_final_operational_n%d.csv", prefix, n_mc)), row.names = FALSE)
  out
}

cumulative_unique_summary <- function(allocations, scenario, n_mc, prefix) {
  out <- summarise_mean_se(
    allocations[allocations$scenario == scenario, ],
    group_cols = c("method", "n_after"),
    metrics = c("unique_doses")
  )
  if ("n_after" %in% names(out)) names(out)[names(out) == "n_after"] <- "n"
  write.csv(out, file.path(DIR_RESULTS, sprintf("%s_cumulative_unique_doses_n%d.csv", prefix, n_mc)), row.names = FALSE)
  out
}

learning_summary <- function(trajectory, scenario, n_mc, prefix) {
  out <- summarise_mean_se(
    trajectory[trajectory$scenario == scenario, ],
    group_cols = c("method", "n"),
    metrics = c("regret", "near_optimal", "dose_distance", "rpsel", "pmad", "unique_doses", "total_cost")
  )
  write.csv(out, file.path(DIR_RESULTS, sprintf("%s_learning_trajectory_n%d.csv", prefix, n_mc)), row.names = FALSE)
  out
}

final_regret_distribution <- function(final_traj, scenario, n_mc, prefix) {
  out <- final_traj[final_traj$scenario == scenario, c("scenario", "replicate", "method", "n", "regret", "rpsel", "unique_doses", "total_cost")]
  write.csv(out, file.path(DIR_RESULTS, sprintf("%s_final_distribution_n%d.csv", prefix, n_mc)), row.names = FALSE)
  out
}

near_optimal_by_stratum_summary <- function(final_recs, scenario, n_mc, prefix) {
  out <- summarise_mean_se(
    final_recs[final_recs$scenario == scenario, ],
    group_cols = c("method", "stratum"),
    metrics = c("near_optimal", "dose_distance", "value_regret", "rpsel")
  )
  write.csv(out, file.path(DIR_RESULTS, sprintf("%s_near_optimal_by_stratum_n%d.csv", prefix, n_mc)), row.names = FALSE)
  out
}

borrowing_final_summary <- function(borrowing, scenario, n_mc, prefix) {
  if (is.null(borrowing)) return(data.frame())
  borrow_s <- borrowing[borrowing$scenario == scenario, ]
  if (nrow(borrow_s) == 0) return(data.frame())
  final_borrow_stratum <- final_rows_by(borrow_s, c("scenario", "method", "replicate", "stratum"))
  final_borrow_rep <- aggregate(rho ~ scenario + method + replicate, final_borrow_stratum, mean, na.rm = TRUE)
  out <- summarise_mean_se(
    final_borrow_rep,
    group_cols = c("method"),
    metrics = c("rho")
  )
  write.csv(out, file.path(DIR_RESULTS, sprintf("%s_borrowing_final_n%d.csv", prefix, n_mc)), row.names = FALSE)
  out
}

novel_fraction_summary <- function(allocations, scenario, n_mc, prefix) {
  alloc_s <- allocations[allocations$scenario == scenario, ]
  out <- summarise_mean_se(
    alloc_s,
    group_cols = c("method"),
    metrics = c("new_dose")
  )
  if (nrow(out) > 0) {
    names(out)[names(out) == "new_dose_mean"] <- "novel_fraction_mean"
    names(out)[names(out) == "new_dose_se"] <- "novel_fraction_se"
    out$repeat_fraction_mean <- 1 - out$novel_fraction_mean
    out$repeat_fraction_se <- out$novel_fraction_se
  }
  write.csv(out, file.path(DIR_RESULTS, sprintf("%s_novel_repeat_fraction_n%d.csv", prefix, n_mc)), row.names = FALSE)
  out
}

summarise_primary_outputs <- function(n_mc = SIM_SETTINGS$n_mc) {
  trajectory <- read.csv(file.path(DIR_RESULTS, sprintf("simulation_trajectory_n%d.csv", n_mc)))
  recommendations <- read.csv(file.path(DIR_RESULTS, sprintf("simulation_recommendations_n%d.csv", n_mc)))
  allocations <- read.csv(file.path(DIR_RESULTS, sprintf("simulation_allocations_n%d.csv", n_mc)))
  borrowing_path <- file.path(DIR_RESULTS, sprintf("simulation_borrowing_n%d.csv", n_mc))
  borrowing <- if (file.exists(borrowing_path)) read.csv(borrowing_path) else NULL

  final_traj <- final_rows_by(trajectory, c("scenario", "method", "replicate"))
  final_recs <- final_rows_by(recommendations, c("scenario", "method", "replicate", "stratum"))

  main_table <- summarise_mean_se(
    final_traj,
    group_cols = c("scenario", "method"),
    metrics = c("n", "regret", "near_optimal", "dose_distance", "rpsel", "pmad", "unique_doses", "total_cost", "budget_fraction")
  )

  if (!is.null(borrowing)) {
    final_borrow_stratum <- final_rows_by(borrowing, c("scenario", "method", "replicate", "stratum"))
    final_borrow_rep <- aggregate(rho ~ scenario + method + replicate, final_borrow_stratum, mean, na.rm = TRUE)
    rho_table <- summarise_mean_se(
      final_borrow_rep,
      group_cols = c("scenario", "method"),
      metrics = c("rho")
    )
    main_table <- merge(main_table, rho_table[, c("scenario", "method", "rho_mean", "rho_se")],
                        by = c("scenario", "method"), all.x = TRUE)
  } else {
    main_table$rho_mean <- NA_real_
    main_table$rho_se <- NA_real_
  }

  main_table <- main_table[order(main_table$scenario, match(main_table$method, METHODS)), ]
  write.csv(main_table, file.path(DIR_RESULTS, sprintf("table_main_results_n%d.csv", n_mc)), row.names = FALSE)

  for (s in 1:4) {
    tab_s <- main_table[main_table$scenario == s, ]
    write.csv(tab_s, file.path(DIR_RESULTS, sprintf("table_scenario%d_results_n%d.csv", s, n_mc)), row.names = FALSE)
  }

  learning_summary(trajectory, 1, n_mc, "figure2_scenario1")
  final_regret_distribution(final_traj, 1, n_mc, "figure2_scenario1")
  borrowing_final_summary(borrowing, 1, n_mc, "figure2_scenario1")
  cumulative_unique_summary(allocations, 1, n_mc, "figure2_scenario1")

  learning_summary(trajectory, 2, n_mc, "figure3_scenario2")
  near_optimal_by_stratum_summary(final_recs, 2, n_mc, "figure3_scenario2")
  allocation_by_stratum_summary(allocations, final_traj, 2, n_mc, "figure3_scenario2")
  final_operational_summary(final_traj, 2, n_mc, "figure3_scenario2")

  final_operational_summary(final_traj, 3, n_mc, "figure4_scenario3")
  cumulative_unique_summary(allocations, 3, n_mc, "figure4_scenario3")
  near_optimal_by_stratum_summary(final_recs, 3, n_mc, "figure4_scenario3")

  learning_summary(trajectory, 4, n_mc, "figure5_scenario4")
  novel_fraction_summary(allocations, 4, n_mc, "figure5_scenario4")
  final_operational_summary(final_traj, 4, n_mc, "figure5_scenario4")
  final_regret_distribution(final_traj, 4, n_mc, "figure5_scenario4")

  sens_all_path <- file.path(DIR_RESULTS, sprintf("sensitivity_trajectory_n%d.csv", n_mc))
  if (file.exists(sens_all_path)) {
    sens <- read.csv(sens_all_path)
    sens_final <- final_rows_by(sens, c("sensitivity", "sensitivity_value", "method", "replicate"))

    fig6_lambda <- summarise_mean_se(
      sens_final[sens_final$sensitivity == "lambda_c", ],
      group_cols = c("sensitivity_value"),
      metrics = c("n", "regret", "unique_doses", "total_cost", "near_optimal")
    )
    names(fig6_lambda)[names(fig6_lambda) == "sensitivity_value"] <- "lambda_c"
    write.csv(fig6_lambda, file.path(DIR_RESULTS, sprintf("figure6_lambda_sensitivity_n%d.csv", n_mc)), row.names = FALSE)

    fig6_cnew <- summarise_mean_se(
      sens_final[sens_final$sensitivity == "c_new", ],
      group_cols = c("sensitivity_value"),
      metrics = c("n", "regret", "unique_doses", "total_cost", "near_optimal")
    )
    names(fig6_cnew)[names(fig6_cnew) == "sensitivity_value"] <- "c_new"
    write.csv(fig6_cnew, file.path(DIR_RESULTS, sprintf("figure6_cnew_sensitivity_n%d.csv", n_mc)), row.names = FALSE)

    fig6_budget <- summarise_mean_se(
      sens_final[sens_final$sensitivity == "budget_multiplier", ],
      group_cols = c("sensitivity_value"),
      metrics = c("n", "regret", "unique_doses", "total_cost", "near_optimal")
    )
    names(fig6_budget)[names(fig6_budget) == "sensitivity_value"] <- "budget_multiplier"
    fig6_budget$b_max <- scenario_budget(SIM_SETTINGS$sensitivity_scenario) * fig6_budget$budget_multiplier
    write.csv(fig6_budget, file.path(DIR_RESULTS, sprintf("figure6_budget_sensitivity_n%d.csv", n_mc)), row.names = FALSE)

    fig6_cohort <- summarise_mean_se(
      sens_final[sens_final$sensitivity == "cohort_size", ],
      group_cols = c("sensitivity_value"),
      metrics = c("n", "regret", "unique_doses", "total_cost", "near_optimal")
    )
    names(fig6_cohort)[names(fig6_cohort) == "sensitivity_value"] <- "cohort_size"
    write.csv(fig6_cohort, file.path(DIR_RESULTS, sprintf("figure6_cohort_sensitivity_n%d.csv", n_mc)), row.names = FALSE)

    sensitivity_block <- function(df, factor, factor_label, setting_label) {
      data.frame(
        factor = factor,
        factor_label = factor_label,
        setting_label = setting_label,
        n_mean = df$n_mean,
        n_se = df$n_se,
        regret_mean = df$regret_mean,
        regret_se = df$regret_se,
        near_optimal_mean = df$near_optimal_mean,
        near_optimal_se = df$near_optimal_se,
        unique_doses_mean = df$unique_doses_mean,
        unique_doses_se = df$unique_doses_se,
        total_cost_mean = df$total_cost_mean,
        total_cost_se = df$total_cost_se,
        n_replicates = df$n_replicates,
        stringsAsFactors = FALSE
      )
    }

    sensitivity_table <- rbind(
      sensitivity_block(
        fig6_lambda,
        "lambda_c",
        "lambda_c",
        format(fig6_lambda$lambda_c, trim = TRUE, scientific = FALSE)
      ),
      sensitivity_block(
        fig6_cnew,
        "c_new",
        "c_new",
        format(fig6_cnew$c_new, trim = TRUE, scientific = FALSE)
      ),
      sensitivity_block(
        fig6_budget,
        "budget",
        "B_max",
        sprintf(
          "%s x primary (%s)",
          format(fig6_budget$budget_multiplier, trim = TRUE, nsmall = 1),
          format(fig6_budget$b_max, trim = TRUE, scientific = FALSE)
        )
      ),
      sensitivity_block(
        fig6_cohort,
        "cohort_size",
        "Cohort size r",
        format(fig6_cohort$cohort_size, trim = TRUE, scientific = FALSE)
      )
    )
    write.csv(
      sensitivity_table,
      file.path(DIR_RESULTS, sprintf("table_sensitivity_results_n%d.csv", n_mc)),
      row.names = FALSE
    )
  }

  invisible(list(main_table = main_table))
}

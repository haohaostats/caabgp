source(file.path("R", "00_settings.R"))
source(file.path("R", "01_surfaces.R"))
source(file.path("R", "02_gp_eb.R"))
source(file.path("R", "03_acquisition.R"))

make_noise_bank <- function(scenario, replicate, max_per_pair = SIM_SETTINGS$n_max + SIM_SETTINGS$n_initial) {
  prevalence <- scenario_prevalence(scenario)
  sigma <- scenario_noise_sd(scenario)
  grid <- candidate_grid()
  set.seed(SIM_SETTINGS$seed_base + scenario * 100000 + replicate)
  out <- array(NA_real_, dim = c(length(prevalence), nrow(grid), max_per_pair))
  for (k in seq_along(prevalence)) {
    out[k, , ] <- stats::rnorm(nrow(grid) * max_per_pair, 0, sigma[k])
  }
  out
}

draw_outcomes <- function(noise_bank, counters, scenario, k, d_idx, n) {
  idx <- counters[k, d_idx] + seq_len(n)
  if (max(idx) > dim(noise_bank)[3]) {
    stop("Noise bank exhausted for scenario ", scenario, ", stratum ", k, ", dose index ", d_idx)
  }
  counters[k, d_idx] <- counters[k, d_idx] + n
  grid <- candidate_grid()
  f <- true_surface(scenario, k, matrix(grid[d_idx, ], nrow = 1))
  list(y = as.numeric(f + noise_bank[k, d_idx, idx]), counters = counters)
}

make_initial_data <- function(scenario, replicate) {
  prevalence <- scenario_prevalence(scenario)
  grid <- candidate_grid()
  init_doses <- SIM_SETTINGS$initial_design
  counts <- initial_counts_by_stratum(scenario)
  noise_bank <- make_noise_bank(scenario, replicate)
  counters <- matrix(0L, nrow = length(prevalence), ncol = nrow(grid))

  rows <- list()
  row_id <- 1
  for (k in seq_along(counts)) {
    dose_order <- rep(seq_len(nrow(init_doses)), length.out = counts[k])
    for (j in dose_order) {
      d <- init_doses[j, ]
      d_idx <- dose_row_index(grid, d)
      draw <- draw_outcomes(noise_bank, counters, scenario, k, d_idx, 1)
      counters <- draw$counters
      rows[[row_id]] <- data.frame(
        d1 = d[1],
        d2 = d[2],
        stratum = k,
        y = draw$y
      )
      row_id <- row_id + 1
    }
  }

  list(
    data = do.call(rbind, rows),
    noise_bank = noise_bank,
    counters = counters
  )
}

initial_cost <- function(data, scenario, c_new = scenario_cost_new(scenario), c_scr = scenario_c_scr(scenario)) {
  prevalence <- scenario_prevalence(scenario)
  manufactured <- unique(dose_key(data$d1, data$d2))
  screening <- sum(vapply(data$stratum, function(k) 1 / max(prevalence[k], SIM_SETTINGS$pi_min), numeric(1)))
  SIM_SETTINGS$c_pt * nrow(data) +
    c_new * length(manufactured) +
    c_scr * screening
}

fit_design <- function(method, data, K, previous = NULL) {
  X <- as.matrix(data[, c("d1", "d2")])
  z <- as.integer(data$stratum)
  y <- data$y

  if (method == "S-GP") {
    previous_par <- if (!is.null(previous)) previous$fit$par else NULL
    return(list(method = method, fit = fit_gp_eb(X, z, y, "S", K, previous_par = previous_par)))
  }

  if (method == "P-GP") {
    previous_par <- if (!is.null(previous)) previous$fit$par else NULL
    return(list(method = method, fit = fit_gp_eb(X, z, y, "P", K, previous_par = previous_par)))
  }

  if (method %in% c("AB-GP", "CA-AB-GP")) {
    previous_par <- if (!is.null(previous)) previous$fit$par else NULL
    return(list(method = method, fit = fit_gp_eb(X, z, y, "AB", K, previous_par = previous_par)))
  }

  if (method == "IND-GP") {
    fits <- vector("list", K)
    for (k in seq_len(K)) {
      idx <- which(z == k)
      prev_par <- if (!is.null(previous) && !is.null(previous$fits[[k]])) previous$fits[[k]]$par else NULL
      fits[[k]] <- fit_gp_eb(X[idx, , drop = FALSE], z[idx], y[idx], "S", K, previous_par = prev_par)
    }
    return(list(method = method, fits = fits))
  }

  stop("Unknown method: ", method)
}

predict_design <- function(design_fit, grid, K) {
  method <- design_fit$method
  rows <- list()

  if (method == "IND-GP") {
    for (k in seq_len(K)) {
      pred <- predict_gp_eb(design_fit$fits[[k]], grid, rep(k, nrow(grid)))
      rows[[k]] <- data.frame(stratum = k, d1 = grid[, 1], d2 = grid[, 2], mu = pred$mu, sd = pred$sd)
    }
    return(do.call(rbind, rows))
  }

  for (k in seq_len(K)) {
    znew <- if (method == "S-GP") rep(1L, nrow(grid)) else rep(k, nrow(grid))
    pred <- predict_gp_eb(design_fit$fit, grid, znew)
    rows[[k]] <- data.frame(stratum = k, d1 = grid[, 1], d2 = grid[, 2], mu = pred$mu, sd = pred$sd)
  }
  do.call(rbind, rows)
}

sigma_by_stratum <- function(design_fit, K) {
  if (design_fit$method == "IND-GP") {
    return(vapply(design_fit$fits, sigma_hat, numeric(1)))
  }
  rep(sigma_hat(design_fit$fit), K)
}

borrowing_rows <- function(design_fit, scenario, replicate, n_current) {
  if (!design_fit$method %in% c("AB-GP", "CA-AB-GP")) return(NULL)
  bi <- borrowing_index(design_fit$fit)
  data.frame(
    scenario = scenario,
    replicate = replicate,
    method = design_fit$method,
    n = n_current,
    stratum = bi$stratum,
    rho = bi$rho
  )
}

recommendation_metrics <- function(pred, scenario, method, replicate, n_current, total_cost, unique_doses, b_max) {
  prevalence <- scenario_prevalence(scenario)
  weights <- scenario_weights(scenario)
  grid <- candidate_grid()
  opt <- true_optima(scenario, grid)
  K <- length(prevalence)

  recs <- vector("list", K)
  for (k in seq_len(K)) {
    pred_k <- pred[pred$stratum == k, ]
    idx <- which.min(pred_k$mu)
    d_hat <- c(pred_k$d1[idx], pred_k$d2[idx])
    truth <- true_surface(scenario, k, matrix(d_hat, nrow = 1))
    dist_grid <- sqrt((d_hat[1] - opt$d1_star[k])^2 + (d_hat[2] - opt$d2_star[k])^2) / 0.25
    value_regret <- as.numeric(truth - opt$f_star[k])
    recs[[k]] <- data.frame(
      scenario = scenario,
      replicate = replicate,
      method = method,
      n = n_current,
      stratum = k,
      d1_hat = d_hat[1],
      d2_hat = d_hat[2],
      mu_hat = pred_k$mu[idx],
      sd_hat = pred_k$sd[idx],
      true_value = truth,
      f_star = opt$f_star[k],
      value_regret = value_regret,
      dose_distance = dist_grid,
      near_optimal = as.numeric(dist_grid <= 1),
      pmad = abs(pred_k$mu[idx] - truth),
      rpsel = sqrt(pred_k$sd[idx]^2 + (pred_k$mu[idx] - truth)^2),
      weight = weights[k],
      prevalence = prevalence[k],
      total_cost = total_cost,
      unique_doses = unique_doses,
      b_max = b_max
    )
  }

  rec_df <- do.call(rbind, recs)
  weighted <- data.frame(
    scenario = scenario,
    replicate = replicate,
    method = method,
    n = n_current,
    regret = sum(rec_df$weight * rec_df$value_regret),
    near_optimal = sum(rec_df$weight * rec_df$near_optimal),
    dose_distance = sum(rec_df$weight * rec_df$dose_distance),
    rpsel = sum(rec_df$weight * rec_df$rpsel),
    pmad = sum(rec_df$weight * rec_df$pmad),
    unique_doses = unique_doses,
    total_cost = total_cost,
    b_max = b_max,
    budget_fraction = total_cost / b_max
  )

  list(weighted = weighted, recommendations = rec_df)
}

run_method_once <- function(
    scenario,
    replicate,
    method,
    initial,
    lambda_c = SIM_SETTINGS$lambda_c,
    c_new = scenario_cost_new(scenario),
    c_scr = scenario_c_scr(scenario),
    b_max = scenario_budget(scenario),
    cohort_size = SIM_SETTINGS$cohort_size,
    n_max = SIM_SETTINGS$n_max,
    stopping_threshold = SIM_SETTINGS$stopping_threshold,
    stopping_patience = SIM_SETTINGS$stopping_patience) {
  prevalence <- scenario_prevalence(scenario)
  weights <- scenario_weights(scenario)
  K <- length(prevalence)
  grid <- candidate_grid()
  data <- initial$data
  noise_bank <- initial$noise_bank
  counters <- initial$counters
  total_cost <- initial_cost(data, scenario, c_new, c_scr)
  manufactured_keys <- unique(dose_key(data$d1, data$d2))

  previous_fit <- NULL
  trajectory <- list()
  recommendation_rows <- list()
  allocation_rows <- list()
  borrowing <- list()
  active_strata <- seq_len(K)
  stopping_counters <- integer(K)
  completed_updates <- 0L
  iter <- 1
  alloc_iter <- 1
  borrow_iter <- 1

  repeat {
    design_fit <- fit_design(method, data, K, previous = previous_fit)
    previous_fit <- design_fit
    pred <- predict_design(design_fit, grid, K)
    metrics <- recommendation_metrics(
      pred = pred,
      scenario = scenario,
      method = method,
      replicate = replicate,
      n_current = nrow(data),
      total_cost = total_cost,
      unique_doses = length(manufactured_keys),
      b_max = b_max
    )
    trajectory[[iter]] <- metrics$weighted
    recommendation_rows[[iter]] <- metrics$recommendations
    bi <- borrowing_rows(design_fit, scenario, replicate, nrow(data))
    if (!is.null(bi)) {
      borrowing[[borrow_iter]] <- bi
      borrow_iter <- borrow_iter + 1
    }
    iter <- iter + 1

    sigma_k <- sigma_by_stratum(design_fit, K)

    if (method != "S-GP" && completed_updates > 0L && length(active_strata)) {
      threshold <- rep_len(as.numeric(stopping_threshold), K)
      acquisition_for_stopping <- personalized_acquisition_table(
        pred = pred,
        grid = grid,
        manufactured_keys = manufactured_keys,
        prevalence = prevalence,
        sigma_by_k = sigma_k,
        c_new = c_new,
        c_scr = c_scr,
        lambda_c = if (method == "CA-AB-GP") lambda_c else 0,
        weights = weights,
        cohort_size = cohort_size,
        remaining_budget = Inf,
        active_strata = active_strata
      )
      for (k in active_strata) {
        scores_k <- acquisition_for_stopping$score[acquisition_for_stopping$stratum == k]
        max_score_k <- if (length(scores_k)) max(scores_k, na.rm = TRUE) else -Inf
        if (is.finite(max_score_k) && max_score_k < threshold[k]) {
          stopping_counters[k] <- stopping_counters[k] + 1L
        } else {
          stopping_counters[k] <- 0L
        }
      }
      active_strata <- active_strata[stopping_counters[active_strata] < stopping_patience]
    }

    if (!length(active_strata)) break
    if (nrow(data) + cohort_size > n_max) break

    remaining_budget <- b_max - total_cost
    if (remaining_budget <= 0) break
    current_cohort <- cohort_size

    if (method == "S-GP") {
      pred_one <- pred[pred$stratum == 1, ]
      action <- select_pooled_feasible_action(
        pred_one = pred_one,
        grid = grid,
        sigma = sigma_k[1],
        manufactured_keys = manufactured_keys,
        prevalence = prevalence,
        c_new = c_new,
        c_scr = c_scr,
        cohort_size = cohort_size,
        remaining_budget = remaining_budget
      )
      if (is.null(action)) break
      feasible_k <- unlist(action$feasible_strata[[1]], use.names = FALSE)
      set.seed(SIM_SETTINGS$seed_base + scenario * 1000000 + replicate * 100 + alloc_iter + cohort_size)
      target_k <- sample(feasible_k, size = 1, prob = weights[feasible_k])
      target_strata <- rep(target_k, cohort_size)
    } else {
      action_lambda <- if (method == "CA-AB-GP") lambda_c else 0
      action <- select_personalized_action(
        pred = pred,
        grid = grid,
        manufactured_keys = manufactured_keys,
        prevalence = prevalence,
        sigma_by_k = sigma_k,
        c_new = c_new,
        c_scr = c_scr,
        lambda_c = action_lambda,
        weights = weights,
        cohort_size = cohort_size,
        remaining_budget = remaining_budget,
        active_strata = active_strata
      )
      if (is.null(action)) break
      target_strata <- rep(action$stratum, cohort_size)
    }

    d <- c(action$d1, action$d2)
    d_idx <- dose_row_index(grid, d)
    new_dose <- !is_manufactured(d, manufactured_keys)
    cost_increment <- 0
    new_rows <- vector("list", length(target_strata))
    for (j in seq_along(target_strata)) {
      k <- target_strata[j]
      draw <- draw_outcomes(noise_bank, counters, scenario, k, d_idx, 1)
      counters <- draw$counters
      new_rows[[j]] <- data.frame(d1 = d[1], d2 = d[2], stratum = k, y = draw$y)
      cost_increment <- cost_increment +
        SIM_SETTINGS$c_pt +
        c_scr / max(prevalence[k], SIM_SETTINGS$pi_min)
    }
    if (new_dose) cost_increment <- cost_increment + c_new
    total_cost <- total_cost + cost_increment
    data <- rbind(data, do.call(rbind, new_rows))
    manufactured_keys <- unique(c(manufactured_keys, dose_key(d[1], d[2])))
    alloc_counts <- tabulate(target_strata, nbins = K)
    alloc_counts4 <- c(alloc_counts, rep(NA_integer_, 4 - K))

    allocation_rows[[alloc_iter]] <- data.frame(
      scenario = scenario,
      replicate = replicate,
      method = method,
      allocation_step = alloc_iter,
      n_before = nrow(data) - length(target_strata),
      n_after = nrow(data),
      stratum = ifelse(method == "S-GP", NA_integer_, action$stratum),
      d1 = d[1],
      d2 = d[2],
      new_dose = as.numeric(new_dose),
      action_score = action$score,
      action_aei = action$aei,
      cost_increment = cost_increment,
      total_cost = total_cost,
      b_max = b_max,
      budget_fraction = total_cost / b_max,
      unique_doses = length(manufactured_keys),
      lambda_c = ifelse(method == "CA-AB-GP", lambda_c, 0),
      c_new = c_new,
      c_scr = c_scr,
      cohort_size = current_cohort,
      planned_cohort_size = cohort_size,
      n_z1 = alloc_counts4[1],
      n_z2 = alloc_counts4[2],
      n_z3 = alloc_counts4[3],
      n_z4 = alloc_counts4[4]
    )
    alloc_iter <- alloc_iter + 1
    completed_updates <- completed_updates + 1L
  }

  list(
    trajectory = do.call(rbind, trajectory),
    recommendations = do.call(rbind, recommendation_rows),
    allocations = if (length(allocation_rows)) do.call(rbind, allocation_rows) else NULL,
    borrowing = if (length(borrowing)) do.call(rbind, borrowing) else NULL
  )
}

run_one_replicate <- function(scenario, replicate, methods = METHODS) {
  initial <- make_initial_data(scenario, replicate)
  out <- lapply(methods, function(m) run_method_once(scenario, replicate, m, initial))
  list(
    trajectory = do.call(rbind, lapply(out, `[[`, "trajectory")),
    recommendations = do.call(rbind, lapply(out, `[[`, "recommendations")),
    allocations = do.call(rbind, lapply(out, `[[`, "allocations")),
    borrowing = do.call(rbind, lapply(out, `[[`, "borrowing"))
  )
}

combine_simulation_outputs <- function(outputs) {
  list(
    trajectory = do.call(rbind, lapply(outputs, `[[`, "trajectory")),
    recommendations = do.call(rbind, lapply(outputs, `[[`, "recommendations")),
    allocations = do.call(rbind, lapply(outputs, `[[`, "allocations")),
    borrowing = do.call(rbind, lapply(outputs, `[[`, "borrowing"))
  )
}

format_progress_time <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) return("--:--:--")
  seconds <- as.integer(round(seconds))
  hours <- seconds %/% 3600L
  minutes <- (seconds %% 3600L) %/% 60L
  seconds <- seconds %% 60L
  sprintf("%02d:%02d:%02d", hours, minutes, seconds)
}

report_chunk_start <- function(stage, chunk_id, total_chunks) {
  cat(sprintf(
    "[%s] %s: starting chunk %d/%d\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stage, chunk_id, total_chunks
  ))
  flush.console()
}

report_chunk_progress <- function(stage, completed, total_chunks, started_at, restored = FALSE) {
  elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
  eta <- if (completed > 0L) elapsed * (total_chunks - completed) / completed else NA_real_
  cat(sprintf(
    "[%s] %s: completed %d/%d (%.1f%%) | elapsed %s | ETA %s%s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stage,
    completed,
    total_chunks,
    100 * completed / total_chunks,
    format_progress_time(elapsed),
    format_progress_time(eta),
    if (isTRUE(restored)) " | restored from checkpoint" else ""
  ))
  flush.console()
}

chunk_checkpoint_path <- function(stage, n_mc, chunk_id) {
  checkpoint_dir <- file.path(
    DIR_CHECKPOINTS,
    sprintf("algorithm_aligned_v1_%s_n%d", stage, n_mc)
  )
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(checkpoint_dir, sprintf("chunk_%04d.rds", chunk_id))
}

save_chunk_checkpoint <- function(value, path) {
  temporary <- paste0(path, ".tmp")
  saveRDS(value, temporary, compress = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) {
    stop("Could not finalize checkpoint: ", path)
  }
  invisible(path)
}

run_primary_simulation <- function(n_mc = SIM_SETTINGS$n_mc, workers = 1, scenarios = 1:4) {
  tasks <- expand.grid(scenario = scenarios, replicate = seq_len(n_mc))
  task_fun <- function(i) {
    run_one_replicate(tasks$scenario[i], tasks$replicate[i])
  }
  chunk_size <- max(1, workers * 4)
  chunks <- split(seq_len(nrow(tasks)), ceiling(seq_along(seq_len(nrow(tasks))) / chunk_size))
  started_at <- Sys.time()
  checkpoint_stage <- sprintf("primary_w%d_s%s", workers, paste(scenarios, collapse = "-"))

  if (workers > 1) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, setwd, ROOT_DIR)
    parallel::clusterEvalQ(cl, {
      source(file.path("R", "04_simulate_methods.R"))
      NULL
    })
    outputs <- list()
    for (chunk_id in seq_along(chunks)) {
      checkpoint <- chunk_checkpoint_path(checkpoint_stage, n_mc, chunk_id)
      if (file.exists(checkpoint)) {
        chunk_output <- readRDS(checkpoint)
        restored <- TRUE
      } else {
        report_chunk_start("Primary simulation", chunk_id, length(chunks))
        chunk_output <- parallel::parLapplyLB(cl, chunks[[chunk_id]], task_fun)
        save_chunk_checkpoint(chunk_output, checkpoint)
        restored <- FALSE
      }
      outputs <- c(outputs, chunk_output)
      report_chunk_progress("Primary simulation", chunk_id, length(chunks), started_at, restored)
    }
  } else {
    outputs <- list()
    for (chunk_id in seq_along(chunks)) {
      checkpoint <- chunk_checkpoint_path(checkpoint_stage, n_mc, chunk_id)
      if (file.exists(checkpoint)) {
        chunk_output <- readRDS(checkpoint)
        restored <- TRUE
      } else {
        report_chunk_start("Primary simulation", chunk_id, length(chunks))
        chunk_output <- lapply(chunks[[chunk_id]], task_fun)
        save_chunk_checkpoint(chunk_output, checkpoint)
        restored <- FALSE
      }
      outputs <- c(outputs, chunk_output)
      report_chunk_progress("Primary simulation", chunk_id, length(chunks), started_at, restored)
    }
  }

  combined <- combine_simulation_outputs(outputs)
  write.csv(combined$trajectory, file.path(DIR_RESULTS, sprintf("simulation_trajectory_n%d.csv", n_mc)), row.names = FALSE)
  write.csv(combined$recommendations, file.path(DIR_RESULTS, sprintf("simulation_recommendations_n%d.csv", n_mc)), row.names = FALSE)
  write.csv(combined$allocations, file.path(DIR_RESULTS, sprintf("simulation_allocations_n%d.csv", n_mc)), row.names = FALSE)
  if (!is.null(combined$borrowing)) {
    write.csv(combined$borrowing, file.path(DIR_RESULTS, sprintf("simulation_borrowing_n%d.csv", n_mc)), row.names = FALSE)
  }
  combined
}

run_lambda_sensitivity <- function(n_mc = SIM_SETTINGS$n_mc, workers = 1, lambda_values = c(0, 0.5, 1)) {
  tasks <- expand.grid(lambda_c = lambda_values, replicate = seq_len(n_mc))
  task_fun <- function(i) {
    scenario <- SIM_SETTINGS$sensitivity_scenario
    initial <- make_initial_data(scenario, tasks$replicate[i])
    out <- run_method_once(
      scenario = scenario,
      replicate = tasks$replicate[i],
      method = "CA-AB-GP",
      initial = initial,
      lambda_c = tasks$lambda_c[i],
      c_new = scenario_cost_new(scenario),
      c_scr = scenario_c_scr(scenario),
      b_max = scenario_budget(scenario)
    )
    out$trajectory$lambda_c <- tasks$lambda_c[i]
    out$recommendations$lambda_c <- tasks$lambda_c[i]
    out$allocations$lambda_c <- tasks$lambda_c[i]
    if (!is.null(out$borrowing)) out$borrowing$lambda_c <- tasks$lambda_c[i]
    out
  }

  if (workers > 1) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, setwd, ROOT_DIR)
    parallel::clusterEvalQ(cl, {
      source(file.path("R", "04_simulate_methods.R"))
      NULL
    })
    outputs <- parallel::parLapplyLB(cl, seq_len(nrow(tasks)), task_fun)
  } else {
    outputs <- lapply(seq_len(nrow(tasks)), task_fun)
  }

  combined <- combine_simulation_outputs(outputs)
  write.csv(combined$trajectory, file.path(DIR_RESULTS, sprintf("lambda_sensitivity_trajectory_n%d.csv", n_mc)), row.names = FALSE)
  write.csv(combined$allocations, file.path(DIR_RESULTS, sprintf("lambda_sensitivity_allocations_n%d.csv", n_mc)), row.names = FALSE)
  combined
}

run_sensitivity_analyses <- function(
    n_mc = SIM_SETTINGS$n_mc,
    workers = 1,
    lambda_values = c(0, 0.5, 1),
    c_new_values = c(5, 10, 15, 20, 30),
    budget_multipliers = c(0.8, 1.0, 1.2),
    cohort_values = c(1, 2, 4)) {
  tasks <- rbind(
    data.frame(setting = "lambda_c", value = lambda_values, replicate = rep(seq_len(n_mc), each = length(lambda_values))),
    data.frame(setting = "c_new", value = c_new_values, replicate = rep(seq_len(n_mc), each = length(c_new_values))),
    data.frame(setting = "budget_multiplier", value = budget_multipliers, replicate = rep(seq_len(n_mc), each = length(budget_multipliers))),
    data.frame(setting = "cohort_size", value = cohort_values, replicate = rep(seq_len(n_mc), each = length(cohort_values)))
  )
  tasks <- tasks[order(tasks$setting, tasks$replicate, tasks$value), ]
  row.names(tasks) <- NULL
  chunk_size <- max(1, workers * 4)
  chunks <- split(seq_len(nrow(tasks)), ceiling(seq_along(seq_len(nrow(tasks))) / chunk_size))
  started_at <- Sys.time()
  checkpoint_stage <- sprintf("sensitivity_w%d", workers)

  task_fun <- function(i) {
    setting <- tasks$setting[i]
    value <- tasks$value[i]
    replicate <- tasks$replicate[i]
    scenario <- SIM_SETTINGS$sensitivity_scenario

    lambda_c <- SIM_SETTINGS$lambda_c
    c_new <- scenario_cost_new(scenario)
    c_scr <- scenario_c_scr(scenario)
    b_max <- scenario_budget(scenario)
    n_max <- SIM_SETTINGS$n_max
    cohort_size <- SIM_SETTINGS$cohort_size

    if (setting == "lambda_c") lambda_c <- value
    if (setting == "c_new") c_new <- value
    if (setting == "budget_multiplier") b_max <- scenario_budget(scenario) * value
    if (setting == "cohort_size") cohort_size <- value

    initial <- make_initial_data(scenario, replicate)
    out <- run_method_once(
      scenario = scenario,
      replicate = replicate,
      method = "CA-AB-GP",
      initial = initial,
      lambda_c = lambda_c,
      c_new = c_new,
      c_scr = c_scr,
      b_max = b_max,
      cohort_size = cohort_size,
      n_max = n_max
    )

    for (slot in c("trajectory", "recommendations", "allocations", "borrowing")) {
      if (!is.null(out[[slot]])) {
        out[[slot]]$sensitivity <- setting
        out[[slot]]$sensitivity_value <- value
        out[[slot]]$lambda_c_setting <- lambda_c
        out[[slot]]$c_new_setting <- c_new
        out[[slot]]$b_max_setting <- b_max
        out[[slot]]$n_max_setting <- n_max
        out[[slot]]$cohort_size_setting <- cohort_size
      }
    }
    out
  }

  if (workers > 1) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, setwd, ROOT_DIR)
    parallel::clusterEvalQ(cl, {
      source(file.path("R", "04_simulate_methods.R"))
      NULL
    })
    outputs <- list()
    for (chunk_id in seq_along(chunks)) {
      checkpoint <- chunk_checkpoint_path(checkpoint_stage, n_mc, chunk_id)
      if (file.exists(checkpoint)) {
        chunk_output <- readRDS(checkpoint)
        restored <- TRUE
      } else {
        report_chunk_start("Sensitivity simulation", chunk_id, length(chunks))
        chunk_output <- parallel::parLapplyLB(cl, chunks[[chunk_id]], task_fun)
        save_chunk_checkpoint(chunk_output, checkpoint)
        restored <- FALSE
      }
      outputs <- c(outputs, chunk_output)
      report_chunk_progress("Sensitivity simulation", chunk_id, length(chunks), started_at, restored)
    }
  } else {
    outputs <- list()
    for (chunk_id in seq_along(chunks)) {
      checkpoint <- chunk_checkpoint_path(checkpoint_stage, n_mc, chunk_id)
      if (file.exists(checkpoint)) {
        chunk_output <- readRDS(checkpoint)
        restored <- TRUE
      } else {
        report_chunk_start("Sensitivity simulation", chunk_id, length(chunks))
        chunk_output <- lapply(chunks[[chunk_id]], task_fun)
        save_chunk_checkpoint(chunk_output, checkpoint)
        restored <- FALSE
      }
      outputs <- c(outputs, chunk_output)
      report_chunk_progress("Sensitivity simulation", chunk_id, length(chunks), started_at, restored)
    }
  }

  combined <- combine_simulation_outputs(outputs)
  write.csv(combined$trajectory, file.path(DIR_RESULTS, sprintf("sensitivity_trajectory_n%d.csv", n_mc)), row.names = FALSE)
  write.csv(combined$recommendations, file.path(DIR_RESULTS, sprintf("sensitivity_recommendations_n%d.csv", n_mc)), row.names = FALSE)
  write.csv(combined$allocations, file.path(DIR_RESULTS, sprintf("sensitivity_allocations_n%d.csv", n_mc)), row.names = FALSE)
  if (!is.null(combined$borrowing)) {
    write.csv(combined$borrowing, file.path(DIR_RESULTS, sprintf("sensitivity_borrowing_n%d.csv", n_mc)), row.names = FALSE)
  }
  combined
}

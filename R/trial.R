trial_cost <- function(data, design, dose_cols = design$dose_cols, stratum_col = "stratum") {
  validate_design(design)
  data <- as.data.frame(data)
  manufactured <- unique_dose_count(data, dose_cols)
  z <- as.integer(data[[stratum_col]])
  screening <- sum(vapply(z, function(k) 1 / max(design$prevalence[k], design$pi_min), numeric(1)))
  design$cost_patient * nrow(data) + design$cost_novel * manufactured + design$cost_screening * screening
}

update_trial_data <- function(data, allocation, outcomes, design, dose_cols = design$dose_cols,
                              stratum_col = "stratum", outcome_col = "y") {
  data <- as.data.frame(data)
  allocation <- as.data.frame(allocation)
  outcomes <- as.numeric(outcomes)
  dose <- allocation[1, dose_cols, drop = FALSE]
  new_rows <- data.frame(dose[rep(1, length(outcomes)), , drop = FALSE])
  new_rows[[stratum_col]] <- as.integer(allocation$stratum[1])
  new_rows[[outcome_col]] <- outcomes
  names(new_rows)[seq_along(dose_cols)] <- dose_cols
  rbind(data, new_rows[, names(data), drop = FALSE])
}

recommend_dose <- function(fit, design = fit$design) {
  pred <- predict_caabgp(fit, design$dose_grid)
  rows <- vector("list", design$n_strata)
  for (k in seq_len(design$n_strata)) {
    pk <- pred[pred$stratum == k, ]
    idx <- which.min(pk$mu)
    rows[[k]] <- data.frame(
      stratum = k,
      pk[idx, design$dose_cols, drop = FALSE],
      mu = pk$mu[idx],
      sd = pk$sd[idx]
    )
  }
  rec <- do.call(rbind, rows)
  row.names(rec) <- NULL
  rec
}

run_caabgp_trial <- function(initial_data, design, outcome_generator,
                             dose_cols = design$dose_cols, stratum_col = "stratum",
                             outcome_col = "y", seed = NULL, maxit = 60) {
  validate_design(design)
  if (!is.null(seed)) set.seed(seed)
  data <- as.data.frame(initial_data)
  total_cost <- trial_cost(data, design, dose_cols, stratum_col)
  manufactured <- manufactured_keys_from_data(data, dose_cols)
  previous <- NULL
  allocations <- list()
  trajectory <- list()
  step <- 1L
  repeat {
    fit <- fit_caabgp(data, design, dose_cols, stratum_col, outcome_col, previous = previous, maxit = maxit)
    previous <- fit
    rec <- recommend_dose(fit, design)
    trajectory[[step]] <- data.frame(
      step = step,
      n = nrow(data),
      total_cost = total_cost,
      unique_doses = length(manufactured),
      borrowing_index = borrowing_index(fit)$rho[1]
    )
    if (nrow(data) >= design$n_max || total_cost >= design$budget_max) break
    remaining <- design$budget_max - total_cost
    cohort <- min(design$cohort_size, design$n_max - nrow(data))
    action <- NULL
    for (r_try in seq(from = cohort, to = 1L, by = -1L)) {
      action <- suggest_next(fit, design, manufactured = manufactured, cohort_size = r_try, remaining_budget = remaining)
      if (!is.null(action)) {
        cohort <- r_try
        break
      }
    }
    if (is.null(action)) break
    dose <- as.numeric(action[1, dose_cols])
    k <- as.integer(action$stratum[1])
    y_new <- outcome_generator(stratum = k, dose = dose, n = cohort)
    data <- update_trial_data(data, action, y_new, design, dose_cols, stratum_col, outcome_col)
    total_cost <- total_cost + action$cost[1]
    manufactured <- unique(c(manufactured, dose_key(matrix(dose, nrow = 1))))
    allocations[[step]] <- data.frame(
      step = step,
      stratum = k,
      action[1, dose_cols, drop = FALSE],
      cohort_size = cohort,
      new_dose = action$new_dose[1],
      acquisition = action$score[1],
      cost_increment = action$cost[1],
      total_cost = total_cost,
      n = nrow(data),
      unique_doses = length(manufactured)
    )
    step <- step + 1L
  }
  final_fit <- fit_caabgp(data, design, dose_cols, stratum_col, outcome_col, previous = previous, maxit = maxit)
  result <- list(
    method = "CA-AB-GP",
    design = design,
    data = data,
    fit = final_fit,
    predictions = predict_caabgp(final_fit, design$dose_grid),
    recommendations = recommend_dose(final_fit, design),
    allocations = if (length(allocations)) do.call(rbind, allocations) else data.frame(),
    trajectory = do.call(rbind, trajectory),
    total_cost = total_cost
  )
  class(result) <- "caabgp_trial"
  result
}

make_initial_data <- function(design, outcome_generator, n_initial = NULL, seed = NULL,
                              dose_cols = design$dose_cols, stratum_col = "stratum",
                              outcome_col = "y") {
  validate_design(design)
  if (!is.null(seed)) set.seed(seed)
  n_initial <- null_coalesce(n_initial, max(2 * design$n_strata, nrow(design$initial_design)))
  counts <- pmax(1L, floor(n_initial * design$weights))
  while (sum(counts) < n_initial) counts[which.max(n_initial * design$weights - counts)] <- counts[which.max(n_initial * design$weights - counts)] + 1L
  rows <- list()
  id <- 1L
  for (k in seq_len(design$n_strata)) {
    dose_idx <- rep(seq_len(nrow(design$initial_design)), length.out = counts[k])
    for (j in dose_idx) {
      dose <- as.numeric(design$initial_design[j, dose_cols])
      y <- outcome_generator(stratum = k, dose = dose, n = 1)
      row <- data.frame(as.list(dose))
      names(row) <- dose_cols
      row[[stratum_col]] <- k
      row[[outcome_col]] <- y
      rows[[id]] <- row
      id <- id + 1L
    }
  }
  do.call(rbind, rows)
}

normal_ei_min <- function(mu, sd, f_best) {
  sd <- pmax(sd, 1e-12)
  z <- (f_best - mu) / sd
  ei <- (f_best - mu) * stats::pnorm(z) + sd * stats::dnorm(z)
  ei[!is.finite(ei)] <- 0
  pmax(ei, 0)
}

aei_min <- function(mu, sd, sigma, xi = 1) {
  sd <- pmax(sd, 1e-12)
  effective_idx <- which.min(mu + xi * sd)
  f_best <- mu[effective_idx]
  ei <- normal_ei_min(mu, sd, f_best)
  penalty <- 1 - sigma / sqrt(sd^2 + sigma^2)
  pmax(ei * pmax(penalty, 0), 0)
}

manufactured_keys_from_data <- function(data, dose_cols) {
  unique(dose_key(as_dose_matrix(data, dose_cols)))
}

is_manufactured <- function(dose, manufactured_keys) {
  dose_key(matrix(as.numeric(dose), nrow = 1)) %in% manufactured_keys
}

evaluation_cost <- function(design, stratum, dose, cohort_size, manufactured_keys) {
  design$cost_patient * cohort_size +
    design$cost_novel * as.numeric(!is_manufactured(dose, manufactured_keys)) +
    design$cost_screening * cohort_size / max(design$prevalence[stratum], design$pi_min)
}

residual_sd_by_stratum <- function(fit) {
  rep(sigma_hat(fit$fit), fit$design$n_strata)
}

suggest_next <- function(fit, design = fit$design, manufactured = NULL, cohort_size = design$cohort_size,
                         remaining_budget = Inf, active_strata = seq_len(design$n_strata)) {
  if (!inherits(fit, "caabgp_fit")) caabgp_stop("fit must be a caabgp_fit object.")
  validate_design(design)
  pred <- predict_caabgp(fit, design$dose_grid)
  grid <- as_dose_matrix(design$dose_grid, design$dose_cols)
  if (is.null(manufactured)) manufactured <- manufactured_keys_from_data(fit$data, fit$dose_cols)
  sigma_k <- residual_sd_by_stratum(fit)
  lambda <- design$lambda_c
  rows <- list()
  id <- 1L
  for (k in active_strata) {
    pk <- pred[pred$stratum == k, ]
    aei <- aei_min(pk$mu, pk$sd, sigma_k[k], xi = design$xi)
    costs <- vapply(seq_len(nrow(grid)), function(i) evaluation_cost(design, k, grid[i, ], cohort_size, manufactured), numeric(1))
    feasible <- costs <= remaining_budget + 1e-8
    score <- design$weights[k] * aei / (costs^lambda)
    score[!feasible] <- -Inf
    rows[[id]] <- data.frame(
      stratum = k,
      design$dose_grid,
      mu = pk$mu,
      sd = pk$sd,
      aei = aei,
      cost = costs,
      score = score,
      new_dose = !vapply(seq_len(nrow(grid)), function(i) is_manufactured(grid[i, ], manufactured), logical(1)),
      feasible = feasible
    )
    id <- id + 1L
  }
  acq <- do.call(rbind, rows)
  if (!any(acq$feasible)) return(NULL)
  if (all(!is.finite(acq$score)) || all(acq$score <= 0, na.rm = TRUE)) {
    acq$score <- -acq$mu
    acq$score[!acq$feasible] <- -Inf
  }
  action <- acq[which.max(acq$score), , drop = FALSE]
  row.names(action) <- NULL
  action
}

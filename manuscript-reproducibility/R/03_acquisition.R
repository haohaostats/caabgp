source(file.path("R", "00_settings.R"))

normal_ei_min <- function(mu, sd, f_best) {
  sd <- pmax(sd, 1e-12)
  z <- (f_best - mu) / sd
  ei <- (f_best - mu) * stats::pnorm(z) + sd * stats::dnorm(z)
  ei[!is.finite(ei)] <- 0
  pmax(ei, 0)
}

aei_min <- function(mu, sd, sigma, xi = SIM_SETTINGS$xi) {
  sd <- pmax(sd, 1e-12)
  effective_idx <- which.min(mu + xi * sd)
  f_best <- mu[effective_idx]
  ei <- normal_ei_min(mu, sd, f_best)
  penalty <- 1 - sigma / sqrt(sd^2 + sigma^2)
  pmax(ei * pmax(penalty, 0), 0)
}

dose_key <- function(d1, d2) {
  paste(sprintf("%.2f", d1), sprintf("%.2f", d2), sep = "_")
}

dose_row_index <- function(grid, d) {
  which(abs(grid[, 1] - d[1]) < 1e-12 & abs(grid[, 2] - d[2]) < 1e-12)[1]
}

is_manufactured <- function(d, manufactured_keys) {
  dose_key(d[1], d[2]) %in% manufactured_keys
}

evaluation_cost <- function(k, d, r, manufactured_keys, prevalence, c_new, c_scr = SIM_SETTINGS$c_scr) {
  SIM_SETTINGS$c_pt * r +
    c_new * as.numeric(!is_manufactured(d, manufactured_keys)) +
    c_scr * r / max(prevalence[k], SIM_SETTINGS$pi_min)
}

personalized_acquisition_table <- function(
    pred,
    grid,
    manufactured_keys,
    prevalence,
    sigma_by_k,
    c_new,
    c_scr = SIM_SETTINGS$c_scr,
    lambda_c = 0,
    weights = prevalence,
    cohort_size = SIM_SETTINGS$cohort_size,
    remaining_budget = Inf,
    active_strata = seq_along(prevalence),
    xi = SIM_SETTINGS$xi) {
  K <- length(prevalence)
  active_strata <- sort(unique(as.integer(active_strata)))
  active_strata <- active_strata[active_strata >= 1L & active_strata <= K]
  if (!length(active_strata)) return(data.frame())
  rows <- list()
  idx <- 1

  for (k in active_strata) {
    pred_k <- pred[pred$stratum == k, ]
    aei <- aei_min(pred_k$mu, pred_k$sd, sigma_by_k[k], xi = xi)
    costs <- vapply(
      seq_len(nrow(grid)),
      function(i) evaluation_cost(k, grid[i, ], cohort_size, manufactured_keys, prevalence, c_new, c_scr),
      numeric(1)
    )
    feasible <- costs <= remaining_budget + 1e-8
    score <- weights[k] * aei / (costs^lambda_c)
    score[!feasible] <- -Inf
    rows[[idx]] <- data.frame(
      stratum = k,
      d1 = grid[, 1],
      d2 = grid[, 2],
      aei = aei,
      cost = costs,
      score = score,
      feasible = feasible
    )
    idx <- idx + 1
  }

  do.call(rbind, rows)
}

select_personalized_action <- function(
    pred,
    grid,
    manufactured_keys,
    prevalence,
    sigma_by_k,
    c_new,
    c_scr = SIM_SETTINGS$c_scr,
    lambda_c = 0,
    weights = prevalence,
    cohort_size = SIM_SETTINGS$cohort_size,
    remaining_budget = Inf,
    active_strata = seq_along(prevalence),
    xi = SIM_SETTINGS$xi) {
  acq <- personalized_acquisition_table(
    pred = pred,
    grid = grid,
    manufactured_keys = manufactured_keys,
    prevalence = prevalence,
    sigma_by_k = sigma_by_k,
    c_new = c_new,
    c_scr = c_scr,
    lambda_c = lambda_c,
    weights = weights,
    cohort_size = cohort_size,
    remaining_budget = remaining_budget,
    active_strata = active_strata,
    xi = xi
  )
  if (!nrow(acq)) return(NULL)
  if (!any(acq$feasible)) return(NULL)
  if (all(acq$score <= 0 | !is.finite(acq$score))) {
    acq$score <- -Inf
    feasible <- acq$feasible
    acq$score <- -pred$mu[match(
      paste(acq$stratum, acq$d1, acq$d2),
      paste(pred$stratum, pred$d1, pred$d2)
    )]
    acq$score[!feasible] <- -Inf
  }
  acq[which.max(acq$score), ]
}

select_pooled_action <- function(pred_one, grid, sigma, xi = SIM_SETTINGS$xi) {
  aei <- aei_min(pred_one$mu, pred_one$sd, sigma, xi = xi)
  if (all(aei <= 0) || all(!is.finite(aei))) {
    idx <- which.min(pred_one$mu)
  } else {
    idx <- which.max(aei)
  }
  data.frame(
    stratum = NA_integer_,
    d1 = grid[idx, 1],
    d2 = grid[idx, 2],
    aei = aei[idx],
    cost = NA_real_,
    score = aei[idx]
  )
}

select_pooled_feasible_action <- function(
    pred_one,
    grid,
    sigma,
    manufactured_keys,
    prevalence,
    c_new,
    c_scr = SIM_SETTINGS$c_scr,
    cohort_size = SIM_SETTINGS$cohort_size,
    remaining_budget = Inf,
    xi = SIM_SETTINGS$xi) {
  K <- length(prevalence)
  aei <- aei_min(pred_one$mu, pred_one$sd, sigma, xi = xi)
  cost_matrix <- vapply(
    seq_len(nrow(grid)),
    function(i) {
      vapply(
        seq_len(K),
        function(k) evaluation_cost(
          k, grid[i, ], cohort_size, manufactured_keys,
          prevalence, c_new, c_scr
        ),
        numeric(1)
      )
    },
    numeric(K)
  )
  feasible_matrix <- cost_matrix <= remaining_budget + 1e-8
  feasible_dose <- colSums(feasible_matrix) > 0L
  if (!any(feasible_dose)) return(NULL)

  score <- aei
  score[!feasible_dose] <- -Inf
  if (all(score <= 0 | !is.finite(score))) {
    score <- -pred_one$mu
    score[!feasible_dose] <- -Inf
  }
  idx <- which.max(score)
  data.frame(
    stratum = NA_integer_,
    d1 = grid[idx, 1],
    d2 = grid[idx, 2],
    aei = aei[idx],
    cost = NA_real_,
    score = score[idx],
    feasible_strata = I(list(which(feasible_matrix[, idx])))
  )
}

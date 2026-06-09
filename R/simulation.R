run_simulation_study <- function(scenarios, design, n_mc = 100, n_initial = NULL, seed = 1,
                                 maxit = 40) {
  validate_design(design)
  if (!is.list(scenarios)) caabgp_stop("scenarios must be a named list of outcome-generator functions.")
  all_trials <- list()
  id <- 1L
  for (sname in names(scenarios)) {
    generator <- scenarios[[sname]]
    for (m in seq_len(n_mc)) {
      init <- make_initial_data(design, generator, n_initial = n_initial, seed = seed + m)
      tr <- run_caabgp_trial(init, design, generator, seed = seed + 10000 * m, maxit = maxit)
      rec <- tr$recommendations
      rec$scenario <- sname
      rec$replicate <- m
      rec$method <- "CA-AB-GP"
      rec$n <- nrow(tr$data)
      rec$total_cost <- tr$total_cost
      rec$unique_doses <- unique_dose_count(tr$data, design$dose_cols)
      all_trials[[id]] <- rec
      id <- id + 1L
    }
  }
  out <- do.call(rbind, all_trials)
  row.names(out) <- NULL
  out
}

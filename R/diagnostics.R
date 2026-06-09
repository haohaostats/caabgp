borrowing_index <- function(x, ...) UseMethod("borrowing_index")

borrowing_index.caabgp_fit <- function(x, ...) {
  p <- x$fit$params
  rho <- p$nu0 / (p$nu0 + p$nub)
  data.frame(stratum = seq_len(x$design$n_strata), rho = rep(rho, x$design$n_strata))
}

borrowing_index.caabgp_trial <- function(x, ...) {
  borrowing_index(x$fit)
}

export_results <- function(result, path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  if (inherits(result, "caabgp_trial")) {
    utils::write.csv(result$data, file.path(path, "trial_data.csv"), row.names = FALSE)
    utils::write.csv(result$predictions, file.path(path, "predictions.csv"), row.names = FALSE)
    utils::write.csv(result$recommendations, file.path(path, "recommendations.csv"), row.names = FALSE)
    utils::write.csv(result$allocations, file.path(path, "allocations.csv"), row.names = FALSE)
    utils::write.csv(result$trajectory, file.path(path, "trajectory.csv"), row.names = FALSE)
  } else if (is.list(result)) {
    for (nm in names(result)) {
      if (is.data.frame(result[[nm]])) utils::write.csv(result[[nm]], file.path(path, paste0(nm, ".csv")), row.names = FALSE)
    }
  } else {
    caabgp_stop("result must be a caabgp_trial or a named list of data frames.")
  }
  invisible(normalizePath(path, winslash = "/", mustWork = FALSE))
}

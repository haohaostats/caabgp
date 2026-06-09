dose_grid <- function(levels = seq(0, 1, by = 0.25), n_agents = 2, names = NULL) {
  check_numeric_vector(levels, "levels")
  if (length(levels) < 2) caabgp_stop("levels must contain at least two dose levels.")
  if (!is.numeric(n_agents) || length(n_agents) != 1 || n_agents < 1) caabgp_stop("n_agents must be a positive integer.")
  n_agents <- as.integer(n_agents)
  if (is.null(names)) names <- paste0("d", seq_len(n_agents))
  if (length(names) != n_agents) caabgp_stop("names must have length n_agents.")
  grid <- expand.grid(rep(list(sort(unique(levels))), n_agents), KEEP.OUT.ATTRS = FALSE)
  names(grid) <- names
  grid
}

caabgp_design <- function(
    dose_grid,
    n_strata,
    prevalence = rep(1 / n_strata, n_strata),
    weights = prevalence,
    cohort_size = 1,
    n_max = 80,
    budget_max = Inf,
    cost_patient = 1,
    cost_novel = 0,
    cost_screening = 0,
    lambda_c = 0.5,
    xi = 1,
    pi_min = 0.02,
    initial_design = NULL,
    dose_cols = NULL) {
  grid <- as.data.frame(dose_grid)
  if (is.null(dose_cols)) dose_cols <- names(grid)
  X <- as_dose_matrix(grid, dose_cols)
  if (!is.numeric(n_strata) || length(n_strata) != 1 || n_strata < 1) caabgp_stop("n_strata must be a positive integer.")
  n_strata <- as.integer(n_strata)
  prevalence <- check_numeric_vector(prevalence, "prevalence", n_strata, positive = TRUE)
  weights <- check_numeric_vector(weights, "weights", n_strata, nonnegative = TRUE)
  if (sum(weights) <= 0) caabgp_stop("At least one stratum weight must be positive.")
  weights <- weights / sum(weights)
  prevalence <- prevalence / sum(prevalence)
  if (is.null(initial_design)) {
    initial_design <- grid[seq_len(min(nrow(grid), max(1, ncol(X) + 3))), dose_cols, drop = FALSE]
  }
  initial_design <- as.data.frame(initial_design)
  names(initial_design) <- dose_cols

  design <- list(
    dose_grid = grid[, dose_cols, drop = FALSE],
    dose_cols = dose_cols,
    n_agents = ncol(X),
    n_strata = n_strata,
    prevalence = prevalence,
    weights = weights,
    cohort_size = as.integer(cohort_size),
    n_max = as.integer(n_max),
    budget_max = budget_max,
    cost_patient = cost_patient,
    cost_novel = cost_novel,
    cost_screening = cost_screening,
    lambda_c = lambda_c,
    xi = xi,
    pi_min = pi_min,
    initial_design = initial_design
  )
  class(design) <- "caabgp_design"
  validate_design(design)
  design
}

validate_design <- function(design) {
  if (!inherits(design, "caabgp_design")) caabgp_stop("design must be a caabgp_design object.")
  if (design$cohort_size < 1) caabgp_stop("cohort_size must be at least 1.")
  if (design$n_max < 1) caabgp_stop("n_max must be positive.")
  check_numeric_vector(c(design$cost_patient, design$cost_novel, design$cost_screening), "costs", nonnegative = TRUE)
  check_numeric_vector(design$lambda_c, "lambda_c", n = 1, nonnegative = TRUE)
  invisible(TRUE)
}

print.caabgp_design <- function(x, ...) {
  cat("Cost-aware adaptive-borrowing GP design\n")
  cat("  dose grid:", nrow(x$dose_grid), "combinations x", x$n_agents, "agents\n")
  cat("  strata:", x$n_strata, "\n")
  cat("  cohort size:", x$cohort_size, "\n")
  cat("  n_max:", x$n_max, " budget_max:", x$budget_max, "\n")
  invisible(x)
}

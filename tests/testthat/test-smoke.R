test_that("minimal CA-AB-GP trial runs", {
  grid <- dose_grid(seq(0, 1, by = 0.5), n_agents = 2)
  design <- caabgp_design(
    dose_grid = grid,
    n_strata = 2,
    prevalence = c(0.5, 0.5),
    weights = c(0.5, 0.5),
    cohort_size = 1,
    n_max = 12,
    budget_max = 100,
    cost_patient = 1,
    cost_novel = 1,
    initial_design = grid[c(1, 3, 5, 7, 9), ]
  )
  truth <- function(stratum, dose, n = 1) {
    center <- if (stratum == 1) c(0, 1) else c(1, 0)
    stats::rnorm(n, -exp(-sum((dose - center)^2) / 0.2), 0.1)
  }
  initial <- make_initial_data(design, truth, n_initial = 6, seed = 1)
  trial <- run_caabgp_trial(initial, design, truth, seed = 2, maxit = 15)
  expect_s3_class(trial, "caabgp_trial")
  expect_equal(nrow(recommend_dose(trial$fit)), 2)
  expect_s3_class(plot_allocation(trial), "ggplot")
})

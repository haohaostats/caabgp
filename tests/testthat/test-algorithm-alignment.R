make_alignment_fixture <- function(
    cohort_size = 2,
    n_max = 20,
    budget_max = 100,
    stopping_threshold = 0,
    stopping_patience = NULL) {
  grid <- dose_grid(seq(0, 1, by = 0.5), n_agents = 2)
  design <- caabgp_design(
    dose_grid = grid,
    n_strata = 2,
    prevalence = c(0.5, 0.5),
    weights = c(0.5, 0.5),
    cohort_size = cohort_size,
    n_max = n_max,
    budget_max = budget_max,
    cost_patient = 1,
    cost_novel = 0,
    cost_screening = 0,
    lambda_c = 0.5,
    stopping_threshold = stopping_threshold,
    stopping_patience = stopping_patience,
    initial_design = grid[c(1, 3, 5, 7, 9), ]
  )
  initial <- data.frame(
    d1 = c(0, 0.5, 1, 0, 0.5, 1),
    d2 = c(0, 0.5, 1, 1, 0.5, 0),
    stratum = c(1, 1, 1, 2, 2, 2),
    y = c(-0.10, -0.40, -0.20, -0.15, -0.35, -0.25)
  )
  generator <- function(stratum, dose, n = 1) rep(-0.3, n)
  list(design = design, initial = initial, generator = generator)
}

test_that("the final cohort is never reduced", {
  fixture <- make_alignment_fixture(budget_max = 7)
  trial <- run_caabgp_trial(
    fixture$initial, fixture$design, fixture$generator, seed = 1, maxit = 8
  )

  expect_equal(nrow(trial$data), nrow(fixture$initial))
  expect_equal(nrow(trial$allocations), 0)
  expect_lte(trial$total_cost, fixture$design$budget_max)
})

test_that("sample-size feasibility is checked for the full cohort", {
  fixture <- make_alignment_fixture(n_max = 7)
  trial <- run_caabgp_trial(
    fixture$initial, fixture$design, fixture$generator, seed = 1, maxit = 8
  )

  expect_equal(nrow(trial$data), nrow(fixture$initial))
  expect_equal(nrow(trial$allocations), 0)
})

test_that("all active strata counters update after every completed cohort", {
  fixture <- make_alignment_fixture(
    stopping_threshold = 1e6,
    stopping_patience = 1L
  )
  trial <- run_caabgp_trial(
    fixture$initial, fixture$design, fixture$generator, seed = 1, maxit = 8
  )

  expect_equal(trial$allocations$cohort_size, fixture$design$cohort_size)
  expect_equal(nrow(trial$allocations), 1)
  expect_equal(trial$stopping_counters, c(1L, 1L))
  expect_length(trial$active_strata, 0)
  expect_equal(sort(unique(trial$recommendations$stratum)), c(1L, 2L))
})

test_that("the default stopping patience is J plus one", {
  fixture <- make_alignment_fixture()
  expect_equal(fixture$design$stopping_patience, 3L)
})

test_that("initial data include at least two patients per stratum", {
  grid <- dose_grid(seq(0, 1, by = 0.5), n_agents = 2)
  design <- caabgp_design(
    grid, n_strata = 3, weights = c(0.8, 0.1, 0.1),
    initial_design = grid[c(1, 3, 5, 7, 9), ]
  )
  generator <- function(stratum, dose, n = 1) rep(stratum, n)

  initial <- make_initial_data(design, generator, n_initial = 7, seed = 1)
  expect_gte(min(tabulate(initial$stratum, nbins = 3)), 2L)
  expect_error(
    make_initial_data(design, generator, n_initial = 5, seed = 1),
    "at least two patients per stratum"
  )
})

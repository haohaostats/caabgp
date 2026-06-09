library(caabgp)

grid <- dose_grid(seq(0, 1, by = 0.25), n_agents = 2)
initial_design <- grid[c(1, 9, 13, 17, 25), ]

design <- caabgp_design(
  dose_grid = grid,
  n_strata = 2,
  prevalence = c(0.55, 0.45),
  weights = c(0.50, 0.50),
  cohort_size = 2,
  n_max = 80,
  budget_max = 220,
  cost_patient = 1,
  cost_novel = 5,
  cost_screening = 0,
  lambda_c = 0.5,
  initial_design = initial_design
)

truth <- function(stratum, dose, n = 1) {
  dose <- as.numeric(dose)
  shared <- -0.95 * exp(-sum((dose - c(0.50, 0.50))^2) / (2 * 0.42^2))
  deviation <- if (stratum == 1) {
    -0.45 * exp(-sum((dose - c(0.25, 0.75))^2) / (2 * 0.16^2))
  } else {
    -0.45 * exp(-sum((dose - c(0.75, 0.25))^2) / (2 * 0.16^2))
  }
  rnorm(n, shared + deviation, 0.05)
}

initial <- make_initial_data(design, truth, n_initial = 20, seed = 21)
trial <- run_caabgp_trial(initial, design, truth, seed = 22)
print(recommend_dose(trial$fit))
print(borrowing_index(trial))

p <- plot_allocation(trial)
save_caabgp_figure(p, tempfile(fileext = ".pdf"), width = 6.5, height = 5)

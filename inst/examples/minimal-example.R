library(caabgp)

grid <- dose_grid(seq(0, 1, by = 0.25), n_agents = 2)
design <- caabgp_design(
  dose_grid = grid,
  n_strata = 2,
  prevalence = c(0.65, 0.35),
  weights = c(0.55, 0.45),
  cohort_size = 2,
  n_max = 40,
  budget_max = 180,
  cost_patient = 1,
  cost_novel = 10,
  cost_screening = 1,
  lambda_c = 0.5,
  initial_design = grid[c(1, 5, 13, 21, 25), ]
)

truth <- function(stratum, dose, n = 1) {
  center <- if (stratum == 1) c(0.25, 0.75) else c(0.75, 0.25)
  mu <- -exp(-sum((dose - center)^2) / (2 * 0.18^2))
  rnorm(n, mu, 0.25)
}

initial <- make_initial_data(design, truth, n_initial = 12, seed = 1)
trial <- run_caabgp_trial(initial, design, truth, seed = 2)
print(recommend_dose(trial$fit))
print(borrowing_index(trial))

p <- plot_allocation(trial)
save_caabgp_figure(p, tempfile(fileext = ".pdf"), width = 6.5, height = 5)

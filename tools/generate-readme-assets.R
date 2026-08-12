## Rebuild the deterministic README showcase from the package quick start.
pkg_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
for (f in list.files(file.path(pkg_root, "R"), pattern = "[.]R$", full.names = TRUE)) {
  sys.source(f, envir = .GlobalEnv)
}

grid <- dose_grid(seq(0, 1, by = 0.25), n_agents = 2)
design <- caabgp_design(
  dose_grid = grid,
  n_strata = 2,
  prevalence = c(0.60, 0.40),
  cohort_size = 2,
  n_max = 40,
  budget_max = 100,
  cost_novel = 5
)

truth <- function(stratum, dose, n = 1) {
  target <- if (stratum == 1) c(0.25, 0.75) else c(0.75, 0.25)
  stats::rnorm(n, -exp(-sum((as.numeric(dose) - target)^2) / 0.15), 0.10)
}

initial <- make_initial_data(design, truth, n_initial = 10, seed = 21)
trial <- run_caabgp_trial(initial, design, truth, seed = 22)

asset_dir <- file.path(pkg_root, "man", "figures")
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

recommendations <- recommend_dose(trial$fit)
rho <- borrowing_index(trial)$rho[1]
showcase <- transform(
  recommendations,
  final_sample_size = nrow(trial$data),
  unique_doses = unique_dose_count(trial$data, design$dose_cols),
  repeat_allocation = if (nrow(trial$allocations)) mean(!trial$allocations$new_dose) else NA_real_,
  total_cost = trial$total_cost,
  borrowing_index = rho
)
utils::write.csv(
  showcase,
  file.path(asset_dir, "quick-start-results.csv"),
  row.names = FALSE
)

p_allocation <- plot_allocation(
  trial,
  main = "A  Sequential allocation",
  base_size = 13
) + ggplot2::theme(legend.position = "bottom")
p_surface_1 <- plot_surface(
  trial,
  stratum = 1,
  main = "B  Stratum 1",
  base_size = 13
) + ggplot2::theme(legend.position = "none")
p_surface_2 <- plot_surface(
  trial,
  stratum = 2,
  main = "C  Stratum 2",
  base_size = 13
) + ggplot2::theme(legend.position = "none")

ggplot2::ggsave(
  file.path(asset_dir, "caabgp-allocation.png"),
  p_allocation,
  width = 6.2,
  height = 5.2,
  dpi = 200,
  bg = "white"
)
ggplot2::ggsave(
  file.path(asset_dir, "caabgp-surface-stratum-1.png"),
  p_surface_1,
  width = 5.4,
  height = 4.8,
  dpi = 200,
  bg = "white"
)
ggplot2::ggsave(
  file.path(asset_dir, "caabgp-surface-stratum-2.png"),
  p_surface_2,
  width = 5.4,
  height = 4.8,
  dpi = 200,
  bg = "white"
)

print(showcase)

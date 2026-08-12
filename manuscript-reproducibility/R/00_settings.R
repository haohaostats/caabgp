## Settings for the empirical-Bayes CA-AB-GP simulation pipeline.
## This file matches the main-text method and Supplementary Table S1:
## prior-free empirical Bayes, squared-exponential kernels, and cost-aware AEI.

ROOT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

DIR_IMAGES <- file.path(ROOT_DIR, "images")
DIR_FIGURES <- file.path(DIR_IMAGES, "Figures")
DIR_RESULTS <- file.path(ROOT_DIR, "results_r")
DIR_CHECKPOINTS <- file.path(ROOT_DIR, "checkpoints")
dir.create(DIR_IMAGES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_CHECKPOINTS, showWarnings = FALSE, recursive = TRUE)

METHODS <- c("S-GP", "P-GP", "IND-GP", "AB-GP", "CA-AB-GP")

SIM_SETTINGS <- list(
  n_mc = 1000,
  n_max = 120,
  n_initial = 20,
  cohort_size = 2,
  xi = 1,
  dose_values = c(0, 0.25, 0.50, 0.75, 1.00),
  initial_design = rbind(
    c(0, 0),
    c(0, 1),
    c(0.5, 0.5),
    c(1, 0),
    c(1, 1)
  ),
  c_pt = 1,
  c_scr = 1,
  pi_min = 0.02,
  lambda_c = 0.5,
  stopping_threshold = 0,
  stopping_patience = 3L,
  sensitivity_scenario = 3,
  seed_base = 20260501
)

candidate_grid <- function() {
  vals <- SIM_SETTINGS$dose_values
  as.matrix(expand.grid(d1 = vals, d2 = vals))
}

scenario_prevalence <- function(scenario) {
  if (scenario == 1) return(c(0.50, 0.50))
  if (scenario == 2) return(c(0.60, 0.30, 0.10))
  if (scenario == 3) return(c(0.50, 0.30, 0.15, 0.05))
  if (scenario == 4) return(c(0.50, 0.50))
  stop("Unknown scenario: ", scenario)
}

scenario_weights <- function(scenario) {
  if (scenario == 1) return(c(0.50, 0.50))
  if (scenario == 2) return(c(0.40, 0.30, 0.30))
  if (scenario == 3) return(c(0.35, 0.25, 0.20, 0.20))
  if (scenario == 4) return(c(0.50, 0.50))
  stop("Unknown scenario: ", scenario)
}

scenario_cost_new <- function(scenario) {
  if (scenario == 1) return(10)
  if (scenario == 2) return(15)
  if (scenario == 3) return(20)
  if (scenario == 4) return(15)
  stop("Unknown scenario: ", scenario)
}

scenario_c_scr <- function(scenario) {
  if (scenario %in% c(1, 4)) return(0)
  if (scenario %in% c(2, 3)) return(SIM_SETTINGS$c_scr)
  stop("Unknown scenario: ", scenario)
}

scenario_budget <- function(scenario) {
  if (scenario == 1) return(450)
  if (scenario == 2) return(500)
  if (scenario == 3) return(550)
  if (scenario == 4) return(450)
  stop("Unknown scenario: ", scenario)
}

scenario_noise_sd <- function(scenario) {
  if (scenario == 1) return(c(0.50, 0.50))
  if (scenario == 2) return(c(0.45, 0.45, 0.65))
  if (scenario == 3) return(c(0.50, 0.50, 0.50, 0.50))
  if (scenario == 4) return(c(0.70, 0.70))
  stop("Unknown scenario: ", scenario)
}

initial_counts_by_stratum <- function(scenario) {
  weights <- scenario_weights(scenario)
  min_count <- 2L
  if (length(weights) * min_count > SIM_SETTINGS$n_initial) {
    stop("Initial sample size is too small to allocate at least two patients per stratum.")
  }

  counts <- rep(min_count, length(weights))
  remaining <- SIM_SETTINGS$n_initial - sum(counts)
  target <- remaining * weights / sum(weights)
  add <- floor(target)
  counts <- counts + add

  while (sum(counts) < SIM_SETTINGS$n_initial) {
    deficit <- target - add
    j <- which.max(deficit)
    counts[j] <- counts[j] + 1L
    add[j] <- add[j] + 1L
  }

  while (sum(counts) > SIM_SETTINGS$n_initial) {
    candidates <- which(counts > min_count)
    j <- candidates[which.max(counts[candidates] - SIM_SETTINGS$n_initial * weights[candidates])]
    counts[j] <- counts[j] - 1L
  }

  as.integer(counts)
}

method_order <- function(x) {
  factor(x, levels = METHODS)
}

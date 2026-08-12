source(file.path("R", "00_settings.R"))
source(file.path("R", "01_surfaces.R"))

stop_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

read_result <- function(stem, n_mc) {
  path <- file.path(DIR_RESULTS, sprintf("%s_n%d.csv", stem, n_mc))
  stop_if(!file.exists(path), paste("Missing result file:", path))
  utils::read.csv(path, check.names = FALSE)
}

validate_final_recommendations <- function(trajectory, recommendations, group_cols, expected_k) {
  final_n <- aggregate(
    trajectory$n,
    trajectory[group_cols],
    max
  )
  names(final_n)[ncol(final_n)] <- "final_n"
  tagged <- merge(recommendations, final_n, by = group_cols, all.x = TRUE)
  final <- tagged[tagged$n == tagged$final_n, , drop = FALSE]
  counts <- aggregate(final$stratum, final[group_cols], function(x) length(unique(x)))
  names(counts)[ncol(counts)] <- "n_final_strata"
  stop_if(any(counts$n_final_strata != expected_k),
          "At least one completed trial is missing an all-stratum final recommendation.")
}

validate_primary <- function(n_mc) {
  allocations <- read_result("simulation_allocations", n_mc)
  trajectory <- read_result("simulation_trajectory", n_mc)
  recommendations <- read_result("simulation_recommendations", n_mc)

  stop_if(any(allocations$cohort_size != allocations$planned_cohort_size),
          "Primary results contain a reduced terminal cohort.")
  stop_if(any(allocations$n_after > SIM_SETTINGS$n_max),
          "Primary results exceed N_max.")
  stop_if(any(allocations$total_cost > allocations$b_max + 1e-8),
          "Primary results exceed B_max.")

  completed <- unique(trajectory[c("scenario", "replicate", "method")])
  counts <- aggregate(completed$replicate, completed[c("scenario", "method")], length)
  names(counts)[3] <- "n_replicates"
  stop_if(any(counts$n_replicates != n_mc),
          "Primary results do not contain the requested number of replicates for every scenario-method pair.")

  for (scenario in sort(unique(trajectory$scenario))) {
    k <- length(scenario_prevalence(scenario))
    subset_traj <- trajectory[trajectory$scenario == scenario, , drop = FALSE]
    subset_rec <- recommendations[recommendations$scenario == scenario, , drop = FALSE]
    validate_final_recommendations(
      subset_traj,
      subset_rec,
      c("scenario", "replicate", "method"),
      k
    )
  }
  invisible(TRUE)
}

validate_sensitivity <- function(n_mc) {
  allocations <- read_result("sensitivity_allocations", n_mc)
  trajectory <- read_result("sensitivity_trajectory", n_mc)
  recommendations <- read_result("sensitivity_recommendations", n_mc)

  stop_if(any(allocations$cohort_size != allocations$planned_cohort_size),
          "Sensitivity results contain a reduced terminal cohort.")
  stop_if(any(allocations$n_after > allocations$n_max_setting),
          "Sensitivity results exceed N_max.")
  stop_if(any(allocations$total_cost > allocations$b_max_setting + 1e-8),
          "Sensitivity results exceed B_max.")

  group_cols <- c(
    "sensitivity", "sensitivity_value", "replicate", "scenario", "method",
    "lambda_c_setting", "c_new_setting", "b_max_setting", "n_max_setting",
    "cohort_size_setting"
  )
  settings <- unique(trajectory[c("sensitivity", "sensitivity_value", "replicate")])
  counts <- aggregate(settings$replicate, settings[c("sensitivity", "sensitivity_value")], length)
  names(counts)[3] <- "n_replicates"
  stop_if(any(counts$n_replicates != n_mc),
          "Sensitivity results do not contain the requested number of replicates for every setting.")

  validate_final_recommendations(
    trajectory,
    recommendations,
    group_cols,
    length(scenario_prevalence(SIM_SETTINGS$sensitivity_scenario))
  )
  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
n_mc <- if (length(args)) as.integer(args[1]) else SIM_SETTINGS$n_mc
validate_sensitivity_outputs <- if (length(args) >= 2) as.logical(args[2]) else TRUE

cat("Validating primary results (n_mc =", n_mc, ")...\n")
validate_primary(n_mc)
if (isTRUE(validate_sensitivity_outputs)) {
  cat("Validating sensitivity results...\n")
  validate_sensitivity(n_mc)
}
cat("Full result validation passed.\n")

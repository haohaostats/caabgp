source(file.path("R", "04_simulate_methods.R"))

validate_output <- function(out, scenario, planned_cohort, n_max, b_max) {
  allocations <- out$allocations
  if (!is.null(allocations) && nrow(allocations)) {
    stopifnot(all(allocations$cohort_size == planned_cohort))
    stopifnot(all(allocations$planned_cohort_size == planned_cohort))
    stopifnot(all(allocations$n_after <= n_max))
    stopifnot(all(allocations$total_cost <= b_max + 1e-8))
  }

  final_n <- max(out$recommendations$n)
  final_recommendations <- out$recommendations[
    out$recommendations$n == final_n,
    ,
    drop = FALSE
  ]
  stopifnot(
    length(unique(final_recommendations$stratum)) ==
      length(scenario_prevalence(scenario))
  )
  invisible(TRUE)
}

scenario <- 2L
replicate <- 1L
initial <- make_initial_data(scenario, replicate)

cat("Checking strict fixed-cohort feasibility...\n")
strict_output <- run_method_once(
  scenario = scenario,
  replicate = replicate,
  method = "CA-AB-GP",
  initial = initial,
  cohort_size = 2L,
  stopping_threshold = 0
)
validate_output(
  strict_output,
  scenario = scenario,
  planned_cohort = 2L,
  n_max = SIM_SETTINGS$n_max,
  b_max = scenario_budget(scenario)
)

cat("Checking active-stratum stopping and final all-stratum recommendations...\n")
stopped_output <- run_method_once(
  scenario = scenario,
  replicate = replicate,
  method = "CA-AB-GP",
  initial = initial,
  cohort_size = 2L,
  stopping_threshold = 1e12,
  stopping_patience = 1L
)
validate_output(
  stopped_output,
  scenario = scenario,
  planned_cohort = 2L,
  n_max = SIM_SETTINGS$n_max,
  b_max = scenario_budget(scenario)
)
stopifnot(nrow(stopped_output$allocations) == 1L)

cat("Algorithm-alignment validation passed.\n")

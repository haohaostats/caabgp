source(file.path("R", "00_settings.R"))
source(file.path("R", "08_table_rendering.R"))

table_dir <- file.path(ROOT_DIR, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

fmt <- function(mean, se) sprintf("%.3f (%.3f)", mean, se)

scenario_titles <- c(
  "Table 3. Fixed-budget operating characteristics in Scenario 1: shared main surface with moderate subgroup deviations.",
  "Table 4. Fixed-budget operating characteristics in Scenario 2: rare, noisy, and clinically prioritized stratum.",
  "Table 5. Fixed-budget operating characteristics in Scenario 3: mixed heterogeneity, rare-stratum prioritization, and high manufacturing cost.",
  "Table 6. Fixed-budget operating characteristics in Scenario 4: high outcome noise and value of replication."
)

for (scenario_id in 1:4) {
  source_csv <- file.path(DIR_RESULTS, sprintf("table_scenario%d_results_n1000.csv", scenario_id))
  dat <- read.csv(source_csv, stringsAsFactors = FALSE)
  display <- data.frame(
    Design = dat$method,
    `Final sample size` = fmt(dat$n_mean, dat$n_se),
    Regret = fmt(dat$regret_mean, dat$regret_se),
    `Near-optimal selection` = fmt(dat$near_optimal_mean, dat$near_optimal_se),
    `Dose distance` = fmt(dat$dose_distance_mean, dat$dose_distance_se),
    RPSEL = fmt(dat$rpsel_mean, dat$rpsel_se),
    `Unique doses` = fmt(dat$unique_doses_mean, dat$unique_doses_se),
    `Total cost` = fmt(dat$total_cost_mean, dat$total_cost_se),
    rho_b = ifelse(is.na(dat$rho_mean), "--", fmt(dat$rho_mean, dat$rho_se)),
    check.names = FALSE
  )
  stem <- sprintf("Table%02d_Scenario%d", scenario_id + 2L, scenario_id)
  render_table_pdf(
    display,
    file.path(table_dir, paste0(stem, ".pdf")),
    scenario_titles[scenario_id],
    "Values are Monte Carlo means, with standard errors in parentheses; 1,000 replicates per design.",
    column_widths = c(1.10, 1.55, 1.05, 1.65, 1.25, 1.05, 1.25, 1.35, 1.05),
    width = 13.0, height = 4.35, header_size = 7.8, body_size = 8.1
  )
  file.copy(source_csv, file.path(table_dir, paste0(stem, ".csv")), overwrite = TRUE)
}

sens_csv <- file.path(DIR_RESULTS, "table_sensitivity_results_n1000.csv")
sens <- read.csv(sens_csv, stringsAsFactors = FALSE, check.names = FALSE)
factor_label <- c(
  lambda_c = "lambda_c",
  c_new = "c_new",
  budget = "B_max",
  cohort_size = "Cohort size r"
)
if (!all(c("factor_label", "setting_label") %in% names(sens))) {
  stop("Sensitivity summary must contain plain-text factor_label and setting_label columns.")
}
setting_label <- sens$setting_label
display_sens <- data.frame(
  `Sensitivity factor` = unname(factor_label[sens$factor]),
  Setting = setting_label,
  `Final sample size` = fmt(sens$n_mean, sens$n_se),
  Regret = fmt(sens$regret_mean, sens$regret_se),
  `Near-optimal selection` = fmt(sens$near_optimal_mean, sens$near_optimal_se),
  `Unique doses` = fmt(sens$unique_doses_mean, sens$unique_doses_se),
  `Total cost` = fmt(sens$total_cost_mean, sens$total_cost_se),
  check.names = FALSE
)
render_table_pdf(
  display_sens,
  file.path(table_dir, "Table07_Sensitivity.pdf"),
  "Table 7. Sensitivity analyses for the proposed CA-AB-GP design under Scenario 3.",
  "Values are Monte Carlo means, with standard errors in parentheses; 1,000 replicates per setting.",
  column_widths = c(1.35, 1.65, 1.55, 1.15, 1.70, 1.35, 1.40),
  width = 12.2, height = 7.0, header_size = 8.0, body_size = 8.1
)
file.copy(sens_csv, file.path(table_dir, "Table07_Sensitivity.csv"), overwrite = TRUE)

message("Generated R-native vector PDF Tables 3--7 from verified n=1000 CSV files.")

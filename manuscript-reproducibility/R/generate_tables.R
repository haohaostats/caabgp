source(file.path("R", "00_settings.R"))

table_dir <- file.path(ROOT_DIR, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

fmt <- function(mean, se) sprintf("%.3f (%.3f)", mean, se)

standalone_table <- function(number, title, column_spec, header, rows, note) {
  c(
    "\\documentclass[border=6pt]{standalone}",
    "\\usepackage{booktabs}",
    "\\usepackage{graphicx}",
    "\\usepackage[scaled=0.92]{helvet}",
    "\\renewcommand{\\familydefault}{\\sfdefault}",
    "\\begin{document}",
    "\\begin{minipage}{19.5cm}",
    sprintf("\\textbf{Table %d. %s}", number, title),
    "",
    "\\vspace{0.35em}\\footnotesize",
    "\\resizebox{\\textwidth}{!}{%",
    sprintf("\\begin{tabular}{%s}", column_spec),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    paste0(rows, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}}",
    "",
    paste0("\\vspace{0.35em}\\scriptsize ", note),
    "\\end{minipage}",
    "\\end{document}"
  )
}

scenario_titles <- c(
  "Fixed-budget operating characteristics in Scenario 1: shared main surface with moderate subgroup deviations.",
  "Fixed-budget operating characteristics in Scenario 2: rare, noisy, and clinically prioritized stratum.",
  "Fixed-budget operating characteristics in Scenario 3: mixed heterogeneity, rare-stratum prioritization, and high manufacturing cost.",
  "Fixed-budget operating characteristics in Scenario 4: high outcome noise and value of replication."
)

for (scenario_id in 1:4) {
  source_csv <- file.path(DIR_RESULTS, sprintf("table_scenario%d_results_n1000.csv", scenario_id))
  dat <- read.csv(source_csv, stringsAsFactors = FALSE)
  rows <- vapply(seq_len(nrow(dat)), function(i) {
    rho <- if (is.na(dat$rho_mean[i])) "--" else fmt(dat$rho_mean[i], dat$rho_se[i])
    paste(dat$method[i], fmt(dat$n_mean[i], dat$n_se[i]),
          fmt(dat$regret_mean[i], dat$regret_se[i]),
          fmt(dat$near_optimal_mean[i], dat$near_optimal_se[i]),
          fmt(dat$dose_distance_mean[i], dat$dose_distance_se[i]),
          fmt(dat$rpsel_mean[i], dat$rpsel_se[i]),
          fmt(dat$unique_doses_mean[i], dat$unique_doses_se[i]),
          fmt(dat$total_cost_mean[i], dat$total_cost_se[i]), rho, sep = " & ")
  }, character(1))
  tex <- standalone_table(
    number = scenario_id + 2L,
    title = scenario_titles[scenario_id],
    column_spec = "lcccccccc",
    header = "Design & Final sample size & Regret & Near-optimal selection & Dose distance & RPSEL & Unique doses & Total cost & $\\widehat\\rho_b$",
    rows = rows,
    note = "Values are Monte Carlo means, with standard errors in parentheses; 1,000 replicates per design."
  )
  stem <- sprintf("Table%02d_Scenario%d", scenario_id + 2L, scenario_id)
  writeLines(tex, file.path(table_dir, paste0(stem, ".tex")), useBytes = TRUE)
  file.copy(source_csv, file.path(table_dir, paste0(stem, ".csv")), overwrite = TRUE)
}

sens_csv <- file.path(DIR_RESULTS, "table_sensitivity_results_n1000.csv")
sens <- read.csv(sens_csv, stringsAsFactors = FALSE, check.names = FALSE)
factor_label <- c(lambda_c = "$\\lambda_c$", c_new = "$c_{\\mathrm{new}}$",
                  budget = "$B_{\\max}$", cohort_size = "Cohort size $r$")
sens_rows <- vapply(seq_len(nrow(sens)), function(i) {
  paste(factor_label[[sens$factor[i]]], sens$setting_latex[i],
        fmt(sens$n_mean[i], sens$n_se[i]),
        fmt(sens$regret_mean[i], sens$regret_se[i]),
        fmt(sens$near_optimal_mean[i], sens$near_optimal_se[i]),
        fmt(sens$unique_doses_mean[i], sens$unique_doses_se[i]),
        fmt(sens$total_cost_mean[i], sens$total_cost_se[i]), sep = " & ")
}, character(1))
sens_tex <- standalone_table(
  number = 7L,
  title = "Sensitivity analyses for CA-AB-GP under Scenario 3.",
  column_spec = "llccccc",
  header = "Sensitivity factor & Setting & Final sample size & Regret & Near-optimal selection & Unique doses & Total cost",
  rows = sens_rows,
  note = "Values are Monte Carlo means, with standard errors in parentheses; 1,000 replicates per setting."
)
writeLines(sens_tex, file.path(table_dir, "Table07_Sensitivity.tex"), useBytes = TRUE)
file.copy(sens_csv, file.path(table_dir, "Table07_Sensitivity.csv"), overwrite = TRUE)

message("Generated standalone LaTeX Tables 3--7 from verified n=1000 CSV files.")

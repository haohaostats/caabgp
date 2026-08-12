root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "R", "07_application_realdata.R"))) stop("Run from the repository root.")

source(file.path("R", "Figure6_7_CaseApplication.R"))

case_results <- file.path(root, "case_application", "results")
case_figures <- file.path(root, "case_application", "figures")
case_tables <- file.path(root, "case_application", "tables")
dir.create(case_results, recursive = TRUE, showWarnings = FALSE)
dir.create(case_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(case_tables, recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(file.path(root, "results_r"), pattern = "^application_.*\\.csv$", full.names = TRUE), case_results, overwrite = TRUE)
file.copy(file.path(root, "figures", c("Figure6.pdf", "Figure7.pdf")), case_figures, overwrite = TRUE)
file.copy(file.path(root, "results_r", "table_application_results.tex"), file.path(case_tables, "Table09_ApplicationResults.tex"), overwrite = TRUE)
file.copy(file.path(root, "results_r", "application_results.csv"), file.path(case_tables, "Table09_ApplicationResults.csv"), overwrite = TRUE)
file.copy(file.path(root, "figures", c("Figure6.pdf", "Figure7.pdf")), file.path(root, "manuscript_outputs", "figures"), overwrite = TRUE)
file.copy(file.path(case_tables, c("Table09_ApplicationResults.tex", "Table09_ApplicationResults.csv")), file.path(root, "manuscript_outputs", "tables"), overwrite = TRUE)
cat("Rebuilt the fixed case application, Figures 6--7, and Table 9.\n")


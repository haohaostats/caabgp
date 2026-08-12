root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "R", "00_settings.R"))) stop("Run from the repository root.")

figure_scripts <- c(
  "Figure1_TrueResponseSurfaces.R",
  "Figure2_Scenario1Diagnostics.R",
  "Figure3_Scenario2Diagnostics.R",
  "Figure4_Scenario3Diagnostics.R",
  "Figure5_Scenario4Diagnostics.R"
)
for (script in figure_scripts) source(file.path("R", script), local = new.env(parent = globalenv()))
source(file.path("R", "generate_tables.R"), local = new.env(parent = globalenv()))

out_fig <- file.path(root, "manuscript_outputs", "figures")
out_tab <- file.path(root, "manuscript_outputs", "tables")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(file.path(root, "figures"), pattern = "^Figure[1-5]\\.pdf$", full.names = TRUE), out_fig, overwrite = TRUE)
file.copy(list.files(file.path(root, "tables"), pattern = "^Table(03|04|05|06|07)_.*\\.(csv|pdf)$", full.names = TRUE), out_tab, overwrite = TRUE)
cat("Rebuilt Figures 1--5 and Tables 3--7 in manuscript_outputs/.\n")

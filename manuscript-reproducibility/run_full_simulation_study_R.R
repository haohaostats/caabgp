## One-command R simulation pipeline for the main-text empirical-Bayes CA-AB-GP study.
##
## Usage from PowerShell:
##   & 'D:\Soft\R-4.5.2\bin\x64\Rscript.exe' run_full_simulation_study_R.R 1000 12 TRUE
##
## Arguments:
##   1. Monte Carlo replicates per scenario, default 1000.
##   2. Parallel workers, default detectCores() - 2.
##   3. Run sensitivity analyses, default TRUE.

source(file.path("R", "00_settings.R"))
source(file.path("R", "04_simulate_methods.R"))
source(file.path("R", "05_summarize_results.R"))
source(file.path("R", "06_plot_figures.R"))

args <- commandArgs(trailingOnly = TRUE)
n_mc <- if (length(args) >= 1) as.integer(args[1]) else SIM_SETTINGS$n_mc
workers <- if (length(args) >= 2) as.integer(args[2]) else max(1, parallel::detectCores() - 2)
run_sensitivity <- if (length(args) >= 3) as.logical(args[3]) else TRUE

cat("Running empirical-Bayes CA-AB-GP simulation pipeline\n")
cat("n_mc =", n_mc, "\n")
cat("workers =", workers, "\n")
cat("sensitivity analyses =", run_sensitivity, "\n")
cat("Chunk checkpoints =", DIR_CHECKPOINTS, "\n")
pipeline_started_at <- Sys.time()

cat("\n[Stage 1/4] Writing true response-surface data\n")
write_true_surface_csv()

cat("\n[Stage 2/4] Running primary simulations\n")
primary_results <- run_primary_simulation(n_mc = n_mc, workers = workers)

if (isTRUE(run_sensitivity)) {
  cat("\n[Stage 3/4] Running sensitivity analyses\n")
  sensitivity_results <- run_sensitivity_analyses(n_mc = n_mc, workers = workers)
} else {
  cat("\n[Stage 3/4] Sensitivity analyses skipped\n")
}

cat("\n[Stage 4/4] Creating numerical summaries and vector figures\n")
summary_results <- summarise_primary_outputs(n_mc = n_mc)
plot_results <- plot_all_figures(n_mc = n_mc)

cat("Done.\n")
cat("Total elapsed time: ", format_progress_time(as.numeric(difftime(Sys.time(), pipeline_started_at, units = "secs"))), "\n", sep = "")
cat("CSV results written to: ", DIR_RESULTS, "\n", sep = "")
cat("PDF figures written to: ", DIR_FIGURES, "\n", sep = "")

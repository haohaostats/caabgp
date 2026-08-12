# CA-AB-GP manuscript reproducibility code

This repository reproduces the simulation-study and real-data-calibrated application outputs reported in the CA-AB-GP manuscript.

## Quick reproduction of manuscript outputs

The repository includes the verified compact summaries from 1,000 Monte Carlo replicates. From the repository root, run:

```r
source("scripts/00_check_environment.R")
source("scripts/02_rebuild_main_text_outputs.R")
source("scripts/03_run_case_application.R")
source("scripts/04_validate_reproduction.R")
```

The regenerated vector PDF figures and LaTeX tables are collected in `manuscript_outputs/`.

## Repository structure

- `R/`: common algorithm, empirical-Bayes GP, acquisition, simulation, summary, and plotting functions.
- `simulation_study/scenario_1/`--`scenario_4/`: verified scenario-specific summaries, manuscript figure, and table.
- `simulation_study/sensitivity_analysis/`: verified Scenario 3 sensitivity results and Table 7.
- `case_application/`: Brock database, fixed application configuration, generated results, Figures 6--7, and Table 9.
- `results_r/`: compact verified inputs used to rebuild the manuscript displays without rerunning 1,000 replicates.
- `scripts/`: environment check, full simulation, display rebuild, case application, and validation entry points.
- `manuscript_outputs/`: publication-ready vector PDFs and LaTeX/CSV tables.

## Full simulation rerun

The full study is computationally intensive. The following command reruns all four scenarios and the sensitivity analyses with deterministic seeds:

```powershell
Rscript scripts/01_run_full_simulation.R 1000 12 TRUE
```

Arguments are the number of Monte Carlo replicates per setting, number of parallel workers, and whether to run sensitivity analyses. Checkpoints are written to `checkpoints/`, so an interrupted run can resume.

After a completed 1,000-replicate rerun, execute `Rscript scripts/02_rebuild_main_text_outputs.R` to rebuild the polished manuscript figures and tables from the newly generated raw results. Figure 3's compact budget trajectory is recomputed automatically whenever the corresponding raw trajectory file is present.

For a short functional test:

```powershell
Rscript scripts/01_run_full_simulation.R 2 1 FALSE
```

Use a clean copy for a short test because raw result filenames include the requested replicate count.

## Application definition

The manuscript application is a **real-data-calibrated heterogeneous replay**, not a direct subgroup analysis of observed patient-level outcomes. Its calibration series is fixed explicitly as `Gandhi2014_1`, with toxicity outcome 1 and efficacy outcome 25. The observed aggregate dose-level results define the common response pattern; two prespecified heterogeneous stratum deviations are then added. All construction constants are listed in `config/application_parameters.csv` and `R/07_application_realdata.R`.

The application entry point directly loads this prespecified study-outcome combination and does not rank candidate datasets. With seed `20260521`, it recreates the reported application table and vector Figures 6--7.

## Software

The analysis uses base R and the standard `parallel` package. It does not require an external R package. A LaTeX installation is optional and is needed only to compile the standalone table `.tex` files to PDF.

## Reproducibility notes

- Random-number seeds are deterministic and defined in the code.
- GP hyperparameters are estimated by empirical Bayes during each analysis; this is part of the prespecified algorithm.
- The four simulation scenarios and sensitivity grids are fixed in `R/00_settings.R` and `R/01_surfaces.R`.
- PDFs are vector graphics. Exact PDF bytes can differ across operating systems because of fonts and PDF metadata, so numerical CSV files are the primary validation targets.

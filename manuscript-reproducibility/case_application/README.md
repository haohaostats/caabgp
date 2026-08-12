# Case application

This directory reproduces the manuscript's real-data-calibrated heterogeneous replay.

The calibration series is fixed to `Gandhi2014_1`, toxicity outcome 1, and efficacy outcome 25. No study or outcome search is performed by the application code. Aggregate observed results establish the common response pattern; prespecified heterogeneous deviations create two replay strata. The complete constants are recorded in `../config/application_parameters.csv` and `../R/07_application_realdata.R`.

Run `source("scripts/03_run_case_application.R")` from the repository root. Outputs are copied into `results/`, `figures/`, and `tables/`.


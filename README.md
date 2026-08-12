# caabgp

[![R-CMD-check](https://github.com/haohaostats/caabgp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haohaostats/caabgp/actions/workflows/R-CMD-check.yaml)
[![GitHub release](https://img.shields.io/github/v/release/haohaostats/caabgp)](https://github.com/haohaostats/caabgp/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`caabgp` implements the empirical-Bayes cost-aware adaptive-borrowing
Gaussian-process Bayesian optimization design for personalized combination
dose-finding trials.

Version 0.2.0 aligns the sequential engine with the current design algorithm:
allocations always use the complete prespecified cohort, stratum-specific
acquisition stopping counters are updated after every cohort, and terminal
recommendations for all strata use the complete trial history.

The package is designed as a reusable trial-design engine.

| Adaptive borrowing | Cost-aware allocation | Operational safeguards |
|:--|:--|:--|
| Learns a shared response surface while retaining stratum-specific deviations. | Balances expected improvement against patient, screening, and novel-dose costs. | Enforces complete cohorts, fixed budgets, sample-size limits, and stratum-specific stopping. |

## Installation

### Local source installation

Download the source package from the latest
[release](https://github.com/haohaostats/caabgp/releases),
then install the local archive:

```r
install.packages("caabgp_0.2.0.tar.gz", repos = NULL, type = "source")
```

After cloning the source repository, install from the parent directory:

```r
install.packages("caabgp", repos = NULL, type = "source")
```

From a shell:

```powershell
R CMD INSTALL caabgp
```

### Remote installation from GitHub

```r
install.packages("remotes")
remotes::install_github("haohaostats/caabgp")
```

## Quick start

```r
library(caabgp)

grid <- dose_grid(seq(0, 1, by = 0.25), n_agents = 2)

design <- caabgp_design(
  dose_grid = grid,
  n_strata = 2,
  prevalence = c(0.60, 0.40),
  cohort_size = 2,
  n_max = 40,
  budget_max = 100,
  cost_novel = 5
)

truth <- function(stratum, dose, n = 1) {
  target <- if (stratum == 1) c(0.25, 0.75) else c(0.75, 0.25)
  rnorm(n, -exp(-sum((as.numeric(dose) - target)^2) / 0.15), 0.10)
}

initial <- make_initial_data(design, truth, n_initial = 10, seed = 21)
trial <- run_caabgp_trial(initial, design, truth, seed = 22)

recommend_dose(trial$fit)
borrowing_index(trial)
```

Smaller outcome values are treated as more favorable. In an actual trial,
use `fit_caabgp()`, `suggest_next()`, and `update_trial_data()` to alternate
between model fitting, allocation, and entry of newly observed outcomes.

### Example output

<p align="center">
  <img src="man/figures/caabgp-allocation.png" width="33%" alt="Sequential allocation map">
  <img src="man/figures/caabgp-surface-stratum-1.png" width="31%" alt="Predictive surface for stratum 1">
  <img src="man/figures/caabgp-surface-stratum-2.png" width="31%" alt="Predictive surface for stratum 2">
</p>

The allocation map shows the number of patients evaluated at each dose
combination. Triangles mark the final stratum-specific recommendations; darker
surface regions indicate lower, more favorable predicted outcomes.

| Stratum | Recommended dose | Predictive mean | Predictive SD |
|:--:|:--:|--:|--:|
| 1 | `(0.25, 0.75)` | -0.914 | 0.064 |
| 2 | `(0.75, 0.25)` | -0.988 | 0.067 |

| Final sample size | Unique doses | Repeat allocation | Total cost | Borrowing index |
|--:|--:|--:|--:|--:|
| 34 | 13 | 33.3% | 99 / 100 | 0.074 |

These deterministic quick-start results use the fixed seeds shown above and
are intended only to illustrate the package workflow.

## Advanced configuration

Additional options for clinical priority weights, patient and screening costs,
cost sensitivity, stratum-specific stopping, initial designs, and acquisition
control are documented in `?caabgp`.

## Manuscript reproducibility

The complete code and verified numerical summaries used to reproduce the
manuscript simulation study and real-data-calibrated application are provided
in [`manuscript-reproducibility/`](manuscript-reproducibility/). This directory
contains the four simulation scenarios, sensitivity analyses, the fixed Brock
calibration application, R-generated vector PDF figures and tables, deterministic seeds,
and validation scripts. It is excluded from the installed R package and does
not change the package API.

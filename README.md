# caabgp

[![R-CMD-check](https://github.com/haohaostats/caabgp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haohaostats/caabgp/actions/workflows/R-CMD-check.yaml)
[![GitHub release](https://img.shields.io/github/v/release/haohaostats/caabgp)](https://github.com/haohaostats/caabgp/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`caabgp` implements the empirical-Bayes cost-aware adaptive-borrowing
Gaussian-process Bayesian optimization design for personalized combination
dose-finding trials.

The package is designed as a reusable trial-design engine.

## Installation

### Local source installation

Download the source package from the latest
[release](https://github.com/haohaostats/caabgp/releases),
then install the local archive:

```r
install.packages("caabgp_0.1.1.tar.gz", repos = NULL, type = "source")
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

## Showcase example

```r
library(caabgp)

grid <- dose_grid(seq(0, 1, by = 0.25), n_agents = 2)
initial_design <- grid[c(1, 9, 13, 17, 25), ]

design <- caabgp_design(
  dose_grid = grid,
  n_strata = 2,
  prevalence = c(0.55, 0.45),
  weights = c(0.50, 0.50),
  cohort_size = 2,
  n_max = 80,
  budget_max = 220,
  cost_patient = 1,
  cost_novel = 5,
  cost_screening = 0,
  lambda_c = 0.5,
  initial_design = initial_design
)

truth <- function(stratum, dose, n = 1) {
  dose <- as.numeric(dose)
  shared <- -0.95 * exp(-sum((dose - c(0.50, 0.50))^2) / (2 * 0.42^2))
  deviation <- if (stratum == 1) {
    -0.45 * exp(-sum((dose - c(0.25, 0.75))^2) / (2 * 0.16^2))
  } else {
    -0.45 * exp(-sum((dose - c(0.75, 0.25))^2) / (2 * 0.16^2))
  }
  rnorm(n, shared + deviation, 0.05)
}

initial <- make_initial_data(design, truth, n_initial = 20, seed = 21)

trial <- run_caabgp_trial(
  initial_data = initial,
  design = design,
  outcome_generator = truth,
  seed = 22
)

recommend_dose(trial$fit)
borrowing_index(trial)

p <- plot_allocation(trial)
save_caabgp_figure(p, "allocation_map.pdf", width = 6.5, height = 5)
```

With the fixed seeds above, CA-AB-GP recommends different dose combinations for
the two strata while retaining substantial cross-stratum borrowing.

The plotting functions return `ggplot2` objects, so they can be further refined
with standard `ggplot2` layers before exporting to PDF or SVG.

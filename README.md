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

Download the source package from the
[v0.1.0 release](https://github.com/haohaostats/caabgp/releases/tag/v0.1.0),
then install the local archive:

```r
install.packages("caabgp_0.1.0.tar.gz", repos = NULL, type = "source")
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

## Minimal example

```r
library(caabgp)

grid <- dose_grid(seq(0, 1, by = 0.25), n_agents = 2)

design <- caabgp_design(
  dose_grid = grid,
  n_strata = 2,
  prevalence = c(0.65, 0.35),
  weights = c(0.55, 0.45),
  cohort_size = 2,
  n_max = 40,
  budget_max = 180,
  cost_patient = 1,
  cost_novel = 10,
  cost_screening = 1,
  lambda_c = 0.5,
  initial_design = grid[c(1, 5, 13, 21, 25), ]
)

truth <- function(stratum, dose, n = 1) {
  center <- if (stratum == 1) c(0.25, 0.75) else c(0.75, 0.25)
  mu <- -exp(-sum((dose - center)^2) / (2 * 0.18^2))
  rnorm(n, mu, 0.25)
}

initial <- make_initial_data(design, truth, n_initial = 12, seed = 1)

trial <- run_caabgp_trial(
  initial_data = initial,
  design = design,
  outcome_generator = truth,
  seed = 2
)

recommend_dose(trial$fit)
borrowing_index(trial)

p <- plot_allocation(trial)
save_caabgp_figure(p, "allocation_map.pdf", width = 6.5, height = 5)
```

The plotting functions return `ggplot2` objects, so they can be further refined
with standard `ggplot2` layers before exporting to PDF or SVG.

\name{caabgp}
\alias{caabgp}
\alias{dose_grid}
\alias{caabgp_design}
\alias{fit_caabgp}
\alias{predict_caabgp}
\alias{suggest_next}
\alias{run_caabgp_trial}
\alias{recommend_dose}
\alias{borrowing_index}
\alias{trial_cost}
\alias{update_trial_data}
\alias{make_initial_data}
\alias{plot_allocation}
\alias{plot_surface}
\alias{save_caabgp_figure}
\alias{theme_caabgp}
\alias{export_results}
\alias{run_simulation_study}
\alias{aei_min}
\title{Cost-aware adaptive-borrowing GP Bayesian optimization}
\description{
Core functions for the empirical-Bayes cost-aware adaptive-borrowing
Gaussian-process Bayesian optimization design for personalized combination
dose-finding trials.
}
\usage{
dose_grid(levels = seq(0, 1, by = 0.25), n_agents = 2, names = NULL)

caabgp_design(
  dose_grid,
  n_strata,
  prevalence = rep(1 / n_strata, n_strata),
  weights = prevalence,
  cohort_size = 1,
  n_max = 80,
  budget_max = Inf,
  cost_patient = 1,
  cost_novel = 0,
  cost_screening = 0,
  lambda_c = 0.5,
  xi = 1,
  pi_min = 0.02,
  initial_design = NULL,
  dose_cols = NULL
)

fit_caabgp(
  data,
  design,
  dose_cols = design$dose_cols,
  stratum_col = "stratum",
  outcome_col = "y",
  previous = NULL,
  maxit = 60
)

predict_caabgp(object, newdata = object$design$dose_grid, ...)

suggest_next(
  fit,
  design = fit$design,
  manufactured = NULL,
  cohort_size = design$cohort_size,
  remaining_budget = Inf,
  active_strata = seq_len(design$n_strata)
)

run_caabgp_trial(
  initial_data,
  design,
  outcome_generator,
  dose_cols = design$dose_cols,
  stratum_col = "stratum",
  outcome_col = "y",
  seed = NULL,
  maxit = 60
)

recommend_dose(fit, design = fit$design)

borrowing_index(x, ...)

trial_cost(data, design, dose_cols = design$dose_cols, stratum_col = "stratum")

update_trial_data(
  data,
  allocation,
  outcomes,
  design,
  dose_cols = design$dose_cols,
  stratum_col = "stratum",
  outcome_col = "y"
)

make_initial_data(
  design,
  outcome_generator,
  n_initial = NULL,
  seed = NULL,
  dose_cols = design$dose_cols,
  stratum_col = "stratum",
  outcome_col = "y"
)

plot_allocation(
  result,
  design = result$design,
  recommendations = NULL,
  main = NULL,
  base_size = 14
)

plot_surface(
  fit,
  stratum = 1,
  design = fit$design,
  recommendations = NULL,
  main = NULL,
  base_size = 14
)

save_caabgp_figure(
  plot,
  filename,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600,
  device = NULL
)

theme_caabgp(base_size = 14, base_family = "sans")

export_results(result, path)

run_simulation_study(
  scenarios,
  design,
  n_mc = 100,
  n_initial = NULL,
  seed = 1,
  maxit = 40
)

aei_min(mu, sd, sigma, xi = 1)
}
\arguments{
\item{levels}{Numeric dose levels used to construct the finite dose grid.}
\item{n_agents}{Number of treatment components.}
\item{names}{Optional dose-column names.}
\item{dose_grid}{Finite manufacturable dose grid as a data frame.}
\item{n_strata}{Number of patient strata.}
\item{prevalence}{Stratum prevalence vector used for screening-cost accounting.}
\item{weights}{Clinical priority weights used in the acquisition function.}
\item{cohort_size}{Default sequential cohort size.}
\item{n_max}{Maximum sample size or safety cap.}
\item{budget_max}{Maximum total operational budget.}
\item{cost_patient}{Per-patient evaluation cost.}
\item{cost_novel}{Cost for evaluating a previously unmanufactured dose combination.}
\item{cost_screening}{Screening or recruitment cost multiplier.}
\item{lambda_c}{Cost-sensitivity parameter for cost-aware AEI.}
\item{xi}{Effective-best conservatism parameter in augmented expected improvement.}
\item{pi_min}{Lower bound for prevalence in screening-cost calculations.}
\item{initial_design}{Optional initial dose combinations.}
\item{dose_cols}{Dose-column names.}
\item{data, initial_data}{Trial data with dose columns, stratum column, and outcome column.}
\item{design}{A \code{caabgp_design} object.}
\item{stratum_col}{Name of the stratum column.}
\item{outcome_col}{Name of the outcome column.}
\item{previous}{Optional previous fit for warm-starting empirical-Bayes optimization.}
\item{maxit}{Maximum number of optimizer iterations for each empirical-Bayes fit.}
\item{object, fit, x}{Fitted model or trial object.}
\item{newdata}{Dose grid for prediction.}
\item{...}{Additional arguments for S3 methods.}
\item{manufactured}{Optional vector of manufactured dose keys.}
\item{remaining_budget}{Remaining operational budget available for the next allocation.}
\item{active_strata}{Strata eligible for the next allocation.}
\item{outcome_generator}{Function \code{function(stratum, dose, n)} returning new outcomes.}
\item{seed}{Optional random seed.}
\item{allocation}{One-row allocation data frame returned by \code{suggest_next}.}
\item{outcomes}{New observed outcomes for the selected allocation.}
\item{n_initial}{Initial sample size.}
\item{result}{A \code{caabgp_trial} object or a compatible allocation data frame.}
\item{stratum}{Stratum index to display in a predictive-surface plot.}
\item{recommendations}{Optional recommendation data frame.}
\item{main}{Optional plot title.}
\item{base_size}{Base font size for publication-style plots.}
\item{base_family}{Base font family for publication-style plots.}
\item{plot}{A \code{ggplot2} plot object.}
\item{filename}{Output path for a saved figure. PDF and SVG outputs are vector graphics.}
\item{width, height}{Figure size.}
\item{units}{Figure size units passed to \code{ggplot2::ggsave}.}
\item{dpi}{Resolution used for raster outputs.}
\item{device}{Optional graphics device passed to \code{ggplot2::ggsave}.}
\item{path}{Directory for exported CSV files.}
\item{scenarios}{Named list of outcome-generator functions.}
\item{n_mc}{Number of Monte Carlo replicates.}
\item{mu, sd}{Predictive mean and standard deviation vectors.}
\item{sigma}{Residual standard deviation used in augmented expected improvement.}
}
\value{
Design, fit, trial, data-frame, \code{ggplot2} plot, or exported-file path
depending on the function.
}
\details{
This first package version is intentionally focused on the proposed CA-AB-GP
method. Comparator designs used in manuscript simulations are not exposed as
package functions. Figures are generated as \code{ggplot2} objects and can be
saved as vector PDF or SVG files using \code{save_caabgp_figure()}.
}

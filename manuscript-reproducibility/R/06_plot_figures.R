source(file.path("R", "00_settings.R"))
source(file.path("R", "01_surfaces.R"))

SCI_COLORS <- c(
  "S-GP" = "#4A5568",
  "P-GP" = "#2B6CB0",
  "IND-GP" = "#D97706",
  "AB-GP" = "#059669",
  "CA-AB-GP" = "#B83280"
)

surface_palette <- grDevices::colorRampPalette(c(
  "#17324D", "#1F6F8B", "#63A7A5", "#D8D7B1", "#F5E5C0", "#B65F42"
))

method_levels_in <- function(x) METHODS[METHODS %in% unique(as.character(x))]

se0 <- function(x) {
  x[!is.finite(x)] <- 0
  x
}

plot_scenario_surface_panel <- function(scenario, main = "") {
  grid <- candidate_grid()
  K <- length(scenario_prevalence(scenario))
  vals <- true_grid_values(scenario, grid)
  xvals <- sort(unique(grid[, 1]))
  yvals <- sort(unique(grid[, 2]))
  z_background <- matrix(apply(vals, 2, min), nrow = length(xvals), ncol = length(yvals))

  graphics::image(
    xvals, yvals, z_background,
    col = surface_palette(90),
    xlab = expression(d[1]),
    ylab = expression(d[2]),
    main = main,
    axes = FALSE,
    useRaster = FALSE
  )
  graphics::axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), cex.axis = 0.8)
  graphics::axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), cex.axis = 0.8, las = 1)
  contour_cols <- c("#FFFFFF", "#1A202C", "#F6AD55", "#68D391")
  for (k in seq_len(K)) {
    zmat <- matrix(vals[k, ], nrow = length(xvals), ncol = length(yvals))
    graphics::contour(xvals, yvals, zmat, add = TRUE, drawlabels = FALSE,
                      nlevels = 5, col = contour_cols[k], lwd = 1.0)
  }
  opt <- true_optima(scenario, grid)
  graphics::points(opt$d1_star, opt$d2_star, pch = 21, bg = "white", col = "#1A202C", cex = 1.25, lwd = 0.7)
  graphics::text(opt$d1_star, opt$d2_star, labels = paste0("Z", opt$stratum), pos = 3, cex = 0.68, col = "white")
  graphics::box(col = "#2D3748")
}

plot_figure1 <- function(n_mc = SIM_SETTINGS$n_mc) {
  write_true_surface_csv()
  grDevices::pdf(file.path(DIR_FIGURES, "Figure1_true_response_surfaces.pdf"), width = 8.6, height = 7.2, useDingbats = FALSE)
  old <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(3.0, 3.2, 2.4, 0.8), oma = c(0.2, 0.2, 1.2, 0.2), family = "sans")
  for (s in 1:4) plot_scenario_surface_panel(s, sprintf("%s  Scenario %d", LETTERS[s], s))
  graphics::mtext("True stratum-specific response surfaces", outer = TRUE, cex = 1.05, font = 2)
}

plot_line_with_se <- function(dat, metric, title, ylab, xlab = "Sample size") {
  dat <- dat[order(dat$n), ]
  methods <- method_levels_in(dat$method)
  mean_col <- paste0(metric, "_mean")
  se_col <- paste0(metric, "_se")
  dat[[se_col]] <- se0(dat[[se_col]])
  ylim <- range(c(dat[[mean_col]] - 1.96 * dat[[se_col]],
                  dat[[mean_col]] + 1.96 * dat[[se_col]]), na.rm = TRUE)
  if (!all(is.finite(ylim)) || diff(ylim) == 0) {
    center <- mean(dat[[mean_col]], na.rm = TRUE)
    ylim <- center + c(-0.05, 0.05)
  }

  graphics::plot(NA, xlim = range(dat$n, na.rm = TRUE), ylim = ylim,
                 xlab = xlab, ylab = ylab, main = title, axes = FALSE)
  graphics::axis(1, cex.axis = 0.8)
  graphics::axis(2, cex.axis = 0.8, las = 1)
  graphics::grid(col = "#E2E8F0", lty = 1)
  for (m in methods) {
    dm <- dat[dat$method == m, ]
    mu <- dm[[mean_col]]
    se <- se0(dm[[se_col]])
    col <- SCI_COLORS[m]
    graphics::polygon(c(dm$n, rev(dm$n)), c(mu - 1.96 * se, rev(mu + 1.96 * se)),
                      border = NA, col = grDevices::adjustcolor(col, 0.12))
    graphics::lines(dm$n, mu, col = col, lwd = ifelse(m == "CA-AB-GP", 2.4, 1.7))
  }
  graphics::box(col = "#2D3748")
}

plot_simple_bar <- function(dat, mean_col, se_col, title, ylab, ylim = NULL, las = 2) {
  dat <- dat[order(match(dat$method, METHODS)), ]
  se <- se0(dat[[se_col]])
  if (is.null(ylim)) {
    ymax <- max(dat[[mean_col]] + 1.96 * se, na.rm = TRUE)
    ylim <- c(0, max(0.05, ymax * 1.20))
  }
  bp <- graphics::barplot(
    dat[[mean_col]],
    names.arg = dat$method,
    col = SCI_COLORS[dat$method],
    border = NA,
    ylim = ylim,
    ylab = ylab,
    main = title,
    las = las,
    cex.names = 0.70
  )
  graphics::arrows(bp, pmax(0, dat[[mean_col]] - 1.96 * se),
                   bp, dat[[mean_col]] + 1.96 * se,
                   angle = 90, code = 3, length = 0.03, col = "#2D3748")
  graphics::box(col = "#2D3748")
}

plot_grouped_bars <- function(dat, group_col, mean_col, se_col, title, ylab, ylim = NULL, legend_pos = "topright") {
  dat$method <- as.character(method_order(dat$method))
  dat <- dat[order(dat[[group_col]], match(dat$method, METHODS)), ]
  mat <- tapply(dat[[mean_col]], list(dat$method, dat[[group_col]]), mean, na.rm = TRUE)
  mat <- mat[METHODS[METHODS %in% rownames(mat)], , drop = FALSE]
  if (is.null(ylim)) ylim <- c(0, max(mat, na.rm = TRUE) * 1.25)
  graphics::barplot(
    mat,
    beside = TRUE,
    col = SCI_COLORS[rownames(mat)],
    border = NA,
    ylim = ylim,
    ylab = ylab,
    main = title,
    las = 1,
    cex.names = 0.80
  )
  graphics::legend(legend_pos, legend = rownames(mat), fill = SCI_COLORS[rownames(mat)], bty = "n", cex = 0.68)
  graphics::box(col = "#2D3748")
}

plot_method_scatter <- function(dat, x_col, y_col, title, xlab, ylab) {
  dat <- dat[order(match(dat$method, METHODS)), ]
  x <- dat[[x_col]]
  y <- dat[[y_col]]
  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)
  if (diff(xr) == 0) xr <- xr + c(-1, 1)
  if (diff(yr) == 0) yr <- yr + c(-0.05, 0.05)
  graphics::plot(x, y, pch = 21, bg = SCI_COLORS[dat$method], col = "white", cex = 1.6,
                 xlim = xr, ylim = yr, xlab = xlab, ylab = ylab, main = title, axes = FALSE)
  graphics::axis(1, cex.axis = 0.8)
  graphics::axis(2, las = 1, cex.axis = 0.8)
  graphics::grid(col = "#E2E8F0", lty = 1)
  graphics::text(x, y, labels = dat$method, pos = 4, cex = 0.68, col = "#1A202C")
  graphics::box(col = "#2D3748")
}

plot_final_distribution <- function(dat, value_col, title, ylab) {
  dat$method <- factor(dat$method, levels = METHODS)
  graphics::boxplot(
    stats::as.formula(paste(value_col, "~ method")),
    data = dat,
    col = grDevices::adjustcolor(SCI_COLORS[levels(dat$method)], 0.42),
    border = "#2D3748",
    outline = FALSE,
    ylab = ylab,
    main = title,
    las = 2,
    cex.axis = 0.7
  )
  graphics::stripchart(
    stats::as.formula(paste(value_col, "~ method")),
    data = dat,
    vertical = TRUE,
    method = "jitter",
    pch = 16,
    cex = 0.35,
    col = grDevices::adjustcolor("#1A202C", 0.28),
    add = TRUE
  )
  graphics::box(col = "#2D3748")
}

plot_stacked_novel_repeat <- function(dat, title) {
  dat <- dat[order(match(dat$method, METHODS)), ]
  mat <- rbind(
    Novel = dat$novel_fraction_mean,
    Repeated = dat$repeat_fraction_mean
  )
  colnames(mat) <- dat$method
  graphics::barplot(
    mat,
    beside = FALSE,
    col = c("#B83280", "#CBD5E0"),
    border = NA,
    ylim = c(0, 1),
    ylab = "Allocation fraction",
    main = title,
    las = 2,
    cex.names = 0.70
  )
  graphics::legend("topright", legend = rownames(mat), fill = c("#B83280", "#CBD5E0"), bty = "n", cex = 0.75)
  graphics::box(col = "#2D3748")
}

plot_scenario1 <- function(n_mc = SIM_SETTINGS$n_mc) {
  learning <- read.csv(file.path(DIR_RESULTS, sprintf("figure2_scenario1_learning_trajectory_n%d.csv", n_mc)))
  final_dist <- read.csv(file.path(DIR_RESULTS, sprintf("figure2_scenario1_final_distribution_n%d.csv", n_mc)))
  borrowing <- read.csv(file.path(DIR_RESULTS, sprintf("figure2_scenario1_borrowing_final_n%d.csv", n_mc)))
  unique <- read.csv(file.path(DIR_RESULTS, sprintf("figure2_scenario1_cumulative_unique_doses_n%d.csv", n_mc)))

  grDevices::pdf(file.path(DIR_FIGURES, "Figure2_scenario1_results.pdf"), width = 9.2, height = 6.6, useDingbats = FALSE)
  old <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(3.6, 3.9, 2.5, 1.1), family = "sans")
  plot_line_with_se(learning, "regret", "A  Weighted regret trajectory", "Weighted regret")
  graphics::legend("topright", legend = METHODS, col = SCI_COLORS[METHODS], lwd = c(1.7, 1.7, 1.7, 1.7, 2.4), bty = "n", cex = 0.68)
  plot_final_distribution(final_dist, "regret", "B  Final weighted regret", "Weighted regret")
  plot_simple_bar(borrowing, "rho_mean", "rho_se", expression("C  Borrowing index " * hat(rho)[b]), expression(hat(rho)[b]), ylim = c(0, 1), las = 1)
  plot_line_with_se(unique, "unique_doses", "D  Cumulative unique doses", "Unique dose combinations")
}

plot_scenario2 <- function(n_mc = SIM_SETTINGS$n_mc) {
  learning <- read.csv(file.path(DIR_RESULTS, sprintf("figure3_scenario2_learning_trajectory_n%d.csv", n_mc)))
  near <- read.csv(file.path(DIR_RESULTS, sprintf("figure3_scenario2_near_optimal_by_stratum_n%d.csv", n_mc)))
  alloc <- read.csv(file.path(DIR_RESULTS, sprintf("figure3_scenario2_allocation_by_stratum_n%d.csv", n_mc)))
  operational <- read.csv(file.path(DIR_RESULTS, sprintf("figure3_scenario2_final_operational_n%d.csv", n_mc)))

  grDevices::pdf(file.path(DIR_FIGURES, "Figure3_scenario2_results.pdf"), width = 9.2, height = 6.6, useDingbats = FALSE)
  old <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(3.6, 3.9, 2.5, 1.1), family = "sans")
  plot_line_with_se(learning, "regret", "A  Weighted regret trajectory", "Weighted regret")
  plot_grouped_bars(near, "stratum", "near_optimal_mean", "near_optimal_se", "B  Stratum-specific accuracy", "Near-optimal probability", ylim = c(0, 1), legend_pos = "bottomright")
  plot_grouped_bars(alloc, "stratum", "n_patients_mean", "n_patients_se", "C  Patient allocation by stratum", "Mean patients", legend_pos = "topright")
  plot_method_scatter(operational, "n_mean", "unique_doses_mean", "D  Budget use", "Final sample size", "Unique dose combinations")
}

plot_scenario3 <- function(n_mc = SIM_SETTINGS$n_mc) {
  operational <- read.csv(file.path(DIR_RESULTS, sprintf("figure4_scenario3_final_operational_n%d.csv", n_mc)))
  unique <- read.csv(file.path(DIR_RESULTS, sprintf("figure4_scenario3_cumulative_unique_doses_n%d.csv", n_mc)))
  near <- read.csv(file.path(DIR_RESULTS, sprintf("figure4_scenario3_near_optimal_by_stratum_n%d.csv", n_mc)))

  grDevices::pdf(file.path(DIR_FIGURES, "Figure4_scenario3_results.pdf"), width = 9.2, height = 6.6, useDingbats = FALSE)
  old <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(3.6, 3.9, 2.5, 1.1), family = "sans")
  plot_method_scatter(operational, "total_cost_mean", "regret_mean", "A  Cost-regret relationship", "Total cost", "Weighted regret")
  plot_line_with_se(unique, "unique_doses", "B  Cumulative unique doses", "Unique dose combinations")
  plot_grouped_bars(near, "stratum", "near_optimal_mean", "near_optimal_se", "C  Stratum-specific accuracy", "Near-optimal probability", ylim = c(0, 1), legend_pos = "bottomright")
  plot_method_scatter(operational, "n_mean", "total_cost_mean", "D  Fixed-budget operating point", "Final sample size", "Total operational cost")
}

plot_scenario4 <- function(n_mc = SIM_SETTINGS$n_mc) {
  learning <- read.csv(file.path(DIR_RESULTS, sprintf("figure5_scenario4_learning_trajectory_n%d.csv", n_mc)))
  novel <- read.csv(file.path(DIR_RESULTS, sprintf("figure5_scenario4_novel_repeat_fraction_n%d.csv", n_mc)))
  operational <- read.csv(file.path(DIR_RESULTS, sprintf("figure5_scenario4_final_operational_n%d.csv", n_mc)))

  grDevices::pdf(file.path(DIR_FIGURES, "Figure5_scenario4_results.pdf"), width = 9.2, height = 6.6, useDingbats = FALSE)
  old <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(3.6, 3.9, 2.5, 1.1), family = "sans")
  plot_line_with_se(learning, "regret", "A  Weighted regret trajectory", "Weighted regret")
  plot_stacked_novel_repeat(novel, "B  Novel versus repeated evaluations")
  plot_simple_bar(operational, "rpsel_mean", "rpsel_se", "C  Root predictive squared error loss", "RPSEL")
  plot_method_scatter(operational, "n_mean", "unique_doses_mean", "D  Replication under fixed budget", "Final sample size", "Unique dose combinations")
}

plot_sensitivity <- function(n_mc = SIM_SETTINGS$n_mc) {
  lambda_path <- file.path(DIR_RESULTS, sprintf("figure6_lambda_sensitivity_n%d.csv", n_mc))
  cnew_path <- file.path(DIR_RESULTS, sprintf("figure6_cnew_sensitivity_n%d.csv", n_mc))
  budget_path <- file.path(DIR_RESULTS, sprintf("figure6_budget_sensitivity_n%d.csv", n_mc))
  cohort_path <- file.path(DIR_RESULTS, sprintf("figure6_cohort_sensitivity_n%d.csv", n_mc))
  if (!all(file.exists(c(lambda_path, cnew_path, budget_path, cohort_path)))) return(invisible(FALSE))

  lambda <- read.csv(lambda_path)
  cnew <- read.csv(cnew_path)
  budget <- read.csv(budget_path)
  cohort <- read.csv(cohort_path)

  grDevices::pdf(file.path(DIR_FIGURES, "Figure6_sensitivity_analyses.pdf"), width = 9.2, height = 6.6, useDingbats = FALSE)
  old <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(3.6, 3.9, 2.5, 1.1), family = "sans")

  graphics::plot(lambda$unique_doses_mean, lambda$regret_mean, type = "b", pch = 21,
                 bg = SCI_COLORS["CA-AB-GP"], col = SCI_COLORS["CA-AB-GP"], lwd = 2,
                 xlab = "Unique dose combinations", ylab = "Weighted regret",
                 main = expression("A  Cost sensitivity " * lambda[c]), axes = FALSE)
  graphics::text(lambda$unique_doses_mean, lambda$regret_mean, labels = lambda$lambda_c, pos = 4, cex = 0.75)
  graphics::axis(1, cex.axis = 0.8); graphics::axis(2, las = 1, cex.axis = 0.8)
  graphics::grid(col = "#E2E8F0", lty = 1); graphics::box(col = "#2D3748")

  graphics::plot(cnew$c_new, cnew$unique_doses_mean, type = "b", pch = 21,
                 bg = SCI_COLORS["CA-AB-GP"], col = SCI_COLORS["CA-AB-GP"], lwd = 2,
                 xlab = expression(c[new]), ylab = "Unique dose combinations",
                 main = expression("B  Manufacturing cost " * c[new]), axes = FALSE)
  graphics::axis(1, at = cnew$c_new, cex.axis = 0.8); graphics::axis(2, las = 1, cex.axis = 0.8)
  graphics::grid(col = "#E2E8F0", lty = 1); graphics::box(col = "#2D3748")

  graphics::plot(budget$b_max, budget$regret_mean, type = "b", pch = 21,
                 bg = SCI_COLORS["CA-AB-GP"], col = SCI_COLORS["CA-AB-GP"], lwd = 2,
                 xlab = expression(B[max]), ylab = "Weighted regret",
                 main = expression("C  Fixed trial budget " * B[max]), axes = FALSE)
  graphics::axis(1, at = budget$b_max, cex.axis = 0.8); graphics::axis(2, las = 1, cex.axis = 0.8)
  graphics::grid(col = "#E2E8F0", lty = 1); graphics::box(col = "#2D3748")

  graphics::plot(cohort$cohort_size, cohort$regret_mean, type = "b", pch = 21,
                 bg = SCI_COLORS["CA-AB-GP"], col = SCI_COLORS["CA-AB-GP"], lwd = 2,
                 xlab = "Cohort size", ylab = "Weighted regret",
                 main = "D  Sequential cohort size", axes = FALSE)
  graphics::axis(1, at = cohort$cohort_size, cex.axis = 0.8); graphics::axis(2, las = 1, cex.axis = 0.8)
  graphics::grid(col = "#E2E8F0", lty = 1); graphics::box(col = "#2D3748")

  invisible(TRUE)
}

plot_all_figures <- function(n_mc = SIM_SETTINGS$n_mc) {
  plot_figure1(n_mc)
  plot_scenario1(n_mc)
  plot_scenario2(n_mc)
  plot_scenario3(n_mc)
  plot_scenario4(n_mc)
  plot_sensitivity(n_mc)
}

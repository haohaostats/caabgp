## Redraw Figure 3 in a saturated color SCI journal style.
##
## Usage from the project root:
##   & 'D:\Soft\R-4.5.2\bin\x64\Rscript.exe' redraw_figure3_sci_color.R
##
## Outputs:
##   images/Figures/Figure3_scenario2_results.pdf
##   images/Figures/Figure3_scenario2_results_SCI_color.pdf

source(file.path("R", "00_settings.R"))

DIR_FIGURES <- file.path(ROOT_DIR, "figures")
DIR_RESULTS <- file.path(ROOT_DIR, "results_r")
dir.create(DIR_FIGURES, recursive = TRUE, showWarnings = FALSE)

n_mc <- 1000

budget_learning <- read.csv(file.path(
  DIR_RESULTS,
  sprintf("figure3_scenario2_budget_regret_trajectory_n%d.csv", n_mc)
))
near_opt <- read.csv(file.path(DIR_RESULTS, sprintf("figure3_scenario2_near_optimal_by_stratum_n%d.csv", n_mc)))
allocation <- read.csv(file.path(DIR_RESULTS, sprintf("figure3_scenario2_allocation_by_stratum_n%d.csv", n_mc)))
operational <- read.csv(file.path(DIR_RESULTS, sprintf("figure3_scenario2_final_operational_n%d.csv", n_mc)))

METHOD_ORDER <- c("S-GP", "P-GP", "IND-GP", "AB-GP", "CA-AB-GP")

method_cols <- c(
  "S-GP" = "#4D4D4D",
  "P-GP" = "#0057B8",
  "IND-GP" = "#F28E2B",
  "AB-GP" = "#2CA02C",
  "CA-AB-GP" = "#D62728"
)

method_fills <- c(
  "S-GP" = "#CFCFCF",
  "P-GP" = "#8AB6F0",
  "IND-GP" = "#FFC078",
  "AB-GP" = "#8ED081",
  "CA-AB-GP" = "#F28B82"
)

line_types <- c(
  "S-GP" = 3,
  "P-GP" = 2,
  "IND-GP" = 4,
  "AB-GP" = 5,
  "CA-AB-GP" = 1
)

line_widths <- c(
  "S-GP" = 2.00,
  "P-GP" = 2.10,
  "IND-GP" = 2.05,
  "AB-GP" = 2.15,
  "CA-AB-GP" = 3.35
)

point_shapes <- c(
  "S-GP" = 21,
  "P-GP" = 24,
  "IND-GP" = 22,
  "AB-GP" = 23,
  "CA-AB-GP" = 21
)

axis_col <- "#222222"
grid_col <- "#ECECEC"
text_col <- "#111111"

se0 <- function(x) {
  x[!is.finite(x)] <- 0
  x
}

methods_present <- function(dat) {
  METHOD_ORDER[METHOD_ORDER %in% unique(as.character(dat$method))]
}

thin_time_grid <- function(dat, max_points = 24) {
  out <- do.call(rbind, lapply(split(dat, dat$method), function(dm) {
    dm <- dm[order(dm$n), , drop = FALSE]
    if (nrow(dm) <= max_points) return(dm)
    keep <- unique(round(seq(1, nrow(dm), length.out = max_points)))
    dm[keep, , drop = FALSE]
  }))
  rownames(out) <- NULL
  out
}

make_budget_regret_trajectory <- function(dat, scenario_id = 2, grid = seq(0.38, 1.00, length.out = 26)) {
  dat <- dat[dat$scenario == scenario_id & dat$method %in% METHOD_ORDER, , drop = FALSE]
  dat <- dat[is.finite(dat$budget_fraction) & is.finite(dat$regret), , drop = FALSE]

  pieces <- lapply(split(dat, interaction(dat$method, dat$replicate, drop = TRUE)), function(dm) {
    dm <- dm[order(dm$budget_fraction, dm$n), , drop = FALSE]
    method <- as.character(dm$method[1])
    replicate <- dm$replicate[1]
    vals <- vapply(grid, function(g) {
      idx <- which(dm$budget_fraction <= g)
      if (length(idx) == 0) {
        dm$regret[1]
      } else {
        dm$regret[max(idx)]
      }
    }, numeric(1))
    data.frame(method = method, replicate = replicate, budget_fraction = grid, regret = vals)
  })

  expanded <- do.call(rbind, pieces)
  rownames(expanded) <- NULL

  stats <- aggregate(
    regret ~ method + budget_fraction,
    data = expanded,
    FUN = function(x) c(mean = mean(x), se = stats::sd(x) / sqrt(length(x)), n = length(x))
  )
  out <- data.frame(
    method = stats$method,
    budget_fraction = stats$budget_fraction,
    regret_mean = stats$regret[, "mean"],
    regret_se = stats$regret[, "se"],
    n_replicates = stats$regret[, "n"]
  )
  out[order(match(out$method, METHOD_ORDER), out$budget_fraction), , drop = FALSE]
}

draw_grid <- function() {
  graphics::axis(1, cex.axis = 1.02, col = axis_col, col.axis = text_col, tck = -0.018)
  graphics::axis(2, cex.axis = 1.02, col = axis_col, col.axis = text_col, las = 1, tck = -0.018)
  graphics::box(col = axis_col, lwd = 1.10)
}

panel_title <- function(label) {
  graphics::mtext(label, side = 3, line = 0.62, adj = 0, cex = 1.08, font = 2, col = text_col)
}

plot_regret_trajectory <- function(dat) {
  dat <- dat[is.finite(dat$regret_mean), , drop = FALSE]
  ylim <- c(0, max(dat$regret_mean + 1.96 * se0(dat$regret_se), na.rm = TRUE) * 1.05)

  graphics::plot(
    NA,
    xlim = c(0.38, 1.00),
    ylim = ylim,
    xlab = "Cumulative budget fraction",
    ylab = "Weighted regret",
    axes = FALSE,
    xaxs = "i",
    yaxs = "i"
  )
  draw_grid()

  for (m in methods_present(dat)) {
    dm <- dat[dat$method == m, , drop = FALSE]
    dm <- dm[order(dm$budget_fraction), , drop = FALSE]
    graphics::lines(dm$budget_fraction, dm$regret_mean, col = method_cols[m],
                    lty = line_types[m], lwd = line_widths[m])
  }

  panel_title("A  Regret over budget")
}

grouped_bar <- function(dat, value_col, se_col, ylab, title, ylim = NULL, rare_label = FALSE) {
  dat$method <- factor(dat$method, levels = METHOD_ORDER)
  dat <- dat[order(dat$stratum, dat$method), , drop = FALSE]
  strata <- sort(unique(dat$stratum))
  mat <- sapply(strata, function(s) {
    vals <- dat[dat$stratum == s, c("method", value_col)]
    setNames(vals[[value_col]], vals$method)[METHOD_ORDER]
  })
  semat <- sapply(strata, function(s) {
    vals <- dat[dat$stratum == s, c("method", se_col)]
    setNames(vals[[se_col]], vals$method)[METHOD_ORDER]
  })
  if (is.null(dim(mat))) {
    mat <- matrix(mat, ncol = length(strata), dimnames = list(METHOD_ORDER, strata))
    semat <- matrix(semat, ncol = length(strata), dimnames = list(METHOD_ORDER, strata))
  }
  mat[!is.finite(mat)] <- NA
  semat[!is.finite(semat)] <- 0

  if (is.null(ylim)) {
    ymax <- max(mat + 1.96 * semat, na.rm = TRUE)
    ylim <- c(0, ymax * 1.18)
  }

  bp <- graphics::barplot(
    mat,
    beside = TRUE,
    col = method_fills[METHOD_ORDER],
    border = method_cols[METHOD_ORDER],
    lwd = 0.80,
    ylim = ylim,
    axes = FALSE,
    xlab = "Stratum",
    ylab = ylab,
    names.arg = paste0("Z", strata),
    cex.names = 1.02
  )
  graphics::axis(2, cex.axis = 1.02, col = axis_col, col.axis = text_col, las = 1, tck = -0.018)
  graphics::box(col = axis_col, lwd = 1.10)

  lower <- pmax(ylim[1], mat - 1.96 * semat)
  upper <- pmin(ylim[2], mat + 1.96 * semat)
  graphics::arrows(bp, lower, bp, upper, angle = 90, code = 3,
                   length = 0.030, lwd = 0.90, col = "#333333")

  if (isTRUE(rare_label)) {
    rare_center <- mean(bp[, ncol(bp)], na.rm = TRUE)
    graphics::text(rare_center, ylim[2] * 0.965, "rare stratum",
                   cex = 0.88, font = 2, col = "#D62728")
  }

  panel_title(title)
}

plot_operational_tradeoff <- function(dat) {
  dat$method <- factor(dat$method, levels = METHOD_ORDER)
  dat <- dat[order(dat$method), , drop = FALSE]

  x <- dat$unique_doses_mean
  y <- dat$n_mean
  xse <- se0(dat$unique_doses_se)
  yse <- se0(dat$n_se)
  xlim <- range(c(x - 1.96 * xse, x + 1.96 * xse), na.rm = TRUE)
  ylim <- range(c(y - 1.96 * yse, y + 1.96 * yse), na.rm = TRUE)
  xpad <- diff(xlim) * 0.22
  ypad <- diff(ylim) * 0.24
  xlim <- xlim + c(-xpad, xpad)
  ylim <- ylim + c(-ypad, ypad)

  graphics::plot(
    NA,
    xlim = xlim,
    ylim = ylim,
    xlab = "Unique dose combinations",
    ylab = "Final sample size",
    axes = FALSE,
    xaxs = "i",
    yaxs = "i"
  )
  draw_grid()

  for (i in seq_len(nrow(dat))) {
    m <- as.character(dat$method[i])
    graphics::segments(x[i] - 1.96 * xse[i], y[i], x[i] + 1.96 * xse[i], y[i],
                       col = method_cols[m], lwd = 1.35)
    graphics::segments(x[i], y[i] - 1.96 * yse[i], x[i], y[i] + 1.96 * yse[i],
                       col = method_cols[m], lwd = 1.35)
    graphics::points(x[i], y[i], pch = 21, bg = method_fills[m],
                     col = method_cols[m], cex = if (m == "CA-AB-GP") 2.30 else 1.95,
                     lwd = if (m == "CA-AB-GP") 1.85 else 1.45)
  }

  label_offsets <- data.frame(
    method = METHOD_ORDER,
    dx = c(-0.18, -0.18, 0.15, 0.12, 0.12),
    dy = c(-1.0, 1.1, -0.9, 1.1, 1.1)
  )
  for (i in seq_len(nrow(dat))) {
    m <- as.character(dat$method[i])
    off <- label_offsets[label_offsets$method == m, , drop = FALSE]
    graphics::text(x[i] + off$dx, y[i] + off$dy, labels = m,
                   cex = 0.88, font = 2, col = method_cols[m])
  }

  panel_title("D  Budget-use trade-off")
}

plot_legend_strip <- function() {
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = METHOD_ORDER,
    col = method_cols[METHOD_ORDER],
    lty = line_types[METHOD_ORDER],
    lwd = line_widths[METHOD_ORDER],
    pch = point_shapes[METHOD_ORDER],
    pt.bg = method_fills[METHOD_ORDER],
    pt.cex = 1.30,
    bty = "n",
    horiz = TRUE,
    x.intersp = 0.80,
    seg.len = 3.10,
    cex = 1.08
  )
}

raw_trajectory_path <- file.path(DIR_RESULTS, sprintf("simulation_trajectory_n%d.csv", n_mc))
if (file.exists(raw_trajectory_path)) {
  raw_trajectory <- read.csv(raw_trajectory_path)
  budget_learning <- make_budget_regret_trajectory(raw_trajectory, scenario_id = 2)
  write.csv(
    budget_learning,
    file.path(DIR_RESULTS, sprintf("figure3_scenario2_budget_regret_trajectory_n%d.csv", n_mc)),
    row.names = FALSE
  )
}

write_figure <- function(path) {
  grDevices::pdf(path, width = 10.2, height = 8.1, family = "Helvetica",
                 useDingbats = FALSE, pointsize = 13)
  old <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old)
    grDevices::dev.off()
  }, add = TRUE)

  layout_matrix <- rbind(c(1, 0, 2), c(3, 0, 4), c(5, 5, 5))
  graphics::layout(layout_matrix, heights = c(1, 1, 0.17), widths = c(1, 0.11, 1))
  graphics::par(mar = c(4.2, 4.5, 2.3, 1.0), mgp = c(2.65, 0.78, 0), tcl = -0.28,
                las = 1, family = "Helvetica", cex.lab = 1.12)

  plot_regret_trajectory(budget_learning)
  grouped_bar(near_opt, "near_optimal_mean", "near_optimal_se",
              "Near-optimal probability", "B  Stratum-specific selection", ylim = c(0, 1.10), rare_label = TRUE)
  grouped_bar(allocation, "n_patients_mean", "n_patients_se",
              "Average patients", "C  Patient allocation by stratum", rare_label = TRUE)
  plot_operational_tradeoff(operational)
  plot_legend_strip()

  invisible(path)
}

out_main <- file.path(DIR_FIGURES, "Figure3.pdf")
out_backup <- out_main

for (path in unique(c(out_backup, out_main))) {
  tryCatch(
    write_figure(path),
    error = function(e) {
      warning(sprintf("Could not write %s: %s", path, conditionMessage(e)))
    }
  )
}

message("Saved color SCI Figure 3 to:")
message("  ", normalizePath(out_backup, winslash = "/", mustWork = FALSE))
message("  ", normalizePath(out_main, winslash = "/", mustWork = FALSE))

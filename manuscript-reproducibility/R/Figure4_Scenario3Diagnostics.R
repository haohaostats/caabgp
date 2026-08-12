## Redraw Figure 4 in a saturated color SCI journal style.
##
## Usage from the project root:
##   & 'D:\Soft\R-4.5.2\bin\x64\Rscript.exe' redraw_figure4_sci_color.R
##
## Outputs:
##   images/Figures/Figure4_scenario3_results.pdf
##   images/Figures/Figure4_scenario3_results_SCI_color.pdf

source(file.path("R", "00_settings.R"))

DIR_FIGURES <- file.path(ROOT_DIR, "figures")
DIR_RESULTS <- file.path(ROOT_DIR, "results_r")
dir.create(DIR_FIGURES, recursive = TRUE, showWarnings = FALSE)

n_mc <- 1000
scenario_id <- 3

cumulative_unique <- read.csv(file.path(DIR_RESULTS, sprintf("figure4_scenario3_cumulative_unique_doses_n%d.csv", n_mc)))
near_opt <- read.csv(file.path(DIR_RESULTS, sprintf("figure4_scenario3_near_optimal_by_stratum_n%d.csv", n_mc)))
operational <- read.csv(file.path(DIR_RESULTS, sprintf("figure4_scenario3_final_operational_n%d.csv", n_mc)))

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

draw_grid <- function() {
  graphics::axis(1, cex.axis = 1.02, col = axis_col, col.axis = text_col, tck = -0.018)
  graphics::axis(2, cex.axis = 1.02, col = axis_col, col.axis = text_col, las = 1, tck = -0.018)
  graphics::box(col = axis_col, lwd = 1.10)
}

panel_title <- function(label) {
  usr <- graphics::par("usr")
  panel <- substr(label, 1, 1)
  title <- sub("^[A-D]  ", "", label)
  y <- usr[4] + 0.055 * diff(usr[3:4])
  graphics::text(usr[1], y, labels = panel, adj = c(0, 0), xpd = NA,
                 cex = 1.08, font = 2, col = text_col)
  graphics::text(usr[1] + 0.065 * diff(usr[1:2]), y, labels = title,
                 adj = c(0, 0), xpd = NA, cex = 1.08, font = 2, col = text_col)
}

draw_callout_labels <- function(x, y, labels, label_x, label_y, cols) {
  label_pos <- c("S-GP" = 3, "P-GP" = 2, "IND-GP" = 4,
                 "AB-GP" = 1, "CA-AB-GP" = 4)
  label_offset <- c("S-GP" = 0.75, "P-GP" = 0.80, "IND-GP" = 0.75,
                    "AB-GP" = 0.75, "CA-AB-GP" = 0.80)
  for (i in seq_along(labels)) {
    graphics::text(x[i], y[i], labels = labels[i], pos = label_pos[[labels[i]]],
                   offset = label_offset[[labels[i]]], col = cols[i],
                   cex = 0.84, font = 2, xpd = NA)
  }
}

make_budget_unique_trajectory <- function(dat, scenario_id = 3, grid = seq(0.48, 1.00, length.out = 27)) {
  dat <- dat[dat$scenario == scenario_id & dat$method %in% METHOD_ORDER, , drop = FALSE]
  dat <- dat[is.finite(dat$budget_fraction) & is.finite(dat$unique_doses), , drop = FALSE]

  pieces <- lapply(split(dat, interaction(dat$method, dat$replicate, drop = TRUE)), function(dm) {
    dm <- dm[order(dm$budget_fraction, dm$n), , drop = FALSE]
    method <- as.character(dm$method[1])
    replicate <- dm$replicate[1]
    vals <- vapply(grid, function(g) {
      idx <- which(dm$budget_fraction <= g)
      if (length(idx) == 0) {
        dm$unique_doses[1]
      } else {
        dm$unique_doses[max(idx)]
      }
    }, numeric(1))
    data.frame(method = method, replicate = replicate, budget_fraction = grid, unique_doses = vals)
  })

  expanded <- do.call(rbind, pieces)
  rownames(expanded) <- NULL

  stats <- aggregate(
    unique_doses ~ method + budget_fraction,
    data = expanded,
    FUN = function(x) c(mean = mean(x), se = stats::sd(x) / sqrt(length(x)), n = length(x))
  )
  out <- data.frame(
    method = stats$method,
    budget_fraction = stats$budget_fraction,
    unique_doses_mean = stats$unique_doses[, "mean"],
    unique_doses_se = stats$unique_doses[, "se"],
    n_replicates = stats$unique_doses[, "n"]
  )
  out[order(match(out$method, METHOD_ORDER), out$budget_fraction), , drop = FALSE]
}

plot_cost_regret <- function(dat) {
  dat$method <- factor(dat$method, levels = METHOD_ORDER)
  dat <- dat[order(dat$method), , drop = FALSE]
  x <- dat$total_cost_mean
  y <- dat$regret_mean
  xse <- se0(dat$total_cost_se)
  yse <- se0(dat$regret_se)
  xlim <- range(c(x - 1.96 * xse, x + 1.96 * xse, 550), na.rm = TRUE)
  ylim <- range(c(0, y - 1.96 * yse, y + 1.96 * yse), na.rm = TRUE)
  xlim <- xlim + c(-0.08, 0.08) * diff(xlim)
  ylim <- ylim + c(-0.04, 0.12) * diff(ylim)

  graphics::plot(NA, xlim = xlim, ylim = ylim,
                 xlab = "Total operational cost",
                 ylab = "Weighted regret",
                 axes = FALSE, xaxs = "i", yaxs = "i")
  draw_grid()
  graphics::abline(v = 550, col = "#777777", lty = 3, lwd = 1.10)

  for (i in seq_len(nrow(dat))) {
    m <- as.character(dat$method[i])
    graphics::segments(x[i] - 1.96 * xse[i], y[i], x[i] + 1.96 * xse[i], y[i],
                       col = method_cols[m], lwd = 1.35)
    graphics::segments(x[i], y[i] - 1.96 * yse[i], x[i], y[i] + 1.96 * yse[i],
                       col = method_cols[m], lwd = 1.35)
    graphics::points(x[i], y[i], pch = point_shapes[m], bg = method_fills[m], col = method_cols[m],
                     cex = if (m == "CA-AB-GP") 2.30 else 1.95,
                     lwd = if (m == "CA-AB-GP") 1.85 else 1.45)
  }

  label_pos <- data.frame(
    method = METHOD_ORDER,
    x = c(544.70, 548.36, 548.55, 549.33, 549.25),
    y = c(0.122, 0.076, 0.236, 0.108, 0.074)
  )
  label_pos <- label_pos[match(as.character(dat$method), label_pos$method), , drop = FALSE]
  draw_callout_labels(
    x, y, as.character(dat$method),
    label_pos$x, label_pos$y,
    method_cols[as.character(dat$method)]
  )

  panel_title("A  Cost-regret relationship")
}

plot_unique_sample_size <- function(dat) {
  dat <- dat[is.finite(dat$unique_doses_mean), , drop = FALSE]
  dat <- dat[dat$n_replicates >= 50, , drop = FALSE]
  dat <- dat[dat$n <= 50, , drop = FALSE]
  ylim <- c(4.8, max(dat$unique_doses_mean + 1.96 * se0(dat$unique_doses_se), na.rm = TRUE) * 1.06)
  graphics::plot(NA, xlim = c(min(dat$n, na.rm = TRUE), 50), ylim = ylim,
                 xlab = "Sample size",
                 ylab = "Unique dose combinations",
                 axes = FALSE, xaxs = "i", yaxs = "i")
  draw_grid()
  for (m in methods_present(dat)) {
    dm <- dat[dat$method == m, , drop = FALSE]
    dm <- dm[order(dm$n), , drop = FALSE]
    graphics::lines(dm$n, dm$unique_doses_mean, col = method_cols[m],
                    lty = line_types[m], lwd = line_widths[m])
  }
  panel_title("B  Cumulative unique doses")
}

plot_near_optimal_by_stratum <- function(dat) {
  dat$method <- factor(dat$method, levels = METHOD_ORDER)
  dat <- dat[order(dat$stratum, dat$method), , drop = FALSE]
  strata <- sort(unique(dat$stratum))
  mat <- sapply(strata, function(s) {
    vals <- dat[dat$stratum == s, c("method", "near_optimal_mean")]
    setNames(vals$near_optimal_mean, vals$method)[METHOD_ORDER]
  })
  semat <- sapply(strata, function(s) {
    vals <- dat[dat$stratum == s, c("method", "near_optimal_se")]
    setNames(vals$near_optimal_se, vals$method)[METHOD_ORDER]
  })
  mat[!is.finite(mat)] <- NA
  semat[!is.finite(semat)] <- 0

  bp <- graphics::barplot(mat, beside = TRUE,
                          col = method_fills[METHOD_ORDER],
                          border = method_cols[METHOD_ORDER],
                          lwd = 0.80,
                          ylim = c(0, 1.10),
                          axes = FALSE,
                          xlab = "Stratum",
                          ylab = "Near-optimal probability",
                          names.arg = paste0("Z", strata),
                          cex.names = 1.02)
  graphics::axis(2, cex.axis = 1.02, col = axis_col, col.axis = text_col, las = 1, tck = -0.018)
  graphics::box(col = axis_col, lwd = 1.10)
  graphics::arrows(bp, pmax(0, mat - 1.96 * semat), bp, pmin(1.10, mat + 1.96 * semat),
                   angle = 90, code = 3, length = 0.030, lwd = 0.90, col = "#333333")
  rare_center <- mean(bp[, ncol(bp)], na.rm = TRUE)
  graphics::text(rare_center, 1.06, "rare stratum",
                 cex = 0.88, font = 2, col = "#D62728")
  panel_title("C  Stratum-specific selection")
}

plot_sample_cost <- function(dat) {
  dat$method <- factor(dat$method, levels = METHOD_ORDER)
  dat <- dat[order(dat$method), , drop = FALSE]
  x <- dat$n_mean
  y <- dat$total_cost_mean
  xse <- se0(dat$n_se)
  yse <- se0(dat$total_cost_se)
  xlim <- range(c(x - 1.96 * xse, x + 1.96 * xse), na.rm = TRUE)
  ylim <- range(c(y - 1.96 * yse, y + 1.96 * yse, 550), na.rm = TRUE)
  xlim <- xlim + c(-0.18, 0.18) * diff(xlim)
  ylim <- ylim + c(-0.12, 0.14) * diff(ylim)

  graphics::plot(NA, xlim = xlim, ylim = ylim,
                 xlab = "Final sample size",
                 ylab = "Total operational cost",
                 axes = FALSE, xaxs = "i", yaxs = "i")
  draw_grid()
  graphics::abline(h = 550, col = "#777777", lty = 3, lwd = 1.10)

  for (i in seq_len(nrow(dat))) {
    m <- as.character(dat$method[i])
    graphics::segments(x[i] - 1.96 * xse[i], y[i], x[i] + 1.96 * xse[i], y[i],
                       col = method_cols[m], lwd = 1.35)
    graphics::segments(x[i], y[i] - 1.96 * yse[i], x[i], y[i] + 1.96 * yse[i],
                       col = method_cols[m], lwd = 1.35)
    graphics::points(x[i], y[i], pch = point_shapes[m], bg = method_fills[m], col = method_cols[m],
                     cex = if (m == "CA-AB-GP") 2.30 else 1.95,
                     lwd = if (m == "CA-AB-GP") 1.85 else 1.45)
  }

  label_pos <- data.frame(
    method = METHOD_ORDER,
    x = c(42.00, 49.40, 41.10, 53.35, 62.70),
    y = c(543.25, 548.35, 549.35, 549.30, 549.35)
  )
  label_pos <- label_pos[match(as.character(dat$method), label_pos$method), , drop = FALSE]
  draw_callout_labels(
    x, y, as.character(dat$method),
    label_pos$x, label_pos$y,
    method_cols[as.character(dat$method)]
  )

  panel_title("D  Sample size and budget use")
}

plot_legend_strip <- function() {
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend("center",
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
                   cex = 1.08)
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

  plot_cost_regret(operational)
  plot_unique_sample_size(cumulative_unique)
  plot_near_optimal_by_stratum(near_opt)
  plot_sample_cost(operational)
  plot_legend_strip()

  invisible(path)
}

out_main <- file.path(DIR_FIGURES, "Figure4.pdf")
out_backup <- out_main

for (path in unique(c(out_backup, out_main))) {
  tryCatch(
    write_figure(path),
    error = function(e) {
      warning(sprintf("Could not write %s: %s", path, conditionMessage(e)))
    }
  )
}

message("Saved color SCI Figure 4 to:")
message("  ", normalizePath(out_backup, winslash = "/", mustWork = FALSE))
message("  ", normalizePath(out_main, winslash = "/", mustWork = FALSE))

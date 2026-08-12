## Redraw Figure 5 in a saturated color SCI journal style.
##
## Usage from the project root:
##   & 'D:\Soft\R-4.5.2\bin\x64\Rscript.exe' redraw_figure5_sci_color.R
##
## Outputs:
##   images/Figures/Figure5_scenario4_results.pdf
##   images/Figures/Figure5_scenario4_results_SCI_color.pdf

source(file.path("R", "00_settings.R"))

DIR_FIGURES <- file.path(ROOT_DIR, "figures")
DIR_RESULTS <- file.path(ROOT_DIR, "results_r")
dir.create(DIR_FIGURES, recursive = TRUE, showWarnings = FALSE)

n_mc <- 1000

learning <- read.csv(file.path(DIR_RESULTS, sprintf("figure5_scenario4_learning_trajectory_n%d.csv", n_mc)))
novel_repeat <- read.csv(file.path(DIR_RESULTS, sprintf("figure5_scenario4_novel_repeat_fraction_n%d.csv", n_mc)))
final_dist <- read.csv(file.path(DIR_RESULTS, sprintf("figure5_scenario4_final_distribution_n%d.csv", n_mc)))
operational <- read.csv(file.path(DIR_RESULTS, sprintf("figure5_scenario4_final_operational_n%d.csv", n_mc)))

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
  graphics::mtext(label, side = 3, line = 0.62, adj = 0, cex = 1.06, font = 2, col = text_col)
}

draw_callout_labels <- function(x, y, labels, label_x, label_y, cols) {
  label_pos <- c("S-GP" = 2, "P-GP" = 3, "IND-GP" = 3,
                 "AB-GP" = 3, "CA-AB-GP" = 1)
  label_offset <- c("S-GP" = 0.75, "P-GP" = 0.75, "IND-GP" = 0.75,
                    "AB-GP" = 0.75, "CA-AB-GP" = 0.75)
  for (i in seq_along(labels)) {
    graphics::text(x[i], y[i], labels = labels[i], pos = label_pos[[labels[i]]],
                   offset = label_offset[[labels[i]]], col = cols[i],
                   cex = 0.84, font = 2, xpd = NA)
  }
}

plot_regret_trajectory <- function(dat) {
  dat <- dat[is.finite(dat$regret_mean) & dat$n_replicates >= 800, , drop = FALSE]
  ylim <- c(0, max(dat$regret_mean + 1.96 * se0(dat$regret_se), na.rm = TRUE) * 1.04)
  graphics::plot(NA,
                 xlim = range(dat$n, na.rm = TRUE),
                 ylim = ylim,
                 xlab = "Sample size",
                 ylab = "Weighted regret",
                 axes = FALSE, xaxs = "i", yaxs = "i")
  draw_grid()
  for (m in methods_present(dat)) {
    dm <- dat[dat$method == m, , drop = FALSE]
    dm <- dm[order(dm$n), , drop = FALSE]
    graphics::lines(dm$n, dm$regret_mean, col = method_cols[m],
                    lty = line_types[m], lwd = line_widths[m])
  }
  panel_title("A  Regret trajectory")
}

plot_novel_repeat <- function(dat) {
  dat$method <- factor(dat$method, levels = METHOD_ORDER)
  dat <- dat[order(dat$method), , drop = FALSE]
  mat <- rbind(
    "Repeat" = dat$repeat_fraction_mean,
    "Novel" = dat$novel_fraction_mean
  )
  colnames(mat) <- as.character(dat$method)
  bp <- graphics::barplot(mat,
                          beside = FALSE,
                          col = c("#D8DEE9", "#3B4252"),
                          border = "#222222",
                          lwd = 0.75,
                          ylim = c(0, 1.18),
                          axes = FALSE,
                          xlab = "Design",
                          ylab = "Allocation fraction",
                          names.arg = METHOD_ORDER,
                          cex.names = 0.93)
  graphics::axis(2, at = seq(0, 1, by = 0.2), cex.axis = 1.02,
                 col = axis_col, col.axis = text_col, las = 1, tck = -0.018)
  graphics::box(col = axis_col, lwd = 1.10)
  graphics::text(bp, dat$repeat_fraction_mean / 2,
                 labels = sprintf("%.0f%%", 100 * dat$repeat_fraction_mean),
                 cex = 0.84, font = 2, col = "#111111")
  graphics::text(bp, dat$repeat_fraction_mean + dat$novel_fraction_mean / 2,
                 labels = sprintf("%.0f%%", 100 * dat$novel_fraction_mean),
                 cex = 0.84, font = 2, col = "#FFFFFF")
  graphics::legend("topright",
                   legend = c("Repeat", "Novel"),
                   fill = c("#D8DEE9", "#3B4252"),
                   border = "#222222",
                   bty = "n",
                   cex = 0.88,
                   inset = c(0.00, 0.01))
  panel_title("B  Exploration versus replication")
}

plot_rpsel_distribution <- function(dat) {
  dat$method <- factor(dat$method, levels = METHOD_ORDER)
  dat <- dat[order(dat$method), , drop = FALSE]
  graphics::boxplot(rpsel ~ method,
                    data = dat,
                    col = method_fills[METHOD_ORDER],
                    border = "#222222",
                    medcol = "#000000",
                    medlwd = 1.60,
                    whisklty = 1,
                    staplewex = 0.45,
                    boxwex = 0.58,
                    xlab = "",
                    ylab = "RPSEL",
                    axes = FALSE,
                    outline = FALSE)
  graphics::axis(1, at = seq_along(METHOD_ORDER), labels = METHOD_ORDER,
                 cex.axis = 0.90, las = 2, col = axis_col, col.axis = text_col, tck = -0.018)
  graphics::axis(2, cex.axis = 1.02, col = axis_col, col.axis = text_col, las = 1, tck = -0.018)
  set.seed(20260521)
  graphics::stripchart(rpsel ~ method,
                       data = dat,
                       vertical = TRUE,
                       method = "jitter",
                       jitter = 0.11,
                       add = TRUE,
                       pch = 16,
                       cex = 0.24,
                       col = grDevices::adjustcolor("#111111", alpha.f = 0.16))
  graphics::box(col = axis_col, lwd = 1.10)
  panel_title("C  RPSEL distribution")
}

plot_sample_unique <- function(dat) {
  dat$method <- factor(dat$method, levels = METHOD_ORDER)
  dat <- dat[order(dat$method), , drop = FALSE]
  x <- dat$n_mean
  y <- dat$unique_doses_mean
  xse <- se0(dat$n_se)
  yse <- se0(dat$unique_doses_se)
  xlim <- range(c(x - 1.96 * xse, x + 1.96 * xse), na.rm = TRUE)
  ylim <- range(c(y - 1.96 * yse, y + 1.96 * yse), na.rm = TRUE)
  xlim <- xlim + c(-0.32, 0.32) * diff(xlim)
  ylim <- ylim + c(-0.20, 0.20) * diff(ylim)

  graphics::plot(NA,
                 xlim = xlim,
                 ylim = ylim,
                 xlab = "Final sample size",
                 ylab = "Unique dose combinations",
                 axes = FALSE, xaxs = "i", yaxs = "i")
  draw_grid()
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
    x = c(116.15, 116.55, 118.05, 115.25, 119.15),
    y = c(15.18, 16.98, 18.95, 18.95, 12.82)
  )
  label_pos <- label_pos[match(as.character(dat$method), label_pos$method), , drop = FALSE]
  draw_callout_labels(
    x, y, as.character(dat$method),
    label_pos$x, label_pos$y,
    method_cols[as.character(dat$method)]
  )
  panel_title("D  Sample size and unique doses")
}

plot_method_legend <- function() {
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

  plot_regret_trajectory(learning)
  plot_novel_repeat(novel_repeat)
  plot_rpsel_distribution(final_dist)
  plot_sample_unique(operational)
  plot_method_legend()

  invisible(path)
}

out_main <- file.path(DIR_FIGURES, "Figure5.pdf")
out_backup <- out_main

for (path in unique(c(out_backup, out_main))) {
  tryCatch(
    write_figure(path),
    error = function(e) {
      warning(sprintf("Could not write %s: %s", path, conditionMessage(e)))
    }
  )
}

message("Saved color SCI Figure 5 to:")
message("  ", normalizePath(out_backup, winslash = "/", mustWork = FALSE))
message("  ", normalizePath(out_main, winslash = "/", mustWork = FALSE))

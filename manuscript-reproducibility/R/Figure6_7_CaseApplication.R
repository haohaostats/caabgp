# Re-run the application from the algorithm-aligned source tree, then redraw
# Figures 6--7 with embedded fonts and direct labels (no leader lines).
source("run_realdata_application_R.R")

# Collect the publication-ready direct-label versions with the other main-text figures.
DIR_FIGURES <- file.path(ROOT_DIR, "figures")
dir.create(DIR_FIGURES, recursive = TRUE, showWarnings = FALSE)

open_vector_pdf <- function(path, width, height) {
  grDevices::cairo_pdf(
    filename = path,
    width = width,
    height = height,
    family = "Arial",
    pointsize = 13,
    onefile = TRUE
  )
}

direct_label_position <- function(x, y) {
  c(
    if (x <= 0.15) x + 0.12 else if (x >= 0.85) x - 0.12 else x + 0.08,
    if (y <= 0.15) y + 0.09 else if (y >= 0.85) y - 0.09 else y + 0.08
  )
}

plot_application_allocation_direct <- function(outputs) {
  out <- file.path(DIR_FIGURES, "Figure6.pdf")
  open_vector_pdf(out, width = 9.4, height = 3.9)
  old <- par(no.readonly = TRUE)
  on.exit({ par(old); dev.off() }, add = TRUE)
  par(mfrow = c(1, 3), mar = c(4.2, 4.1, 2.4, 1.0),
      mgp = c(2.45, 0.75, 0), tcl = -0.28, cex.lab = 1.1)
  grid <- candidate_grid()
  for (m in METHOD_ORDER_APP) {
    dat <- outputs[[m]]$data
    rec <- outputs[[m]]$recommendations
    counts <- aggregate(stratum ~ d1 + d2, dat, length)
    names(counts)[3] <- "n"
    plot(NA, xlim = c(-0.06, 1.06), ylim = c(-0.06, 1.06),
         xlab = "Component 1 level", ylab = "Component 2 level",
         axes = FALSE, xaxs = "i", yaxs = "i")
    draw_axis()
    points(grid[, 1], grid[, 2], pch = 16, col = "#CFCFCF", cex = 0.82)
    points(counts$d1, counts$d2, pch = 21, bg = method_fills_app[m],
           col = method_cols_app[m], cex = 0.75 + 0.20 * sqrt(counts$n), lwd = 1.35)
    for (i in seq_len(nrow(rec))) {
      off <- recommendation_offset(rec$stratum[i])
      x_rec <- rec$d1_hat[i] + off[1]
      y_rec <- rec$d2_hat[i] + off[2]
      points(x_rec, y_rec, pch = 24,
             bg = stratum_cols_app[as.character(rec$stratum[i])],
             col = "#111111", cex = 1.65, lwd = 1.1)
      lab <- direct_label_position(x_rec, y_rec)
      draw_callout_label(lab[1], lab[2], paste0("Z", rec$stratum[i]),
                         fill = "#111111", text_col = "#FFFFFF", cex = 0.82)
    }
    mtext(m, side = 3, line = 0.65, adj = 0, cex = 1.18, font = 2)
  }
  out
}

plot_application_surfaces_direct <- function(pred, recs) {
  out <- file.path(DIR_FIGURES, "Figure7.pdf")
  pred <- pred[pred$method == "CA-AB-GP", ]
  zlim <- range(pred$mu, finite = TRUE)
  open_vector_pdf(out, width = 8.8, height = 4.2)
  old <- par(no.readonly = TRUE)
  on.exit({ par(old); dev.off() }, add = TRUE)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.12))
  par(mar = c(4.2, 4.1, 2.4, 0.9),
      mgp = c(2.45, 0.75, 0), tcl = -0.28, cex.lab = 1.1)
  vals <- SIM_SETTINGS$dose_values
  for (k in 1:2) {
    pk <- pred[pred$stratum == k, ]
    z <- matrix(pk$mu, nrow = length(vals), ncol = length(vals))
    image(vals, vals, z, col = surface_cols_app, zlim = zlim, axes = FALSE,
          xlab = "Component 1 level", ylab = "Component 2 level", useRaster = FALSE)
    contour(vals, vals, z, add = TRUE, drawlabels = FALSE, nlevels = 7,
            col = grDevices::adjustcolor("#FFFFFF", 0.75), lwd = 0.75)
    draw_axis()
    rk <- recs[recs$method == "CA-AB-GP" & recs$stratum == k, ]
    points(rk$d1_hat, rk$d2_hat, pch = 24,
           bg = stratum_cols_app[as.character(k)], col = "#FFFFFF",
           cex = 1.95, lwd = 2.1)
    points(rk$d1_hat, rk$d2_hat, pch = 24, bg = NA,
           col = "#111111", cex = 2.08, lwd = 0.9)
    lab <- direct_label_position(rk$d1_hat, rk$d2_hat)
    draw_callout_label(lab[1], lab[2], paste0("Z", k),
                       fill = "#111111", text_col = "#FFFFFF", cex = 0.92)
    mtext(paste0("Stratum ", k), side = 3, line = 0.65,
          adj = 0, cex = 1.18, font = 2)
  }
  par(mar = c(4.2, 0.4, 2.4, 3.0), mgp = c(1.6, 0.45, 0))
  zmat <- matrix(seq(zlim[1], zlim[2], length.out = 256), nrow = 1)
  image(1, seq(zlim[1], zlim[2], length.out = 256), zmat,
        col = surface_cols_app, axes = FALSE, xlab = "", ylab = "", useRaster = FALSE)
  axis(4, las = 1, cex.axis = 0.9, tck = -0.18)
  mtext("Predicted outcome", side = 4, line = 1.85, cex = 0.86, font = 2)
  box(lwd = 1.0)
  out
}

fig6_direct <- plot_application_allocation_direct(app$outputs)
fig7_direct <- plot_application_surfaces_direct(app$predictions, app$recommendations)
cat("Redrawn direct-label vector figures:\n")
cat(normalizePath(fig6_direct, winslash = "/", mustWork = FALSE), "\n")
cat(normalizePath(fig7_direct, winslash = "/", mustWork = FALSE), "\n")

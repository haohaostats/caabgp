## Redraw Figure 1 in a polished journal style.
##
## Usage from the project root:
##   & 'D:\Soft\R-4.5.2\bin\x64\Rscript.exe' redraw_figure1_sci.R
##
## Outputs:
##   images/Figures/Figure1_true_response_surfaces.pdf
##   images/Figures/Figure1_true_response_surfaces_SCI.pdf
##   results_r/figure1_true_surfaces_dense_grid.csv

source(file.path("R", "00_settings.R"))
source(file.path("R", "01_surfaces.R"))

DIR_FIGURES <- file.path(ROOT_DIR, "figures")
DIR_RESULTS <- file.path(ROOT_DIR, "results_r")
dir.create(DIR_FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_RESULTS, recursive = TRUE, showWarnings = FALSE)

dense_candidate_grid <- function(n = 151) {
  x <- seq(0, 1, length.out = n)
  as.matrix(expand.grid(d1 = x, d2 = x))
}

scenario_labels <- c(
  "A  Scenario 1",
  "B  Scenario 2",
  "C  Scenario 3",
  "D  Scenario 4"
)

stratum_cols <- c("#FFFFFF", "#FFD166", "#06D6A0", "#EF476F")
stratum_ltys <- c(1, 1, 1, 1)

surface_cols <- grDevices::colorRampPalette(c(
  "#06152B", "#0B2C5D", "#155A8A", "#14889A",
  "#51B89F", "#B8D99A", "#F3E7B1", "#FFF7D6"
))(256)

panel_frame_col <- "#22313F"
text_col <- "#1F2933"
grid_col <- grDevices::adjustcolor("#FFFFFF", alpha.f = 0.26)

scale_to_cols <- function(z, zlim) {
  idx <- round((z - zlim[1]) / diff(zlim) * (length(surface_cols) - 1)) + 1
  idx <- pmax(1, pmin(length(surface_cols), idx))
  idx
}

draw_panel <- function(scenario, label, zlim) {
  dense <- dense_candidate_grid()
  K <- length(scenario_prevalence(scenario))
  vals <- true_grid_values(scenario, dense)
  xvals <- sort(unique(dense[, 1]))
  yvals <- sort(unique(dense[, 2]))

  ## Show the best achievable subgroup-specific response at each dose location.
  z_background <- matrix(apply(vals, 2, min), nrow = length(xvals), ncol = length(yvals))

  graphics::image(
    xvals, yvals, z_background,
    zlim = zlim,
    col = surface_cols,
    axes = FALSE,
    xlab = expression(d[1]),
    ylab = expression(d[2]),
    useRaster = FALSE
  )

  ## A very light manufacturable grid anchors the continuous surface to the trial design.
  grid <- candidate_grid()
  graphics::points(grid[, 1], grid[, 2], pch = 16, cex = 0.16, col = grid_col)

  for (k in seq_len(K)) {
    zk <- matrix(vals[k, ], nrow = length(xvals), ncol = length(yvals))
    graphics::contour(
      xvals, yvals, zk,
      add = TRUE,
      drawlabels = FALSE,
      nlevels = 6,
      col = grDevices::adjustcolor(stratum_cols[k], alpha.f = 0.86),
      lwd = 0.62,
      lty = stratum_ltys[k]
    )
  }

  opt <- true_optima(scenario, candidate_grid())
  for (k in seq_len(nrow(opt))) {
    graphics::points(opt$d1_star[k], opt$d2_star[k], pch = 21, bg = stratum_cols[k],
                     col = "#111827", lwd = 0.85, cex = 1.45)
    graphics::text(opt$d1_star[k], opt$d2_star[k], labels = paste0("Z", opt$stratum[k]),
                   pos = 3, offset = 0.34, cex = 0.86, col = "#111827", font = 2)
  }

  graphics::axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                 tck = -0.018, cex.axis = 1.02, col = panel_frame_col, col.axis = text_col)
  graphics::axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                 tck = -0.018, cex.axis = 1.02, las = 1, col = panel_frame_col, col.axis = text_col)
  graphics::box(col = panel_frame_col, lwd = 1.10)
  graphics::mtext(label, side = 3, line = 0.72, adj = 0, cex = 1.18, font = 2, col = text_col)
}

draw_colorbar <- function(zlim) {
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  graphics::par(mar = c(3.35, 0.25, 2.10, 3.00), mgp = c(1.60, 0.45, 0))
  z <- matrix(seq(zlim[1], zlim[2], length.out = 256), nrow = 1)
  graphics::image(
    x = 1,
    y = seq(zlim[1], zlim[2], length.out = 256),
    z = z,
    col = surface_cols,
    axes = FALSE,
    xlab = "",
    ylab = "",
    useRaster = FALSE
  )
  ticks <- pretty(zlim, n = 5)
  graphics::axis(4, at = ticks, labels = sprintf("%.2f", ticks), las = 1,
                 cex.axis = 0.95, tck = -0.18, col = panel_frame_col, col.axis = text_col)
  graphics::mtext("Outcome", side = 4, line = 1.75,
                  cex = 0.88, font = 2, col = text_col)
  graphics::mtext("lower is better", side = 1, line = 1.55,
                  cex = 0.78, col = text_col)
  graphics::box(col = panel_frame_col, lwd = 1.00)
}

write_dense_surface_csv <- function() {
  dense <- dense_candidate_grid()
  rows <- list()
  idx <- 1
  for (scenario in 1:4) {
    K <- length(scenario_prevalence(scenario))
    for (k in seq_len(K)) {
      rows[[idx]] <- data.frame(
        scenario = scenario,
        stratum = k,
        d1 = dense[, 1],
        d2 = dense[, 2],
        f = true_surface(scenario, k, dense)
      )
      idx <- idx + 1
    }
  }
  out <- do.call(rbind, rows)
  write.csv(out, file.path(DIR_RESULTS, "figure1_true_surfaces_dense_grid.csv"), row.names = FALSE)
  invisible(out)
}

draw_figure1_sci <- function() {
  dense <- dense_candidate_grid()
  all_background <- unlist(lapply(1:4, function(s) {
    vals <- true_grid_values(s, dense)
    apply(vals, 2, min)
  }))
  zlim <- range(all_background, finite = TRUE)

  target <- file.path(DIR_FIGURES, "Figure1.pdf")

  for (out_file in target) {
    grDevices::pdf(out_file, width = 9.40, height = 7.05, useDingbats = FALSE,
                   family = "Helvetica", pointsize = 13)
    old <- graphics::par(no.readonly = TRUE)
    layout_matrix <- matrix(c(1, 2, 5,
                              3, 4, 5), nrow = 2, byrow = TRUE)
    graphics::layout(layout_matrix, widths = c(1, 1, 0.20), heights = c(1, 1))
    graphics::par(oma = c(0.65, 0.65, 1.05, 0.30), mar = c(3.25, 3.60, 2.45, 0.85),
                  mgp = c(2.20, 0.62, 0), tcl = -0.26, family = "Helvetica",
                  cex.lab = 1.05)

    for (scenario in 1:4) draw_panel(scenario, scenario_labels[scenario], zlim)
    draw_colorbar(zlim)

    graphics::par(old)
    grDevices::dev.off()
  }

  write_dense_surface_csv()
  message("Saved polished Figure 1 to:")
  message("  ", target)
  invisible(target)
}

draw_figure1_sci()

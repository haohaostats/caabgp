## Redraw Figure 2 in a polished colorblind-friendly SCI journal style.
##
## Usage from the project root:
##   & 'D:\Soft\R-4.5.2\bin\x64\Rscript.exe' redraw_figure2_sci_color.R
##
## Outputs:
##   images/Figures/Figure2_scenario1_results.pdf
##   images/Figures/Figure2_scenario1_results_SCI_color.pdf

source(file.path("R", "00_settings.R"))

DIR_FIGURES <- file.path(ROOT_DIR, "figures")
DIR_RESULTS <- file.path(ROOT_DIR, "results_r")
dir.create(DIR_FIGURES, recursive = TRUE, showWarnings = FALSE)

n_mc <- 1000

learning <- read.csv(file.path(DIR_RESULTS, sprintf("figure2_scenario1_learning_trajectory_n%d.csv", n_mc)))
final_dist <- read.csv(file.path(DIR_RESULTS, sprintf("figure2_scenario1_final_distribution_n%d.csv", n_mc)))
borrowing <- read.csv(file.path(DIR_RESULTS, sprintf("figure2_scenario1_borrowing_final_n%d.csv", n_mc)))
unique_doses <- read.csv(file.path(DIR_RESULTS, sprintf("figure2_scenario1_cumulative_unique_doses_n%d.csv", n_mc)))

METHOD_ORDER <- c("S-GP", "P-GP", "IND-GP", "AB-GP", "CA-AB-GP")

## High-saturation journal palette, with CA-AB-GP highlighted.
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
  n_vals <- sort(unique(dat$n))
  if (length(n_vals) <= max_points) return(dat)
  keep <- unique(round(seq(1, length(n_vals), length.out = max_points)))
  dat[dat$n %in% n_vals[keep], , drop = FALSE]
}

theme_panel <- function() {
  graphics::axis(1, cex.axis = 1.02, col = axis_col, col.axis = text_col, tck = -0.018)
  graphics::axis(2, cex.axis = 1.02, col = axis_col, col.axis = text_col, las = 1, tck = -0.018)
  graphics::box(col = axis_col, lwd = 1.10)
}

panel_title <- function(label) {
  graphics::mtext(label, side = 3, line = 0.70, adj = 0, cex = 1.22, font = 2, col = text_col)
}

place_endpoint_labels <- function(y, labels, ylim, min_gap_frac = 0.045) {
  ord <- order(y)
  y_adj <- y
  min_gap <- diff(ylim) * min_gap_frac
  if (length(y_adj) > 1) {
    for (i in 2:length(ord)) {
      prev <- ord[i - 1]
      cur <- ord[i]
      if (y_adj[cur] - y_adj[prev] < min_gap) y_adj[cur] <- y_adj[prev] + min_gap
    }
  }
  y_adj <- pmin(y_adj, ylim[2] - diff(ylim) * 0.035)
  data.frame(y = y_adj, label = labels, original_y = y)
}

plot_trajectory <- function(dat, metric, ylab, title, ylim = NULL, direct_label = FALSE) {
  mean_col <- paste0(metric, "_mean")
  dat <- dat[is.finite(dat[[mean_col]]), , drop = FALSE]
  dat <- thin_time_grid(dat, max_points = 24)

  if (is.null(ylim)) {
    yr <- range(dat[[mean_col]], na.rm = TRUE)
    pad <- diff(yr) * 0.08
    if (!is.finite(pad) || pad == 0) pad <- 0.05
    ylim <- yr + c(-pad, pad)
  }

  xlim <- if (isTRUE(direct_label)) {
    c(min(dat$n, na.rm = TRUE), max(dat$n, na.rm = TRUE) + 27)
  } else {
    range(dat$n, na.rm = TRUE)
  }

  graphics::plot(
    NA,
    xlim = xlim,
    ylim = ylim,
    xlab = "Sample size",
    ylab = ylab,
    axes = FALSE,
    xaxs = "i",
    yaxs = "i"
  )
  theme_panel()

  for (m in methods_present(dat)) {
    dm <- dat[dat$method == m, , drop = FALSE]
    dm <- dm[order(dm$n), , drop = FALSE]
    graphics::lines(dm$n, dm[[mean_col]], col = method_cols[m],
                    lty = line_types[m], lwd = line_widths[m])
  }

  if (isTRUE(direct_label)) {
    endpoints <- do.call(rbind, lapply(methods_present(dat), function(m) {
      dm <- dat[dat$method == m, , drop = FALSE]
      dm <- dm[order(dm$n), , drop = FALSE]
      data.frame(method = m, n = dm$n[nrow(dm)], y = dm[[mean_col]][nrow(dm)])
    }))
    labs <- place_endpoint_labels(endpoints$y, endpoints$method, ylim)
    label_x <- max(dat$n, na.rm = TRUE) + 4.0
    for (i in seq_len(nrow(endpoints))) {
      graphics::segments(endpoints$n[i], endpoints$y[i], label_x - 0.9, labs$y[i],
                         col = method_cols[endpoints$method[i]], lwd = 0.70)
      graphics::text(label_x, labs$y[i], labels = labs$label[i],
                     adj = 0, cex = 0.60, col = method_cols[endpoints$method[i]], font = 2)
    }
  }

  panel_title(title)
}

plot_final_regret_box <- function(dat) {
  dat$method <- factor(dat$method, levels = METHOD_ORDER)
  graphics::boxplot(
    regret ~ method,
    data = dat,
    outline = FALSE,
    axes = FALSE,
    col = method_fills[METHOD_ORDER],
    border = "#222222",
    medcol = "#000000",
    medlwd = 1.35,
    whisklty = 1,
    staplewex = 0.45,
    boxwex = 0.55,
    xlab = "",
    ylab = "Final weighted regret",
    ylim = range(dat$regret, na.rm = TRUE) * c(0.96, 1.04)
  )
  graphics::axis(1, at = seq_along(METHOD_ORDER), labels = METHOD_ORDER,
                 cex.axis = 0.92, las = 2, col = axis_col, col.axis = text_col, tck = -0.018)
  graphics::axis(2, cex.axis = 1.02, col = axis_col, col.axis = text_col, las = 1, tck = -0.018)
  graphics::stripchart(
    regret ~ method,
    data = dat,
    vertical = TRUE,
    method = "jitter",
    jitter = 0.12,
    add = TRUE,
    pch = 16,
    cex = 0.28,
    col = grDevices::adjustcolor("#111111", alpha.f = 0.18)
  )
  graphics::box(col = axis_col, lwd = 1.10)
  panel_title("B  Final regret distribution")
}

plot_borrowing <- function(dat) {
  dat <- dat[dat$method %in% c("AB-GP", "CA-AB-GP"), , drop = FALSE]
  dat$method <- factor(dat$method, levels = c("AB-GP", "CA-AB-GP"))
  dat <- dat[order(dat$method), , drop = FALSE]
  x <- seq_len(nrow(dat))
  se <- se0(dat$rho_se)

  graphics::plot(
    NA,
    xlim = c(0.55, length(x) + 0.45),
    ylim = c(0, 1),
    xlab = "",
    ylab = expression("Borrowing index " * widehat(rho)[b]),
    axes = FALSE,
    xaxs = "i",
    yaxs = "i"
  )
  graphics::axis(1, at = x, labels = as.character(dat$method),
                 cex.axis = 1.02, col = axis_col, col.axis = text_col, tck = -0.018)
  graphics::axis(2, cex.axis = 1.02, col = axis_col, col.axis = text_col, las = 1, tck = -0.018)
  graphics::segments(x, pmax(0, dat$rho_mean - 1.96 * se),
                     x, pmin(1, dat$rho_mean + 1.96 * se),
                     col = method_cols[as.character(dat$method)], lwd = 1.55)
  graphics::segments(x - 0.08, pmax(0, dat$rho_mean - 1.96 * se),
                     x + 0.08, pmax(0, dat$rho_mean - 1.96 * se),
                     col = method_cols[as.character(dat$method)], lwd = 1.55)
  graphics::segments(x - 0.08, pmin(1, dat$rho_mean + 1.96 * se),
                     x + 0.08, pmin(1, dat$rho_mean + 1.96 * se),
                     col = method_cols[as.character(dat$method)], lwd = 1.55)
  graphics::points(x, dat$rho_mean, pch = 21, bg = method_fills[as.character(dat$method)],
                   col = method_cols[as.character(dat$method)], cex = 1.90, lwd = 1.35)
  graphics::box(col = axis_col, lwd = 1.10)
  panel_title(expression("C  Estimated borrowing index " * widehat(rho)[b]))
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
    bty = "n",
    horiz = TRUE,
    x.intersp = 0.80,
    seg.len = 3.10,
    cex = 1.08
  )
}

draw_figure2_color <- function() {
  target <- file.path(DIR_FIGURES, "Figure2.pdf")

  regret_ylim <- range(learning$regret_mean, na.rm = TRUE)
  regret_ylim <- c(0, regret_ylim[2] * 1.08)
  unique_ylim <- range(unique_doses$unique_doses_mean, na.rm = TRUE)
  unique_ylim <- c(max(0, unique_ylim[1] * 0.92), unique_ylim[2] * 1.08)

  written <- character(0)
  for (out_file in target) {
    ok <- tryCatch({
      grDevices::pdf(out_file, width = 9.40, height = 8.10, useDingbats = FALSE,
                     family = "Helvetica", pointsize = 13)
      TRUE
    }, error = function(e) {
      warning("Could not open output PDF. Close it if it is open in a viewer: ", out_file)
      FALSE
    })
    if (!ok) next

    old <- graphics::par(no.readonly = TRUE)
    layout_matrix <- matrix(c(1, 2,
                              3, 4,
                              5, 5), nrow = 3, byrow = TRUE)
    graphics::layout(layout_matrix, heights = c(1, 1, 0.17), widths = c(1, 1))
    graphics::par(oma = c(0.25, 0.25, 0.35, 0.25),
                  mar = c(4.20, 4.50, 2.30, 1.25),
                  mgp = c(2.65, 0.78, 0),
                  tcl = -0.28,
                  family = "Helvetica",
                  cex.lab = 1.12)

    plot_trajectory(learning, "regret", "Weighted regret", "A  Regret trajectory", ylim = regret_ylim)
    plot_final_regret_box(final_dist)
    plot_borrowing(borrowing)
    plot_trajectory(unique_doses, "unique_doses", "Unique dose combinations",
                    "D  Cumulative unique doses", ylim = unique_ylim, direct_label = FALSE)
    plot_legend_strip()

    graphics::par(old)
    grDevices::dev.off()
    written <- c(written, out_file)
  }

  message("Saved color SCI Figure 2 to:")
  for (f in written) message("  ", f)
  invisible(written)
}

draw_figure2_color()

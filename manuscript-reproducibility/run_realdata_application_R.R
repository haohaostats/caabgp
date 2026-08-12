source(file.path("R", "07_application_realdata.R"))

METHOD_ORDER_APP <- APP_SETTINGS$methods
method_cols_app <- c("P-GP" = "#0057B8", "AB-GP" = "#2CA02C", "CA-AB-GP" = "#D62728")
method_fills_app <- c("P-GP" = "#8AB6F0", "AB-GP" = "#8ED081", "CA-AB-GP" = "#F28B82")
stratum_cols_app <- c("1" = "#0057B8", "2" = "#D62728")

fmt <- function(x, digits = 3) sprintf(paste0("%.", digits, "f"), x)
fmt_pct <- function(x) paste0(sprintf("%.1f", 100 * x), "\\%")
dose_level_label <- function(x) paste0("L", as.integer(round(x / 0.25)) + 1L)
dose_pair_label <- function(d1, d2) paste0("(", dose_level_label(d1), ",", dose_level_label(d2), ")")

write_application_table <- function(recs) {
  recs$method <- factor(recs$method, levels = METHOD_ORDER_APP)
  recs <- recs[order(recs$method, recs$stratum), ]
  lines <- c(
    "\\begin{table*}[t]",
    "\t\\centering",
    "\t\\caption{Final recommendations and operating summaries in the real-data-calibrated application. Values are based on the final analysis of the sequential replay study.}",
    "\t\\label{tab:application_results}",
    "\t\\footnotesize",
    "\t\\setlength{\\tabcolsep}{4pt}",
    "\t\\begin{tabular}{llcccccccc}",
    "\t\t\\toprule",
    "\t\t\\textbf{Design} & \\textbf{Stratum} & $\\boldsymbol{\\widehat{\\bm d}_k}$ & $\\boldsymbol{\\mu_{T,k}(\\widehat{\\bm d}_k)}$ & $\\boldsymbol{s_{T,k}(\\widehat{\\bm d}_k)}$ & \\textbf{Final sample size} & \\textbf{Unique doses} & \\textbf{Repeat allocation} & \\textbf{Total cost} & $\\boldsymbol{\\widehat\\rho_b}$ \\\\",
    "\t\t\\midrule"
  )
  for (m in METHOD_ORDER_APP) {
    dm <- recs[recs$method == m, ]
    rho <- if (all(is.na(dm$rho))) "--" else fmt(dm$rho[1])
    for (i in seq_len(nrow(dm))) {
      d <- dose_pair_label(dm$d1_hat[i], dm$d2_hat[i])
      common <- if (i == 1) {
        paste0(
          "\\multirow{2}{*}{", m, "} & Stratum ", dm$stratum[i], " & $", d, "$ & ",
          fmt(dm$mu_hat[i]), " & ", fmt(dm$sd_hat[i]), " & \\multirow{2}{*}{", dm$n[i],
          "} & \\multirow{2}{*}{", dm$unique_doses[i],
          "} & \\multirow{2}{*}{", fmt_pct(dm$repeat_fraction[i]),
          "} & \\multirow{2}{*}{", fmt(dm$total_cost[i]),
          "} & ", if (rho == "--") "--" else paste0("\\multirow{2}{*}{", rho, "}"), " \\\\"
        )
      } else {
        paste0("& Stratum ", dm$stratum[i], " & $", d, "$ & ",
               fmt(dm$mu_hat[i]), " & ", fmt(dm$sd_hat[i]), " & & & & & ",
               if (rho == "--") "--" else "", " \\\\")
      }
      lines <- c(lines, paste0("\t\t", common))
    }
    if (m != tail(METHOD_ORDER_APP, 1)) lines <- c(lines, "\t\t\\addlinespace")
  }
  lines <- c(
    lines,
    "\t\t\\bottomrule",
    "\t\\end{tabular}",
    "\t\\vspace{0.5em}",
    "\t\\begin{minipage}{0.95\\textwidth}",
    "\t\t\\footnotesize",
    "\t\t\\textit{Note:} $\\widehat{\\bm d}_k$ denotes the final recommended dose-combination level for stratum $k$. Dose levels L1--L5 correspond to the ordered prespecified feasible levels of each treatment component from lowest to highest. $\\mu_{T,k}(\\widehat{\\bm d}_k)$ and $s_{T,k}(\\widehat{\\bm d}_k)$ are the empirical-Bayes predictive mean and standard deviation at the recommended dose. Repeat allocation is the percentage of post-initial allocations assigned to dose combinations that had already been evaluated. The borrowing index $\\widehat\\rho_b$ is reported only for adaptive-borrowing designs.",
    "\t\\end{minipage}",
    "\\end{table*}"
  )
  out <- file.path(DIR_RESULTS, "table_application_results.tex")
  writeLines(lines, out)
  out
}

draw_axis <- function() {
  axis(1, at = SIM_SETTINGS$dose_values, labels = paste0("L", seq_along(SIM_SETTINGS$dose_values)),
       cex.axis = 1.0, tck = -0.018, col = "#222222")
  axis(2, at = SIM_SETTINGS$dose_values, labels = paste0("L", seq_along(SIM_SETTINGS$dose_values)),
       cex.axis = 1.0, las = 1, tck = -0.018, col = "#222222")
  box(col = "#222222", lwd = 1.1)
}

draw_callout_label <- function(x, y, label, fill = "#111111", text_col = "#FFFFFF", cex = 0.92) {
  pad_x <- 0.035
  pad_y <- 0.028
  w <- strwidth(label, cex = cex, font = 2) + pad_x
  h <- strheight(label, cex = cex, font = 2) + pad_y
  rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2,
       col = fill, border = "#FFFFFF", lwd = 0.75)
  text(x, y, labels = label, cex = cex, font = 2, col = text_col)
}

recommendation_offset <- function(stratum) {
  off <- 0.022
  if (stratum == 1L) return(c(-off, -off))
  if (stratum == 2L) return(c(off, off))
  c(0, 0)
}

recommendation_label_position <- function(x, y, stratum) {
  dx <- if (x < 0.14) 0.12 else if (x > 0.86) -0.12 else if (stratum == 1L) -0.10 else 0.10
  dy <- if (y < 0.14) 0.12 else if (y > 0.86) -0.10 else if (stratum == 1L) -0.085 else 0.095
  c(
    pmin(1.02, pmax(-0.02, x + dx)),
    pmin(1.02, pmax(-0.02, y + dy))
  )
}

plot_application_allocation <- function(outputs) {
  out <- file.path(DIR_FIGURES, "Figure6.pdf")
  pdf(out, width = 9.4, height = 3.9, family = "Helvetica", useDingbats = FALSE, pointsize = 13)
  old <- par(no.readonly = TRUE)
  on.exit({ par(old); dev.off() }, add = TRUE)
  par(mfrow = c(1, 3), mar = c(4.2, 4.1, 2.4, 1.0), mgp = c(2.45, 0.75, 0), tcl = -0.28, cex.lab = 1.1)
  grid <- candidate_grid()
  for (m in METHOD_ORDER_APP) {
    dat <- outputs[[m]]$data
    rec <- outputs[[m]]$recommendations
    counts <- aggregate(stratum ~ d1 + d2, dat, length)
    names(counts)[3] <- "n"
    plot(NA, xlim = c(-0.06, 1.06), ylim = c(-0.06, 1.06), xlab = "Component 1 level", ylab = "Component 2 level",
         axes = FALSE, xaxs = "i", yaxs = "i")
    draw_axis()
    points(grid[, 1], grid[, 2], pch = 16, col = "#CFCFCF", cex = 0.82)
    points(counts$d1, counts$d2, pch = 21, bg = method_fills_app[m], col = method_cols_app[m],
           cex = 0.75 + 0.20 * sqrt(counts$n), lwd = 1.35)
    for (i in seq_len(nrow(rec))) {
      off <- recommendation_offset(rec$stratum[i])
      x_rec <- rec$d1_hat[i] + off[1]
      y_rec <- rec$d2_hat[i] + off[2]
      points(x_rec, y_rec, pch = 24, bg = stratum_cols_app[as.character(rec$stratum[i])],
             col = "#111111", cex = 1.65, lwd = 1.1)
      lab <- recommendation_label_position(x_rec, y_rec, rec$stratum[i])
      arrows(lab[1], lab[2], x_rec, y_rec, col = "#111111",
             lwd = 0.95, length = 0.075, angle = 22)
      draw_callout_label(lab[1], lab[2], paste0("Z", rec$stratum[i]),
                         fill = "#111111", text_col = "#FFFFFF", cex = 0.82)
    }
    mtext(m, side = 3, line = 0.65, adj = 0, cex = 1.18, font = 2)
  }
  out
}

surface_cols_app <- grDevices::colorRampPalette(c(
  "#06152B", "#0B2C5D", "#155A8A", "#14889A",
  "#51B89F", "#B8D99A", "#F3E7B1", "#FFF7D6"
))(256)

plot_application_surfaces <- function(pred, recs) {
  out <- file.path(DIR_FIGURES, "Figure7.pdf")
  pred <- pred[pred$method == "CA-AB-GP", ]
  zlim <- range(pred$mu, finite = TRUE)
  pdf(out, width = 8.8, height = 4.2, family = "Helvetica", useDingbats = FALSE, pointsize = 13)
  old <- par(no.readonly = TRUE)
  on.exit({ par(old); dev.off() }, add = TRUE)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.12))
  par(mar = c(4.2, 4.1, 2.4, 0.9), mgp = c(2.45, 0.75, 0), tcl = -0.28, cex.lab = 1.1)
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
    points(rk$d1_hat, rk$d2_hat, pch = 24, bg = stratum_cols_app[as.character(k)],
           col = "#FFFFFF", cex = 1.95, lwd = 2.1)
    points(rk$d1_hat, rk$d2_hat, pch = 24, bg = NA,
           col = "#111111", cex = 2.08, lwd = 0.9)
    lab_x <- pmin(1.00, pmax(0.08, rk$d1_hat + ifelse(rk$d1_hat < 0.80, 0.16, -0.16)))
    lab_y <- pmin(1.00, pmax(0.08, rk$d2_hat + ifelse(rk$d2_hat < 0.80, 0.12, -0.12)))
    arrows(lab_x, lab_y, rk$d1_hat, rk$d2_hat, col = "#FFFFFF",
           lwd = 3.2, length = 0.08, angle = 22)
    arrows(lab_x, lab_y, rk$d1_hat, rk$d2_hat, col = "#111111",
           lwd = 1.0, length = 0.08, angle = 22)
    draw_callout_label(lab_x, lab_y, paste0("Z", k), fill = "#111111",
                       text_col = "#FFFFFF", cex = 0.92)
    mtext(paste0("Stratum ", k), side = 3, line = 0.65, adj = 0, cex = 1.18, font = 2)
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

cat("Running real-data-calibrated application...\n")
app <- run_realdata_application()
table_path <- write_application_table(app$recommendations)
fig6 <- plot_application_allocation(app$outputs)
fig7 <- plot_application_surfaces(app$predictions, app$recommendations)
cat("Application results written to:\n")
cat("  ", normalizePath(table_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(fig6, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(fig7, winslash = "/", mustWork = FALSE), "\n", sep = "")

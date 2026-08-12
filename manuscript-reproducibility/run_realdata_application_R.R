source(file.path("R", "07_application_realdata.R"))
source(file.path("R", "08_table_rendering.R"))

METHOD_ORDER_APP <- APP_SETTINGS$methods
method_cols_app <- c("P-GP" = "#0057B8", "AB-GP" = "#2CA02C", "CA-AB-GP" = "#D62728")
method_fills_app <- c("P-GP" = "#8AB6F0", "AB-GP" = "#8ED081", "CA-AB-GP" = "#F28B82")
stratum_cols_app <- c("1" = "#0057B8", "2" = "#D62728")
surface_cols_app <- grDevices::colorRampPalette(c(
  "#06152B", "#0B2C5D", "#155A8A", "#14889A",
  "#51B89F", "#B8D99A", "#F3E7B1", "#FFF7D6"
))(256)

fmt <- function(x, digits = 3) sprintf(paste0("%.", digits, "f"), x)
fmt_pct <- function(x) paste0(sprintf("%.1f", 100 * x), "%")
dose_level_label <- function(x) paste0("L", as.integer(round(x / 0.25)) + 1L)
dose_pair_label <- function(d1, d2) paste0("(", dose_level_label(d1), ",", dose_level_label(d2), ")")

draw_axis <- function() {
  axis(1, at = SIM_SETTINGS$dose_values,
       labels = paste0("L", seq_along(SIM_SETTINGS$dose_values)),
       cex.axis = 1.0, tck = -0.018, col = "#222222")
  axis(2, at = SIM_SETTINGS$dose_values,
       labels = paste0("L", seq_along(SIM_SETTINGS$dose_values)),
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

write_application_table <- function(recs) {
  recs$method <- factor(recs$method, levels = METHOD_ORDER_APP)
  recs <- recs[order(recs$method, recs$stratum), ]
  display <- data.frame(
    Design = as.character(recs$method),
    Stratum = paste0("Z", recs$stratum),
    `Recommended dose` = dose_pair_label(recs$d1_hat, recs$d2_hat),
    `Predictive mean` = fmt(recs$mu_hat),
    `Predictive SD` = fmt(recs$sd_hat),
    `Final sample size` = recs$n,
    `Unique doses` = recs$unique_doses,
    `Repeat allocation` = fmt_pct(recs$repeat_fraction),
    `Total cost` = fmt(recs$total_cost),
    rho_b = ifelse(is.na(recs$rho), "--", fmt(recs$rho)),
    check.names = FALSE
  )
  out <- file.path(DIR_RESULTS, "table_application_results.pdf")
  render_table_pdf(
    display, out,
    "Table 9. Final recommendations and operating summaries in the real-data-calibrated application.",
    paste(
      "Values are based on the final analysis of the sequential replay study.",
      "Dose levels L1--L5 run from the lowest to highest prespecified level.",
      "Repeat allocation is the percentage of post-initial allocations assigned to previously evaluated combinations."
    ),
    column_widths = c(1.15, 0.85, 1.45, 1.25, 1.10, 1.35, 1.10, 1.35, 1.15, 0.95),
    width = 13.0, height = 5.0, header_size = 7.8, body_size = 8.3
  )
  out
}

cat("Running real-data-calibrated application...\n")
app <- run_realdata_application()
table_path <- write_application_table(app$recommendations)
cat("Application results written to:\n")
cat("  ", normalizePath(table_path, winslash = "/", mustWork = FALSE), "\n", sep = "")

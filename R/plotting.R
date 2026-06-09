theme_caabgp <- function(base_size = 14, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.25, margin = ggplot2::margin(b = 8)),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "#111111"),
      axis.line = ggplot2::element_line(linewidth = 0.55, color = "#111111"),
      axis.ticks = ggplot2::element_line(linewidth = 0.45, color = "#111111"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

plot_allocation <- function(result, design = result$design, recommendations = NULL,
                            main = NULL, base_size = 14) {
  validate_design(design)
  if (inherits(result, "caabgp_trial")) {
    dat <- result$data
    recommendations <- null_coalesce(recommendations, result$recommendations)
  } else {
    dat <- as.data.frame(result)
  }
  grid <- as.data.frame(design$dose_grid)
  X <- as_dose_matrix(grid, design$dose_cols)
  if (ncol(X) != 2) caabgp_stop("plot_allocation currently supports two dose dimensions.")

  counts <- stats::aggregate(dat[[design$dose_cols[1]]], dat[, design$dose_cols, drop = FALSE], length)
  names(counts)[ncol(counts)] <- "n"
  names(grid) <- design$dose_cols

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = grid,
      ggplot2::aes(x = .data[[design$dose_cols[1]]], y = .data[[design$dose_cols[2]]]),
      color = "#C7C7C7",
      size = 2.3
    ) +
    ggplot2::geom_point(
      data = counts,
      ggplot2::aes(
        x = .data[[design$dose_cols[1]]],
        y = .data[[design$dose_cols[2]]],
        size = .data$n
      ),
      shape = 21,
      fill = "#D62728",
      color = "#7F1D1D",
      stroke = 0.9,
      alpha = 0.82
    ) +
    ggplot2::scale_size_area(max_size = 9, name = "Patients") +
    ggplot2::coord_equal(
      xlim = range(X[, 1]),
      ylim = range(X[, 2]),
      expand = TRUE
    ) +
    ggplot2::labs(
      title = null_coalesce(main, "Sequential allocation map"),
      x = design$dose_cols[1],
      y = design$dose_cols[2]
    ) +
    theme_caabgp(base_size = base_size)

  if (!is.null(recommendations)) {
    rec <- as.data.frame(recommendations)
    p <- p +
      ggplot2::geom_point(
        data = rec,
        ggplot2::aes(x = .data[[design$dose_cols[1]]], y = .data[[design$dose_cols[2]]]),
        shape = 24,
        fill = "#0057B8",
        color = "#111111",
        stroke = 0.9,
        size = 4.2
      ) +
      ggplot2::geom_label(
        data = rec,
        ggplot2::aes(
          x = .data[[design$dose_cols[1]]],
          y = .data[[design$dose_cols[2]]],
          label = paste0("Z", .data$stratum)
        ),
        nudge_x = 0.035,
        nudge_y = 0.035,
        fill = "white",
        color = "#111111",
        label.size = 0.25,
        size = base_size / 4
      )
  }
  p
}

plot_surface <- function(fit, stratum = 1, design = fit$design, recommendations = NULL,
                         main = NULL, base_size = 14) {
  validate_design(design)
  pred <- if (inherits(fit, "caabgp_trial")) fit$predictions else predict_caabgp(fit, design$dose_grid)
  grid <- as_dose_matrix(design$dose_grid, design$dose_cols)
  if (ncol(grid) != 2) caabgp_stop("plot_surface currently supports two dose dimensions.")
  pk <- pred[pred$stratum == stratum, , drop = FALSE]
  if (!nrow(pk)) caabgp_stop("No predictions found for stratum ", stratum, ".")

  p <- ggplot2::ggplot(
    pk,
    ggplot2::aes(
      x = .data[[design$dose_cols[1]]],
      y = .data[[design$dose_cols[2]]],
      fill = .data$mu
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_contour(
      ggplot2::aes(z = .data$mu),
      color = "white",
      linewidth = 0.42,
      alpha = 0.85,
      bins = 8
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#062B4F", "#145C7F", "#2A9D8F", "#BBD7A2", "#FFF3C4"),
      name = "Predicted\noutcome"
    ) +
    ggplot2::coord_equal(
      xlim = range(grid[, 1]),
      ylim = range(grid[, 2]),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = null_coalesce(main, paste("Predicted surface, stratum", stratum)),
      x = design$dose_cols[1],
      y = design$dose_cols[2]
    ) +
    theme_caabgp(base_size = base_size)

  rec <- recommendations
  if (is.null(rec) && inherits(fit, "caabgp_trial")) rec <- fit$recommendations
  if (is.null(rec) && inherits(fit, "caabgp_fit")) rec <- recommend_dose(fit, design)
  if (!is.null(rec)) {
    rec <- as.data.frame(rec)
    rec <- rec[rec$stratum == stratum, , drop = FALSE]
    if (nrow(rec)) {
      p <- p +
        ggplot2::geom_point(
          data = rec,
          ggplot2::aes(x = .data[[design$dose_cols[1]]], y = .data[[design$dose_cols[2]]]),
          inherit.aes = FALSE,
          shape = 24,
          fill = "#E31A1C",
          color = "#111111",
          stroke = 0.9,
          size = 4.5
        ) +
        ggplot2::geom_label(
          data = rec,
          ggplot2::aes(
            x = .data[[design$dose_cols[1]]],
            y = .data[[design$dose_cols[2]]],
            label = paste0("Z", .data$stratum)
          ),
          inherit.aes = FALSE,
          nudge_x = 0.04,
          nudge_y = 0.04,
          fill = "white",
          color = "#111111",
          label.size = 0.25,
          size = base_size / 4
        )
    }
  }
  p
}

save_caabgp_figure <- function(plot, filename, width = 7, height = 5,
                               units = "in", dpi = 600, device = NULL) {
  if (is.null(device)) {
    ext <- tolower(tools::file_ext(filename))
    device <- if (ext %in% c("pdf", "svg", "eps", "ps")) ext else NULL
  }
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    device = device,
    bg = "white",
    limitsize = FALSE
  )
  invisible(normalizePath(filename, winslash = "/", mustWork = FALSE))
}

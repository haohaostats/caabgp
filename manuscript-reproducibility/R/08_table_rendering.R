## R-native vector table rendering for manuscript outputs.
## Uses only base R and the recommended grid package.

open_table_pdf <- function(path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(path, width = width, height = height,
                         family = "sans", pointsize = 10, onefile = TRUE)
  } else {
    grDevices::pdf(path, width = width, height = height,
                   family = "Helvetica", pointsize = 10,
                   useDingbats = FALSE, onefile = TRUE)
  }
}

wrap_table_text <- function(x, width) {
  paste(strwrap(x, width = width), collapse = "\n")
}

render_table_pdf <- function(data, path, title, note,
                             column_widths = NULL,
                             width = 12.5, height = 4.5,
                             header_size = 8.2, body_size = 8.0) {
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  data[] <- lapply(data, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "--"
    x
  })
  nr <- nrow(data)
  nc <- ncol(data)
  if (!nr || !nc) stop("Cannot render an empty table.")

  if (is.null(column_widths)) column_widths <- rep(1, nc)
  if (length(column_widths) != nc || any(column_widths <= 0)) {
    stop("column_widths must contain one positive value per column.")
  }
  column_widths <- column_widths / sum(column_widths)
  x_left <- c(0, cumsum(column_widths)[-nc])
  x_mid <- x_left + column_widths / 2

  open_table_pdf(path, width, height)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()

  navy <- "#17365D"
  navy_light <- "#DCE6F1"
  stripe <- "#F5F7FA"
  border <- "#AEB8C4"
  text <- "#17212B"

  grid::grid.text(
    wrap_table_text(title, 125),
    x = grid::unit(0.015, "npc"), y = grid::unit(0.965, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontfamily = "sans", fontface = "bold", fontsize = 11.5, col = text)
  )

  table_top <- 0.845
  table_bottom <- if (nzchar(note)) 0.165 else 0.08
  row_height <- (table_top - table_bottom) / (nr + 1)

  grid::grid.rect(
    x = 0.5, y = table_top - row_height / 2,
    width = 0.98, height = row_height,
    gp = grid::gpar(fill = navy, col = navy, lwd = 0.8)
  )
  for (j in seq_len(nc)) {
    grid::grid.text(
      names(data)[j],
      x = grid::unit(0.01 + 0.98 * x_mid[j], "npc"),
      y = grid::unit(table_top - row_height / 2, "npc"),
      gp = grid::gpar(fontfamily = "sans", fontface = "bold",
                     fontsize = header_size, col = "white")
    )
  }

  for (i in seq_len(nr)) {
    y <- table_top - (i + 0.5) * row_height
    fill <- if (i %% 2L) "white" else stripe
    grid::grid.rect(
      x = 0.5, y = y, width = 0.98, height = row_height,
      gp = grid::gpar(fill = fill, col = border, lwd = 0.35)
    )
    for (j in seq_len(nc)) {
      face <- if (j == 1L) "bold" else "plain"
      col <- if (j == 1L && data[i, j] == "CA-AB-GP") "#C62828" else text
      grid::grid.text(
        data[i, j],
        x = grid::unit(0.01 + 0.98 * x_mid[j], "npc"), y = grid::unit(y, "npc"),
        gp = grid::gpar(fontfamily = "sans", fontface = face,
                       fontsize = body_size, col = col)
      )
    }
  }

  grid::grid.lines(
    x = grid::unit(c(0.01, 0.99), "npc"),
    y = grid::unit(c(table_top, table_top), "npc"),
    gp = grid::gpar(col = navy, lwd = 1.3)
  )
  grid::grid.lines(
    x = grid::unit(c(0.01, 0.99), "npc"),
    y = grid::unit(c(table_bottom, table_bottom), "npc"),
    gp = grid::gpar(col = navy, lwd = 1.1)
  )

  if (nzchar(note)) {
    grid::grid.text(
      wrap_table_text(note, 180),
      x = grid::unit(0.015, "npc"), y = grid::unit(0.115, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontfamily = "sans", fontsize = 7.8, col = "#3E4A56")
    )
  }
  invisible(path)
}

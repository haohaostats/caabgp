caabgp_stop <- function(...) stop(..., call. = FALSE)

utils::globalVariables(c(".data"))

check_numeric_vector <- function(x, name, n = NULL, positive = FALSE, nonnegative = FALSE) {
  if (!is.numeric(x) || any(!is.finite(x))) caabgp_stop(name, " must be a finite numeric vector.")
  if (!is.null(n) && length(x) != n) caabgp_stop(name, " must have length ", n, ".")
  if (positive && any(x <= 0)) caabgp_stop(name, " must be positive.")
  if (nonnegative && any(x < 0)) caabgp_stop(name, " must be nonnegative.")
  x
}

as_dose_matrix <- function(x, dose_cols = NULL) {
  if (is.matrix(x)) {
    X <- x
  } else {
    x <- as.data.frame(x)
    if (is.null(dose_cols)) {
      dose_cols <- grep("^d[0-9]+$", names(x), value = TRUE)
      if (!length(dose_cols)) dose_cols <- names(x)
    }
    X <- as.matrix(x[, dose_cols, drop = FALSE])
  }
  storage.mode(X) <- "double"
  if (!nrow(X) || !ncol(X) || any(!is.finite(X))) caabgp_stop("Dose grid/data must contain finite dose values.")
  X
}

dose_key <- function(dose, digits = 6) {
  dose <- as.matrix(dose)
  apply(dose, 1, function(z) paste(format(round(z, digits), nsmall = 2, trim = TRUE), collapse = "_"))
}

unique_dose_count <- function(data, dose_cols) {
  length(unique(dose_key(as_dose_matrix(data, dose_cols))))
}

match_dose_row <- function(grid, dose) {
  grid <- as.matrix(grid)
  dose <- as.numeric(dose)
  hit <- which(rowSums(abs(sweep(grid, 2, dose, "-"))) < 1e-10)
  if (!length(hit)) caabgp_stop("Dose is not present in the design dose grid.")
  hit[1]
}

null_coalesce <- function(x, y) if (is.null(x)) y else x

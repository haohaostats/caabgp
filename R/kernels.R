scaled_sqdist <- function(X1, X2, ell) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  ell <- as.numeric(ell)
  out <- matrix(0, nrow(X1), nrow(X2))
  for (j in seq_len(ncol(X1))) {
    out <- out + outer(X1[, j], X2[, j], "-")^2 / (2 * ell[j]^2)
  }
  out
}

se_kernel <- function(X1, X2, ell) {
  exp(-scaled_sqdist(X1, X2, ell))
}

safe_chol <- function(Sigma) {
  jitter <- 1e-8
  for (attempt in seq_len(8)) {
    L <- tryCatch(chol(Sigma + diag(jitter, nrow(Sigma))), error = function(e) NULL)
    if (!is.null(L)) return(list(L = L, jitter = jitter))
    jitter <- jitter * 10
  }
  caabgp_stop("Covariance matrix is not positive definite after jitter escalation.")
}

chol_solve <- function(L, B) {
  backsolve(L, forwardsolve(t(L), B))
}

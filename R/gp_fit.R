unpack_params <- function(par, model, J) {
  if (model == "S") {
    list(nu = exp(par[1]), ell = exp(par[seq_len(J) + 1]), sigma = exp(par[J + 2]))
  } else if (model == "P") {
    list(
      nu = exp(par[1]),
      ell = exp(par[seq_len(J) + 1]),
      ell_z = exp(par[J + 2]),
      sigma = exp(par[J + 3])
    )
  } else if (model == "AB") {
    list(
      nu0 = exp(par[1]),
      nub = exp(par[2]),
      ell0 = exp(par[seq_len(J) + 2]),
      ellb = exp(par[seq_len(J) + J + 2]),
      sigma = exp(par[2 * J + 3])
    )
  } else {
    caabgp_stop("Unknown GP model: ", model)
  }
}

covariance_matrix <- function(X, z, par, model, K, add_noise = TRUE) {
  X <- as.matrix(X)
  J <- ncol(X)
  p <- unpack_params(par, model, J)
  n <- nrow(X)
  if (model == "S") {
    Sigma <- p$nu * se_kernel(X, X, p$ell)
  } else if (model == "P") {
    Kd <- se_kernel(X, X, p$ell)
    Dz <- outer(z, z, "!=") * 2
    Kz <- exp(-Dz / (2 * p$ell_z^2))
    Sigma <- p$nu * Kd * Kz
  } else if (model == "AB") {
    Sigma <- p$nu0 * se_kernel(X, X, p$ell0)
    Kb <- se_kernel(X, X, p$ellb)
    for (k in seq_len(K)) {
      idx <- which(z == k)
      if (length(idx)) Sigma[idx, idx] <- Sigma[idx, idx] + p$nub * Kb[idx, idx]
    }
  } else {
    caabgp_stop("Unknown GP model: ", model)
  }
  if (add_noise) Sigma <- Sigma + diag(p$sigma^2, n)
  Sigma
}

cross_covariance <- function(fit, Xnew, znew) {
  X <- fit$X
  z <- fit$z
  p <- fit$params
  model <- fit$model
  Xnew <- as.matrix(Xnew)
  if (model == "S") return(p$nu * se_kernel(X, Xnew, p$ell))
  if (model == "P") {
    Kd <- se_kernel(X, Xnew, p$ell)
    Dz <- outer(z, znew, "!=") * 2
    Kz <- exp(-Dz / (2 * p$ell_z^2))
    return(p$nu * Kd * Kz)
  }
  if (model == "AB") {
    C <- p$nu0 * se_kernel(X, Xnew, p$ell0)
    Kb <- se_kernel(X, Xnew, p$ellb)
    for (k in seq_len(fit$K)) {
      rows <- which(z == k)
      cols <- which(znew == k)
      if (length(rows) && length(cols)) C[rows, cols] <- C[rows, cols] + p$nub * Kb[rows, cols]
    }
    return(C)
  }
  caabgp_stop("Unknown GP model: ", model)
}

latent_variance_diag <- function(fit, znew) {
  p <- fit$params
  if (fit$model == "S") return(rep(p$nu, length(znew)))
  if (fit$model == "P") return(rep(p$nu, length(znew)))
  if (fit$model == "AB") return(rep(p$nu0 + p$nub, length(znew)))
  caabgp_stop("Unknown GP model: ", fit$model)
}

gp_nll <- function(par, X, z, y, model, K) {
  Sigma <- tryCatch(covariance_matrix(X, z, par, model, K, add_noise = TRUE), error = function(e) NULL)
  if (is.null(Sigma)) return(1e50)
  chol_info <- tryCatch(safe_chol(Sigma), error = function(e) NULL)
  if (is.null(chol_info)) return(1e50)
  L <- chol_info$L
  one <- rep(1, length(y))
  inv_y <- chol_solve(L, y)
  inv_one <- chol_solve(L, one)
  beta <- as.numeric(sum(inv_one * y) / sum(inv_one))
  resid <- y - beta
  alpha <- chol_solve(L, resid)
  logdet <- 2 * sum(log(diag(L)))
  0.5 * sum(resid * alpha) + 0.5 * logdet + 0.5 * length(y) * log(2 * pi)
}

default_starts <- function(y, model, J) {
  raw_var <- if (length(y) > 1) stats::var(y) else NA_real_
  vy <- if (is.finite(raw_var)) max(raw_var, 0.05^2) else 0.10
  ell_a <- rep(0.35, J)
  ell_b <- rep(0.55, J)
  ell_c <- rep(0.22, J)
  if (model == "S") {
    starts <- rbind(c(0.70 * vy, ell_a, 0.20), c(0.50 * vy, ell_b, 0.15), c(1.00 * vy, ell_c, 0.25))
  } else if (model == "P") {
    starts <- rbind(c(0.70 * vy, ell_a, 0.75, 0.20), c(0.50 * vy, ell_b, 0.45, 0.15), c(1.00 * vy, ell_c, 1.25, 0.25))
  } else if (model == "AB") {
    starts <- rbind(
      c(0.50 * vy, 0.20 * vy, ell_a, ell_a, 0.20),
      c(0.75 * vy, 0.10 * vy, ell_b, ell_c, 0.15),
      c(0.25 * vy, 0.50 * vy, ell_c, ell_b, 0.25)
    )
  } else {
    caabgp_stop("Unknown GP model: ", model)
  }
  log(starts)
}

parameter_bounds <- function(model, J) {
  if (model == "S") {
    lower <- c(1e-5, rep(0.04, J), 0.02)
    upper <- c(20.0, rep(2.00, J), 2.00)
  } else if (model == "P") {
    lower <- c(1e-5, rep(0.04, J), 0.08, 0.02)
    upper <- c(20.0, rep(2.00, J), 5.00, 2.00)
  } else if (model == "AB") {
    lower <- c(1e-5, 1e-5, rep(0.04, J), rep(0.04, J), 0.02)
    upper <- c(20.0, 20.0, rep(2.00, J), rep(2.00, J), 2.00)
  } else {
    caabgp_stop("Unknown GP model: ", model)
  }
  list(lower = log(lower), upper = log(upper))
}

fit_gp_eb <- function(X, z, y, model, K, previous_par = NULL, maxit = 60) {
  X <- as.matrix(X)
  z <- as.integer(z)
  y <- as.numeric(y)
  J <- ncol(X)
  bounds <- parameter_bounds(model, J)
  starts <- default_starts(y, model, J)
  if (!is.null(previous_par) && length(previous_par) == ncol(starts)) starts <- matrix(previous_par, nrow = 1)
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    fit_i <- tryCatch(
      stats::optim(
        par = starts[i, ],
        fn = gp_nll,
        X = X,
        z = z,
        y = y,
        model = model,
        K = K,
        method = "L-BFGS-B",
        lower = bounds$lower,
        upper = bounds$upper,
        control = list(maxit = maxit, factr = 1e7)
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_i) && (is.null(best) || fit_i$value < best$value)) best <- fit_i
  }
  if (is.null(best)) caabgp_stop("Empirical-Bayes GP optimization failed for model ", model)
  Sigma <- covariance_matrix(X, z, best$par, model, K, add_noise = TRUE)
  chol_info <- safe_chol(Sigma)
  L <- chol_info$L
  one <- rep(1, length(y))
  inv_one <- chol_solve(L, one)
  beta <- as.numeric(sum(inv_one * y) / sum(inv_one))
  alpha <- chol_solve(L, y - beta)
  structure(
    list(
      model = model,
      K = K,
      X = X,
      z = z,
      y = y,
      par = best$par,
      params = unpack_params(best$par, model, J),
      beta = beta,
      L = L,
      alpha = alpha,
      inv_one = inv_one,
      one_inv_one = sum(inv_one),
      nll = best$value,
      convergence = best$convergence,
      jitter = chol_info$jitter
    ),
    class = "caabgp_gpfit"
  )
}

predict_gp_eb <- function(fit, Xnew, znew) {
  Xnew <- as.matrix(Xnew)
  znew <- as.integer(znew)
  C <- cross_covariance(fit, Xnew, znew)
  mu <- as.numeric(fit$beta + t(C) %*% fit$alpha)
  V <- chol_solve(fit$L, C)
  base_var <- latent_variance_diag(fit, znew)
  mean_unc <- (1 - as.numeric(t(C) %*% fit$inv_one))^2 / fit$one_inv_one
  sd <- sqrt(pmax(base_var - colSums(C * V) + mean_unc, 1e-12))
  data.frame(mu = mu, sd = sd)
}

sigma_hat <- function(fit) fit$params$sigma

fit_caabgp <- function(data, design, dose_cols = design$dose_cols,
                       stratum_col = "stratum", outcome_col = "y", previous = NULL, maxit = 60) {
  validate_design(design)
  data <- as.data.frame(data)
  X <- as_dose_matrix(data, dose_cols)
  z <- as.integer(data[[stratum_col]])
  y <- as.numeric(data[[outcome_col]])
  if (any(!z %in% seq_len(design$n_strata))) caabgp_stop("stratum values must be integers from 1 to n_strata.")
  if (any(!is.finite(y))) caabgp_stop("outcomes must be finite.")
  prev <- if (!is.null(previous) && !is.null(previous$fit)) previous$fit$par else NULL
  fit <- list(method = "CA-AB-GP", fit = fit_gp_eb(X, z, y, "AB", design$n_strata, previous_par = prev, maxit = maxit))
  fit$data <- data
  fit$dose_cols <- dose_cols
  fit$stratum_col <- stratum_col
  fit$outcome_col <- outcome_col
  fit$design <- design
  class(fit) <- "caabgp_fit"
  fit
}

predict_caabgp <- function(object, newdata = object$design$dose_grid, ...) {
  if (!inherits(object, "caabgp_fit")) caabgp_stop("object must be a caabgp_fit.")
  grid <- as.data.frame(newdata)
  Xnew <- as_dose_matrix(grid, object$dose_cols)
  rows <- vector("list", object$design$n_strata)
  for (k in seq_len(object$design$n_strata)) {
    pred <- predict_gp_eb(object$fit, Xnew, rep(k, nrow(Xnew)))
    rows[[k]] <- data.frame(stratum = k, grid[, object$dose_cols, drop = FALSE], pred)
  }
  do.call(rbind, rows)
}

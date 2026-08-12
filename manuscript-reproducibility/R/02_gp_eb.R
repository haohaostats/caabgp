source(file.path("R", "00_settings.R"))

clamp_positive <- function(x, eps = 1e-10) {
  pmax(x, eps)
}

scaled_sqdist <- function(X1, X2, ell) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
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
  stop("Covariance matrix is not positive definite after jitter escalation.")
}

chol_solve <- function(L, B) {
  backsolve(L, forwardsolve(t(L), B))
}

unpack_params <- function(par, model, K) {
  if (model == "S") {
    list(
      nu = exp(par[1]),
      ell = exp(par[2:3]),
      sigma = exp(par[4])
    )
  } else if (model == "P") {
    list(
      nu = exp(par[1]),
      ell = exp(par[2:3]),
      ell_z = exp(par[4]),
      sigma = exp(par[5])
    )
  } else if (model == "AB") {
    list(
      nu0 = exp(par[1]),
      nub = exp(par[2]),
      ell0 = exp(par[3:4]),
      ellb = exp(par[5:6]),
      sigma = exp(par[7])
    )
  } else {
    stop("Unknown GP model: ", model)
  }
}

covariance_matrix <- function(X, z, par, model, K, add_noise = TRUE) {
  p <- unpack_params(par, model, K)
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
      if (length(idx) > 0) {
        Sigma[idx, idx] <- Sigma[idx, idx] + p$nub * Kb[idx, idx]
      }
    }
  }

  if (add_noise) Sigma <- Sigma + diag(p$sigma^2, n)
  Sigma
}

cross_covariance <- function(fit, Xnew, znew) {
  X <- fit$X
  z <- fit$z
  p <- fit$params
  model <- fit$model
  K <- fit$K
  Xnew <- as.matrix(Xnew)

  if (model == "S") {
    return(p$nu * se_kernel(X, Xnew, p$ell))
  }

  if (model == "P") {
    Kd <- se_kernel(X, Xnew, p$ell)
    Dz <- outer(z, znew, "!=") * 2
    Kz <- exp(-Dz / (2 * p$ell_z^2))
    return(p$nu * Kd * Kz)
  }

  if (model == "AB") {
    C <- p$nu0 * se_kernel(X, Xnew, p$ell0)
    Kb <- se_kernel(X, Xnew, p$ellb)
    for (k in seq_len(K)) {
      rows <- which(z == k)
      cols <- which(znew == k)
      if (length(rows) > 0 && length(cols) > 0) {
        C[rows, cols] <- C[rows, cols] + p$nub * Kb[rows, cols]
      }
    }
    return(C)
  }

  stop("Unknown GP model: ", model)
}

latent_variance_diag <- function(fit, znew) {
  p <- fit$params
  if (fit$model == "S") return(rep(p$nu, length(znew)))
  if (fit$model == "P") return(rep(p$nu, length(znew)))
  if (fit$model == "AB") return(rep(p$nu0 + p$nub, length(znew)))
  stop("Unknown GP model: ", fit$model)
}

profile_beta <- function(Sigma, y) {
  L <- safe_chol(Sigma)$L
  one <- rep(1, length(y))
  inv_y <- chol_solve(L, y)
  inv_one <- chol_solve(L, one)
  beta <- as.numeric(sum(inv_one * y) / sum(inv_one))
  list(beta = beta, L = L, inv_one = inv_one, one_inv_one = sum(inv_one))
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

default_starts <- function(y, model, K) {
  raw_var <- if (length(y) > 1) stats::var(y) else NA_real_
  vy <- if (is.finite(raw_var)) max(raw_var, 0.05^2) else 0.10
  if (model == "S") {
    starts <- rbind(
      log(c(0.70 * vy, 0.35, 0.35, 0.20)),
      log(c(0.50 * vy, 0.20, 0.45, 0.15)),
      log(c(1.00 * vy, 0.55, 0.25, 0.25))
    )
  } else if (model == "P") {
    starts <- rbind(
      log(c(0.70 * vy, 0.35, 0.35, 0.75, 0.20)),
      log(c(0.50 * vy, 0.20, 0.45, 0.45, 0.15)),
      log(c(1.00 * vy, 0.55, 0.25, 1.25, 0.25))
    )
  } else if (model == "AB") {
    starts <- rbind(
      log(c(0.50 * vy, 0.20 * vy, 0.35, 0.35, 0.35, 0.35, 0.20)),
      log(c(0.75 * vy, 0.10 * vy, 0.50, 0.50, 0.25, 0.25, 0.15)),
      log(c(0.25 * vy, 0.50 * vy, 0.25, 0.25, 0.50, 0.50, 0.25))
    )
  } else {
    stop("Unknown GP model: ", model)
  }
  starts
}

parameter_bounds <- function(model, K) {
  if (model == "S") {
    lower <- log(c(1e-5, 0.04, 0.04, 0.02))
    upper <- log(c(20.0, 2.00, 2.00, 2.00))
  } else if (model == "P") {
    lower <- log(c(1e-5, 0.04, 0.04, 0.08, 0.02))
    upper <- log(c(20.0, 2.00, 2.00, 5.00, 2.00))
  } else if (model == "AB") {
    lower <- log(c(1e-5, 1e-5, 0.04, 0.04, 0.04, 0.04, 0.02))
    upper <- log(c(20.0, 20.0, 2.00, 2.00, 2.00, 2.00, 2.00))
  } else {
    stop("Unknown GP model: ", model)
  }
  list(lower = lower, upper = upper)
}

fit_gp_eb <- function(X, z, y, model, K, previous_par = NULL, maxit = 40) {
  X <- as.matrix(X)
  z <- as.integer(z)
  y <- as.numeric(y)
  bounds <- parameter_bounds(model, K)
  starts <- default_starts(y, model, K)
  if (!is.null(previous_par) && length(previous_par) == ncol(starts)) {
    starts <- matrix(previous_par, nrow = 1)
  }

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

  if (is.null(best)) {
    stop("Empirical-Bayes GP optimization failed for model ", model)
  }

  Sigma <- covariance_matrix(X, z, best$par, model, K, add_noise = TRUE)
  chol_info <- safe_chol(Sigma)
  L <- chol_info$L
  one <- rep(1, length(y))
  inv_one <- chol_solve(L, one)
  beta <- as.numeric(sum(inv_one * y) / sum(inv_one))
  alpha <- chol_solve(L, y - beta)

  list(
    model = model,
    K = K,
    X = X,
    z = z,
    y = y,
    par = best$par,
    params = unpack_params(best$par, model, K),
    beta = beta,
    L = L,
    alpha = alpha,
    inv_one = inv_one,
    one_inv_one = sum(inv_one),
    nll = best$value,
    convergence = best$convergence,
    jitter = chol_info$jitter
  )
}

predict_gp_eb <- function(fit, Xnew, znew) {
  Xnew <- as.matrix(Xnew)
  znew <- as.integer(znew)
  C <- cross_covariance(fit, Xnew, znew)
  inv_C <- chol_solve(fit$L, C)
  mu <- as.numeric(fit$beta + t(C) %*% fit$alpha)

  trend_adjust <- (1 - as.numeric(t(fit$inv_one) %*% C))^2 / fit$one_inv_one
  v <- latent_variance_diag(fit, znew) - colSums(C * inv_C) + trend_adjust
  data.frame(mu = mu, sd = sqrt(clamp_positive(v)))
}

borrowing_index <- function(fit) {
  if (fit$model != "AB") return(NULL)
  rho <- fit$params$nu0 / (fit$params$nu0 + fit$params$nub)
  data.frame(stratum = seq_len(fit$K), rho = rep(rho, fit$K))
}

sigma_hat <- function(fit) {
  fit$params$sigma
}

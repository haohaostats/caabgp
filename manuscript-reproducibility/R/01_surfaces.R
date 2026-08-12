source(file.path("R", "00_settings.R"))

mvnorm_density_2d <- function(d, mu, Sigma) {
  d <- as.matrix(d)
  centered <- sweep(d, 2, mu, "-")
  invS <- solve(Sigma)
  q <- rowSums((centered %*% invS) * centered)
  denom <- 2 * pi * sqrt(det(Sigma))
  exp(-0.5 * q) / denom
}

efficacy_basin <- function(d, mu, Sigma, amplitude) {
  -amplitude * mvnorm_density_2d(d, mu, Sigma)
}

scenario_covariances <- function() {
  list(
    Sigma0 = matrix(c(0.08, 0.00, 0.00, 0.08), 2, 2, byrow = TRUE),
    SigmaM = matrix(c(0.12, 0.03, 0.03, 0.12), 2, 2, byrow = TRUE),
    SigmaH = matrix(c(0.10, 0.04, 0.04, 0.06), 2, 2, byrow = TRUE)
  )
}

surface_basis_values <- function(d) {
  d <- as.matrix(d)
  S <- scenario_covariances()
  list(
    G0 = efficacy_basin(d, c(0.60, 0.60), S$Sigma0, 1.00),
    GL = efficacy_basin(d, c(0.35, 0.70), S$SigmaM, 1.00),
    GR = efficacy_basin(d, c(0.70, 0.35), S$SigmaM, 1.00),
    GQ = efficacy_basin(d, c(0.25, 0.25), S$SigmaH, 1.10),
    G0_swap = efficacy_basin(d[, c(2, 1), drop = FALSE], c(0.60, 0.60), S$Sigma0, 1.00)
  )
}

true_surface <- function(scenario, stratum, d) {
  d <- as.matrix(d)
  B <- surface_basis_values(d)

  if (scenario == 1) {
    if (stratum == 1) return(0.60 * B$G0 + 0.40 * B$GL)
    if (stratum == 2) return(0.60 * B$G0 + 0.40 * B$GR)
  }

  if (scenario == 2) {
    if (stratum == 1) return(0.65 * B$G0 + 0.35 * B$GL)
    if (stratum == 2) return(0.65 * B$G0 + 0.35 * B$GR)
    if (stratum == 3) return(0.55 * B$G0 + 0.45 * B$GQ)
  }

  if (scenario == 3) {
    if (stratum == 1) return(0.70 * B$G0 + 0.30 * B$G0_swap)
    if (stratum == 2) return(0.60 * B$G0 + 0.40 * B$GL)
    if (stratum == 3) return(0.60 * B$G0 + 0.40 * B$GR)
    if (stratum == 4) return(0.50 * B$G0 + 0.50 * B$GQ)
  }

  if (scenario == 4) {
    if (stratum == 1) return(0.55 * B$G0 + 0.45 * B$GL)
    if (stratum == 2) return(0.55 * B$G0 + 0.45 * B$GR)
  }

  stop("Unknown scenario/stratum combination: scenario=", scenario, ", stratum=", stratum)
}

true_grid_values <- function(scenario, grid = candidate_grid()) {
  K <- length(scenario_prevalence(scenario))
  out <- matrix(NA_real_, nrow = K, ncol = nrow(grid))
  for (k in seq_len(K)) out[k, ] <- true_surface(scenario, k, grid)
  out
}

true_optima <- function(scenario, grid = candidate_grid()) {
  vals <- true_grid_values(scenario, grid)
  idx <- apply(vals, 1, which.min)
  data.frame(
    scenario = scenario,
    stratum = seq_len(nrow(vals)),
    d1_star = grid[idx, 1],
    d2_star = grid[idx, 2],
    f_star = vals[cbind(seq_len(nrow(vals)), idx)]
  )
}

write_true_surface_csv <- function() {
  grid <- candidate_grid()
  rows <- list()
  opt_rows <- list()
  i <- 1
  j <- 1
  for (s in 1:4) {
    K <- length(scenario_prevalence(s))
    opt_rows[[j]] <- true_optima(s, grid)
    j <- j + 1
    for (k in seq_len(K)) {
      rows[[i]] <- data.frame(
        scenario = s,
        stratum = k,
        d1 = grid[, 1],
        d2 = grid[, 2],
        f = true_surface(s, k, grid)
      )
      i <- i + 1
    }
  }
  write.csv(do.call(rbind, rows), file.path(DIR_RESULTS, "figure1_true_surfaces_grid.csv"), row.names = FALSE)
  write.csv(do.call(rbind, opt_rows), file.path(DIR_RESULTS, "figure1_true_surface_optima.csv"), row.names = FALSE)
}

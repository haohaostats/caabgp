source(file.path("R", "00_settings.R"))
source(file.path("R", "02_gp_eb.R"))
source(file.path("R", "03_acquisition.R"))

APP_SETTINGS <- list(
  K = 2,
  prevalence = c(0.65, 0.35),
  weights = c(0.55, 0.45),
  omega = 0.5,
  n_initial = 20,
  n_max = 80,
  b_max = 350,
  cohort_size = 2,
  c_pt = 1,
  c_new = 15,
  c_scr = 1,
  lambda_c = 0.5,
  stopping_threshold = 0,
  stopping_patience = 3L,
  seed = 20260521,
  methods = c("P-GP", "AB-GP", "CA-AB-GP"),
  calibration_study = "Gandhi2014_1",
  toxicity_outcome = 1L,
  efficacy_outcome = 25L,
  smoothing_bandwidth = 0.30,
  deviation_scale_sd = 2.60,
  deviation_scale_range = 0.75,
  opposite_penalty_fraction = 0.18,
  interaction_scale = 0.08,
  primary_bump_scale = 0.17,
  opposite_bump_scale = 0.22
)

xml_unescape <- function(x) {
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  x <- gsub("&lt;", "<", x, fixed = TRUE)
  x <- gsub("&gt;", ">", x, fixed = TRUE)
  x <- gsub("&quot;", "\"", x, fixed = TRUE)
  x <- gsub("&apos;", "'", x, fixed = TRUE)
  x
}

extract_matches <- function(x, pattern) {
  m <- gregexpr(pattern, x, perl = TRUE)[[1]]
  if (m[1] < 0) return(character(0))
  regmatches(x, list(m))[[1]]
}

extract_first <- function(x, pattern) {
  out <- extract_matches(x, pattern)
  if (length(out) == 0) return(NA_character_)
  sub(pattern, "\\1", out[1], perl = TRUE)
}

extract_attr <- function(x, attr) {
  pattern <- paste0(".*\\s", attr, "=\"([^\"]*)\".*")
  if (!grepl(pattern, x, perl = TRUE)) return(NA_character_)
  sub(pattern, "\\1", x, perl = TRUE)
}

excel_col_index <- function(cell_ref) {
  letters <- gsub("[0-9]", "", cell_ref)
  chars <- strsplit(letters, "", fixed = TRUE)[[1]]
  out <- 0L
  for (ch in chars) out <- out * 26L + match(ch, LETTERS)
  out
}

read_xlsx_base <- function(path, sheet_name) {
  tmp <- tempfile("xlsx_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  utils::unzip(path, exdir = tmp)

  workbook <- paste(readLines(file.path(tmp, "xl", "workbook.xml"), warn = FALSE), collapse = "")
  rels <- paste(readLines(file.path(tmp, "xl", "_rels", "workbook.xml.rels"), warn = FALSE), collapse = "")

  sheet_blocks <- extract_matches(workbook, "<sheet[^>]*/>")
  sheet_names <- vapply(sheet_blocks, extract_attr, character(1), attr = "name")
  sheet_ids <- vapply(sheet_blocks, extract_attr, character(1), attr = "r:id")
  idx <- match(sheet_name, sheet_names)
  if (is.na(idx)) stop("Sheet not found: ", sheet_name)

  rel_blocks <- extract_matches(rels, "<Relationship[^>]*/>")
  rel_id <- vapply(rel_blocks, extract_attr, character(1), attr = "Id")
  rel_target <- vapply(rel_blocks, extract_attr, character(1), attr = "Target")
  target <- rel_target[match(sheet_ids[idx], rel_id)]
  sheet_path <- file.path(tmp, "xl", target)

  shared <- character(0)
  ss_path <- file.path(tmp, "xl", "sharedStrings.xml")
  if (file.exists(ss_path)) {
    ss <- paste(readLines(ss_path, warn = FALSE), collapse = "")
    si <- extract_matches(ss, "<si>.*?</si>")
    shared <- vapply(si, function(block) {
      txt <- extract_matches(block, "<t[^>]*>.*?</t>")
      txt <- gsub("<t[^>]*>|</t>", "", txt)
      xml_unescape(paste(txt, collapse = ""))
    }, character(1))
  }

  sx <- paste(readLines(sheet_path, warn = FALSE), collapse = "")
  row_blocks <- extract_matches(sx, "<row[^>]*>.*?</row>")
  rows <- list()
  max_col <- 0L
  for (i in seq_along(row_blocks)) {
    cells <- extract_matches(row_blocks[i], "<c[^>]*>.*?</c>")
    if (length(cells) == 0) next
    vals <- list()
    for (cell in cells) {
      cell_start <- sub("^(<c[^>]*>).*$", "\\1", cell, perl = TRUE)
      ref <- extract_attr(cell_start, "r")
      type <- extract_attr(cell_start, "t")
      col <- excel_col_index(ref)
      max_col <- max(max_col, col)
      val <- extract_first(cell, ".*<v>(.*?)</v>.*")
      if (is.na(val) && grepl("<is>", cell, fixed = TRUE)) {
        txt <- extract_matches(cell, "<t[^>]*>.*?</t>")
        txt <- gsub("<t[^>]*>|</t>", "", txt)
        val <- paste(txt, collapse = "")
      }
      if (is.na(val)) val <- ""
      if (!is.na(type) && type == "s" && nzchar(val)) {
        val <- shared[as.integer(val) + 1L]
      }
      vals[[as.character(col)]] <- xml_unescape(val)
    }
    row_vec <- rep("", max_col)
    for (nm in names(vals)) row_vec[as.integer(nm)] <- vals[[nm]]
    rows[[length(rows) + 1L]] <- row_vec
  }

  max_len <- max(vapply(rows, length, integer(1)))
  mat <- do.call(rbind, lapply(rows, function(x) {
    length(x) <- max_len
    x[is.na(x)] <- ""
    x
  }))
  header <- make.names(mat[1, ], unique = TRUE)
  df <- as.data.frame(mat[-1, , drop = FALSE], stringsAsFactors = FALSE)
  names(df) <- header
  df
}

as_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.+-]", "", x)))

parse_dose_components <- function(x) {
  nums <- regmatches(x, gregexpr("[0-9]+\\.?[0-9]*", x, perl = TRUE))[[1]]
  nums <- suppressWarnings(as.numeric(nums))
  nums <- nums[is.finite(nums)]
  if (length(nums) >= 2L) return(nums[1:2])
  if (length(nums) == 1L) return(c(nums[1], NA_real_))
  c(NA_real_, NA_real_)
}

standardize_range <- function(x) {
  rx <- range(x, na.rm = TRUE)
  if (!all(is.finite(rx)) || diff(rx) <= 0) return(rep(0.5, length(x)))
  (x - rx[1]) / diff(rx)
}

kernel_smooth_2d <- function(grid, x_obs, y_obs, bandwidth = 0.30) {
  out <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    dist2 <- rowSums((sweep(x_obs, 2, grid[i, ], "-"))^2)
    w <- exp(-dist2 / (2 * bandwidth^2))
    if (sum(w) <= 1e-12) {
      out[i] <- y_obs[which.min(dist2)]
    } else {
      out[i] <- sum(w * y_obs) / sum(w)
    }
  }
  out
}

load_brock_calibration <- function(path = file.path(ROOT_DIR, "case_application", "data", "Brock_Database.xlsx")) {
  if (!file.exists(path)) stop("Missing Brock database: ", path)
  outcomes <- read_xlsx_base(path, "Outcomes")
  events <- read_xlsx_base(path, "BinaryOutcomeEvents")
  series <- read_xlsx_base(path, "BinaryOutcomeAnalysisSeries")

  names(outcomes) <- tolower(names(outcomes))
  names(events) <- tolower(names(events))
  names(series) <- tolower(names(series))

  events$studyid <- trimws(events$study)
  events$outcomeid <- as.integer(as_num(events$outcomeid))
  events$n <- as_num(events$n)
  events$events <- as_num(events$events)
  events$dose <- trimws(events$dose)

  outcomes$outcomeid <- as.integer(as_num(outcomes$outcomeid))
  outcomes$outcome_text <- paste(outcomes$outcometext, outcomes$outcomeclass, outcomes$note, sep = " ")
  outcomes$outcome_text <- tolower(outcomes$outcome_text)
  outcomes$is_toxicity <- grepl("dlt|tox|safety|adverse", outcomes$outcome_text)
  outcomes$is_efficacy <- grepl("response|efficacy|recist|disease control|activity", outcomes$outcome_text)

  ev <- merge(
    events,
    outcomes[, c("outcomeid", "outcome_text", "is_toxicity", "is_efficacy")],
    by = "outcomeid",
    all.x = TRUE
  )
  ev <- ev[nzchar(ev$studyid) & is.finite(ev$n) & ev$n > 0 & nzchar(ev$dose), ]

  evs <- ev[ev$studyid == APP_SETTINGS$calibration_study, ]
  tox <- evs[evs$outcomeid == APP_SETTINGS$toxicity_outcome,
             c("dose", "n", "events", "outcome_text")]
  eff <- evs[evs$outcomeid == APP_SETTINGS$efficacy_outcome,
             c("dose", "n", "events", "outcome_text")]
  names(tox) <- c("dose", "n_t", "x_t", "tox_text")
  names(eff) <- c("dose", "n_e", "x_e", "eff_text")
  dat <- merge(tox, eff, by = "dose")
  dat <- dat[is.finite(dat$n_t) & is.finite(dat$n_e), ]
  if (nrow(dat) < 4L) {
    stop("The prespecified Brock calibration series is missing or incomplete. Check APP_SETTINGS and Brock_Database.xlsx.")
  }
  best <- list(
    studyid = APP_SETTINGS$calibration_study,
    toxicity_outcome = APP_SETTINGS$toxicity_outcome,
    efficacy_outcome = APP_SETTINGS$efficacy_outcome,
    n_levels = nrow(dat),
    total_n = sum(dat$n_t + dat$n_e, na.rm = TRUE),
    data = dat
  )
  dat <- best$data

  dose_components <- t(vapply(dat$dose, parse_dose_components, numeric(2)))
  dat$dose1_raw <- dose_components[, 1]
  dat$dose2_raw <- dose_components[, 2]
  dat$has_two_dimensional_dose <- is.finite(dat$dose1_raw) & is.finite(dat$dose2_raw)
  if (all(dat$has_two_dimensional_dose)) {
    dat$d1_std <- standardize_range(dat$dose1_raw)
    dat$d2_std <- standardize_range(dat$dose2_raw)
    dat <- dat[order(dat$dose1_raw, dat$dose2_raw), ]
  } else {
    dat$dose_number <- as_num(dat$dose)
    dat <- dat[order(dat$dose), ]
  }
  row.names(dat) <- NULL

  dat$p_e <- (dat$x_e + 0.5) / (dat$n_e + 1)
  dat$p_t <- (dat$x_t + 0.5) / (dat$n_t + 1)
  dat$utility <- dat$p_e - APP_SETTINGS$omega * dat$p_t
  dat$y <- -dat$utility

  if (all(dat$has_two_dimensional_dose)) {
    obs_dose <- as.matrix(dat[, c("d1_std", "d2_std")])
    y_grid <- kernel_smooth_2d(candidate_grid(), obs_dose, dat$y,
                               bandwidth = APP_SETTINGS$smoothing_bandwidth)
    grid_t <- NULL
  } else {
    dat$t <- if (nrow(dat) == 1) 0 else seq(0, 1, length.out = nrow(dat))
    grid_t <- seq(0, 1, by = 0.25)
    y_grid <- stats::approx(dat$t, dat$y, xout = grid_t, rule = 2)$y
    obs_dose <- NULL
  }

  scale_y <- stats::sd(y_grid)
  if (!is.finite(scale_y) || scale_y < 0.03) scale_y <- 0.08

  list(
    metadata = data.frame(
      studyid = best$studyid,
      toxicity_outcome = best$toxicity_outcome,
      efficacy_outcome = best$efficacy_outcome,
      n_levels = best$n_levels,
      total_n = best$total_n,
      omega = APP_SETTINGS$omega
    ),
    dose_level = dat,
    obs_dose = obs_dose,
    grid_t = grid_t,
    y_grid = y_grid,
    scale_y = scale_y
  )
}

gaussian_bump <- function(grid, center, scale = 0.12) {
  exp(-rowSums((sweep(grid, 2, center, "-"))^2) / (2 * scale^2))
}

calibrated_app_surface <- function(calib, grid = candidate_grid()) {
  if (!is.null(calib$obs_dose)) {
    common <- kernel_smooth_2d(grid, calib$obs_dose, calib$dose_level$y,
                               bandwidth = APP_SETTINGS$smoothing_bandwidth)
  } else {
    h <- pmin(1, pmax(0, 0.55 * grid[, 1] + 0.45 * grid[, 2]))
    common <- stats::approx(calib$grid_t, calib$y_grid, xout = h, rule = 2)$y
  }
  s <- calib$scale_y
  observed_range <- diff(range(common, finite = TRUE))
  deviation_amplitude <- max(APP_SETTINGS$deviation_scale_sd * s,
                             APP_SETTINGS$deviation_scale_range * observed_range)
  opposite_penalty <- APP_SETTINGS$opposite_penalty_fraction * deviation_amplitude
  interaction <- APP_SETTINGS$interaction_scale * s * (grid[, 1] - 0.5) * (grid[, 2] - 0.5)
  y1 <- common + interaction -
    deviation_amplitude * gaussian_bump(grid, c(0.00, 1.00), APP_SETTINGS$primary_bump_scale) +
    opposite_penalty * gaussian_bump(grid, c(1.00, 0.00), APP_SETTINGS$opposite_bump_scale)
  y2 <- common - interaction -
    deviation_amplitude * gaussian_bump(grid, c(1.00, 0.00), APP_SETTINGS$primary_bump_scale) +
    opposite_penalty * gaussian_bump(grid, c(0.00, 1.00), APP_SETTINGS$opposite_bump_scale)
  data.frame(
    d1 = rep(grid[, 1], APP_SETTINGS$K),
    d2 = rep(grid[, 2], APP_SETTINGS$K),
    stratum = rep(seq_len(APP_SETTINGS$K), each = nrow(grid)),
    y_mean = c(y1, y2)
  )
}

app_surface_value <- function(surface, k, d) {
  idx <- which(surface$stratum == k & abs(surface$d1 - d[1]) < 1e-12 & abs(surface$d2 - d[2]) < 1e-12)
  surface$y_mean[idx[1]]
}

app_initial_counts <- function() {
  weights <- APP_SETTINGS$weights
  min_count <- 2L
  counts <- rep(min_count, length(weights))
  remaining <- APP_SETTINGS$n_initial - sum(counts)
  add <- floor(remaining * weights / sum(weights))
  counts <- counts + add
  while (sum(counts) < APP_SETTINGS$n_initial) {
    target <- APP_SETTINGS$n_initial * weights / sum(weights)
    j <- which.max(target - counts)
    counts[j] <- counts[j] + 1L
  }
  as.integer(counts)
}

make_application_initial <- function(surface, residual_sd) {
  set.seed(APP_SETTINGS$seed)
  init_doses <- SIM_SETTINGS$initial_design
  counts <- app_initial_counts()
  rows <- list()
  id <- 1L
  for (k in seq_along(counts)) {
    dose_order <- rep(seq_len(nrow(init_doses)), length.out = counts[k])
    for (j in dose_order) {
      d <- init_doses[j, ]
      mu <- app_surface_value(surface, k, d)
      rows[[id]] <- data.frame(d1 = d[1], d2 = d[2], stratum = k, y = stats::rnorm(1, mu, residual_sd[k]))
      id <- id + 1L
    }
  }
  do.call(rbind, rows)
}

application_cost <- function(data) {
  manufactured <- unique(dose_key(data$d1, data$d2))
  screening <- sum(vapply(data$stratum, function(k) 1 / max(APP_SETTINGS$prevalence[k], SIM_SETTINGS$pi_min), numeric(1)))
  APP_SETTINGS$c_pt * nrow(data) +
    APP_SETTINGS$c_new * length(manufactured) +
    APP_SETTINGS$c_scr * screening
}

fit_application_design <- function(method, data, previous = NULL) {
  X <- as.matrix(data[, c("d1", "d2")])
  z <- as.integer(data$stratum)
  y <- data$y
  if (method == "P-GP") {
    prev <- if (!is.null(previous)) previous$fit$par else NULL
    return(list(method = method, fit = fit_gp_eb(X, z, y, "P", APP_SETTINGS$K, previous_par = prev)))
  }
  if (method %in% c("AB-GP", "CA-AB-GP")) {
    prev <- if (!is.null(previous)) previous$fit$par else NULL
    return(list(method = method, fit = fit_gp_eb(X, z, y, "AB", APP_SETTINGS$K, previous_par = prev)))
  }
  stop("Unsupported application method: ", method)
}

predict_application_design <- function(design_fit, grid = candidate_grid()) {
  rows <- list()
  for (k in seq_len(APP_SETTINGS$K)) {
    pred <- predict_gp_eb(design_fit$fit, grid, rep(k, nrow(grid)))
    rows[[k]] <- data.frame(stratum = k, d1 = grid[, 1], d2 = grid[, 2], mu = pred$mu, sd = pred$sd)
  }
  do.call(rbind, rows)
}

application_recommendations <- function(pred, method, n_current, unique_doses, total_cost, repeat_fraction, rho = NA_real_) {
  rows <- lapply(seq_len(APP_SETTINGS$K), function(k) {
    pk <- pred[pred$stratum == k, ]
    idx <- which.min(pk$mu)
    data.frame(
      method = method,
      stratum = k,
      d1_hat = pk$d1[idx],
      d2_hat = pk$d2[idx],
      mu_hat = pk$mu[idx],
      sd_hat = pk$sd[idx],
      n = n_current,
      unique_doses = unique_doses,
      repeat_fraction = repeat_fraction,
      total_cost = total_cost,
      rho = rho
    )
  })
  do.call(rbind, rows)
}

run_application_method <- function(method, initial_data, surface, residual_sd) {
  grid <- candidate_grid()
  data <- initial_data
  total_cost <- application_cost(data)
  manufactured <- unique(dose_key(data$d1, data$d2))
  previous <- NULL
  allocations <- list()
  active_strata <- seq_len(APP_SETTINGS$K)
  stopping_counters <- integer(APP_SETTINGS$K)
  completed_updates <- 0L
  step <- 1L

  repeat {
    fit <- fit_application_design(method, data, previous)
    previous <- fit
    pred <- predict_application_design(fit, grid)

    sigma_k <- rep(sigma_hat(fit$fit), APP_SETTINGS$K)
    if (completed_updates > 0L && length(active_strata)) {
      acquisition_for_stopping <- personalized_acquisition_table(
        pred = pred,
        grid = grid,
        manufactured_keys = manufactured,
        prevalence = APP_SETTINGS$prevalence,
        sigma_by_k = sigma_k,
        c_new = APP_SETTINGS$c_new,
        c_scr = APP_SETTINGS$c_scr,
        lambda_c = if (method == "CA-AB-GP") APP_SETTINGS$lambda_c else 0,
        weights = APP_SETTINGS$weights,
        cohort_size = APP_SETTINGS$cohort_size,
        remaining_budget = Inf,
        active_strata = active_strata
      )
      for (k in active_strata) {
        scores_k <- acquisition_for_stopping$score[acquisition_for_stopping$stratum == k]
        max_score_k <- if (length(scores_k)) max(scores_k, na.rm = TRUE) else -Inf
        if (is.finite(max_score_k) && max_score_k < APP_SETTINGS$stopping_threshold) {
          stopping_counters[k] <- stopping_counters[k] + 1L
        } else {
          stopping_counters[k] <- 0L
        }
      }
      active_strata <- active_strata[
        stopping_counters[active_strata] < APP_SETTINGS$stopping_patience
      ]
    }

    if (!length(active_strata)) break
    if (nrow(data) + APP_SETTINGS$cohort_size > APP_SETTINGS$n_max) break

    remaining <- APP_SETTINGS$b_max - total_cost
    if (remaining <= 0) break
    cohort <- APP_SETTINGS$cohort_size
    lambda <- if (method == "CA-AB-GP") APP_SETTINGS$lambda_c else 0
    action <- select_personalized_action(
      pred = pred,
      grid = grid,
      manufactured_keys = manufactured,
      prevalence = APP_SETTINGS$prevalence,
      sigma_by_k = sigma_k,
      c_new = APP_SETTINGS$c_new,
      c_scr = APP_SETTINGS$c_scr,
      lambda_c = lambda,
      weights = APP_SETTINGS$weights,
      cohort_size = cohort,
      remaining_budget = remaining,
      active_strata = active_strata
    )
    if (is.null(action)) break

    d <- c(action$d1, action$d2)
    k <- action$stratum
    new_dose <- !is_manufactured(d, manufactured)
    cost_inc <- APP_SETTINGS$c_pt * cohort +
      APP_SETTINGS$c_scr * cohort / max(APP_SETTINGS$prevalence[k], SIM_SETTINGS$pi_min) +
      APP_SETTINGS$c_new * as.numeric(new_dose)
    if (cost_inc > remaining + 1e-8) break

    y_new <- stats::rnorm(cohort, app_surface_value(surface, k, d), residual_sd[k])
    data <- rbind(data, data.frame(d1 = d[1], d2 = d[2], stratum = k, y = y_new))
    total_cost <- total_cost + cost_inc
    manufactured <- unique(c(manufactured, dose_key(d[1], d[2])))
    allocations[[step]] <- data.frame(
      method = method,
      step = step,
      stratum = k,
      d1 = d[1],
      d2 = d[2],
      cohort_size = cohort,
      new_dose = as.numeric(new_dose),
      cost_increment = cost_inc,
      total_cost = total_cost,
      n = nrow(data),
      unique_doses = length(manufactured)
    )
    step <- step + 1L
    completed_updates <- completed_updates + 1L
  }

  fit <- fit_application_design(method, data, previous)
  pred <- predict_application_design(fit, grid)
  alloc <- if (length(allocations)) do.call(rbind, allocations) else data.frame()
  post_initial_alloc <- if (nrow(alloc)) sum(alloc$cohort_size) else 0
  repeat_alloc <- if (nrow(alloc) && post_initial_alloc > 0) {
    sum(alloc$cohort_size[alloc$new_dose == 0]) / post_initial_alloc
  } else {
    NA_real_
  }
  rho <- if (method %in% c("AB-GP", "CA-AB-GP")) borrowing_index(fit$fit)$rho[1] else NA_real_
  rec <- application_recommendations(pred, method, nrow(data), length(manufactured), total_cost, repeat_alloc, rho)
  list(data = data, fit = fit, pred = pred, allocations = alloc, recommendations = rec)
}

run_realdata_application <- function() {
  calib <- load_brock_calibration()
  surface <- calibrated_app_surface(calib)
  residual_sd <- rep(max(0.04, calib$scale_y * 0.85), APP_SETTINGS$K)
  initial <- make_application_initial(surface, residual_sd)
  outputs <- lapply(APP_SETTINGS$methods, run_application_method, initial_data = initial, surface = surface, residual_sd = residual_sd)
  names(outputs) <- APP_SETTINGS$methods

  all_recs <- do.call(rbind, lapply(outputs, `[[`, "recommendations"))
  all_alloc <- do.call(rbind, lapply(outputs, `[[`, "allocations"))
  all_pred <- do.call(rbind, lapply(names(outputs), function(m) {
    cbind(method = m, outputs[[m]]$pred)
  }))

  write.csv(calib$metadata, file.path(DIR_RESULTS, "application_calibration_metadata.csv"), row.names = FALSE)
  write.csv(calib$dose_level, file.path(DIR_RESULTS, "application_calibration_dose_levels.csv"), row.names = FALSE)
  write.csv(surface, file.path(DIR_RESULTS, "application_calibrated_surface.csv"), row.names = FALSE)
  write.csv(initial, file.path(DIR_RESULTS, "application_initial_data.csv"), row.names = FALSE)
  write.csv(all_recs, file.path(DIR_RESULTS, "application_results.csv"), row.names = FALSE)
  write.csv(all_alloc, file.path(DIR_RESULTS, "application_allocations.csv"), row.names = FALSE)
  write.csv(all_pred, file.path(DIR_RESULTS, "application_predictions.csv"), row.names = FALSE)

  list(calibration = calib, surface = surface, initial = initial, outputs = outputs,
       recommendations = all_recs, allocations = all_alloc, predictions = all_pred)
}

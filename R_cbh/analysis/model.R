# Exploratory joint age-band model. Source after R_cbh/load_pipeline.R.
cbh_trial_spec <- function(hiv = c("incidence", "prevalence")) {
  hiv <- match.arg(hiv)
  spec <- list(id = "age_band_complete_case_v1", weighting = "unweighted",
    covariates = c("sex", "multiple_birth", "birth_order", "maternal_age_birth",
      "maternal_education_years", "wealth_quintile", "urban", "log_hiv_prev",
      "log_gdp_pc", "log_health_expenditure_pc", "political_stability"),
    scaled = c("birth_order", "maternal_age_birth", "maternal_education_years",
      "log_hiv_prev", "log_gdp_pc", "log_health_expenditure_pc", "political_stability"),
    pfpr_k = 5L, time_k = 6L, seed = 20260907L, nthreads = 2L)
  if (hiv == "incidence") {
    spec$id <- "age_band_hiv_incidence_v2"
    spec$covariates[spec$covariates == "log_hiv_prev"] <- "log_hiv_incidence"
    spec$scaled[spec$scaled == "log_hiv_prev"] <- "log_hiv_incidence"
    spec$incidence_panel <- "data/derived_cbh/hiv_incidence/child_incidence_country_year.csv"
    spec$incidence_draws <- "data/derived_cbh/hiv_incidence/child_incidence_draws.rds"
  }
  spec
}

cbh_attach_incidence <- function(d, panel) {
  cbh_require(d, c("country", "entry_year"), "Child-band incidence join")
  cbh_require(panel, c("iso3", "year", "hiv_incidence_per1000", "hiv_incidence_status"), "Child incidence panel")
  cbh_unique(panel, c("iso3", "year"), "Child incidence panel")
  j <- match(paste(d$country, d$entry_year), paste(panel$iso3, panel$year))
  v <- panel$hiv_incidence_per1000[j]
  valid <- is.finite(v) & v > 0
  d$hiv_incidence_per1000 <- ifelse(valid, v, NA_real_)
  d$hiv_incidence_status <- ifelse(is.na(j), "missing_country_year", panel$hiv_incidence_status[j])
  d$hiv_incidence_status[!is.na(j) & !valid] <- "missing_or_invalid_value"
  d$log_hiv_incidence <- log(d$hiv_incidence_per1000)
  d
}

cbh_trial_prepare <- function(d, age_levels, spec) {
  required <- c("death", "age_band", "pfpr_pct", "calendar_year", "band_years",
                "survey", "country", "region", spec$covariates)
  cbh_require(d, required, "Complete-case model input")
  if (!nrow(d) || any(!complete.cases(d[required]))) stop("Incomplete or empty model input.")
  for (v in required) if (is.numeric(d[[v]]) && any(!is.finite(d[[v]]))) stop("Nonfinite model input.")
  d$age_band <- factor(d$age_band, levels = age_levels, ordered = FALSE)
  if (anyNA(d$age_band) || any(table(d$age_band) == 0)) stop("All seven age bands must be represented.")
  for (v in c("survey", "country", "region")) d[[v]] <- factor(d[[v]])
  d$sex <- factor(d$sex, levels = c("female", "male"))
  d$wealth_quintile <- factor(d$wealth_quintile, levels = 1:5, ordered = FALSE)
  if (anyNA(d$sex) || anyNA(d$wealth_quintile)) stop("Unexpected categorical covariate value.")
  d$country_age <- interaction(d$country, d$age_band, drop = TRUE)
  stopifnot(all(d$death %in% 0:1), all(d$band_years > 0))
  scaling <- data.frame(variable = spec$scaled,
    mean = vapply(d[spec$scaled], mean, numeric(1)),
    sd = vapply(d[spec$scaled], stats::sd, numeric(1)), row.names = NULL)
  if (any(!is.finite(scaling$sd) | scaling$sd <= 0)) stop("A continuous covariate has no variation.")
  for (i in seq_len(nrow(scaling))) {
    v <- scaling$variable[i]
    d[[paste0("z_", v)]] <- (d[[v]] - scaling$mean[i]) / scaling$sd[i]
  }
  list(data = d, scaling = scaling)
}

cbh_trial_data <- function(input_dir, spec) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Installed data.table is required.")
  meta <- readRDS(file.path(input_dir, "manifest.rds"))
  if (!isTRUE(meta$complete)) stop("The dataset build is incomplete.")
  m <- meta$manifest
  if (any(m$status == "failed")) stop("Resolve failed dataset builds before fitting.")
  skipped <- m[!m$status %in% c("built", "cached"), c("survey", "status")]
  m <- m[m$status %in% c("built", "cached"), ]
  required <- c("death", "age_band", "pfpr_pct", "calendar_year", "band_years",
                "survey", "country", "region", spec$covariates)
  retained <- unique(c(required, "child_id", "mother_id", "psu", "stratum",
    "stratum_variable", "survey_weight", "entry_year", "band_entry_cmc", "band_end_cmc",
    "interview_cmc", "age_band_index", "hiv_prev_pct", "death_band_b6_b7_disagree"))
  incidence <- NULL
  if (!is.null(spec$incidence_panel)) {
    incidence <- cbh_read_csv(spec$incidence_panel)
    if (!is.null(spec$imputation_draw)) {
      draws <- readRDS(spec$incidence_draws)
      stopifnot(length(unique(incidence$imputation_signature)) == 1L,
                identical(unique(incidence$imputation_signature), draws$signature))
      j <- match(paste(incidence$iso3, incidence$year), paste(draws$keys$iso3, draws$keys$year))
      stopifnot(!anyNA(j), spec$imputation_draw <= nrow(draws$log_incidence))
      incidence$hiv_incidence_per1000 <- exp(draws$log_incidence[spec$imputation_draw, j])
      rm(draws)
    }
    retained <- unique(c(retained, "hiv_incidence_per1000", "hiv_incidence_status"))
  }
  rows <- selection <- missing <- list()
  for (i in seq_len(nrow(m))) {
    object <- readRDS(file.path(input_dir, m$file[i]))
    if (!identical(object$signature, m$signature[i])) stop("Dataset shard differs from manifest.")
    d <- object$data
    if (!is.null(incidence)) d <- cbh_attach_incidence(d, incidence)
    cbh_require(d, retained, "Dataset shard")
    core <- d$model_ready
    keep <- core & complete.cases(d[required])
    for (v in required) if (is.numeric(d[[v]])) keep <- keep & is.finite(d[[v]])
    selection[[i]] <- data.frame(survey = m$survey[i], country = m$country[i],
      eligible_rows = nrow(d), pfpr_available_rows = sum(core),
      excluded_incomplete = sum(core & !keep), complete_case_rows = sum(keep),
      complete_case_deaths = sum(d$death[keep]))
    missing[[i]] <- data.frame(survey = m$survey[i], country = m$country[i], variable = spec$covariates,
      pfpr_available_rows = sum(core), missing_rows = vapply(spec$covariates, function(v) {
        bad <- is.na(d[[v]])
        if (is.numeric(d[[v]])) bad <- bad | !is.finite(d[[v]])
        sum(core & bad)
      }, integer(1)))
    d <- d[keep, retained, drop = FALSE]
    d$analysis_weight <- if (nrow(d)) d$survey_weight / mean(d$survey_weight) else numeric()
    rows[[i]] <- d
    rm(object, d)
    if (i %% 20L == 0L) message("Complete-case selection: ", i, "/", nrow(m), " surveys")
  }
  d <- as.data.frame(data.table::rbindlist(rows, use.names = TRUE))
  rm(rows); invisible(gc(FALSE))
  prepared <- cbh_trial_prepare(d, meta$config$age_bands$age_band, spec)
  list(data = prepared$data, scaling = prepared$scaling, selection = cbh_bind(selection),
       missing = cbh_bind(missing), skipped = skipped, input_manifest = meta)
}

cbh_trial_formula <- function(spec) {
  x <- ifelse(spec$covariates %in% spec$scaled, paste0("z_", spec$covariates), spec$covariates)
  stats::as.formula(paste0("death ~ 0 + age_band + ",
    "s(pfpr_pct, by=age_band, bs='cr', k=", spec$pfpr_k, ") + ",
    "s(calendar_year, by=age_band, bs='cr', k=", spec$time_k, ") + ",
    "age_band:(", paste(x, collapse = " + "), ") + ",
    "s(survey, bs='re') + s(country_age, bs='re') + s(region, bs='re') + offset(log(band_years))"))
}

cbh_trial_fit <- function(d, spec, trace = TRUE, start = NULL, smoothing = NULL) {
  if (!identical(spec$weighting, "unweighted")) stop("This trial implements the unweighted conditional likelihood.")
  warnings <- character()
  set.seed(spec$seed)
  timing <- system.time({
    fit <- withCallingHandlers(mgcv::bam(cbh_trial_formula(spec), data = d,
      family = stats::binomial(link = "cloglog"), method = "fREML", discrete = TRUE,
      nthreads = spec$nthreads, gc.level = 1, na.action = stats::na.fail,
      coef = start, sp = smoothing,
      control = mgcv::gam.control(trace = trace, maxit = 100)),
      warning = function(w) {
        warnings <<- unique(c(warnings, conditionMessage(w)))
        invokeRestart("muffleWarning")
      })
  })
  if (length(fit$fitted.values) != nrow(d)) stop("Model fitting changed the number of records.")
  list(fit = fit, elapsed_seconds = unname(timing["elapsed"]), warnings = warnings)
}

cbh_trial_contrasts <- function(fit, d, from = 40, to = 20) {
  age_levels <- levels(d$age_band)
  index <- match(age_levels, d$age_band)
  a <- b <- d[index, , drop = FALSE]
  a$pfpr_pct <- from; b$pfpr_pct <- to
  # Random effects and fixed confounders cancel in the within-band contrast.
  la <- stats::predict(fit, newdata = a, type = "lpmatrix")
  lb <- stats::predict(fit, newdata = b, type = "lpmatrix")
  L <- lb - la
  delta <- drop(L %*% stats::coef(fit))
  variance <- rowSums((L %*% fit$Vp) * L)
  if (any(!is.finite(variance) | variance < -1e-8)) stop("Invalid contrast variance.")
  se <- sqrt(pmax(0, variance))
  check <- drop(stats::predict(fit, b, type = "link") - stats::predict(fit, a, type = "link"))
  stopifnot(max(abs(delta - check)) < 1e-7)
  support <- t(vapply(age_levels, function(age) {
    stats::quantile(d$pfpr_pct[d$age_band == age], c(.025, .975), names = FALSE)
  }, numeric(2)))
  data.frame(age_band = age_levels, pfpr_from_pct = from, pfpr_to_pct = to,
    hazard_ratio = exp(delta), lower_95 = exp(delta - 1.96 * se), upper_95 = exp(delta + 1.96 * se),
    percent_change = 100 * expm1(delta), pfpr_p025 = support[, 1], pfpr_p975 = support[, 2],
    contrast_within_central_95pct = pmin(from, to) >= support[, 1] & pmax(from, to) <= support[, 2])
}

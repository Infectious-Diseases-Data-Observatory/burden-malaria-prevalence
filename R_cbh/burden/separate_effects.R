# Extract only within-age PfPR effects. Cross-age covariance is not assumed zero
# or estimated here, and these objects cannot produce country-total intervals.
cbh_separate_burden_components <- function(source_dir, cache_dir, ages) {
  paths <- file.path(source_dir, "sensitivity_single_imputation", paste0("age_", seq_along(ages), ".rds"))
  stopifnot(all(file.exists(paths)))
  signature <- cbh_hash(list(vapply(paths, cbh_file_hash, ""),
    cbh_file_hash("R_cbh/burden/separate_effects.R"), as.character(packageVersion("mgcv"))))
  path <- file.path(cache_dir, "burden_pfpr_components.rds")
  cache <- if (file.exists(path)) readRDS(path) else NULL
  if (is.null(cache) || !identical(cache$signature, signature)) {
    pieces <- lapply(seq_along(ages), function(i) {
      message("Extracting separate-age PfPR curve: ", ages[i])
      obj <- readRDS(paths[i]); f <- obj$fit
      stopifnot(isTRUE(f$converged), all(is.finite(coef(f))), all(is.finite(f$Vp)),
        max(abs(exp(f$model[["offset(log(band_years))"]]) - c(1, 5, 6, 12, 12, 12, 12)[i]/12)) < 1e-12,
        identical(f$family$link, "cloglog"))
      ss <- Filter(function(s) identical(s$term, "pfpr_pct"), f$smooth)
      stopifnot(length(ss) == 1L, identical(ss[[1]]$by, "NA"))
      ix <- ss[[1]]$first.para:ss[[1]]$last.para
      pp <- f$model$pfpr_pct
      # Check the compact representation against the full fitted prediction matrix
      # at zero as well as nonzero exposures, including the burden extrapolation.
      nd <- f$model[rep(1L, 5L), , drop = FALSE]
      nd$band_years <- exp(nd[["offset(log(band_years))"]])
      nd$pfpr_pct <- c(0, 20, 40, 60, max(pp)); ref <- nd; ref$pfpr_pct <- 20
      L <- predict(f, nd, type = "lpmatrix", discrete = FALSE) -
        predict(f, ref, type = "lpmatrix", discrete = FALSE)
      C <- mgcv::PredictMat(ss[[1]], nd) - mgcv::PredictMat(ss[[1]], ref)
      stopifnot(max(abs(drop(L %*% coef(f)) - drop(C %*% coef(f)[ix]))) < 1e-8,
        max(abs(rowSums((L %*% f$Vp) * L) - rowSums((C %*% f$Vp[ix, ix]) * C))) < 1e-8)
      list(age_band = ages[i], coef = coef(f)[ix], covariance = f$Vp[ix, ix], smooth = ss[[1]],
        input_signature = obj$input_signature, fit_signature = obj$signature,
        support = data.frame(age_band = ages[i], minimum = min(pp), maximum = max(pp),
          p025 = unname(quantile(pp, .025)), p975 = unname(quantile(pp, .975))))
    })
    stopifnot(length(unique(vapply(pieces, function(x) x$input_signature, ""))) == 1L)
    cache <- list(signature = signature, pieces = pieces,
      imputation = "posterior_median", cross_age_covariance_available = FALSE)
    cbh_atomic_rds(cache, path)
  }
  stopifnot(identical(vapply(cache$pieces, function(x) x$age_band, ""), ages))
  list(cache = cache, paths = paths)
}

cbh_separate_burden_contrasts <- function(components, valid, saved_contrasts) {
  do.call(rbind, lapply(components$cache$pieces, function(piece) {
    d <- valid[valid$age_band == piece$age_band, ]; ref <- d; ref$pfpr_pct <- 0
    L <- mgcv::PredictMat(piece$smooth, ref) - mgcv::PredictMat(piece$smooth, d)
    variance <- rowSums((L %*% piece$covariance) * L)
    stopifnot(all(is.finite(variance)), min(variance) >= -1e-9)
    # Also match the published single-imputation sensitivity contrast.
    C <- mgcv::PredictMat(piece$smooth, data.frame(pfpr_pct = 20)) -
      mgcv::PredictMat(piece$smooth, data.frame(pfpr_pct = 40))
    previous <- saved_contrasts[saved_contrasts$model == "separate" & saved_contrasts$age_band == piece$age_band, ]
    stopifnot(nrow(previous) == 1L,
      abs(exp(drop(C %*% piece$coef)) - previous$hazard_ratio_40_to_20) < 1e-8)
    data.frame(iso3 = d$iso3, age_band = d$age_band, draw = NA_integer_,
      log_hr_zero_vs_current = drop(L %*% piece$coef), variance = pmax(variance, 0))
  }))
}

# Shared helpers for single-imputation subgroup sensitivity analyses.
cbh_sensitivity_formula <- function(spec, single_age = FALSE) {
  if (!single_age) return(cbh_trial_formula(spec))
  x <- ifelse(spec$covariates %in% spec$scaled,
              paste0("z_", spec$covariates), spec$covariates)
  stats::as.formula(paste0(
    "death ~ s(pfpr_pct, bs='cr', k=", spec$pfpr_k,
    ") + s(calendar_year, bs='cr', k=", spec$time_k, ") + ",
    paste(x, collapse = " + "),
    " + s(survey, bs='re') + s(country, bs='re') + s(region, bs='re')",
    " + offset(log(band_years))"))
}

cbh_sensitivity_curves <- function(fit, d, id) {
  do.call(rbind, lapply(levels(d$age_band), function(age) {
    use <- which(d$age_band == age)
    support <- quantile(d$pfpr_pct[use], c(0, .025, .975, 1), names = FALSE)
    # A common grid permits direct comparison; unsupported predictions are flagged.
    grid <- sort(unique(c(seq(0, 100, by = .5), support)))
    nd <- d[rep(use[1], length(grid)), , drop = FALSE]
    ref <- nd
    nd$pfpr_pct <- grid; ref$pfpr_pct <- 20
    L <- predict(fit, nd, type = "lpmatrix", discrete = FALSE) -
      predict(fit, ref, type = "lpmatrix", discrete = FALSE)
    est <- drop(L %*% coef(fit))
    variance <- rowSums((L %*% fit$Vp) * L)
    direct <- predict(fit, nd, type = "link", discrete = FALSE) -
      predict(fit, ref, type = "link", discrete = FALSE)
    stopifnot(all(is.finite(variance)), min(variance) > -1e-8,
      max(abs(est - direct)) < 1e-7, abs(est[grid == 20]) < 1e-10,
      abs(variance[grid == 20]) < 1e-10)
    se <- sqrt(pmax(variance, 0))
    data.frame(fit_id = id, age_band = age, pfpr_pct = grid,
      log_hazard_ratio = est, standard_error = se,
      lower_95 = est - 1.96 * se, upper_95 = est + 1.96 * se,
      pfpr_min = support[1], pfpr_p025 = support[2],
      pfpr_p975 = support[3], pfpr_max = support[4],
      within_observed_support = grid >= support[1] & grid <= support[4],
      within_central_support = grid >= support[2] & grid <= support[3],
      reference_within_observed_support = 20 >= support[1] & 20 <= support[4])
  }))
}

cbh_sensitivity_diagnostics <- function(saved, d, id) {
  f <- saved$fit
  h <- f$outer.info$hess
  data.frame(fit_id = id, rows = nrow(d), deaths = sum(d$death),
    countries = nlevels(d$country), surveys = nlevels(d$survey),
    regions = nlevels(d$region), converged = isTRUE(f$converged),
    iterations = f$iter, rank = f$rank, coefficients = length(coef(f)),
    finite_coefficients = all(is.finite(coef(f))), finite_covariance = all(is.finite(f$Vp)),
    elapsed_seconds = saved$elapsed_seconds,
    min_smoothing_hessian_eigenvalue = if (length(h)) min(eigen(h, symmetric = TRUE, only.values = TRUE)$values) else NA_real_,
    max_smoothing_gradient = if (length(f$outer.info$grad)) max(abs(f$outer.info$grad)) else NA_real_,
    warnings = paste(saved$warnings, collapse = " | "))
}

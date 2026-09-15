#!/usr/bin/env Rscript
# Read-only audit of saved primary inputs/fits; exports aggregate checks only.
# No fitting, data rebuilding, source retrieval or modification of prior results.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/sensitivity/model.R")
library(mgcv)
library(data.table)
out <- "results/cbh/code_audit_2026_09_15"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
root <- "results/cbh/map_snow_gamma2_v1"
manifest <- cbh_read_csv(file.path(root, "fit_manifest.csv"))
manifest <- manifest[manifest$series == "map_full", ]
ages <- cbh_config()$age_bands$age_band
manifest <- manifest[match(ages, manifest$age_band), ]
stopifnot(nrow(manifest) == 7L, !anyNA(manifest), !anyDuplicated(manifest$fit_id))
data_path <- "data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/complete_case_dataset.rds"
prepared <- readRDS(data_path)
input <- prepared$data
stopifnot(nrow(input) == 5885022L, sum(input$death) == 82415L,
  !anyNA(input$child_id), data.table::uniqueN(input, by = c("child_id", "age_band")) == nrow(input))
spec <- cbh_trial_spec()
expected <- cbh_sensitivity_formula(spec, single_age = TRUE)
rows <- list()
for (i in seq_along(ages)) {
  m <- manifest[i, ]
  stopifnot(identical(cbh_file_hash(m$model_file), m$md5),
    identical(cbh_file_hash(m$source_path), m$source_md5))
  saved <- readRDS(m$model_file)
  fit <- saved$fit
  dd <- input[input$age_band == ages[i], , drop = FALSE]
  stopifnot(identical(saved$input_signature, prepared$signature), saved$gamma == 2,
    identical(fit$family$family, "binomial"), identical(fit$family$link, "cloglog"),
    identical(all.vars(formula(fit)), all.vars(expected)),
    identical(attr(terms(formula(fit)), "term.labels"), attr(terms(expected), "term.labels")),
    isTRUE(fit$converged), fit$rank == length(coef(fit)),
    all(is.finite(coef(fit))), all(is.finite(fit$Vp)),
    all(fit$prior.weights == 1), nrow(dd) == nrow(fit$model),
    inherits(fit$smooth[[1]], "cr.smooth"), inherits(fit$smooth[[2]], "cr.smooth"),
    length(saved$knots$pfpr_pct) == 5L, length(saved$knots$calendar_year) == 6L,
    identical(saved$knots$pfpr_pct, fit$smooth[[1]]$xp),
    identical(saved$knots$calendar_year, fit$smooth[[2]]$xp),
    all(vapply(fit$smooth[3:5], inherits, logical(1), "random.effect")),
    identical(vapply(fit$smooth[3:5], function(s) s$term, ""), c("survey", "country", "region")))
  # Compare every outcome, predictor and offset with the saved modelling data,
  # not just row/death counts. Factor labels are compared rather than integer codes.
  for (v in names(fit$model)) {
    want <- if (v == "offset(log(band_years))") log(dd$band_years) else dd[[v]]
    got <- fit$model[[v]]
    if (is.factor(got)) got <- as.character(got)
    if (is.factor(want)) want <- as.character(want)
    stopifnot(!is.null(want), isTRUE(all.equal(got, want, check.attributes = FALSE)))
  }
  nd <- dd[rep(1L, 5L), ]; nd$pfpr_pct <- c(0, 10, 20, 40, 50)
  ref <- nd; ref$pfpr_pct <- 20
  L <- predict(fit, nd, type = "lpmatrix", discrete = FALSE) -
    predict(fit, ref, type = "lpmatrix", discrete = FALSE)
  s <- fit$smooth[[1]]; ix <- s$first.para:s$last.para
  C <- PredictMat(s, nd) - PredictMat(s, ref)
  stopifnot(max(abs(drop(L %*% coef(fit)) - drop(C %*% coef(fit)[ix]))) < 1e-8,
    max(abs(rowSums((L %*% fit$Vp) * L) - rowSums((C %*% fit$Vp[ix, ix]) * C))) < 1e-8)
  hmin <- min(eigen(fit$outer.info$hess, symmetric = TRUE, only.values = TRUE)$values)
  stopifnot(is.finite(hmin), hmin > 0)
  rows[[i]] <- data.frame(age_band = ages[i], records = nrow(dd), deaths = sum(dd$death),
    distinct_children = uniqueN(dd$child_id), gamma = saved$gamma,
    input_columns_verified = ncol(fit$model), model_md5_verified = TRUE,
    reference_md5_verified = TRUE, input_signature_verified = TRUE,
    formula_verified = TRUE, unweighted = TRUE, knots_verified = TRUE,
    compact_contrasts_verified = TRUE, min_hessian_eigenvalue = hmin,
    max_smoothing_gradient = max(abs(fit$outer.info$grad)),
    sum_fitted_death_probabilities = sum(fit$fitted.values))
  message("Verified primary age ", ages[i], ": ", nrow(dd), " records")
  rm(saved, fit, dd); invisible(gc(FALSE))
}
cbh_atomic_csv(do.call(rbind, rows), file.path(out, "primary_fit_checks.csv"))
cbh_atomic_csv(data.frame(distinct_children = uniqueN(input$child_id), records = nrow(input),
  deaths = sum(input$death), surveys = uniqueN(input$survey), countries = uniqueN(input$country)),
  file.path(out, "primary_sample.csv"))
rm(prepared, input); invisible(gc(FALSE))

# Verify the saved primary national calculations, without a new raster extraction.
burden <- cbh_read_csv(file.path(root, "burden/country_age_estimates.csv"))
burden <- burden[burden$series == "map_full", ]
totals <- cbh_read_csv(file.path(root, "burden/country_totals.csv"))
totals <- totals[totals$series == "map_full", ]
stopifnot(nrow(burden) == 3L * 45L * 7L, nrow(totals) == 3L * 45L,
  setequal(unique(burden$year), c(2005, 2015, 2024)),
  all(burden$population_weight_year == 2020), all(burden$exposure_source == "MAP"))
for (unit in c("deaths", "rate_per100000")) {
  baseline <- burden[[paste0("ihme_", unit)]]
  counter <- burden[[paste0("counterfactual_", unit)]]
  attr <- burden[[paste0("attributable_", unit)]]
  stopifnot(max(abs(baseline - counter - attr), na.rm = TRUE) < 1e-7,
    max(abs(attr - baseline * (1 - burden$hr_zero_vs_current)), na.rm = TRUE) < 1e-7)
}
annual <- lapply(sort(unique(totals$year)), function(y) {
  z <- totals[totals$year == y, ]
  for (iso in z$iso3) {
    b <- burden[burden$year == y & burden$iso3 == iso, ]
    old <- b[b$age_band %in% ages[5:7], ]
    stopifnot(nrow(b) == 7L, diff(range(old$ihme_deaths)) < 1e-7,
      diff(range(old$ihme_rate_per100000)) < 1e-7,
      isTRUE(all.equal(sum(b$attributable_deaths), z$attributable_under5_deaths[z$iso3 == iso])))
  }
  data.frame(year = y, country_rows = nrow(z), estimable = sum(is.finite(z$attributable_under5_deaths)),
    attributable_deaths = sum(z$attributable_under5_deaths, na.rm = TRUE),
    accounting_verified = TRUE, equal_ages_2_to_4_allocation_verified = TRUE)
})
cbh_atomic_csv(do.call(rbind, annual), file.path(out, "burden_checks.csv"))
inputs <- c(data_path, manifest$model_file, manifest$source_path,
  file.path(root, c("fit_manifest.csv", "burden/country_age_estimates.csv", "burden/country_totals.csv")),
  "R_cbh/audit/01_verify_saved_primary.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")),
  file.path(out, "verification_provenance.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(out, "session_info.txt"))
message("Saved-primary audit passed; aggregate checks in ", out)

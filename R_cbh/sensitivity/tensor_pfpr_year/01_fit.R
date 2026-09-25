#!/usr/bin/env Rscript
# Sensitivity (added 25 September 2026, at the user's request): the primary age-band models
# without the survey-region random intercept, with PfPR[2-10] and calendar time modelled jointly
# by a tensor product instead of two separate splines. Fitted to the collapsed survey-region x
# entry-month cells (binomial deaths out of children entering the band), which gives the same
# likelihood as the child-level Bernoulli model (R_cbh/sensitivity/bam_vs_gam). Per band:
#   separate: s(pfpr_pct, cr, k=5) + s(calendar_year, cr, k=6)             (comparator)
#   tensor  : te(pfpr_pct, calendar_year, cr x cr, k = c(5, 6))
#   ti      : s(pfpr_pct) + s(calendar_year) + ti(pfpr_pct, calendar_year)  (interaction test)
# all with the 17 covariates, survey and country random intercepts, offset(log(band_years)),
# the primary reference knots, gamma = 2, bam fREML without discretisation. A fit whose smoothing
# Hessian is not positive definite takes the primary's strict restart (kept only if it converges, is
# positive definite, lowers the smoothing gradient and does not worsen fREML); a fit that still fails
# stops the run. Cached by a signature that includes this script, so any change to the fit settings refits.
source("R_cbh/load_pipeline.R")
source("R_cbh/sensitivity/model.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
suppressPackageStartupMessages({library(data.table); library(mgcv)})
args <- commandArgs(trailingOnly = TRUE); stopifnot(all(args %in% "--force"))
base <- cbh_primary_settings("regional_mics")
id <- "tensor_pfpr_year_dhsmics_map_gamma2_v1"; out <- file.path("results/cbh", id); private <- file.path("data/derived_cbh/models", id)
for (p in c(out, private)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ages <- cbh_config()$age_bands$age_band

covs <- paste0("z_", cbh_primary_regional_spec(base)$covariates)
rest <- paste(c(covs, "s(survey, bs = \"re\")", "s(country, bs = \"re\")", "offset(log(band_years))"), collapse = " + ")
forms <- list(
  separate = as.formula(paste("cbind(deaths, survivors) ~ s(pfpr_pct, bs = \"cr\", k = 5) + s(calendar_year, bs = \"cr\", k = 6) +", rest)),
  tensor = as.formula(paste("cbind(deaths, survivors) ~ te(pfpr_pct, calendar_year, bs = c(\"cr\", \"cr\"), k = c(5, 6)) +", rest)),
  ti = as.formula(paste("cbind(deaths, survivors) ~ s(pfpr_pct, bs = \"cr\", k = 5) + s(calendar_year, bs = \"cr\", k = 6) + ti(pfpr_pct, calendar_year, bs = c(\"cr\", \"cr\"), k = c(5, 6)) +", rest)))
writeLines(unlist(lapply(names(forms), function(m) c(paste0("## ", m), deparse(forms[[m]]), ""))), file.path(out, "model_formulas.txt"))
vars <- unique(c("region", "calendar_year", "pfpr_pct", "band_years", "survey", "country", covs))

data_hash <- cbh_file_hash(base$data)
stopifnot(identical(unique(fread(file.path(base$out, "fit_manifest.csv"))$prepared_data_md5), data_hash))
dat <- as.data.table(readRDS(base$data)$data)[, c("death", "age_band", vars), with = FALSE]
stopifnot(nrow(dat) == base$expected_records, sum(dat$death) == base$expected_deaths)
knot_table <- fread(base$knots)
code <- c("R_cbh/sensitivity/tensor_pfpr_year/01_fit.R", base$knots, "R_cbh/primary/specification.R", "R_cbh/primary/settings.R", "R_cbh/sensitivity/model.R")
script_hash <- cbh_file_hash("R_cbh/sensitivity/tensor_pfpr_year/01_fit.R")
unlink(file.path(out, "restarts.csv"))
cbh_atomic_csv(data.frame(file = c(base$data, code), md5 = c(data_hash, vapply(code, cbh_file_hash, ""))), file.path(out, "fit_input_provenance.csv"))

fit_bam <- function(form, cells, knots, strict = FALSE, start = NULL) {
  warnings <- character(); set.seed(20260907L)
  control <- if (strict) gam.control(epsilon = 1e-9, mgcv.tol = 1e-9, efs.tol = .001, maxit = 200) else gam.control(maxit = 100)
  timing <- system.time(f <- withCallingHandlers(bam(form, data = cells, knots = knots, family = binomial(link = "cloglog"),
    method = "fREML", discrete = FALSE, gamma = 2, select = FALSE, nthreads = 2, na.action = na.fail, coef = start, control = control),
    warning = function(w) { warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleWarning") }))
  list(fit = f, elapsed_seconds = unname(timing["elapsed"]), warnings = warnings)
}
diags <- restarts <- summaries <- list(); components <- list(); cellsum <- list()
for (a in ages) {
  i <- match(a, ages); d <- dat[age_band == a]
  cells <- d[, c(lapply(.SD, function(x) x[1]), .(deaths = sum(death), children = .N)), by = .(region, calendar_year), .SDcols = setdiff(vars, c("region", "calendar_year"))]
  chk <- d[, lapply(.SD, uniqueN), by = .(region, calendar_year), .SDcols = setdiff(vars, c("region", "calendar_year"))]
  stopifnot(all(as.matrix(chk[, -(1:2)]) == 1L), sum(cells$children) == nrow(d), sum(cells$deaths) == sum(d$death))
  cells[, survivors := children - deaths]; for (v in c("survey", "country")) cells[[v]] <- factor(cells[[v]])
  # PfPR support by period (children-weighted), to show where a time interaction is identified.
  cells[, period := cut(calendar_year, c(2000, 2005, 2010, 2015, 2020, 2025), right = FALSE, labels = c("2000-04", "2005-09", "2010-14", "2015-19", "2020-23"))]
  wq <- function(x, w, p) { o <- order(x); x[o][findInterval(p * sum(w), cumsum(w[o])) + 1] }
  cellsum[[a]] <- cells[, .(children = sum(children), deaths = sum(deaths), pfpr_p025 = wq(pfpr_pct, children, .025), pfpr_median = wq(pfpr_pct, children, .5),
    pfpr_p975 = wq(pfpr_pct, children, .975)), by = period][order(period)][, age_band := a]
  kt <- knot_table[age_band == a]
  knots <- lapply(c(pfpr_pct = "pfpr_pct", calendar_year = "calendar_year"), function(v) kt[variable == v][order(index)]$value)
  for (m in names(forms)) {
    fit_id <- sprintf("%s_age_%d", m, i); path <- file.path(private, paste0(fit_id, ".rds"))
    sig <- cbh_hash(list(fit_id, knots, deparse(forms[[m]]), data_hash, script_hash, R.version.string, as.character(packageVersion("mgcv"))))
    saved <- if (file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
    if (is.null(saved) || !identical(saved$signature, sig)) {
      message(a, " ", m, ": fitting on ", nrow(cells), " cells")
      saved <- fit_bam(forms[[m]], cells, knots); saved$signature <- sig; saved$version <- "original"
      diag <- cbh_sensitivity_diagnostics(saved, cells, fit_id)
      if (!diag$converged || !is.finite(diag$min_smoothing_hessian_eigenvalue) || diag$min_smoothing_hessian_eigenvalue <= 0) {
        message(a, " ", m, ": strict restart")
        strict <- fit_bam(forms[[m]], cells, knots, strict = TRUE, start = coef(saved$fit))
        sd <- cbh_sensitivity_diagnostics(strict, cells, fit_id)
        take <- sd$converged && sd$finite_covariance && sd$rank == sd$coefficients && is.finite(sd$min_smoothing_hessian_eigenvalue) &&
          sd$min_smoothing_hessian_eigenvalue > 0 && sd$max_smoothing_gradient < diag$max_smoothing_gradient && strict$fit$gcv.ubre <= saved$fit$gcv.ubre
        restarts[[fit_id]] <- data.frame(fit_id, original_min_hessian = diag$min_smoothing_hessian_eigenvalue, strict_min_hessian = sd$min_smoothing_hessian_eigenvalue,
          original_fREML = unname(saved$fit$gcv.ubre), strict_fREML = unname(strict$fit$gcv.ubre), strict_selected = take)
        if (take) { strict$signature <- sig; strict$version <- "strict_restart"; saved <- strict }
      }
      cbh_atomic_rds(saved, path)
    }
    f <- saved$fit; diag <- cbh_sensitivity_diagnostics(saved, cells, fit_id)
    stopifnot(diag$converged, diag$finite_covariance, diag$rank == diag$coefficients, is.finite(diag$min_smoothing_hessian_eigenvalue), diag$min_smoothing_hessian_eigenvalue > 0)
    # Kernel log-likelihood (excludes the binomial coefficient); AIC() adds it and uses the corrected df (edf2).
    ll <- { eta <- predict(f, cells, type = "link"); p <- -expm1(-exp(eta)); sum(cells$deaths * log(p) + cells$survivors * log1p(-p)) }
    diags[[fit_id]] <- cbind(diag[, c("fit_id", "converged", "rank", "coefficients", "min_smoothing_hessian_eigenvalue", "max_smoothing_gradient", "elapsed_seconds")],
      data.frame(model = m, age_band = a, cells = nrow(cells), version = saved$version, loglik = ll, aic = AIC(f), aic_df = attr(logLik(f), "df"), edf_total = sum(f$edf),
        deviance_explained = summary(f)$dev.expl, fREML = unname(f$gcv.ubre)))
    st <- as.data.frame(summary(f)$s.table); st$term <- rownames(st); st$model <- m; st$age_band <- a; summaries[[fit_id]] <- st
    keep <- which(vapply(f$smooth, function(s) any(c("pfpr_pct", "calendar_year") %in% s$term), NA))
    components[[fit_id]] <- list(model = m, age_band = a, smooths = f$smooth[keep],
      index = lapply(f$smooth[keep], function(s) s$first.para:s$last.para), coef = coef(f), Vp = f$Vp)
    cbh_atomic_csv(do.call(rbind, diags), file.path(out, "fit_diagnostics.csv"))
    cbh_atomic_csv(do.call(rbind, summaries), file.path(out, "smooth_summaries.csv"))
    if (length(restarts)) cbh_atomic_csv(do.call(rbind, restarts), file.path(out, "restarts.csv"))
    rm(f, saved); gc(FALSE)
  }
  rm(d, cells, chk); gc(FALSE)
}
# Compact components (smooth objects, coefficients and covariance) for the report.
components <- lapply(components, function(z) { ix <- sort(unique(unlist(z$index))); z$coef_all_names <- names(z$coef); z$Vp <- z$Vp[ix, ix, drop = FALSE]; z$coef <- z$coef[ix]; z$ix <- ix; z })
cbh_atomic_rds(components, file.path(private, "pfpr_year_components.rds"))
cbh_atomic_csv(as.data.frame(rbindlist(cellsum)), file.path(out, "pfpr_support_by_period.csv"))
message("Tensor-product fits complete: ", out)

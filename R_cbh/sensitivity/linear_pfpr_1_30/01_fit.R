#!/usr/bin/env Rscript
# Sensitivity (added 25 September 2026, at the user's request): PfPR[2-10] forced to be linear on the
# 1-30% range. Records are restricted to band entries with PfPR 1-30% (inclusive) and each age band is
# fitted twice on that subset:
#   linear : pfpr_pct as a linear term (log hazard linear in PfPR)
#   spline : s(pfpr_pct, cr, k = 5), knots equally spaced on 1-30% (comparator on the same records)
# Otherwise both are the v7 primary: s(calendar_year, cr, k = 6) on the primary reference knots, the 17
# covariates (scaled over the full v7 sample), survey, country and survey-region random intercepts,
# offset(log(band_years)), gamma = 2, bam fREML. Fitted to the collapsed survey-region x entry-month
# cells (binomial deaths out of children entering the band), which gives the same likelihood as the
# child-level Bernoulli model (R_cbh/sensitivity/bam_vs_gam), with covariate discretisation as in the v7
# primary (discrete = TRUE; the non-discretised fit with the region random intercept took 7-18 minutes per
# model, and discretisation changed hazard ratios by <0.0005 in bam_vs_gam). A fit whose
# smoothing Hessian is not positive definite takes the strict restart used in tensor_pfpr_year/01_fit.R;
# a fit that still fails stops the run. Cached by a signature that includes this script.
source("R_cbh/load_pipeline.R")
source("R_cbh/sensitivity/model.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
suppressPackageStartupMessages({library(data.table); library(mgcv)})
args <- commandArgs(trailingOnly = TRUE); stopifnot(all(args %in% "--force"))
base <- cbh_primary_settings("regional_mics")
id <- "linear_pfpr_1_30_dhsmics_map_gamma2_v1"; out <- file.path("results/cbh", id); private <- file.path("data/derived_cbh/models", id)
for (p in c(out, private)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ages <- cbh_config()$age_bands$age_band
range_pct <- c(1, 30)

rhs <- gsub("\\s+", " ", deparse1(cbh_primary_regional_formula(base)[[3]]))
spl <- "s(pfpr_pct, bs = \"cr\", k = 5)"
stopifnot(grepl(spl, rhs, fixed = TRUE), grepl("s(region, bs = \"re\")", rhs, fixed = TRUE))
forms <- list(linear = as.formula(paste("cbind(deaths, survivors) ~", sub(spl, "pfpr_pct", rhs, fixed = TRUE))),
              spline = as.formula(paste("cbind(deaths, survivors) ~", rhs)))
writeLines(unlist(lapply(names(forms), function(m) c(paste0("## ", m), deparse(forms[[m]]), ""))), file.path(out, "model_formulas.txt"))
covs <- paste0("z_", cbh_primary_regional_spec(base)$covariates)
vars <- unique(c("region", "calendar_year", "pfpr_pct", "band_years", "survey", "country", covs))

data_hash <- cbh_file_hash(base$data)
stopifnot(identical(unique(fread(file.path(base$out, "fit_manifest.csv"))$prepared_data_md5), data_hash))
dat <- as.data.table(readRDS(base$data)$data)[, c("death", "age_band", vars), with = FALSE]
stopifnot(nrow(dat) == base$expected_records, sum(dat$death) == base$expected_deaths)
dat[, in_range := pfpr_pct >= range_pct[1] & pfpr_pct <= range_pct[2]]
selection <- rbind(dat[, .(subset = "v7 primary", records = .N, deaths = sum(death), surveys = uniqueN(survey), countries = uniqueN(country), regions = uniqueN(region)), by = age_band],
                   dat[in_range == TRUE, .(subset = "PfPR 1-30%", records = .N, deaths = sum(death), surveys = uniqueN(survey), countries = uniqueN(country), regions = uniqueN(region)), by = age_band])
cbh_atomic_csv(as.data.frame(selection), file.path(out, "selection_by_age.csv"))
dat <- dat[in_range == TRUE][, in_range := NULL]
knot_table <- fread(base$knots)
script_hash <- cbh_file_hash("R_cbh/sensitivity/linear_pfpr_1_30/01_fit.R")
code <- c("R_cbh/sensitivity/linear_pfpr_1_30/01_fit.R", base$knots, "R_cbh/primary/specification.R", "R_cbh/primary/settings.R", "R_cbh/sensitivity/model.R")
unlink(file.path(out, "restarts.csv"))
cbh_atomic_csv(data.frame(file = c(base$data, code), md5 = c(data_hash, vapply(code, cbh_file_hash, ""))), file.path(out, "fit_input_provenance.csv"))

fit_bam <- function(form, cells, knots, strict = FALSE, start = NULL) {
  warnings <- character(); set.seed(20260907L)
  control <- if (strict) gam.control(epsilon = 1e-9, mgcv.tol = 1e-9, efs.tol = .001, maxit = 200) else gam.control(maxit = 100)
  timing <- system.time(f <- withCallingHandlers(bam(form, data = cells, knots = knots, family = binomial(link = "cloglog"),
    method = "fREML", discrete = TRUE, gamma = 2, select = FALSE, nthreads = 2, na.action = na.fail, coef = start, control = control),
    warning = function(w) { warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleWarning") }))
  list(fit = f, elapsed_seconds = unname(timing["elapsed"]), warnings = warnings)
}
diags <- restarts <- summaries <- list(); components <- list()
for (a in ages) {
  i <- match(a, ages); d <- dat[age_band == a]
  cells <- d[, c(lapply(.SD, function(x) x[1]), .(deaths = sum(death), children = .N)), by = .(region, calendar_year), .SDcols = setdiff(vars, c("region", "calendar_year"))]
  chk <- d[, lapply(.SD, uniqueN), by = .(region, calendar_year), .SDcols = setdiff(vars, c("region", "calendar_year"))]
  stopifnot(all(as.matrix(chk[, -(1:2)]) == 1L), sum(cells$children) == nrow(d), sum(cells$deaths) == sum(d$death))
  cells[, survivors := children - deaths]; for (v in c("survey", "country", "region")) cells[[v]] <- factor(cells[[v]])
  kt <- knot_table[age_band == a]
  cal_knots <- kt[variable == "calendar_year"][order(index)]$value
  stopifnot(min(cells$calendar_year) >= min(cal_knots), max(cells$calendar_year) <= max(cal_knots))
  knot_sets <- list(linear = list(calendar_year = cal_knots), spline = list(pfpr_pct = seq(range_pct[1], range_pct[2], length.out = 5), calendar_year = cal_knots))
  for (m in names(forms)) {
    fit_id <- sprintf("%s_age_%d", m, i); path <- file.path(private, paste0(fit_id, ".rds")); knots <- knot_sets[[m]]
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
    diag$deaths <- sum(cells$deaths)
    diags[[fit_id]] <- cbind(diag[, c("fit_id", "rows", "deaths", "surveys", "countries", "regions", "converged", "rank", "coefficients",
      "min_smoothing_hessian_eigenvalue", "max_smoothing_gradient", "elapsed_seconds")],
      data.frame(model = m, age_band = a, children = sum(cells$children), version = saved$version, loglik = ll, aic = AIC(f), aic_df = attr(logLik(f), "df"),
        edf_total = sum(f$edf), deviance_explained = summary(f)$dev.expl, fREML = unname(f$gcv.ubre)))
    st <- as.data.frame(summary(f)$s.table); st$term <- rownames(st); st$model <- m; st$age_band <- a; summaries[[fit_id]] <- st
    if (m == "linear") {
      j <- match("pfpr_pct", names(coef(f)))
      components[[fit_id]] <- list(model = m, age_band = a, smooth = NULL, coef = coef(f)[j], covariance = f$Vp[j, j, drop = FALSE],
        p_value = summary(f)$p.table["pfpr_pct", "Pr(>|z|)"])
    } else {
      s <- f$smooth[[1]]; stopifnot(identical(s$term, "pfpr_pct")); ix <- s$first.para:s$last.para
      components[[fit_id]] <- list(model = m, age_band = a, smooth = s, coef = coef(f)[ix], covariance = f$Vp[ix, ix])
    }
    components[[fit_id]]$support <- quantile(rep(cells$pfpr_pct, cells$children), c(0, .025, .5, .975, 1), names = FALSE)
    cbh_atomic_csv(do.call(rbind, diags), file.path(out, "fit_diagnostics.csv"))
    cbh_atomic_csv(do.call(rbind, summaries), file.path(out, "smooth_summaries.csv"))
    if (length(restarts)) cbh_atomic_csv(do.call(rbind, restarts), file.path(out, "restarts.csv"))
    rm(f, saved); gc(FALSE)
  }
  rm(d, cells, chk); gc(FALSE)
}
cbh_atomic_rds(components, file.path(private, "pfpr_components.rds"))
message("Linear-PfPR (1-30%) fits complete: ", out)

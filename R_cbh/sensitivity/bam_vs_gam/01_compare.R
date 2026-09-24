#!/usr/bin/env Rscript
# Does fitting with bam (fREML, discretised covariates) instead of gam (REML) change the primary
# fit? (added 24 September 2026). The exposure and every covariate are constant within a survey-
# region x entry-month cell, and band width is constant within a band, so the child-level
# Bernoulli likelihood equals a binomial likelihood on cells (deaths out of children entering):
# gam can then fit exactly the same model. For each requested band this script fits
#   gam  : gam(method = "REML") on cells
#   bam_c: bam(method = "fREML", discrete = FALSE) on cells (fREML without discretisation)
# and compares both with the saved v7 primary fit (bam, fREML, discrete = TRUE, child-level).
# Usage: Rscript 01_compare.R [--no-region] [age bands, default all]. --no-region compares the model
# without the survey-region random intercept (R_cbh/sensitivity/no_region_re) instead of v7, which
# gam fits far faster (about 230 instead of 1,450 coefficients). Aggregates only are written.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
suppressPackageStartupMessages({library(data.table); library(mgcv)})
base <- cbh_primary_settings("regional_mics")
id <- "bam_vs_gam_dhsmics_map_gamma2_v1"; out <- file.path("results/cbh", id); private <- file.path("data/derived_cbh/models", id)
for (p in c(out, private)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ages <- cbh_config()$age_bands$age_band
args <- commandArgs(trailingOnly = TRUE); no_region <- "--no-region" %in% args; args <- setdiff(args, "--no-region")
run_ages <- if (length(args)) args else ages; stopifnot(all(run_ages %in% ages))
suffix <- if (no_region) "_no_region" else ""

form_child <- cbh_primary_regional_formula(base)
if (no_region) form_child <- as.formula(gsub(" + s(region, bs = \"re\")", "", gsub("\\s+", " ", deparse1(form_child)), fixed = TRUE))
rhs <- sub("^death ~ ", "", gsub("\\s+", " ", deparse1(form_child)))
form_cell <- as.formula(paste("cbind(deaths, survivors) ~", rhs))
vars <- union(setdiff(all.vars(form_child), "death"), "region")
manifest <- if (no_region) fread("results/cbh/no_region_re_dhsmics_map_gamma2_v1/fit_manifest.csv") else fread(file.path(base$out, "fit_manifest.csv"))[series == "map_full"]
knot_table <- fread(base$knots)
dat <- as.data.table(readRDS(base$data)$data)[, c("death", vars, "age_band"), with = FALSE]
stopifnot(nrow(dat) == base$expected_records)

contrasts <- function(f) {
  s <- f$smooth[[1]]; ix <- s$first.para:s$last.para
  L <- PredictMat(s, data.frame(pfpr_pct = c(20, 0))) - PredictMat(s, data.frame(pfpr_pct = c(40, 20)))
  est <- drop(L %*% coef(f)[ix]); se <- sqrt(rowSums((L %*% f$Vp[ix, ix]) * L))
  data.table(contrast = c("40% to 20%", "20% to 0%"), hr = exp(est), lower_95 = exp(est - 1.96 * se), upper_95 = exp(est + 1.96 * se))
}
edf_of <- function(f) vapply(f$smooth, function(s) sum(f$edf[s$first.para:s$last.para]), 0)
rows <- list()
for (a in run_ages) {
  d <- dat[age_band == a]
  # Collapse to cells; every model variable must be constant within (region, calendar_year).
  cells <- d[, c(lapply(.SD, function(x) x[1]), .(deaths = sum(death), children = .N)), by = .(region, calendar_year), .SDcols = setdiff(vars, c("region", "calendar_year"))]
  chk <- d[, lapply(.SD, uniqueN), by = .(region, calendar_year), .SDcols = setdiff(vars, c("region", "calendar_year"))]
  stopifnot(all(as.matrix(chk[, -(1:2)]) == 1L), sum(cells$children) == nrow(d), sum(cells$deaths) == sum(d$death))
  cells[, survivors := children - deaths]
  for (v in c("survey", "country", "region")) cells[[v]] <- factor(cells[[v]])
  kt <- knot_table[age_band == a]
  knots <- lapply(c(pfpr_pct = "pfpr_pct", calendar_year = "calendar_year"), function(v) kt[variable == v][order(index)]$value)
  saved <- readRDS(manifest[age_band == a]$model_file)$fit
  fits <- list(bam_discrete_child = saved)
  message(a, ": ", nrow(d), " children in ", nrow(cells), " cells; fitting gam (REML)")
  t_gam <- system.time(fits$gam_reml_cell <- gam(form_cell, data = cells, knots = knots, family = binomial(link = "cloglog"),
    method = "REML", gamma = 2, select = FALSE, na.action = na.fail, control = gam.control(nthreads = 2)))["elapsed"]
  message(a, ": fitting bam (fREML, no discretisation) on cells")
  t_bam <- system.time(fits$bam_fREML_cell <- bam(form_cell, data = cells, knots = knots, family = binomial(link = "cloglog"),
    method = "fREML", discrete = FALSE, gamma = 2, select = FALSE, nthreads = 2, na.action = na.fail))["elapsed"]
  elapsed <- c(bam_discrete_child = NA_real_, gam_reml_cell = unname(t_gam), bam_fREML_cell = unname(t_bam))
  saveRDS(lapply(fits[-1], function(f) list(coef = coef(f), Vp = f$Vp, sp = f$sp, edf = edf_of(f))), file.path(private, sprintf("fits%s_%d.rds", suffix, match(a, ages))))
  # Log-likelihood of the child-level Bernoulli model at each fit's coefficients, on the same scale.
  ll <- function(f) { eta <- predict(f, cells, type = "link"); p <- -expm1(-exp(eta)); sum(cells$deaths * log(p) + cells$survivors * log1p(-p)) }
  for (m in names(fits)) {
    f <- fits[[m]]; e <- edf_of(f); cz <- grep("^z_", names(coef(f)), value = TRUE)
    rows[[paste(a, m)]] <- cbind(data.table(age_band = a, fit = m, cells = nrow(cells), elapsed_seconds = elapsed[[m]], converged = isTRUE(f$converged),
      loglik = ll(f), edf_total = sum(f$edf), edf_pfpr = e[1],
      sp_pfpr = f$sp[1], max_abs_covariate_logHR_diff = max(abs(coef(f)[cz] - coef(saved)[cz]))), as.data.table(as.list(setNames(e[-1], sub("^s\\((.*)\\)$", "edf_\\1", vapply(f$smooth[-1], `[[`, "", "label"))))), dcast(contrasts(f)[, .(k = 1, contrast, hr)], k ~ contrast, value.var = "hr")[, -1])
  }
  res <- rbindlist(rows, fill = TRUE); cbh_atomic_csv(as.data.frame(res), file.path(out, paste0("comparison", suffix, ".csv")))
  print(res[age_band == a]); rm(fits, saved, d, cells); gc(FALSE)
}
message("bam vs gam comparison written: ", out)

#!/usr/bin/env Rscript
# Support and identification checks for the PfPR x calendar-time interaction (added 25 September
# 2026, following an independent review of the tensor-product sensitivity). For the 24-35, 36-47
# and 48-59 month bands, on the same collapsed cells, without the region random intercept:
#   support  : children and deaths by PfPR class x period
#   binned   : PfPR class x period indicators (reference class 20-40%) replacing s(pfpr_pct)
#   svy_fe   : separate splines + ti(PfPR, year) with survey fixed effects (interaction within surveys)
#   fixed_re : separate splines + ti with the survey and country RE smoothing parameters held at the
#              separate-spline fit's values (does the interaction add fit beyond a survey RE of that size?)
#   window   : separate, tensor and ti refitted on band entries 2005-2019 only
# Aggregates only are written.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
suppressPackageStartupMessages({library(data.table); library(mgcv)})
base <- cbh_primary_settings("regional_mics")
out <- "results/cbh/tensor_pfpr_year_dhsmics_map_gamma2_v2"
ages <- cbh_config()$age_bands$age_band; bands <- c("24-35", "36-47", "48-59")
covs <- paste0("z_", cbh_primary_regional_spec(base)$covariates)
vars <- unique(c("region", "calendar_year", "pfpr_pct", "band_years", "survey", "country", covs))
dat <- as.data.table(readRDS(base$data)$data)[age_band %in% bands, c("death", "age_band", vars), with = FALSE]
knot_table <- fread(base$knots)
S <- "s(pfpr_pct, bs = \"cr\", k = 5) + s(calendar_year, bs = \"cr\", k = 6)"
TI <- "ti(pfpr_pct, calendar_year, bs = c(\"cr\", \"cr\"), k = c(5, 6))"
TE <- "te(pfpr_pct, calendar_year, bs = c(\"cr\", \"cr\"), k = c(5, 6))"
RE <- "s(survey, bs = \"re\") + s(country, bs = \"re\")"
fm <- function(...) as.formula(paste("cbind(deaths, survivors) ~", paste(c(..., covs, "offset(log(band_years))"), collapse = " + ")))
fitb <- function(form, d, knots, sp = NULL) { set.seed(20260907L)
  suppressWarnings(bam(form, data = d, knots = knots, family = binomial(link = "cloglog"), method = "fREML", discrete = FALSE,
    gamma = 2, nthreads = 2, na.action = na.fail, sp = sp, control = gam.control(maxit = 100))) }
term_edf <- function(f, pat) sum(f$edf[unlist(lapply(f$smooth[grepl(pat, vapply(f$smooth, `[[`, "", "label"))], function(s) s$first.para:s$last.para))])
ti_p <- function(f) { st <- summary(f)$s.table; r <- grep("^ti\\(", rownames(st)); if (length(r)) st[r, "p-value"] else NA_real_ }
hr_by_year <- function(f, d, yrs = c(2005.5, 2020.5)) {
  nd <- d[rep(1L, length(yrs))]; nd$calendar_year <- yrs; a <- copy(nd); a$pfpr_pct <- 0; b <- copy(nd); b$pfpr_pct <- 20
  L <- predict(f, a, type = "lpmatrix") - predict(f, b, type = "lpmatrix"); exp(drop(L %*% coef(f)))
}
cls_lab <- c("<2", "2-5", "5-10", "10-20", "20-40", "40+"); per_lab <- c("2000-04", "2005-09", "2010-14", "2015-19", "2020-23")
support <- binned <- tests <- list()
for (a in bands) {
  d <- dat[age_band == a]
  cells <- d[, c(lapply(.SD, function(x) x[1]), .(deaths = sum(death), children = .N)), by = .(region, calendar_year), .SDcols = setdiff(vars, c("region", "calendar_year"))]
  cells[, survivors := children - deaths]; for (v in c("survey", "country")) cells[[v]] <- factor(cells[[v]])
  cells[, pclass := cut(pfpr_pct, c(0, 2, 5, 10, 20, 40, 101), right = FALSE, labels = cls_lab)]
  cells[, period := cut(calendar_year, c(2000, 2005, 2010, 2015, 2020, 2025), right = FALSE, labels = per_lab)]
  support[[a]] <- cells[, .(children = sum(children), deaths = sum(deaths), countries = uniqueN(country)), by = .(pclass, period)][, age_band := a]
  kt <- knot_table[age_band == a]
  knots <- lapply(c(pfpr_pct = "pfpr_pct", calendar_year = "calendar_year"), function(v) kt[variable == v][order(index)]$value)
  message(a, ": binned PfPR class x period")
  cells[, cls4 := relevel(factor(cut(pfpr_pct, c(0, 5, 10, 20, 40, 101), right = FALSE, labels = c("<5", "5-10", "10-20", "20-40", "40+"))), ref = "20-40")]
  cells[, per4 := factor(cut(calendar_year, c(2000, 2006, 2011, 2016, 2025), right = FALSE, labels = c("2000-05", "2006-10", "2011-15", "2016-23")))]
  cells[, cp := interaction(cls4, per4, drop = TRUE)]
  cells[, cp := relevel(cp, ref = "20-40.2011-15")]
  fb <- fitb(as.formula(paste("cbind(deaths, survivors) ~ cp + s(calendar_year, bs = \"cr\", k = 6) +", RE, "+", paste(c(covs, "offset(log(band_years))"), collapse = " + "))), cells, knots)
  pt <- summary(fb)$p.table; b <- setNames(c(0, pt[grep("^cp", rownames(pt)), 1]), c("20-40.2011-15", sub("^cp", "", grep("^cp", rownames(pt), value = TRUE))))
  V <- vcov(fb)[grep("^cp", names(coef(fb))), grep("^cp", names(coef(fb)))]
  for (per in levels(cells$per4)) for (cl in c("<5", "5-10", "10-20")) {
    num <- paste(cl, per, sep = "."); den <- paste("20-40", per, sep = ".")
    if (!num %in% names(b) || !den %in% names(b)) next
    w <- setNames(numeric(nrow(V)), rownames(V)); if (paste0("cp", num) %in% names(w)) w[paste0("cp", num)] <- 1; if (paste0("cp", den) %in% names(w)) w[paste0("cp", den)] <- -1
    binned[[paste(a, per, cl)]] <- data.table(age_band = a, period = per, contrast = sprintf("%s vs 20-40%%", cl), log_hr = unname(b[num] - b[den]),
      se = sqrt(drop(t(w) %*% V %*% w)), deaths_in_class = cells[cls4 == cl & per4 == per, sum(deaths)])
  }
  message(a, ": separate, ti with fixed RE variances, ti with survey fixed effects")
  fs <- fitb(fm(S, RE), cells, knots)
  fti <- fitb(fm(S, TI, RE), cells, knots)
  ffix <- fitb(fm(S, TI, RE), cells, knots, sp = c(-1, -1, -1, -1, fs$sp[c("s(survey)", "s(country)")]))
  fsfe <- fitb(fm(S, "survey"), cells, knots); ftfe <- fitb(fm(S, TI, "survey"), cells, knots)
  w0 <- cells[calendar_year >= 2005 & calendar_year < 2020]; w0[, survey := droplevels(survey)]; w0[, country := droplevels(country)]
  message(a, ": 2005-2019 window")
  ws <- fitb(fm(S, RE), w0, knots); wt <- fitb(fm(TE, RE), w0, knots); wi <- fitb(fm(S, TI, RE), w0, knots)
  row <- function(check, f, ref, d) { h <- hr_by_year(f, d)
    data.table(age_band = a, check, delta_aic = AIC(f) - AIC(ref), ti_edf = term_edf(f, "^ti\\("), ti_p = ti_p(f), survey_edf = term_edf(f, "^s\\(survey\\)"), hr_2005 = h[1], hr_2020 = h[2]) }
  tests[[a]] <- rbind(row("ti, free (as in 01_fit.R)", fti, fs, cells), row("ti, survey and country RE variances fixed at separate fit", ffix, fs, cells),
    row("ti, survey fixed effects (versus separate with survey fixed effects)", ftfe, fsfe, cells),
    row("tensor, entries 2005-2019 (versus separate, same window)", wt, ws, w0), row("ti, entries 2005-2019 (versus separate, same window)", wi, ws, w0))
  cbh_atomic_csv(as.data.frame(rbindlist(tests)), file.path(out, "interaction_checks.csv"))
  cbh_atomic_csv(as.data.frame(rbindlist(binned)), file.path(out, "binned_pfpr_by_period.csv"))
  cbh_atomic_csv(as.data.frame(rbindlist(support)), file.path(out, "support_pfpr_class_by_period.csv"))
  rm(d, cells, w0, fb, fs, fti, ffix, fsfe, ftfe, ws, wt, wi); gc(FALSE)
}
print(rbindlist(tests)); message("Interaction checks written: ", out)

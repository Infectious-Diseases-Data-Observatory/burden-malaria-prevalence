#!/usr/bin/env Rscript
# Multiple imputation of the remaining survey-region covariate gaps on the combined
# DHS and MICS overlay, with the specification of the DHS-only version
# (R_cbh/covariates/16_impute_regional.R): chained equations with predictive mean
# matching over the 13 regional covariates and regional MAP PfPR at the survey year,
# with survey year, the national series at survey year and a country factor as
# predictors. DHS and MICS survey-regions are imputed jointly; observed values are
# never changed.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/specification.R")
source("R_cbh/sensitivity/imputation_dhsmics/settings.R")
suppressMessages(library(mice))
st <- cbh_imputed_dhsmics_settings()
out_private <- st$imputation_private; out <- st$imputation_results
spec <- cbh_primary_regional_spec(st); regional <- spec$regional
stopifnot(length(regional) == 13L)
dhs <- cbh_read_csv(st$dhs_overlay)[c("survey", "regkey", regional)]; dhs$programme <- "DHS"
mics <- cbh_read_csv(st$mics_overlay)[c("survey", "regkey", regional)]; mics$programme <- "MICS"
wide <- rbind(dhs, mics); cbh_unique(wide, c("survey", "regkey"), "Combined regional overlay")
stopifnot(all(startsWith(mics$survey, "MC_")), !any(startsWith(dhs$survey, "MC_")))
reg <- do.call(rbind, lapply(st$registries, function(x) cbh_read_csv(x)[c("svkey", "iso3", "year")]))
cbh_unique(reg, "svkey", "Combined registry")
wide$country <- reg$iso3[match(wide$survey, reg$svkey)]
wide$survey_year <- reg$year[match(wide$survey, reg$svkey)]
stopifnot(!anyNA(wide$country), !anyNA(wide$survey_year))
map <- do.call(rbind, lapply(st$annual_maps, function(x) cbh_read_csv(x)[c("svkey", "regkey", "year", "pfpr2_10")]))
wide$pfpr_survey_year <- map$pfpr2_10[match(paste(wide$survey, wide$regkey, wide$survey_year), paste(map$svkey, map$regkey, map$year))]
# Survey-regions without MAP prevalence have no analysis records (checked when the
# analysis dataset is prepared); they are dropped rather than imputed.
no_map <- wide[!is.finite(wide$pfpr_survey_year), c("programme", "survey", "regkey")]
stopifnot(all(no_map$programme == "MICS"))
wide <- wide[is.finite(wide$pfpr_survey_year), ]
nat <- cbh_read_csv(st$imputed_national)
k <- match(paste(wide$country, wide$survey_year), paste(nat$iso3, nat$year)); stopifnot(!anyNA(k))
wide$log_gdp_pc <- nat$log_gdp_pc[k]; wide$log_health_expenditure_pc <- nat$log_health_expenditure_pc[k]
wide$political_stability <- nat$political_stability[k]

## ---- imputation model -----------------------------------------------------------------------
vars <- c(regional, "pfpr_survey_year")
X <- wide[c(vars, "survey_year", "log_gdp_pc", "log_health_expenditure_pc", "political_stability", "country")]
X$country <- factor(X$country)
before <- X[regional]
meth <- make.method(X); meth[vars] <- "pmm"; meth[setdiff(names(X), vars)] <- ""
pred <- make.predictorMatrix(X)
imp <- mice(X, m = st$m, maxit = st$maxit, method = meth, predictorMatrix = pred, seed = st$imputation_seed, printFlag = FALSE)
stopifnot(!any(is.na(complete(imp, 1)[vars])))
completed <- lapply(seq_len(st$m), function(i) complete(imp, i)[regional])
point <- Reduce(`+`, completed) / st$m
for (v in regional) stopifnot(all(point[[v]] >= min(before[[v]], na.rm = TRUE) - 1e-9),
                              all(point[[v]] <= max(before[[v]], na.rm = TRUE) + 1e-9))

## ---- outputs ----------------------------------------------------------------------------------
res <- wide[c("programme", "survey", "country", "survey_year", "regkey")]
for (v in regional) {
  res[[v]] <- point[[v]]
  res[[paste0(v, "_before_model_imputation")]] <- before[[v]]
  res[[paste0(v, "_model_imputed")]] <- !is.finite(before[[v]])
}
stopifnot(all(vapply(regional, function(v) { o <- !res[[paste0(v, "_model_imputed")]]
  all(abs(res[[v]][o] - before[[v]][o]) < 1e-12) }, logical(1))))
cbh_atomic_csv(res, file.path(out_private, "regional_covariates_wide.csv"))
cbh_atomic_rds(list(imputations = completed, keys = wide[c("survey", "regkey")],
  imputed = as.data.frame(lapply(before, function(z) !is.finite(z))), m = st$m, maxit = st$maxit,
  seed = st$imputation_seed, method = "pmm", predictors = names(X), loggedEvents = imp$loggedEvents,
  dropped_without_map = no_map), file.path(out_private, "regional_imputations.rds"))
by_var <- do.call(rbind, lapply(regional, function(v) { ii <- res[[paste0(v, "_model_imputed")]]
  data.frame(variable = v, survey_regions = nrow(res), imputed_regions = sum(ii),
    imputed_regions_dhs = sum(ii & res$programme == "DHS"), imputed_regions_mics = sum(ii & res$programme == "MICS"),
    imputed_surveys = length(unique(res$survey[ii])), observed_mean = mean(before[[v]], na.rm = TRUE),
    imputed_mean = if (any(ii)) mean(point[[v]][ii]) else NA_real_,
    between_imputation_sd = if (any(ii)) mean(apply(sapply(completed, function(d) d[[v]][ii]), 1, sd)) else NA_real_) }))
cbh_atomic_csv(by_var, file.path(out, "regional_imputation_by_variable.csv"))
by_survey <- do.call(rbind, lapply(regional, function(v) { ii <- res[[paste0(v, "_model_imputed")]]
  if (!any(ii)) return(NULL)
  t <- as.data.frame(table(survey = res$survey[ii]), stringsAsFactors = FALSE); names(t)[2] <- "regions_imputed"
  t$programme <- res$programme[match(t$survey, res$survey)]; t$country <- res$country[match(t$survey, res$survey)]
  t$variable <- v; t[c("programme", "survey", "country", "variable", "regions_imputed")] }))
cbh_atomic_csv(by_survey, file.path(out, "regional_imputation_by_survey.csv"))
cbh_atomic_csv(no_map, file.path(out, "regions_without_map_dropped.csv"))
ragg::agg_png(file.path(out, "mice_convergence.png"), width = 1800, height = 2400, res = 200)
print(plot(imp, layout = c(4, 7))); dev.off()
ragg::agg_png(file.path(out, "mice_density.png"), width = 2200, height = 1600, res = 200)
plotted <- regional[vapply(regional, function(v) sum(!is.finite(before[[v]])) >= 2L, logical(1))]
print(densityplot(imp, as.formula(paste("~", paste(plotted, collapse = "+"))))); dev.off()
fmt_na <- function(x) ifelse(is.na(x), "—", sprintf("%.2f", x))
writeLines(c("# Regional covariate imputation (DHS and MICS, imputed_dhsmics_v6)", "",
  sprintf("Chained-equation multiple imputation (`mice`, predictive mean matching, m = %d, maxit = %d, seed %d) at the survey-region level on %d survey-regions (%d DHS, %d MICS) from %d surveys, imputed jointly. Targets: the 13 regional covariates and regional MAP PfPR at the survey year; predictors: all of these plus survey year, the national series at survey year (log GDP, log health expenditure, political stability, from `01_impute_national.R`) and a country factor. The specification is that of the DHS-only version (`R_cbh/covariates/16_impute_regional.R`); with the larger donor pool and a different random stream, DHS imputed values differ from that version. Observed values are unchanged. The point overlay is the mean of the %d imputations; all are kept for propagation.",
    st$m, st$maxit, st$imputation_seed, nrow(res), sum(res$programme == "DHS"), sum(res$programme == "MICS"), length(unique(res$survey)), st$m), "",
  sprintf("%d MICS survey-regions without MAP prevalence at the survey year and without analysis records were dropped: %s.", nrow(no_map), paste(no_map$survey, no_map$regkey, collapse = ", ")), "",
  "| Variable | Survey-regions imputed (DHS / MICS) | Surveys | Observed mean | Imputed mean | Between-imputation SD |", "|---|---:|---:|---:|---:|---:|",
  sprintf("| %s | %d (%d / %d) | %d | %.2f | %s | %s |", by_var$variable, by_var$imputed_regions, by_var$imputed_regions_dhs,
    by_var$imputed_regions_mics, by_var$imputed_surveys, by_var$observed_mean, fmt_na(by_var$imputed_mean), fmt_na(by_var$between_imputation_sd)), "",
  "MICS gaps are whole-survey wasting and stunting gaps (Nigeria 2021, Madagascar 2012 South, Somalia 2011 North-East and Somaliland 2011). Somalia's only other survey with anthropometry (MICS 2006) reports one national wasting and stunting value, so the country factor carries little Somalia-specific information.", "",
  "[Convergence traces](mice_convergence.png) · [Observed versus imputed densities](mice_density.png) · [By variable](regional_imputation_by_variable.csv) · [By survey](regional_imputation_by_survey.csv)"),
  file.path(out, "REGIONAL_IMPUTATION.md"))
paths <- c(st$dhs_overlay, st$mics_overlay, st$registries, st$annual_maps, st$imputed_national,
  "R_cbh/sensitivity/imputation_dhsmics/02_impute_regional.R", "R_cbh/sensitivity/imputation_dhsmics/settings.R",
  "R_cbh/primary/specification.R")
cbh_atomic_csv(data.frame(file = paths, md5 = vapply(paths, cbh_file_hash, "")), file.path(out, "regional_imputation_provenance.csv"))
print(by_var[, c("variable", "imputed_regions", "imputed_regions_dhs", "imputed_regions_mics", "observed_mean", "imputed_mean")])
message("Regional covariate imputation (DHS and MICS) complete")

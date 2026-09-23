#!/usr/bin/env Rscript
# National annual covariate gaps for the DHS+MICS imputed-covariate sensitivity (v6).
# Same methods as the DHS-only version (R_cbh/covariates/15_impute_national.R) on the
# 40 countries of the combined sample, with the extra gaps the MICS countries bring:
#   political stability : no WGI round in 2001 -> mean of 2000 and 2002 (deterministic);
#                         South Sudan before independence (2000-2010) -> GAM, with draws
#   GDP per capita (log): South Sudan outside 2008-2015 -> GAM, with draws
#   health expenditure  : Zimbabwe 2000-2009, Somalia 2000-2012, South Sudan outside
#                         2017-2023, every country's 2024 -> the v4 GAM (year smooth,
#                         log GDP, country intercept and slope), with draws
# GAM draws: coefficients from the posterior (Vp) plus residual noise. Health-expenditure
# draw k uses GDP draw k where GDP is imputed. Aggregate outputs only.
source("R_cbh/load_pipeline.R")
source("R_cbh/sensitivity/imputation_dhsmics/settings.R")
library(mgcv)
st <- cbh_imputed_dhsmics_settings()
out_private <- st$imputation_private; out <- st$imputation_results
for (p in c(out_private, out)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
cfg <- cbh_config()
years <- cfg$first_entry_year:cfg$last_entry_year
manifests <- rbind(cbh_read_csv(st$dhs_manifest)[c("survey", "country", "status")],
                   cbh_read_csv(st$mics_manifest)[c("survey", "country", "status")])
countries <- sort(unique(manifests$country[manifests$status %in% c("built", "cached")]))
stopifnot(length(countries) == st$expected_countries)
polstab <- cbh_read_csv(cfg$annual_panels$governance$path)
hexp <- cbh_read_csv(cfg$annual_panels$health_spending$path)
gdp <- cbh_read_csv(cfg$annual_panels$gdp$path)
grid <- expand.grid(iso3 = countries, year = years, stringsAsFactors = FALSE)
pick <- function(panel, col) panel[[col]][match(paste(grid$iso3, grid$year), paste(panel$iso3, panel$year))]
grid$political_stability <- pick(polstab, "polstab")
grid$gdp_pc <- pick(gdp, "gdp_pc"); grid$health_expenditure_pc <- pick(hexp, "hexp_pc")
pos_log <- function(x) ifelse(is.finite(x) & x > 0, log(x), NA_real_)
grid$log_gdp_pc <- pos_log(grid$gdp_pc)
grid$log_health_expenditure_pc <- pos_log(grid$health_expenditure_pc)
grid$political_stability_source <- ifelse(is.finite(grid$political_stability), "wgi_observed", NA)
grid$log_gdp_pc_source <- ifelse(is.finite(grid$log_gdp_pc), "wdi_observed", NA)
grid$log_health_expenditure_pc_source <- ifelse(is.finite(grid$log_health_expenditure_pc), "ghed_observed", NA)
grid$iso3_f <- factor(grid$iso3, levels = countries)
# Where a model fill sits relative to the country's own observed series.
position <- function(v, i) vapply(i, function(r) {
  y <- grid$year[grid$iso3 == grid$iso3[r] & is.finite(grid[[v]]) & grepl("observed$", grid[[paste0(v, "_source")]])]
  if (!length(y)) "no_series" else if (grid$year[r] < min(y)) "backcast_before_series"
  else if (grid$year[r] > max(y)) "forecast_after_series" else "within_series_gap"
}, "")

## ---- political stability: interpolate 2001 where 2000 and 2002 exist -------------------
ps0 <- function(iso, yr) grid$political_stability[grid$iso3 == iso & grid$year == yr]
miss_ps <- which(!is.finite(grid$political_stability))
interp <- miss_ps[grid$year[miss_ps] == 2001L & vapply(miss_ps, function(i)
  is.finite(ps0(grid$iso3[i], 2000L)) && is.finite(ps0(grid$iso3[i], 2002L)), TRUE)]
for (i in interp) grid$political_stability[i] <- (ps0(grid$iso3[i], 2000L) + ps0(grid$iso3[i], 2002L)) / 2
grid$political_stability_source[interp] <- "wgi_2001_interpolated"

## ---- GAM fills with draws ------------------------------------------------------------------
n_draws <- st$m; set.seed(st$imputation_seed)
gam_fill <- function(v, rhs, rows_fit) {
  obs <- grid[rows_fit, ]
  fit <- gam(reformulate(rhs, v), data = transform(obs, iso3 = iso3_f), method = "REML")
  list(fit = fit, residual_sd = sqrt(fit$sig2), n_observed = nrow(obs))
}
draw_cells <- function(g, miss, nd_list) {
  beta <- MASS::mvrnorm(n_draws, coef(g$fit), g$fit$Vp)
  noise <- matrix(rnorm(length(miss) * n_draws, 0, g$residual_sd), ncol = n_draws)
  vapply(seq_len(n_draws), function(k) drop(predict(g$fit, nd_list[[k]], type = "lpmatrix") %*% beta[k, ]), numeric(length(miss))) + noise
}
re_terms <- c("s(year,k=6)", "s(iso3,bs='re')", "s(iso3,year,bs='re')")
# Political stability: countries with a gap beyond the 2001 round (South Sudan).
miss_ps <- which(!is.finite(grid$political_stability))
stopifnot(all(grid$iso3[miss_ps] == "SSD"))
g_ps <- gam_fill("political_stability", re_terms, which(grid$political_stability_source %in% "wgi_observed"))
nd <- data.frame(year = grid$year[miss_ps], iso3 = grid$iso3_f[miss_ps])
ps_draws <- draw_cells(g_ps, miss_ps, rep(list(nd), n_draws))
grid$political_stability_source[miss_ps] <- paste0("gam_", position("political_stability", miss_ps))
grid$political_stability[miss_ps] <- drop(predict(g_ps$fit, nd))
# GDP (log): South Sudan outside 2008-2015.
miss_gdp <- which(!is.finite(grid$log_gdp_pc))
stopifnot(all(grid$iso3[miss_gdp] == "SSD"))
g_gdp <- gam_fill("log_gdp_pc", re_terms, which(is.finite(grid$log_gdp_pc)))
nd <- data.frame(year = grid$year[miss_gdp], iso3 = grid$iso3_f[miss_gdp])
gdp_draws <- draw_cells(g_gdp, miss_gdp, rep(list(nd), n_draws))
grid$log_gdp_pc_source[miss_gdp] <- paste0("gam_", position("log_gdp_pc", miss_gdp))
grid$log_gdp_pc[miss_gdp] <- drop(predict(g_gdp$fit, nd))
grid$gdp_pc <- exp(grid$log_gdp_pc)
# Health expenditure (log): the v4 model, fitted to country-years with observed GDP, so
# South Sudan's 2017-2023 values (no GDP) inform nothing and its fills are population-level
# predictions with the prior country-effect uncertainty.
miss_he <- which(!is.finite(grid$log_health_expenditure_pc))
g_he <- gam_fill("log_health_expenditure_pc", c("s(year,k=6)", "log_gdp_pc", "s(iso3,bs='re')", "s(iso3,year,bs='re')"),
  which(is.finite(grid$log_health_expenditure_pc) & grid$log_gdp_pc_source == "wdi_observed"))
gdp_draw_for <- function(rows, k) { v <- grid$log_gdp_pc[rows]; j <- match(rows, miss_gdp); v[!is.na(j)] <- gdp_draws[j[!is.na(j)], k]; v }
he_nd <- lapply(seq_len(n_draws), function(k)
  data.frame(year = grid$year[miss_he], log_gdp_pc = gdp_draw_for(miss_he, k), iso3 = grid$iso3_f[miss_he]))
he_draws <- draw_cells(g_he, miss_he, he_nd)
nd <- data.frame(year = grid$year[miss_he], log_gdp_pc = grid$log_gdp_pc[miss_he], iso3 = grid$iso3_f[miss_he])
grid$log_health_expenditure_pc_source[miss_he] <- paste0("gam_", position("log_health_expenditure_pc", miss_he))
grid$log_health_expenditure_pc[miss_he] <- drop(predict(g_he$fit, nd))
grid$health_expenditure_pc <- exp(grid$log_health_expenditure_pc)
grid$iso3_f <- NULL
stopifnot(all(is.finite(grid$political_stability)), all(is.finite(grid$log_gdp_pc)),
  all(is.finite(grid$log_health_expenditure_pc)), all(is.finite(ps_draws)), all(is.finite(gdp_draws)),
  all(is.finite(he_draws)), !anyNA(grid[grep("_source$", names(grid))]))

## ---- outputs ---------------------------------------------------------------------------------
cbh_atomic_csv(grid, file.path(out_private, "national_covariates_imputed.csv"))
cbh_atomic_csv(grid, file.path(out, "national_covariates_imputed.csv"))
block <- function(rows, draws, g) list(keys = grid[rows, c("iso3", "year")], draws = draws,
  model = paste(deparse(formula(g$fit)), collapse = ""))
cbh_atomic_rds(list(log_health_expenditure_pc = block(miss_he, he_draws, g_he),
  log_gdp_pc = block(miss_gdp, gdp_draws, g_gdp), political_stability = block(miss_ps, ps_draws, g_ps),
  political_stability_2001_note = "deterministic interpolation; no draws",
  n_draws = n_draws, seed = st$imputation_seed), file.path(out_private, "national_imputations.rds"))
summ <- function(v, g, miss) data.frame(variable = v, model = paste(deparse(formula(g$fit)), collapse = ""),
  n_observed_in_fit = g$n_observed, n_imputed = length(miss),
  imputed_countries = paste(sort(unique(grid$iso3[miss])), collapse = " "),
  residual_sd = g$residual_sd, deviance_explained = summary(g$fit)$dev.expl)
check <- rbind(summ("political_stability", g_ps, miss_ps), summ("log_gdp_pc", g_gdp, miss_gdp),
               summ("log_health_expenditure_pc", g_he, miss_he))
check$interpolated_2001 <- c(length(interp), 0L, 0L)
cbh_atomic_csv(check, file.path(out, "national_imputation_summary.csv"))
sources <- do.call(rbind, lapply(c("political_stability", "log_gdp_pc", "log_health_expenditure_pc"), function(v) {
  s <- paste0(v, "_source"); t <- aggregate(list(country_years = grid$year), list(source = grid[[s]], iso3 = grid$iso3), length)
  t <- t[!grepl("observed$", t$source), ]; if (nrow(t)) cbind(variable = v, t) else NULL }))
cbh_atomic_csv(sources, file.path(out, "national_imputation_by_country.csv"))
focus <- grid[grid$iso3 %in% c("SOM", "SSD", "ZWE"), c("iso3", "year", "political_stability", "political_stability_source",
  "gdp_pc", "log_gdp_pc_source", "health_expenditure_pc", "log_health_expenditure_pc_source")]
cbh_atomic_csv(focus, file.path(out, "somalia_south_sudan_zimbabwe_series.csv"))
# The years the analysis actually uses (band entry years of the affected MICS surveys).
used <- function(iso, yrs) sprintf("%s %d-%d", iso, min(yrs), max(yrs))
writeLines(c("# National covariate imputation (DHS and MICS, imputed_dhsmics_v6)", "",
  sprintf("Countries: the %d countries with built DHS or MICS surveys; years %d-%d. Methods follow the DHS-only version (`R_cbh/covariates/15_impute_national.R`); the MICS countries add Somalia, South Sudan, the Central African Republic and Guinea-Bissau.",
    length(countries), min(years), max(years)), "",
  sprintf("**Political stability (WGI):** %d country-years in 2001 (no WGI round) filled by the mean of each country's 2000 and 2002 estimates (deterministic, `wgi_2001_interpolated`). South Sudan has no WGI series before independence: its %d country-years 2000-2010 are filled from `%s` fitted by REML to %d observed country-years (residual SD %.3f), with %d draws saved.",
    length(interp), length(miss_ps), check$model[1], check$n_observed_in_fit[1], check$residual_sd[1], n_draws), "",
  sprintf("**GDP per capita (log, current US$):** South Sudan is observed only for 2008-2015; its %d missing country-years are filled from `%s` (%d observed country-years, residual SD %.3f), with draws. Only 2005-2007 enter the analysis (MC_SSD2010 band entry years).",
    length(miss_gdp), check$model[2], check$n_observed_in_fit[2], check$residual_sd[2]), "",
  sprintf("**Health expenditure per capita (log, current US$):** %d country-years filled from the DHS-only model `%s`, fitted to the %d observed country-years with observed GDP (residual SD %.3f, deviance explained %.1f%%): Zimbabwe 2000-2009, Somalia 2000-2012, South Sudan outside 2017-2023 and every country in 2024. South Sudan's observed 2017-2023 values have no GDP and are not used in the fit, so its fills are population-level predictions with the country-effect prior uncertainty. Health-expenditure draw k uses GDP draw k where GDP is imputed.",
    length(miss_he), check$model[3], check$n_observed_in_fit[3], check$residual_sd[3], 100 * check$deviance_explained[3]), "",
  "Fill labels give each fill's position relative to the country's own observed series (`gam_backcast_before_series`, `gam_within_series_gap`, `gam_forecast_after_series`). The Somalia (13 years) and South Sudan (up to 17 years) health-expenditure fills and the South Sudan pre-independence political stability are long extrapolations; they affect MC_SOM2006, MC_SOM2011NE, MC_SOM2011SL and MC_SSD2010 only. Their values are listed in the country table below.", "",
  "Child HIV incidence for São Tomé and Príncipe (MICS 2014 and 2019) comes from the extended incidence panel (latent adolescent series); the other MICS countries have their own UNAIDS series.", "",
  "[Imputed panel](national_covariates_imputed.csv) · [Summary](national_imputation_summary.csv) · [Fills by country](national_imputation_by_country.csv) · [Somalia, South Sudan and Zimbabwe series](somalia_south_sudan_zimbabwe_series.csv)"),
  file.path(out, "NATIONAL_IMPUTATION.md"))
paths <- c(cfg$annual_panels$governance$path, cfg$annual_panels$health_spending$path, cfg$annual_panels$gdp$path,
  st$dhs_manifest, st$mics_manifest, "R_cbh/sensitivity/imputation_dhsmics/01_impute_national.R",
  "R_cbh/sensitivity/imputation_dhsmics/settings.R")
cbh_atomic_csv(data.frame(file = paths, md5 = vapply(paths, cbh_file_hash, "")), file.path(out, "national_imputation_provenance.csv"))
print(check[, c("variable", "n_observed_in_fit", "n_imputed", "imputed_countries", "residual_sd")])
message("National covariate imputation (DHS and MICS) complete")

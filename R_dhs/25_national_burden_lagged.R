# =============================================================================
# 25_national_burden_lagged.R — national malaria-attributable death estimates from
# the primary model fitted with LAGGED prevalence.
#
# SPECIFICATION, as requested:
#   outcome    all-cause post-neonatal mortality from chmort with the default
#              Period = 60 months (the manuscript's own outcome);
#   exposure   regional PfPR2-10 averaged over the mortality window rather than
#              taken at the survey year (lags 0-4, the approximately centre-matched
#              window of script 17);
#   structure  additive penalised spline in PfPR2-10, NO interaction with calendar
#              time, ridge covariate block, country random intercept and random
#              PfPR slope, log-exposure offset;
#   prediction country random effects set to zero (population average), and the
#              national prevalence input is the SAME 5-year trailing mean the model
#              was fitted on - predicting a windowed-exposure model with a
#              single-year input would mismatch the fitted scale.
#
# TARGET YEAR. 2025 was requested but cannot be produced: the MAP
# 202508_Global_Pf_Parasite_Rate product ends at 2024 (verified - getRaster returns
# nothing for 2025 and script 02 deliberately refuses to substitute 2024), and the
# local IGME/World Bank mortality series also ends at 2024. The table is therefore
# for 2024, the latest year with both inputs complete, using the 2020-2024 mean
# prevalence and 2024 all-cause post-neonatal deaths.
#
# Outputs (results/dhs_rebuild/):
#   national_burden_lagged_2024.csv
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("mgcv", "MASS", "countrycode"))
set.seed(20260729)

TARGET_YEAR <- 2024L
WINDOW <- 4L          # trailing mean over TARGET_YEAR-4 .. TARGET_YEAR, matching lags 0-4
NSIM <- 4000L

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
PFPR_LONG_CSV <- file.path(DERIVED_DIR, "map_pfpr_lagged_long.csv")
if (!file.exists(PFPR_LONG_CSV)) {
  message("Run R_dhs/17_sensitivity_timing.R first to build the lagged prevalence panel.")
  quit(save = "no", status = 0)
}
bundle <- readRDS(MODEL_BUNDLE_RDS)
catalog <- bundle$catalog
catalog$included_in_main <- as.logical(catalog$included_in_main)
year_center <- unique(bundle$year_center)[1]

## ---- fit the primary model on WINDOWED regional prevalence ------------------
analysis <- read_analysis_data()
analysis$k <- paste(analysis$svkey, analysis$regkey, sep = "|")
pfpr_long <- read.csv(PFPR_LONG_CSV, stringsAsFactors = FALSE)
pfpr_long$k <- paste(pfpr_long$svkey, pfpr_long$regkey, sep = "|")
window_mean <- local({
  sub <- pfpr_long[pfpr_long$lag <= WINDOW, ]
  tapply(sub$pfpr2_10, sub$k, mean)
})
fit_data <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
prev <- as.numeric(window_mean[fit_data$k])
fit_data$pfpr2_10 <- prev
fit_data$pfpr10 <- prev / 10
fit_data <- fit_data[is.finite(fit_data$postneonatal_mortality) &
                       fit_data$postneonatal_mortality > 0 &
                       is.finite(fit_data$exposure) & fit_data$exposure > 0 &
                       is.finite(fit_data$pfpr10) & fit_data$pfpr10 > 0, , drop = FALSE]
fit <- fit_ridge_gam(fit_data, "postneonatal_mortality", catalog,
                     "spline_no_interaction", method = "REML",
                     preprocessing = bundle$preprocessing)
model <- fit$model
st <- summary(model)$s.table
cat(sprintf("model fitted on %d region-years, %d countries; s(PfPR) edf %.2f, p %.3g\n",
            nrow(fit_data), length(unique(fit_data$iso3)),
            st[grep("^s\\(pfpr10\\)", rownames(st)), "edf"],
            st[grep("^s\\(pfpr10\\)", rownames(st)), "p-value"]))
linear <- fit_ridge_gam(fit_data, "postneonatal_mortality", catalog,
                        "linear_no_interaction", method = "REML",
                        preprocessing = bundle$preprocessing)
lin <- model_summary_row(linear)
cat(sprintf("linear summary: %+.1f%% per +10 PfPR points [%.1f, %.1f]\n",
            lin$pct_change_per_10, lin$pct_change_lo, lin$pct_change_hi))

## ---- national inputs: trailing-mean prevalence and IGME all-cause deaths ----
mortality <- read.csv(file.path(DATA_DIR, "wb_mortality_timeseries.csv"), stringsAsFactors = FALSE)
prevalence <- read.csv(file.path(DATA_DIR, "pfpr_by_country_year.csv"), stringsAsFactors = FALSE)
if (!TARGET_YEAR %in% mortality$year || !TARGET_YEAR %in% prevalence$year) {
  stop("Target year ", TARGET_YEAR, " is absent from the national inputs.")
}
latest <- mortality[mortality$year == TARGET_YEAR,
                    c("iso3", "u5mr", "nnmr", "births", "allcause_1mo5y")]
latest$region <- countrycode::countrycode(latest$iso3, "iso3c", "region", warn = FALSE)
latest <- latest[!is.na(latest$region) & latest$region == "Sub-Saharan Africa" &
                   is.finite(latest$allcause_1mo5y) & latest$allcause_1mo5y > 0 &
                   is.finite(latest$births), , drop = FALSE]
trailing <- do.call(rbind, lapply(unique(latest$iso3), function(is3) {
  v <- prevalence$pfpr_pct[prevalence$iso3 == is3 &
                             prevalence$year >= TARGET_YEAR - WINDOW &
                             prevalence$year <= TARGET_YEAR]
  if (!length(v)) return(NULL)
  data.frame(iso3 = is3, pfpr_window = mean(v), pfpr_years = length(v),
             pfpr_target_year = prevalence$pfpr_pct[prevalence$iso3 == is3 &
                                                      prevalence$year == TARGET_YEAR][1])
}))
national <- merge(latest, trailing, by = "iso3")
national <- national[is.finite(national$pfpr_window) & national$pfpr_window > 0, , drop = FALSE]
cat(sprintf("national inputs: %d sub-Saharan countries for %d (window %d-%d)\n",
            nrow(national), TARGET_YEAR, TARGET_YEAR - WINDOW, TARGET_YEAR))

## ---- attributable fraction and deaths, with parameter uncertainty -----------
high <- newdata_at_mean(model, national$pfpr_window / 10,
                        year_c = TARGET_YEAR - year_center)
low <- newdata_at_mean(model, rep(AF_REFERENCE / 10, nrow(national)),
                       year_c = TARGET_YEAR - year_center)
dX <- population_lpmatrix(model, high) - population_lpmatrix(model, low)
af <- pmax(1 - exp(-as.numeric(dX %*% coef(model))), 0)
draws <- MASS::mvrnorm(NSIM, mu = coef(model), Sigma = vcov(model))
af_sims <- pmax(1 - exp(-(dX %*% t(draws))), 0)
deaths <- af * national$allcause_1mo5y
deaths_sims <- af_sims * national$allcause_1mo5y

national$attributable_fraction_pct <- 100 * af
national$malaria_deaths <- deaths
national$malaria_deaths_lo <- apply(deaths_sims, 1, quantile, 0.025)
national$malaria_deaths_hi <- apply(deaths_sims, 1, quantile, 0.975)
national$deaths_per_1000_births <- 1000 * deaths / national$births
national$country <- countrycode::countrycode(national$iso3, "iso3c", "country.name", warn = FALSE)
national <- national[order(-national$malaria_deaths), , drop = FALSE]

totals <- colSums(deaths_sims)
out <- national[, c("iso3", "country", "pfpr_window", "pfpr_target_year",
                    "allcause_1mo5y", "attributable_fraction_pct",
                    "malaria_deaths", "malaria_deaths_lo", "malaria_deaths_hi",
                    "deaths_per_1000_births", "births", "pfpr_years")]
write.csv(out, file.path(RESULTS_DIR, sprintf("national_burden_lagged_%d.csv", TARGET_YEAR)),
          row.names = FALSE)

cat(sprintf("\n=== Malaria-attributable post-neonatal deaths, %d ===\n", TARGET_YEAR))
cat(sprintf("    Period-60 outcome, prevalence = %d-%d mean, no time interaction\n\n",
            TARGET_YEAR - WINDOW, TARGET_YEAR))
show <- out
show$pfpr_window <- round(show$pfpr_window, 1)
show$allcause_1mo5y <- round(show$allcause_1mo5y)
show$attributable_fraction_pct <- round(show$attributable_fraction_pct, 1)
show$malaria_deaths <- round(show$malaria_deaths)
show$malaria_deaths_lo <- round(show$malaria_deaths_lo)
show$malaria_deaths_hi <- round(show$malaria_deaths_hi)
show$deaths_per_1000_births <- round(show$deaths_per_1000_births, 1)
print(show[, c("country", "pfpr_window", "allcause_1mo5y", "attributable_fraction_pct",
               "malaria_deaths", "malaria_deaths_lo", "malaria_deaths_hi",
               "deaths_per_1000_births")], row.names = FALSE)
cat(sprintf("\nTOTAL across %d countries: %.0f malaria-attributable post-neonatal deaths\n",
            nrow(out), sum(deaths)))
cat(sprintf("  95%% interval from model parameter uncertainty: %.0f to %.0f\n",
            quantile(totals, 0.025), quantile(totals, 0.975)))
cat(sprintf("  all-cause post-neonatal deaths in these countries: %.0f\n",
            sum(national$allcause_1mo5y)))
cat(sprintf("  implied overall attributable fraction: %.1f%%\n",
            100 * sum(deaths) / sum(national$allcause_1mo5y)))

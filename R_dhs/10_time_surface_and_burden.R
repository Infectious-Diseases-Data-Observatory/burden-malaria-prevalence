# =============================================================================
# 10_time_surface_and_burden.R
# Compare three time/prevalence structures on the common DHS/MIS sample:
#   * additive: s(PfPR) + s(year)
#   * ANOVA tensor interaction: s(PfPR) + s(year) + ti(PfPR, year)
#   * full tensor surface: te(PfPR, year)
#
# Then apply population-average attributable fractions (country random effects
# set to zero) to national all-cause post-neonatal deaths. The latest country
# comparison evaluates the mortality surface at 2025 but necessarily uses the
# latest available MAP PfPR, IGME mortality, and live-birth inputs from 2024.
#
# WHO supplies all-age malaria deaths through 2024. The comparison reports both
# those estimates and a clearly labelled under-five proxy using WHO's statement
# that approximately 75% of African-region malaria deaths occur under age five.
#
# Optional environment variable for an already-downloaded WHO API response:
#   WHO_JSON_PATH=/path/to/who_malaria_est_deaths_2024.json
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c(
  "mgcv", "ggplot2", "countrycode", "jsonlite", "httr", "MASS"
))
set.seed(20260727)

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
catalog <- bundle$catalog
catalog$included_in_main <- as.logical(catalog$included_in_main)

prepare_model_data <- function(data) {
  keep <- is.finite(data$postneonatal_mortality) &
    data$postneonatal_mortality > 0 &
    is.finite(data$exposure) & data$exposure > 0 &
    is.finite(data$pfpr10) & is.finite(data$year_c) &
    !is.na(data$iso3)
  dd <- data[keep, , drop = FALSE]
  dd$country <- factor(dd$iso3)
  dd$deaths <- round(dd$postneonatal_mortality / 1000 * dd$exposure)
  ridge <- make_ridge_matrix(dd, catalog, bundle$preprocessing)
  dd$G <- ridge$matrix
  dd
}

fit_full_te <- function(data, method) {
  penalty <- list(G = list(diag(ncol(data$G))))
  model <- mgcv::gam(
    deaths ~
      te(pfpr10, year_c, k = c(6, 6)) +
      G +
      s(country, bs = "re") +
      s(country, pfpr10, bs = "re") +
      offset(log(exposure)),
    family = mgcv::nb(),
    method = method,
    paraPen = penalty,
    data = data
  )
  model
}

model_data <- prepare_model_data(analysis)
message("Fitting full te(PfPR, year) model using ML and REML.")
te_ml <- fit_full_te(model_data, "ML")
te_reml <- fit_full_te(model_data, "REML")

additive_reml <- fit_ridge_gam(
  analysis,
  outcome = "postneonatal_mortality",
  catalog = catalog,
  specification = "spline_no_interaction",
  method = "REML",
  preprocessing = bundle$preprocessing
)$model
ti_reml <- fit_ridge_gam(
  analysis,
  outcome = "postneonatal_mortality",
  catalog = catalog,
  specification = "spline_time_interaction",
  method = "REML",
  preprocessing = bundle$preprocessing
)$model

te_table <- summary(te_ml)$s.table
te_row <- grep(
  "^te\\(pfpr10,year_c\\)",
  rownames(te_table),
  value = TRUE
)
if (!length(te_row)) stop("Could not find the te(PfPR, year) summary row.")

surface_comparison <- rbind(
  data.frame(
    model = bundle$comparison$specification,
    AIC = bundle$comparison$AIC,
    prevalence_edf = bundle$comparison$pfpr_smooth_edf,
    prevalence_p = bundle$comparison$pfpr_p,
    interaction_edf = bundle$comparison$time_interaction_edf,
    interaction_p = bundle$comparison$time_interaction_p,
    full_surface_edf = NA_real_,
    full_surface_p = NA_real_
  ),
  data.frame(
    model = "full_te_surface",
    AIC = AIC(te_ml),
    prevalence_edf = NA_real_,
    prevalence_p = NA_real_,
    interaction_edf = NA_real_,
    interaction_p = NA_real_,
    full_surface_edf = te_table[te_row[1], "edf"],
    full_surface_p = te_table[te_row[1], "p-value"]
  )
)
surface_comparison$dAIC <- surface_comparison$AIC -
  min(surface_comparison$AIC)
surface_comparison <- surface_comparison[
  order(surface_comparison$AIC),
  ,
  drop = FALSE
]
write.csv(
  surface_comparison,
  file.path(RESULTS_DIR, "time_surface_model_comparison.csv"),
  row.names = FALSE
)

surface_models <- list(
  spline_no_time_interaction = additive_reml,
  spline_ti_interaction = ti_reml,
  full_te_surface = te_reml
)
saveRDS(
  surface_models,
  file.path(DERIVED_DIR, "time_surface_models.rds")
)

anchor_years <- c(2000, 2010, 2014, 2020, 2024, 2025)
anchor_prevalence <- c(5, 10, 20, 30, 50)
anchor_grid <- expand.grid(
  model = names(surface_models),
  year = anchor_years,
  prevalence = anchor_prevalence,
  stringsAsFactors = FALSE
)
anchor_rows <- lapply(seq_len(nrow(anchor_grid)), function(i) {
  row <- anchor_grid[i, , drop = FALSE]
  af <- af_from_model(
    surface_models[[row$model]],
    prevalence = row$prevalence,
    year_c = row$year - bundle$year_center
  )
  data.frame(
    model = row$model,
    year = row$year,
    prevalence = row$prevalence,
    af = af$af,
    lo = af$lo,
    hi = af$hi
  )
})
anchor_results <- do.call(rbind, anchor_rows)
write.csv(
  anchor_results,
  file.path(RESULTS_DIR, "time_surface_af_anchors.csv"),
  row.names = FALSE
)

prevalence_grid <- seq(1, 70, length.out = 240)
curve_years <- c(2000, 2014, 2024, 2025)
curve_rows <- list()
for (model_name in names(surface_models)) {
  for (year in curve_years) {
    af <- af_from_model(
      surface_models[[model_name]],
      prevalence = prevalence_grid,
      year_c = year - bundle$year_center
    )
    af$model <- model_name
    af$year <- year
    curve_rows[[paste(model_name, year)]] <- af
  }
}
surface_curves <- do.call(rbind, curve_rows)
rownames(surface_curves) <- NULL
surface_labels <- c(
  spline_no_time_interaction = "No interaction: s(PfPR) + s(year)",
  spline_ti_interaction = "ANOVA interaction: s(PfPR) + s(year) + ti",
  full_te_surface = "Full te(PfPR, year) surface"
)
write.csv(
  surface_curves,
  file.path(RESULTS_DIR, "time_surface_af_curves.csv"),
  row.names = FALSE
)

surface_plot <- ggplot2::ggplot(
  surface_curves,
  ggplot2::aes(prevalence, 100 * af, colour = model)
) +
  ggplot2::geom_line(linewidth = 0.95) +
  ggplot2::facet_wrap(~year, ncol = 2) +
  ggplot2::scale_x_log10(
    breaks = c(1, 2, 5, 10, 20, 50, 70)
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      spline_no_time_interaction = "#00A83B",
      spline_ti_interaction = "#5B8FF9",
      full_te_surface = "#F8766D"
    ),
    labels = surface_labels
  ) +
  ggplot2::labs(
    x = expression(italic(Pf) * "PR"[2-10] * " (%), log scale"),
    y = "Population-average attributable fraction\n(%, versus 1% PfPR)",
    colour = "Time structure"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(RESULTS_DIR, "time_surface_af_comparison.png"),
  surface_plot,
  width = 10,
  height = 7.5,
  dpi = 300
)

national_time_series <- merge(
  read.csv(
    file.path(DATA_DIR, "wb_mortality_timeseries.csv"),
    stringsAsFactors = FALSE
  )[, c("iso3", "year", "allcause_1mo5y")],
  read.csv(
    file.path(DATA_DIR, "pfpr_by_country_year.csv"),
    stringsAsFactors = FALSE
  ),
  by = c("iso3", "year")
)
national_time_series$region <- countrycode::countrycode(
  national_time_series$iso3,
  "iso3c",
  "region",
  warn = FALSE
)
national_time_series <- national_time_series[
  national_time_series$region == "Sub-Saharan Africa" &
    is.finite(national_time_series$allcause_1mo5y) &
    is.finite(national_time_series$pfpr_pct),
  ,
  drop = FALSE
]

burden_time_rows <- lapply(names(surface_models), function(model_name) {
  model <- surface_models[[model_name]]
  af <- af_from_model(
    model,
    prevalence = national_time_series$pfpr_pct,
    year_c = national_time_series$year - bundle$year_center
  )$af
  country_rows <- data.frame(
    model = model_name,
    iso3 = national_time_series$iso3,
    year = national_time_series$year,
    pfpr_pct = national_time_series$pfpr_pct,
    allcause_postneonatal = national_time_series$allcause_1mo5y,
    attributable_fraction = af,
    malaria_deaths = af * national_time_series$allcause_1mo5y
  )
  country_rows
})
burden_time_country <- do.call(rbind, burden_time_rows)
burden_time_total <- aggregate(
  malaria_deaths ~ model + year,
  burden_time_country,
  sum
)
write.csv(
  burden_time_total,
  file.path(RESULTS_DIR, "time_surface_burden_timeseries.csv"),
  row.names = FALSE
)

burden_time_plot <- ggplot2::ggplot(
  burden_time_total,
  ggplot2::aes(year, malaria_deaths / 1000, colour = model)
) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_colour_manual(
    values = c(
      spline_no_time_interaction = "#00A83B",
      spline_ti_interaction = "#5B8FF9",
      full_te_surface = "#F8766D"
    ),
    labels = surface_labels
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Model-attributable post-neonatal deaths (thousands)",
    colour = "Time structure"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(RESULTS_DIR, "time_surface_burden_timeseries.png"),
  burden_time_plot,
  width = 10,
  height = 6,
  dpi = 300
)

read_who_country_estimates <- function() {
  output <- file.path(
    RESULTS_DIR,
    "who_wmr2025_country_malaria_deaths_2024.csv"
  )
  if (file.exists(output)) {
    cached <- read.csv(output, stringsAsFactors = FALSE)
    if ("who_region" %in% names(cached)) return(cached)
  }

  api_url <- paste0(
    "https://ghoapi.azureedge.net/api/MALARIA_EST_DEATHS",
    "?$filter=TimeDim%20eq%202024&$format=json"
  )
  json_path <- Sys.getenv("WHO_JSON_PATH", unset = "")
  if (nzchar(json_path) && file.exists(json_path)) {
    payload <- jsonlite::fromJSON(json_path)
  } else {
    message("Fetching WHO country malaria-death estimates from the GHO API.")
    response <- httr::GET(api_url, httr::timeout(180))
    httr::stop_for_status(response)
    payload <- jsonlite::fromJSON(
      httr::content(response, "text", encoding = "UTF-8")
    )
  }
  x <- payload$value
  x <- x[
    x$SpatialDimType == "COUNTRY" &
      x$TimeDim == 2024 &
      is.finite(x$NumericValue),
    ,
    drop = FALSE
  ]
  out <- data.frame(
    iso3 = x$SpatialDim,
    who_region = x$ParentLocationCode,
    who_all_age_2024 = x$NumericValue,
    who_all_age_lo_2024 = x$Low,
    who_all_age_hi_2024 = x$High,
    source_url = api_url
  )
  write.csv(out, output, row.names = FALSE)
  out
}

pfpr_latest <- read.csv(
  file.path(DATA_DIR, "pfpr_by_country_2024.csv"),
  stringsAsFactors = FALSE
)
mortality_latest <- read.csv(
  file.path(DATA_DIR, "igme_mortality_by_country.csv"),
  stringsAsFactors = FALSE
)
births_latest <- read.csv(
  file.path(DATA_DIR, "wb_livebirths_by_country.csv"),
  stringsAsFactors = FALSE
)[, c("iso3", "births")]

latest <- merge(
  pfpr_latest,
  mortality_latest[, c("iso3", "u5mr_year", "m_1mo_5y")],
  by = "iso3"
)
latest <- merge(latest, births_latest, by = "iso3")
latest$region <- countrycode::countrycode(
  latest$iso3,
  "iso3c",
  "region",
  warn = FALSE
)
latest$country <- countrycode::countrycode(
  latest$iso3,
  "iso3c",
  "country.name",
  warn = FALSE
)
latest <- latest[
  latest$region == "Sub-Saharan Africa" &
    is.finite(latest$pfpr_pct) &
    is.finite(latest$m_1mo_5y) &
    is.finite(latest$births),
  ,
  drop = FALSE
]
latest$allcause_postneonatal_2024 <-
  latest$m_1mo_5y / 1000 * latest$births

for (model_name in names(surface_models)) {
  af <- af_from_model(
    surface_models[[model_name]],
    prevalence = latest$pfpr_pct,
    year_c = 2025 - bundle$year_center
  )$af
  latest[[paste0("af_", model_name, "_2025")]] <- af
  latest[[paste0("deaths_", model_name, "_2025")]] <-
    af * latest$allcause_postneonatal_2024
}

ihme <- read.csv(
  file.path(DATA_DIR, "ihme_malaria_u5_deaths_by_country.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
ihme <- ihme[
  ihme$Year == 2025 &
    ihme$Measure == "Deaths" &
    ihme$Age == "Under 5" &
    ihme$Unit == "Number",
  ,
  drop = FALSE
]
ihme$iso3 <- countrycode::countrycode(
  ihme$Location,
  "country.name",
  "iso3c",
  warn = FALSE
)
ihme <- ihme[!is.na(ihme$iso3), , drop = FALSE]
ihme <- ihme[, c("iso3", "Value", "Lower", "Upper")]
names(ihme)[2:4] <- c(
  "ihme_u5_2025",
  "ihme_u5_lo_2025",
  "ihme_u5_hi_2025"
)

who <- read_who_country_estimates()
country_comparison <- merge(latest, ihme, by = "iso3", all.x = TRUE)
country_comparison <- merge(
  country_comparison,
  who[, setdiff(names(who), "source_url")],
  by = "iso3",
  all.x = TRUE
)
country_comparison$who_u5_proxy_2024 <-
  0.75 * country_comparison$who_all_age_2024
country_comparison$who_u5_proxy_lo_2024 <-
  0.75 * country_comparison$who_all_age_lo_2024
country_comparison$who_u5_proxy_hi_2024 <-
  0.75 * country_comparison$who_all_age_hi_2024

best_surface <- surface_comparison$model[1]
best_column <- switch(
  best_surface,
  spline_no_interaction = "deaths_spline_no_time_interaction_2025",
  spline_time_interaction = "deaths_spline_ti_interaction_2025",
  full_te_surface = "deaths_full_te_surface_2025",
  stop("The AIC-best model is not a spline surface model: ", best_surface)
)
country_comparison$best_model <- best_surface
country_comparison$best_model_deaths <- country_comparison[[best_column]]
country_comparison$ratio_vs_ihme <-
  country_comparison$best_model_deaths /
  country_comparison$ihme_u5_2025
country_comparison$ratio_vs_who_u5_proxy <-
  country_comparison$best_model_deaths /
  country_comparison$who_u5_proxy_2024
country_comparison$absolute_difference_vs_ihme <-
  country_comparison$best_model_deaths -
  country_comparison$ihme_u5_2025
country_comparison$absolute_difference_vs_who_u5_proxy <-
  country_comparison$best_model_deaths -
  country_comparison$who_u5_proxy_2024
country_comparison$very_different_vs_ihme <-
  is.finite(country_comparison$ratio_vs_ihme) &
  is.finite(country_comparison$absolute_difference_vs_ihme) &
  abs(country_comparison$absolute_difference_vs_ihme) >= 1000 &
  (
    country_comparison$ratio_vs_ihme >= 2 |
      country_comparison$ratio_vs_ihme <= 0.5
  )
country_comparison$very_different_vs_who <-
  is.finite(country_comparison$ratio_vs_who_u5_proxy) &
  is.finite(country_comparison$absolute_difference_vs_who_u5_proxy) &
  abs(country_comparison$absolute_difference_vs_who_u5_proxy) >= 1000 &
  (
    country_comparison$ratio_vs_who_u5_proxy >= 2 |
      country_comparison$ratio_vs_who_u5_proxy <= 0.5
  )
country_comparison <- country_comparison[
  order(-abs(country_comparison$absolute_difference_vs_ihme)),
  ,
  drop = FALSE
]
write.csv(
  country_comparison,
  file.path(RESULTS_DIR, "latest_country_burden_comparison.csv"),
  row.names = FALSE
)

simulate_total <- function(model, data, simulations = 4000) {
  high <- newdata_at_mean(
    model,
    data$pfpr_pct / 10,
    year_c = 2025 - bundle$year_center
  )
  low <- newdata_at_mean(
    model,
    rep(AF_REFERENCE / 10, nrow(data)),
    year_c = 2025 - bundle$year_center
  )
  dX <- population_lpmatrix(model, high) -
    population_lpmatrix(model, low)
  draws <- MASS::mvrnorm(
    simulations,
    mu = coef(model),
    Sigma = vcov(model)
  )
  eta <- dX %*% t(draws)
  af <- pmax(1 - exp(-eta), 0)
  totals <- colSums(
    af * data$allcause_postneonatal_2024
  )
  c(
    point = sum(
      pmax(1 - exp(-as.numeric(dX %*% coef(model))), 0) *
        data$allcause_postneonatal_2024
    ),
    lo = unname(quantile(totals, 0.025)),
    hi = unname(quantile(totals, 0.975))
  )
}

model_total_rows <- lapply(names(surface_models), function(model_name) {
  estimate <- simulate_total(surface_models[[model_name]], latest)
  data.frame(
    source = model_name,
    estimate_year = 2025,
    exposure_denominator_year = 2024,
    deaths = estimate["point"],
    lo = estimate["lo"],
    hi = estimate["hi"],
    note = paste(
      "Model surface evaluated at 2025; MAP PfPR, IGME mortality",
      "and live births are latest available 2024 values"
    )
  )
})
model_totals <- do.call(rbind, model_total_rows)
ihme_ssa <- read.csv(
  file.path(REPO_ROOT, "results", "ihme_u5_deaths_ssa_timeseries.csv"),
  stringsAsFactors = FALSE
)
ihme_ssa <- ihme_ssa[ihme_ssa$year == 2025, , drop = FALSE]
who_africa <- read.csv(
  file.path(REPO_ROOT, "results", "who_wmr2025_africa_deaths.csv"),
  stringsAsFactors = FALSE
)
who_africa <- who_africa[who_africa$year == 2024, , drop = FALSE]
reference_totals <- rbind(
  data.frame(
    source = "IHME_GBD_SSA_aggregate_under5",
    estimate_year = 2025,
    exposure_denominator_year = 2025,
    deaths = ihme_ssa$point,
    lo = ihme_ssa$lo,
    hi = ihme_ssa$hi,
    note = "IHME/GBD supplied SSA aggregate; under-five deaths"
  ),
  data.frame(
    source = "IHME_GBD_matched_countries_under5",
    estimate_year = 2025,
    exposure_denominator_year = 2025,
    deaths = sum(country_comparison$ihme_u5_2025, na.rm = TRUE),
    lo = sum(country_comparison$ihme_u5_lo_2025, na.rm = TRUE),
    hi = sum(country_comparison$ihme_u5_hi_2025, na.rm = TRUE),
    note = "Sum of IHME/GBD under-five estimates for the 43 model countries"
  ),
  data.frame(
    source = "WHO_African_region_all_age",
    estimate_year = 2024,
    exposure_denominator_year = 2024,
    deaths = who_africa$point,
    lo = who_africa$lo,
    hi = who_africa$hi,
    note = "WHO WMR 2025 African-region aggregate; all ages, year 2024"
  ),
  data.frame(
    source = "WHO_African_region_under5_proxy_75pct",
    estimate_year = 2024,
    exposure_denominator_year = 2024,
    deaths = 0.75 * who_africa$point,
    lo = 0.75 * who_africa$lo,
    hi = 0.75 * who_africa$hi,
    note = paste(
      "WHO African-region all-age total multiplied by 0.75;",
      "approximate under-five proxy"
    )
  ),
  data.frame(
    source = "WHO_matched_countries_all_age",
    estimate_year = 2024,
    exposure_denominator_year = 2024,
    deaths = sum(country_comparison$who_all_age_2024, na.rm = TRUE),
    lo = sum(country_comparison$who_all_age_lo_2024, na.rm = TRUE),
    hi = sum(country_comparison$who_all_age_hi_2024, na.rm = TRUE),
    note = "Sum of WHO all-age estimates for the 43 model countries"
  ),
  data.frame(
    source = "WHO_matched_countries_under5_proxy_75pct",
    estimate_year = 2024,
    exposure_denominator_year = 2024,
    deaths = sum(country_comparison$who_u5_proxy_2024, na.rm = TRUE),
    lo = sum(country_comparison$who_u5_proxy_lo_2024, na.rm = TRUE),
    hi = sum(country_comparison$who_u5_proxy_hi_2024, na.rm = TRUE),
    note = paste(
      "WHO all-age country estimates multiplied by 0.75;",
      "approximate, not a WHO country-level under-five series"
    )
  )
)
latest_totals <- rbind(model_totals, reference_totals)
write.csv(
  latest_totals,
  file.path(RESULTS_DIR, "latest_burden_totals.csv"),
  row.names = FALSE
)

message("\nTime-surface model comparison:")
print(surface_comparison, row.names = FALSE)
message("\nLatest burden comparison over ", nrow(latest), " matched SSA countries:")
print(latest_totals, row.names = FALSE)
message("\nCountries meeting the 'very different' rule versus IHME:")
print(
  country_comparison[
    country_comparison$very_different_vs_ihme,
    c(
      "iso3", "country", "pfpr_pct", "best_model_deaths",
      "ihme_u5_2025", "ratio_vs_ihme",
      "absolute_difference_vs_ihme"
    )
  ],
  row.names = FALSE
)

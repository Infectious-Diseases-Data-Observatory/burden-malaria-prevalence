# Versioned regional specifications; preserve the 18-variable benchmark.
source("R_cbh/covariates/regional.R")
cbh_primary_regional_spec <- function(settings) {
  if(!isTRUE(settings$nutrition))return(cbh_regional_spec())
  regional <- c("mean_maternal_age_first_birth","mean_maternal_education_years","mean_wealth_quintile",
    "urban_pct","dtp3_pct","measles_pct","facility_delivery_pct","short_birth_interval_pct",
    "improved_water_pct","improved_sanitation_pct","electricity_pct","wasting_pct","stunting_pct")
  annual <- c("log_hiv_incidence","log_gdp_pc","log_health_expenditure_pc","political_stability")
  list(id="regional17_nutrition",regional=regional,annual=annual,covariates=c(regional,annual))
}
cbh_primary_regional_formula <- function(settings) {
  as.formula(paste("death ~ s(pfpr_pct,bs='cr',k=5) + s(calendar_year,bs='cr',k=6) +",
    paste(paste0("z_",cbh_primary_regional_spec(settings)$covariates),collapse=" + "),
    "+ s(survey,bs='re') + s(country,bs='re') + s(region,bs='re') + offset(log(band_years))"))
}

# National covariate imputation (DHS and MICS, imputed_dhsmics_v6)

Countries: the 40 countries with built DHS or MICS surveys; years 2000-2024. Methods follow the DHS-only version (`R_cbh/covariates/15_impute_national.R`); the MICS countries add Somalia, South Sudan, the Central African Republic and Guinea-Bissau.

**Political stability (WGI):** 39 country-years in 2001 (no WGI round) filled by the mean of each country's 2000 and 2002 estimates (deterministic, `wgi_2001_interpolated`). South Sudan has no WGI series before independence: its 11 country-years 2000-2010 are filled from `political_stability ~ s(year, k = 6) + s(iso3, bs = "re") + s(iso3,     year, bs = "re")` fitted by REML to 950 observed country-years (residual SD 0.399), with 10 draws saved.

**GDP per capita (log, current US$):** South Sudan is observed only for 2008-2015; its 17 missing country-years are filled from `log_gdp_pc ~ s(year, k = 6) + s(iso3, bs = "re") + s(iso3, year,     bs = "re")` (983 observed country-years, residual SD 0.232), with draws. Only 2005-2007 enter the analysis (MC_SSD2010 band entry years).

**Health expenditure per capita (log, current US$):** 80 country-years filled from the DHS-only model `log_health_expenditure_pc ~ s(year, k = 6) + log_gdp_pc + s(iso3,     bs = "re") + s(iso3, year, bs = "re")`, fitted to the 913 observed country-years with observed GDP (residual SD 0.230, deviance explained 94.2%): Zimbabwe 2000-2009, Somalia 2000-2012, South Sudan outside 2017-2023 and every country in 2024. South Sudan's observed 2017-2023 values have no GDP and are not used in the fit, so its fills are population-level predictions with the country-effect prior uncertainty. Health-expenditure draw k uses GDP draw k where GDP is imputed.

Fill labels give each fill's position relative to the country's own observed series (`gam_backcast_before_series`, `gam_within_series_gap`, `gam_forecast_after_series`). The Somalia (13 years) and South Sudan (up to 17 years) health-expenditure fills and the South Sudan pre-independence political stability are long extrapolations; they affect MC_SOM2006, MC_SOM2011NE, MC_SOM2011SL and MC_SSD2010 only. Their values are listed in the country table below.

Child HIV incidence for São Tomé and Príncipe (MICS 2014 and 2019) comes from the extended incidence panel (latent adolescent series); the other MICS countries have their own UNAIDS series.

[Imputed panel](national_covariates_imputed.csv) · [Summary](national_imputation_summary.csv) · [Fills by country](national_imputation_by_country.csv) · [Somalia, South Sudan and Zimbabwe series](somalia_south_sudan_zimbabwe_series.csv)

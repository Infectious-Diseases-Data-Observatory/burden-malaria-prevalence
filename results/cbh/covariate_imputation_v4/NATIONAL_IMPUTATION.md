# National covariate imputation (imputed_v4)

Countries: 36 built survey countries; years 2000-2024.

**Political stability (WGI):** 36 country-years missing, all in 2001 (no WGI round). Filled by the mean of each country's 2000 and 2002 estimates; deterministic, flagged `wgi_2001_interpolated`.

**Health expenditure per capita (log, current US$):** 46 country-years missing: Zimbabwe 2000-2009 (11) and 2024 for all 36 countries. Filled from `log_health_expenditure_pc ~ s(year, k = 6) + log_gdp_pc + s(iso3,     bs = "re") + s(iso3, year, bs = "re")` fitted by REML on 854 observed country-years (residual SD 0.225 on the log scale, deviance explained 94.5%). Point values are fitted means; 10 draws from the coefficient posterior plus residual noise are saved for multiple-imputation propagation. Flags: `gam_within_series_model` (Zimbabwe) and `gam_one_year_extrapolation` (2024).

GDP per capita has no gaps and is unchanged. Child HIV incidence for Liberia and Sao Tome and Principe is imputed separately by the extended incidence model (`R_cbh/hiv/01_fit_incidence.R --extended`).

[Imputed panel](national_covariates_imputed.csv) · [Summary](national_imputation_summary.csv) · [Zimbabwe series](zimbabwe_health_expenditure.csv)

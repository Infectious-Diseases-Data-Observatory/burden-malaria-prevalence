# Regional covariate imputation (DHS and MICS, imputed_dhsmics_v6)

Chained-equation multiple imputation (`mice`, predictive mean matching, m = 10, maxit = 20, seed 20260918) at the survey-region level on 1470 survey-regions (1126 DHS, 344 MICS) from 166 surveys, imputed jointly. Targets: the 13 regional covariates and regional MAP PfPR at the survey year; predictors: all of these plus survey year, the national series at survey year (log GDP, log health expenditure, political stability, from `01_impute_national.R`) and a country factor. The specification is that of the DHS-only version (`R_cbh/covariates/16_impute_regional.R`); with the larger donor pool and a different random stream, DHS imputed values differ from that version. Observed values are unchanged. The point overlay is the mean of the 10 imputations; all are kept for propagation.

5 MICS survey-regions without MAP prevalence at the survey year and without analysis records were dropped: MC_MRT2011 dakhletnouadhibou, MC_MRT2015 dakhletnouadhibou, MC_SWZ2010 manzini, MC_SWZ2014 manzini, MC_SWZ2021 manzini.

| Variable | Survey-regions imputed (DHS / MICS) | Surveys | Observed mean | Imputed mean | Between-imputation SD |
|---|---:|---:|---:|---:|---:|
| mean_maternal_age_first_birth | 0 (0 / 0) | 0 | 19.36 | — | — |
| mean_maternal_education_years | 0 (0 / 0) | 0 | 4.43 | — | — |
| mean_wealth_quintile | 66 (66 / 0) | 9 | 2.83 | 2.94 | 0.27 |
| urban_pct | 0 (0 / 0) | 0 | 32.38 | — | — |
| dtp3_pct | 0 (0 / 0) | 0 | 68.72 | — | — |
| measles_pct | 0 (0 / 0) | 0 | 70.37 | — | — |
| facility_delivery_pct | 36 (36 / 0) | 5 | 62.03 | 53.50 | 9.71 |
| short_birth_interval_pct | 0 (0 / 0) | 0 | 18.93 | — | — |
| improved_water_pct | 0 (0 / 0) | 0 | 67.77 | — | — |
| improved_sanitation_pct | 0 (0 / 0) | 0 | 40.94 | — | — |
| electricity_pct | 12 (12 / 0) | 1 | 33.82 | 19.37 | 7.01 |
| wasting_pct | 111 (93 / 18) | 14 | 7.62 | 9.05 | 3.36 |
| stunting_pct | 89 (71 / 18) | 13 | 33.25 | 36.76 | 6.46 |

MICS gaps are whole-survey wasting and stunting gaps (Nigeria 2021, Madagascar 2012 South, Somalia 2011 North-East and Somaliland 2011). Somalia's only other survey with anthropometry (MICS 2006) reports one national wasting and stunting value, so the country factor carries little Somalia-specific information.

[Convergence traces](mice_convergence.png) · [Observed versus imputed densities](mice_density.png) · [By variable](regional_imputation_by_variable.csv) · [By survey](regional_imputation_by_survey.csv)

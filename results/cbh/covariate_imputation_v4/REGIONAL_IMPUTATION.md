# Regional covariate imputation (imputed_v4)

Chained-equation multiple imputation (`mice`, predictive mean matching, m = 10, maxit = 20, seed 20260918) at the survey-region level on 1126 survey-regions from 120 surveys. Targets: the 13 regional covariates and regional MAP PfPR at the survey year; predictors: all of these plus survey year, the national series at survey year (log GDP, log health expenditure, political stability; 2001/2024 gaps already filled by `15_impute_national.R`) and a country factor. Observed values are unchanged. The point overlay is the mean of the 10 imputations; all 10 are kept for propagation.

| Variable | Survey-regions imputed | Surveys | Observed mean | Imputed mean | Between-imputation SD |
|---|---:|---:|---:|---:|---:|
| mean_maternal_age_first_birth | 0 | 0 | 19.32 | — | — |
| mean_maternal_education_years | 0 | 0 | 4.43 | — | — |
| mean_wealth_quintile | 66 | 9 | 2.84 | 2.92 | 0.25 |
| urban_pct | 0 | 0 | 32.22 | — | — |
| dtp3_pct | 0 | 0 | 70.00 | — | — |
| measles_pct | 0 | 0 | 71.84 | — | — |
| facility_delivery_pct | 36 | 5 | 62.03 | 52.84 | 10.05 |
| short_birth_interval_pct | 0 | 0 | 18.73 | — | — |
| improved_water_pct | 0 | 0 | 67.44 | — | — |
| improved_sanitation_pct | 0 | 0 | 40.73 | — | — |
| electricity_pct | 12 | 1 | 33.48 | 17.41 | 7.15 |
| wasting_pct | 93 | 10 | 7.64 | 8.79 | 3.32 |
| stunting_pct | 71 | 9 | 33.78 | 35.03 | 6.28 |

Whole-survey gaps borrow from other surveys of the same country (country factor) and from the covariate relationships; predictive mean matching donates observed values, so imputations stay within observed ranges. Including PfPR as an auxiliary variable follows the usual multiple-imputation advice to include the analysis model's variables; the imputed covariates are adjustment variables, not the exposure or outcome. Imputation uncertainty is not reflected in the point-imputation fit; it is propagated in the separate multiple-imputation check.

[Convergence traces](mice_convergence.png) · [Observed versus imputed densities](mice_density.png) · [By variable](regional_imputation_by_variable.csv) · [By survey](regional_imputation_by_survey.csv)

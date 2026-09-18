# Revised 17-covariate missingness audit

Denominator: 6,357,802 MAP-eligible child–age-band records, 1,113 survey-regions, 120 surveys and 36 countries, before complete-case selection. Percentages below are the unweighted share of records lacking their assigned regional/annual predictor, not missing individual responses. Regional missingness counts treat each survey-region once; national annual measures may be missing for only some entry years.

Before filling uses original regional values and observed numeric child HIV incidence. After filling applies the existing HIV posterior estimates, exact survey-year UNICEF DTP3/measles fallback and arithmetic mean of available regions in the same survey. Censored HIV values lack a usable numeric point estimate before filling and are counted separately in hiv_original_value_status.csv. No new HIV model, cross-survey borrowing, national nutrition fallback or mortality refit is performed.

| Covariate | Records missing before filling | Records missing after filling | Regions entirely missing after filling |
|---|---:|---:|---:|
| Maternal age at first birth | 0.00% | 0.00% | 0 |
| Maternal education (years) | 0.00% | 0.00% | 0 |
| Household wealth-quintile score | 0.56% | 0.56% | 62 |
| Urban residence (%) | 0.00% | 0.00% | 0 |
| DTP3 coverage (%) | 3.75% | 0.00% | 0 |
| Measles coverage (%) | 3.75% | 0.00% | 0 |
| Facility delivery (%) | 3.74% | 3.26% | 36 |
| Birth interval <24 months (%) | 0.75% | 0.00% | 0 |
| Improved drinking water (%) | 0.75% | 0.00% | 0 |
| Improved sanitation (%) | 0.75% | 0.00% | 0 |
| Household electricity (%) | 2.15% | 1.40% | 12 |
| Wasting prevalence (%) | 7.31% | 6.83% | 93 |
| Stunting prevalence (%) | 6.16% | 5.68% | 71 |
| Child HIV incidence (log) | 14.25% | 2.24% | 26 |
| GDP per capita (log) | 0.00% | 0.00% | 0 |
| Health expenditure per capita (log) | 0.99% | 0.99% | 10 |
| Political stability | 2.83% | 2.83% | 0 |

Complete cases after substitution retain **5,465,305 records, 1,686,004 children and 75,726 deaths, in 916 survey-regions, 95 surveys and 34 countries**. 197 survey-regions have no remaining eligible records; 198 additional regions retain only part of their records. Covariate-specific losses overlap and must not be summed.

Against the fitted 18-variable sample (5,680,117 records), the new specification excludes 214,812 records and recovers 0. These are availability results; no 17-variable mortality model has been fitted.

Age at first birth uses v212, valid completed ages 8–49, weighted by v005 among distinct interviewed mothers with a birth in the last 60 months. It is not age at each recent birth. Wasting and stunting use the official WHO-standard DHS indicators CN_NUTS_C_WH2 and CN_NUTS_C_HA2 (below −2 SD, including severe cases). Published regional denominators and reviewed geography rules are preserved; ambiguous or unmatched estimates remain unavailable before the same-survey fallback. Regional nutritional summaries describe measured surviving children at survey, not anthropometry of children who died.

[CSV table](covariate_missingness.csv) · [Complete-case counts](complete_case_summary.csv) · [Regional missingness](missingness_by_survey_region.csv) · [Regional selection](selection_by_survey_region.csv) · [Nutrition selections](nutrition_selected.csv) · [Unmatched nutrition labels](nutrition_unmatched.csv) · [Maternal summary denominators](first_birth_age_summary.csv)

Reproduce: Rscript R_cbh/covariates/11_extract_first_birth_age.R; python3 R_cbh/covariates/12_fetch_nutrition.py; Rscript R_cbh/covariates/13_audit_planned17.R. Public API retrieval is explicit and cached. New audit files are separate from all existing modelling data and results. No TeX files are written.

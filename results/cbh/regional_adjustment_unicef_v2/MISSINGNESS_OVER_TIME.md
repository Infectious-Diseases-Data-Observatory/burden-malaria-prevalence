# Missingness after UNICEF national fallback

Updated 17 September 2026. Denominator: the fixed previous primary sample of 5,885,022 child–age-band records, before applying the expanded regional adjustment requirements. These are percentages of records, not distinct children or survey-regions. Years refer to band entry; the DTP3/measles fallback itself uses country and survey year. This denominator allows direct comparison with the earlier temporal audit and is not the full pool from which the revised model will be selected.

| Band-entry period | Records | Hib3 missing | PCV missing | Rotavirus missing | Any required covariate missing |
|---|---:|---:|---:|---:|---:|
| 2000–2004 | 782,851 | 82.5% | 100.0% | 100.0% | 100.0% |
| 2005–2009 | 1,526,224 | 37.0% | 100.0% | 100.0% | 100.0% |
| 2010–2014 | 1,614,251 | 5.1% | 56.7% | 84.0% | 84.0% |
| 2015–2019 | 1,225,515 | 0.0% | 5.1% | 23.5% | 26.9% |
| 2020–2023 | 736,181 | 0.0% | 0.6% | 10.7% | 14.3% |

DTP3 and measles are now complete in every represented year. National substitution filled 207,937 records for each. Hib3, PCV and rotavirus missingness is unchanged because those inputs already used the same national UNICEF estimates. No unverified pre-introduction zeros or interpolation were applied.

Requiring all 22 covariates retains 1,784,215 of these previous-primary records (30.3%); 927 more than before substitution. The substantial early-period exclusion remains. Across all MAP-eligible records, the separate full availability audit retains 1,784,353 records in 419 survey-regions; that is the candidate revised sample, not a newly fitted mortality model.

Country/survey composition changes by year. No 2001 or 2024 band-entry records occur in this comparison. Complete vaccine data in represented 2022–2023 records do not establish coverage availability for every country.

Reproduce using `Rscript R_cbh/covariates/08_missingness_time_after_imputation.R`. The script reconciles all 22 individual covariate totals and joint exclusions with the before/after availability audits. Detailed outputs: [annual](missingness_by_entry_year.csv), [period](missingness_by_period.csv).

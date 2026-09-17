# Reduced primary regional adjustment — 17 September 2026

The primary adjustment set now has 18 scalar predictors. Hib3, PCV, rotavirus and exclusive breastfeeding are excluded from the formula and complete-case filter. Urban/rural residence is retained as urban_pct: the survey-weighted percentage urban among distinct mothers with a recent birth in each survey-region. It has no missing assigned values in any of the audited samples.

DTP3 and measles remain, with the existing exact country/survey-year UNICEF fallback. Other missing regional summaries use the arithmetic mean of finite available regions in the same survey, independently for each variable. Donor regions are counted once, observed values are preserved, and whole-survey gaps remain missing. The four retained national annual predictors (HIV incidence, GDP, health expenditure and political stability) retain their band-entry-year assignment and existing HIV imputation. No pre-introduction vaccine zero-fill is needed by the reduced specification.

Starting from **6,357,802 MAP-eligible child-band records in 120 surveys**, complete cases retain **5,680,117 records, 1,755,838 children, 78,634 deaths, 973 survey-regions, 100 surveys and 34 countries**. Remaining missingness is **10.66%** of starting records. These are data-availability counts; the revised mortality models have not been fitted.

## Missingness in the previous fitted sample

Denominator: 5,885,022 child-band records and 1,015 survey-regions. Percentages count records lacking their assigned aggregate covariate, not missing individual questionnaire responses.

| Covariate | Records missing | Entirely unavailable survey-regions |
|---|---:|---:|---:|
| Male births | 0.00% | 0 |
| Multiple births | 0.00% | 0 |
| Mean birth order | 0.00% | 0 |
| Mean maternal age at birth | 0.00% | 0 |
| Mean maternal education | 0.00% | 0 |
| Mean wealth-quintile score | 0.00% | 0 |
| Urban residence | 0.00% | 0 |
| DTP3 coverage | 0.00% | 0 |
| Measles coverage | 0.00% | 0 |
| Facility delivery | 3.12% | 30 |
| Birth interval <24 months | 0.00% | 0 |
| Improved drinking water | 0.00% | 0 |
| Improved sanitation | 0.00% | 0 |
| Electricity | 1.42% | 12 |
| Child HIV incidence (fixed imputation) | 0.00% | 0 |
| GDP per capita | 0.00% | 0 |
| Health expenditure | 0.00% | 0 |
| Political stability | 0.00% | 0 |

## Alternative denominators

| Starting sample | Records | Survey-regions | Regions lost | Regions lost (%) |
|---|---:|---:|---:|---:|
| previous_primary | 5,885,022 | 1,015 | 42 | 4.1% |
| eligible_MAP | 6,357,802 | 1,113 | 140 | 12.6% |
| regional_baseline | 5,948,359 | 1,015 | 42 | 4.1% |

Missingness across covariates overlaps. A survey-region is lost only if it retains zero records; partially affected regions are listed separately in the ledger. Source extraction tables can retain excluded variables for historical reproduction, but those variables are absent from the current wide overlay, formula and selection rule.

[Complete-case counts](complete_case_summary.csv), [missingness](missingness_summary.csv), [survey-region ledger](selection_by_survey_region.csv), [regional imputation flags and donor counts](regional_mean_imputation.csv), [formula](planned_formula.txt), [validation](validation.txt).

Current private overlay: `data/derived_cbh/regional_adjustment/reduced_v3/regional_covariates_wide.csv`. Run `Rscript R_cbh/covariates/run.R` to reproduce. Previous 22-variable results are preserved under `regional_adjustment_unicef_v2`; use `--legacy-expanded` to reproduce them.

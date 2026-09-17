# Regional adjustment availability audit — 17 September 2026

**Primary policy: national UNICEF vaccination fallback.** Preserve observed regional DTP3/measles values; substitute the exact country/survey-year WUENIC estimate only when the regional value is missing. The three existing national vaccines use band-entry year. Pre-series/no-series assumptions remain missing unless independently verified; no interpolation, nearest-year carry or automatic zero filling is applied.

| Vaccine | Previously included records imputed | Records still missing |
|---|---:|---:|
| DTP3 coverage | 207,937 | 0 |
| Measles coverage | 207,937 | 0 |
| Hib3 coverage | 0 | 1,292,421 |
| PCV coverage | 0 | 3,291,804 |
| Rotavirus coverage | 0 | 4,032,597 |

This is deterministic national substitution, not a statistical multiple-imputation model. Regional heterogeneity and uncertainty in the substituted estimates are not propagated. Preserve the [pre-imputation audit](../regional_adjustment_v1/README.md) as the comparison. The regional covariate long table describes source estimates before fallback; `regional_vaccine_imputation.csv` and the versioned wide overlay contain the final vaccination values and flags.

**Complete cases retain 419 of the previous 1,015 survey-regions: 596 are lost (58.7%).** A further 217 retained regions lose some records. Region counts distinguish surveys/boundary versions.

Starting from all eligible MAP records, the revised regional adjustment set yields **1,784,353 child-band records, 598,386 distinct children and 17,822 deaths, from 41 surveys in 24 countries and 419 survey-regions**. This is an availability result, not a fitted model.

Within the previous fitted sample, 1,784,215 records remain and 69.68% are excluded. Starting before its individual-level complete-case filter recovers 138 additional records (48 children) using regional means.

## Remaining missingness after UNICEF substitution

Denominator: **5,885,022 child-band records and 1,015 survey-regions**. Percentages are the unweighted share of mortality records lacking the assigned regional or annual covariate, not the percentage of respondents with an unanswered question. A regional mean may remain available despite some missing responses. The final column counts regions with no usable value for that covariate in any of their previously included records. Annual covariates can also remove part of a region's records.

| Covariate | Records missing | Entirely unavailable survey-regions |
|---|---:|---:|
| Male births | 0.00% | 0 |
| Multiple births | 0.00% | 0 |
| Mean birth order | 0.00% | 0 |
| Mean maternal age at birth | 0.00% | 0 |
| Mean maternal education | 0.00% | 0 |
| Mean wealth-quintile score | 0.00% | 0 |
| Urban residence | 0.00% | 0 |
| DTP3 coverage | 0.00% | 0 |
| Measles coverage | 0.00% | 0 |
| Facility delivery | 3.52% | 36 |
| Exclusive breastfeeding | 15.63% | 230 |
| Birth interval <24 months | 0.69% | 9 |
| Improved drinking water | 0.69% | 9 |
| Improved sanitation | 0.69% | 9 |
| Electricity | 2.11% | 21 |
| Child HIV incidence (fixed imputation) | 0.00% | 0 |
| GDP per capita | 0.00% | 0 |
| Health expenditure | 0.00% | 0 |
| Political stability | 0.00% | 0 |
| Hib3 coverage | 21.96% | 132 |
| PCV coverage | 55.94% | 466 |
| Rotavirus coverage | 68.52% | 563 |

Losses overlap; neither percentages nor regional losses should be added across covariates. All eleven earlier adjustment concepts have 0% missingness within the previous fitted sample by its original selection; the regional versions also have complete coverage on that denominator.

## Alternative denominators

| Starting sample | Records | Survey-regions | Regions lost | Regions lost (%) |
|---|---:|---:|---:|---:|
| previous_primary | 5,885,022 | 1,015 | 596 | 58.7% |
| eligible_MAP | 6,357,802 | 1,113 | 694 | 62.4% |
| regional_baseline | 5,948,359 | 1,015 | 596 | 58.7% |

`eligible_MAP` is every valid complete-band record with geography and exact entry-year MAP exposure. `regional_baseline` requires the regional versions of the previous adjustment concepts and the four existing annual predictors, before adding vaccines/delivery/feeding/interval/WASH/electricity. These definitions explain why regional averaging can increase initial record availability without adding a region to the former fitted sample.

## Extraction findings and unresolved issues

- All five vaccine measures are selected. National Hib3/PCV/rotavirus missingness is predominantly pre-series placeholders previously labelled assumed not introduced, plus countries with no series. These are not independently verified zero-coverage estimates. They remain missing under the existing pipeline rule. Resolving introduction dates/source coverage could substantially change these losses.
- Published DHS regional vaccination, delivery, water, sanitation and electricity replace fragile birth-recode label matching. Historical boundary alternatives and nested source regions are resolved explicitly; individual household coverage is not estimated using women's weights.
- Birth interval is a regional percentage among non-first births. First births are not excluded from the mortality model because B11 is structurally unavailable.
- Breastfeeding remains unavailable for 230 previously included regions: 159 have no usable recode estimate and 71 have recode estimates held out for questionnaire review. Published regional estimates take precedence where available. Ten surveys trigger the >5-percentage-point national discrepancy screen; their source labels, active items and comparison values are in `breastfeeding_validation.csv`. This is a conservative implementation screen, not a calibration or proof that the remaining questionnaires are equivalent.
- Among retained regions, 3 breastfeeding estimates have fewer than 25 observations and 60 have fewer than 50. They remain included in this availability audit; small-denominator uncertainty and a declared sensitivity need assessment before final inference.
- Some absent regional values still need source/geography review. The full mapping ledger includes unmatched alternative subregions outside the model geography; these are not all lost analysis regions. Use `selection_by_survey_region.csv` and the wide covariate table to identify the actual losses.
- Facility-delivery recall can be five, three or two years. DHS covariates are survey-time summaries; national annual covariates retain band-entry-year assignment. Regional averages change the adjustment interpretation and do not control individual-level confounding in the same way as the earlier model.

## Files and reproduction

- [Complete-case summary](complete_case_summary.csv), [missingness by variable](missingness_summary.csv), [survey-region ledger](selection_by_survey_region.csv), [country losses](selection_by_country.csv), and [age/survey selection](selection_by_survey_age.csv).
- [Regional values and denominators](regional_covariates.csv), [breastfeeding validation](breastfeeding_validation.csv), [published regional mapping](published_region_mapping.csv), [source choices](published_selections.csv), and [vaccine source statuses](vaccine_source_status.csv).
- [Pipeline and definitions](../../../R_cbh/covariates/README.md). The wide overlay is `data/derived_cbh/regional_adjustment/regional_covariates_wide.csv`. Raw recodes and respondent identifiers remain under ignored `data/`.
- The audit reproduces all old sample totals exactly and checks unique joins, equal denominators for interval components, valid percentages and non-overlapping complete-case counts. No mortality fits, existing figures or TeX files were changed.

Sources: [DHS API](https://api.dhsprogram.com/), the [DHS reference feeding code](https://github.com/DHSProgram/DHS-Indicators-R/blob/main/Chap11_NT/NT_IYCF.R), and [DHS youngest-child selection](https://github.com/DHSProgram/DHS-Indicators-R/blob/main/Chap11_NT/!NTmain.R). Public source URLs and hashes are retained in the data snapshot's `source_manifest.csv`; local inputs are fingerprinted in `audit_input_manifest.csv`.

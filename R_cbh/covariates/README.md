# Expanded regional adjustment

The 17 September 2026 primary specification uses **22 regional/annual scalar confounders**. All DHS confounders are survey-region means or percentages, not individual predictors. National annual inputs remain constant across regions within a country-year. The child-band outcome, MAP exposure at entry, time spline, gamma=2 and seven separate age-band models are unchanged.

This stage builds and audits a covariate overlay. It does **not** refit mortality models or modify existing base shards, fitted objects, manuscript figures or TeX. `run_all.R --primary` still reproduces the preceding individual-adjustment benchmark. A new fitting version and prepared dataset are needed before revised results can be promoted.

**Current primary policy:** retain usable regional vaccination estimates and replace missing DTP3/measles coverage with exact country/survey-year UNICEF WUENIC estimates. Hib3/PCV/rotavirus continue to use national band-entry-year estimates. The current overlay/audit version is `regional_adjustment_unicef_v2`; the preceding no-substitution v1 audit remains preserved. No absent year is automatically converted to zero.

Run the local pipeline from the project root with `Rscript R_cbh/covariates/run.R`. It extracts/caches recode summaries, reads the existing UNICEF snapshot, applies national substitution and regenerates availability reports and validation. Public DHS retrieval remains a separate explicit command. `--complete-case-only` reproduces the earlier no-substitution policy.

Individual stages:

```sh
# Explicit online retrieval of public aggregate data; cached pages reused.
Rscript R_cbh/covariates/01_fetch_published.R
# Add --refresh only for an intentional new source snapshot.

# Offline weighted regional extraction from authorized local birth recodes:
Rscript R_cbh/covariates/02_extract_regional.R
# Extract five vaccine series from the existing UNICEF snapshot:
Rscript R_cbh/covariates/07_prepare_unicef.R
# Offline assembly, national fallback, join validation and availability:
Rscript R_cbh/covariates/03_audit_missingness.R
Rscript R_cbh/covariates/04_report.R
Rscript R_cbh/covariates/05_validate.R
# Before/after missingness by band-entry year on the fixed previous sample:
Rscript R_cbh/covariates/08_missingness_time_after_imputation.R
```

No packages are installed. The data retrieval sends only public survey/indicator IDs. It retains original API records, indicator definitions, denominators, recall windows, geographic labels, URLs and hashes under `data/derived_cbh/regional_adjustment/`. The DHS API's `DataId` is not globally unique; semantic/source keys are checked instead. Published estimates can be repeated at nested geographic levels; equal values with equal denominators can be deduplicated, conflicting values remain missing unless a reviewed source-version rule resolves them.

## Definitions and joins

See [analysis plan Section 2.4](../../docs/ANALYSIS_PLAN.md) for the complete covariate dictionary. `regional.R` defines the selected variables and the new formula through `cbh_regional_spec()` and `cbh_regional_formula()`.

- Sex, multiple birth, birth order and maternal age: women's-weighted summaries over live births in the 60 months before interview, counting births once. This population is independent of survival, complete-band eligibility and exposure availability.
- Education, mean wealth-quintile score and urban residence: women's-weighted summaries over distinct interviewed mothers with a birth in that interval. These are **not all-household averages**. Wealth is an ordinal 1–5 mean, replacing the former individual categorical predictor.
- Regional DTP3/measles, delivery, WASH and electricity: published DHS estimates. Household denominators/weights come from the published source, avoiding the legacy practice of applying women's weights to deduplicated BR households. Delivery uses five-year estimates where available, otherwise an explicitly flagged three-/two-year window. The vaccination age window follows the published indicator, normally 12–23 months; source-specific windows require checking when comparing releases.
- Short birth interval: sum published 7–17- and 18–23-month percentages only when their denominators agree. The denominator excludes first births; first births remain eligible mortality observations.
- Exclusive breastfeeding: prefer published regional `CN_BFSS_C_EBF`; otherwise compute a survey-weighted estimate among youngest children aged 0–5 months living with their mother. Use B19 where available, otherwise V008−B3, youngest-child selection by BIDX, and child-specific M4. Never repeat maternal V404 across children. Require observed feeding domains and distinguish missing/don't-know answers from “no”. Survey-wide all-NA placeholders are not active questions. The implemented conservative screen flags a survey for questionnaire review when its recode national estimate differs from published `CN_IYCB_C_EXB` by more than 5 percentage points; it does not calibrate estimates to match. Published regional estimates remain usable in those surveys. A difference below this threshold is **not proof of full questionnaire equivalence**. Conditional skips and subsampling still need review where flagged, and unsupported estimates remain unavailable in this audit.
- HIV incidence, GDP, health expenditure, political stability, Hib3, PCV and rotavirus: retain exact country/band-entry-year assignment from the base pipeline (with the same fixed child-incidence imputation). Assumed pre-series/no-series vaccine zeros remain missing; verified observed coverage of zero is valid.

Published source names retain exclusion/inclusion clauses. Only date qualifiers are normalized, and `published_region_versions.csv` explicitly selects old/new geographies for DRC, Ghana, Senegal, Sierra Leone and Uganda. `published_region_overrides.csv` documents capital-city inclusion rules. Reviewed aliases are explicit in `published_region_aliases.csv`. No fuzzy matching or averaging of overlapping regions is used. After source assembly, missing DTP3/measles values receive the declared national UNICEF fallback; this is the only newly authorized regional substitution. A directly reported broad-region estimate takes precedence over fine units mapped to that broad region.

The current wide overlay is `data/derived_cbh/regional_adjustment/unicef_v2/regional_covariates_wide.csv` and has one row per `(survey, regkey)`. It includes `*_before_imputation`, `*_imputed`, `*_source_year` and `*_imputation_source` columns for DTP3/measles. `cbh_unicef_fill()` in `unicef.R` performs exact-year substitution and preserves observed values; use survey year for these regional measures and entry year for Hib3/PCV/rotavirus. The versioned manifest records the policy, source hashes and base dataset manifest. National fallback uncertainty is not propagated; retain an observed-regional-only sensitivity.

The wide overlay has one row per `(survey, regkey)`. Join it to **all eligible base shards**, then attach child HIV incidence and apply complete-case selection on `cbh_regional_spec()$covariates`. Do not start by restricting to the previous complete-case sample: regional means can recover records with missing individual answers. All 22 predictors are then centered/scaled using newly saved preprocessing. No child-level confounder terms belong in the new formula.

## Audit outputs and limits

`08_missingness_time_after_imputation.R` writes `missingness_by_entry_year.csv` and `missingness_by_period.csv` in the current results directory. It compares all 22 covariates and joint missingness before/after UNICEF substitution on the same previous-primary sample, reconciling totals with both availability audits. Percentages count child–age-band records, not distinct children or regions; national vaccine exposure years are band-entry years, while DTP3/measles substitution uses survey year. Years without eligible records are absent, not assigned zero missingness. Country composition changes over time.

[Current results](../../results/cbh/regional_adjustment_unicef_v2/README.md) ([pre-substitution audit](../../results/cbh/regional_adjustment_v1/README.md)) include three distinct denominators: the previous primary sample, all MAP-eligible records, and the regional version of the previous adjustment set before the added covariates. “Lost region” means **zero** surviving child-band rows, while a partially reduced region retains at least one. Region identifiers are survey-specific, not unique areas across time. Multiple missing covariates overlap; do not sum their loss counts.

`recode_regional.csv` retains eligible/observed response counts, weighted denominators and flags for n<25/n<50. No automatic small-denominator exclusion was added. Published suppression may instead appear as an absent source estimate. The availability audit counts a mean as observed if its chosen source is usable; it does not assume every contributing respondent had a complete answer. It checks that the old selection reproduces 5,885,022 records, 82,415 deaths, 1,817,912 children and 1,015 survey-regions, without exporting respondent identifiers.

Before refitting, resolve or explicitly retain the questionnaire/geography exclusions, review source-window comparability, and assess the temporal/country selection caused by national vaccine gaps. The regional averages also have sampling error; their current point estimates do not propagate that uncertainty.

Reference definitions and feeding selection were checked against the [DHS API](https://api.dhsprogram.com/) and the DHS Program's [IYCF R code](https://github.com/DHSProgram/DHS-Indicators-R/blob/main/Chap11_NT/NT_IYCF.R) and [denominator selection](https://github.com/DHSProgram/DHS-Indicators-R/blob/main/Chap11_NT/!NTmain.R). Source snapshots are dated 17 September 2026.

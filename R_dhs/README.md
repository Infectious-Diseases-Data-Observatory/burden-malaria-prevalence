# Rebuilt DHS/MIS analysis

This folder is a parallel, DHS/MIS-only replacement for the evolving scripts in
`R/`. The legacy scripts must remain in place until the rebuilt outputs have
been reviewed and accepted.

See `MANUSCRIPT_ALIGNMENT.md` for the DHS/MIS statements and numerical results
that should be updated in the current manuscript after the rebuild is accepted.

## Analytical contract

- Unit: DHS/MIS survey-region-year.
- Survey years: 2000–2025; MAP is extracted for each survey's own year. A
  missing 2025 MAP raster is flagged and is never silently replaced with 2024.
- Outcomes: all-under-5 mortality, neonatal mortality, and post-neonatal
  mortality, with the latter as the primary outcome.
- Exposure: annual population-weighted MAP PfPR2-10 for the matching survey
  region and year.
- Primary row flag: regional PfPR2-10 at least 1%, positive outcomes/exposure,
  and country mean PfPR2-10 greater than 1%. Both component flags are retained.
- Covariates: DHS survey-region, exact-year national WUENIC vaccine coverage
  and nearest-year national predictors.
- Direct vaccine candidates: regional pentavalent-dose-3, PCV-dose-3 and
  completion of the survey-specific rotavirus schedule. See
  `COVARIATE_SOURCES.md` for their recodes, availability and recommended
  longitudinal sources.
- National vaccine candidates: exact-year WUENIC Hib3, PCV-completion and
  final-dose rotavirus coverage extracted from the local UNICEF global
  dataflow file. Structural pre-introduction years are set to zero and
  provenance flags are retained.
- Missingness: variables with at most 5% missingness receive single imputation
  by country median with overall-median fallback; variables above 5% are flagged
  and excluded.
- Covariate adjustment: proportion variables are logit transformed, continuous
  variables remain on their declared scale, all are standardised, and the full
  eligible block is ridge penalised.
- Model comparison: linear versus spline PfPR, each with versus without a
  PfPR-by-calendar-time interaction, plus a full `te(PfPR, year)` surface. All
  include country random intercepts, random linear PfPR slopes and a
  log-exposure offset. The four decomposed specifications include a separate
  smooth calendar-year term; the full tensor surface contains the year main
  effect.
- Model selection: ML/AIC on one shared sample, followed by REML refitting.
- Outcomes: the five-model AIC comparison is run separately for post-neonatal
  and neonatal mortality. The post-neonatal winner remains the prespecified
  structure for the primary negative-control comparison; the independently
  selected neonatal winner is retained as a diagnostic.

## Script order

1. `00_config.R` — configuration and shared functions.
2. `01_access_dhs_data.R` — survey registry and optional authorised downloads.
3. `02_extract_map_pfpr.R` — annual MAP PfPR by survey-region.
4. `02b_extract_unicef_immunisation.R` — compact country-year WUENIC Hib3,
   PCV-completion and rotavirus-completion panel.
5. `03_build_analysis_dataset.R` — outcomes, covariates, missingness and flags.
6. `04_fit_main_models.R` — five-model post-neonatal and neonatal comparisons,
   outcome-specific AIC table and outcome refits.
7. `05_make_main_plots.R` — figures from saved models only, including a forest
   plot of ridge-standardized conditional covariate associations.
8. `06_sensitivity_samples.R` — prevalence and complete-case samples.
9. `07_sensitivity_model_structure.R` — random slopes, covariate levels and
   spline basis dimensions.
10. `08_sensitivity_likelihood.R` — log-Gaussian rate sensitivity retaining the
   ridge block and primary random-effects/PfPR structure.
11. `09_validate_reproduction.R` — aggregate migration checks against the saved
   legacy headline models and results.
12. `10_time_surface_and_burden.R` — comparison of additive, `ti`, and full
   `te` time/prevalence surfaces, followed by latest-year national burden
   comparisons with IHME and WHO. It is downstream of the DHS/MIS-only model
   rebuild but produces the burden/temporal/IHME-WHO numbers, so `run_all.R`
   runs it whenever its national-burden inputs are present (and it can also be
   run standalone).
13. `11_study_flow.R` — analysis-flow diagram (survey-inclusion funnel plus the
   external MAP/UNICEF/UNAIDS/World Bank data inflows). Self-contained.
14. `12_triangulation_rct.R` — RCT triangulation of the prevalence→mortality
   link: predicts each ITN/curtain trial's mortality reduction from its
   age-standardised prevalence drop using the fitted model bundle, and compares
   with the observed trial effects (manuscript figure fig3_rct).
15. `13_ihme_share.R` — Method-1 comparison figure: IHME/GBD malaria share of
   post-neonatal deaths versus national PfPR2-10 (manuscript figure
   fig_ihme_share_spline). Reads the Method-1 country panel from
   `R/02_component1_country.R`; skips itself if that panel is absent.

## Commands

```sh
# Full rebuild. Does not download DHS recodes unless script 01 is separately
# invoked with --download, but MAP extraction may request uncached annual data.
Rscript R_dhs/run_all.R

# Inputs from scripts 01-02 already exist.
Rscript R_dhs/run_all.R --analysis-only

# Migration check using only the existing aggregate panel.
Rscript R_dhs/run_all.R --from-legacy-aggregate
```

Raw DHS records remain under `data/` and must never be committed. Aggregate model
outputs and plots are written under `results/dhs_rebuild/`. Nothing is copied to
Overleaf automatically.

The UNICEF source file
`data/fusion_GLOBAL_DATAFLOW_UNICEF_1.0_all.csv` and its compact derived panel
remain ignored under `data/`. `run_all.R` refreshes the compact panel when the
source file is newer. `COVARIATE_SOURCES.md` documents the implemented vaccine
and paediatric-HIV variables and candidate SMC sources.

Child (0-14) HIV prevalence is now included. It is derived in
`03_build_analysis_dataset.R` from the UNAIDS 2025 estimates workbook
`data/HIV_Epidemiology_Children_Adolescents_2025.xlsx` (estimated number of
children 0-14 living with HIV) divided by the World Bank 0-14 population
(`SP.POP.0014.TO`), joined by ISO3 and exact survey year and entered on the log
scale as `log_hiv_prev`. Nigeria and Comoros publish only a 15-19 UNAIDS series,
so they stay missing and fall under the standard <=5% imputation rule with a
status flag. SMC is still not included; see `COVARIATE_SOURCES.md`.

When the two saved legacy headline artifacts are present, `run_all.R` also writes
`reproduction_check.csv` and `reproduction_summary.txt`. These test substantive
reproduction using sample counts, outcome effect direction/significance,
ridge-linear effect estimates, and nonlinear attributable-fraction anchors.

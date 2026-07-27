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
- Covariates: DHS survey-region and nearest-year national predictors.
- Missingness: variables with at most 5% missingness receive single imputation
  by country median with overall-median fallback; variables above 5% are flagged
  and excluded.
- Covariate adjustment: proportion variables are logit transformed, continuous
  variables remain on their declared scale, all are standardised, and the full
  eligible block is ridge penalised.
- Model comparison: linear versus spline PfPR, each with versus without a
  PfPR-by-calendar-time interaction. All include country random intercepts and
  random linear PfPR slopes, a smooth calendar-year term, and a log-exposure
  offset.
- Model selection: ML/AIC on one shared sample, followed by REML refitting.

## Script order

1. `00_config.R` — configuration and shared functions.
2. `01_access_dhs_data.R` — survey registry and optional authorised downloads.
3. `02_extract_map_pfpr.R` — annual MAP PfPR by survey-region.
4. `03_build_analysis_dataset.R` — outcomes, covariates, missingness and flags.
5. `04_fit_main_models.R` — 2×2 model comparison and outcome refits.
6. `05_make_main_plots.R` — figures from saved models only.
7. `06_sensitivity_samples.R` — prevalence and complete-case samples.
8. `07_sensitivity_model_structure.R` — random slopes, covariate levels and
   spline basis dimensions.
9. `08_sensitivity_likelihood.R` — log-Gaussian rate sensitivity retaining the
   ridge block and primary random-effects/PfPR structure.
10. `09_validate_reproduction.R` — aggregate migration checks against the saved
   legacy headline models and results.

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

When the two saved legacy headline artifacts are present, `run_all.R` also writes
`reproduction_check.csv` and `reproduction_summary.txt`. These test substantive
reproduction using sample counts, outcome effect direction/significance,
ridge-linear effect estimates, and nonlinear attributable-fraction anchors.

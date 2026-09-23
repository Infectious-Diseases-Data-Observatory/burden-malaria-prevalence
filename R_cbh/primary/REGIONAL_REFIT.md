# Current primary analysis: regional adjustment, DHS and MICS surveys

Updated 23 September 2026. Run from the project root:

```sh
Rscript run_all.R --primary
# The same run called directly (fresh fits, default version regional_mics):
Rscript R_cbh/primary/run_regional.R
# Resume with verified fit caches:
Rscript R_cbh/primary/run_regional.R --resume
# Effects, diagnostics, comparison, figures and tables from saved fits:
Rscript R_cbh/primary/run_regional.R --report-only
```

The current output version is `results/cbh/primary_map_regional17_dhsmics_gamma2_v5/`. It combines DHS and UNICEF MICS
birth histories. Its [report](../../results/cbh/primary_map_regional17_dhsmics_gamma2_v5/REPORT.md) compares the seven
MAP gamma=2 fits with the DHS-only primary `results/cbh/primary_map_regional17_gamma2_v3/`, which is preserved unchanged.

| Sample | DHS part (= v3) | Added by MICS | v5 |
|---|---:|---:|---:|
| Child-band records | 5,465,305 | 2,033,154 | 7,498,459 |
| Children | 1,686,004 | 606,085 | 2,292,089 |
| Deaths | 75,726 | 26,556 | 102,282 |
| Survey-regions | 916 | 295 | 1,211 |
| Surveys | 95 DHS | 37 MICS | 132 |
| Countries | 34 | 2 (CAF, GNB) | 36 |

None of the five MIS surveys passes complete-case selection: they do not publish the facility-delivery and
anthropometry indicators. All seven v5 fits converged at the default tolerance, with full rank, finite covariance and
positive smoothing-Hessian eigenvalues, so no tighter-tolerance restart was selected (`fit_manifest.csv`,
`fit_diagnostics.csv`).

The 17-variable model uses 13 survey-region summaries and four national annual covariates, as defined by
`cbh_primary_regional_spec()` in `primary/specification.R`. Sex, multiple births and birth order are removed; maternal age
at first birth replaces age at each birth, and regional wasting and stunting are included. Urban percentage is included.
Hib3, PCV, rotavirus and exclusive breastfeeding are excluded from both the model and its complete-case requirements. For
MICS surveys all 13 regional summaries are computed from the microdata with the DHS definitions and the same fallbacks
(`R_mics/09_regional_covariates.R`). MICS education is converted from level and grade to years, and MICS facility delivery
refers to the last birth in the two years before interview. Guinea 2016, Comoros 2022 and Chad 2019 use the national
WUENIC DTP3 and measles estimates because their recall doses are unusable.

## Prerequisites

The runner does not rebuild any of these inputs; the preparation and reporting stages record their hashes.

1. DHS child-band shards: `Rscript R_cbh/01_make_analysis_data.R` writes `data/derived_cbh/manifest.rds` and
   `child_bands/` ([dataset pipeline](../README.md)).
2. Fixed child HIV incidence panel: `R_cbh/hiv/01_fit_incidence.R` writes
   `data/derived_cbh/hiv_incidence/child_incidence_country_year.csv`.
3. DHS regional overlay `planned17_audit`: the `R_cbh/covariates/` stages ([covariates](../covariates/README.md)) write
   `data/derived_cbh/regional_adjustment/planned17_audit/regional_covariates_wide.csv` and the audit files in
   `results/cbh/planned17_covariate_missingness/`. Their audits (03, 13, 14) read the HIV panel of item 2, so run it first.
4. MICS build: `Rscript run_all.R --mics-build` runs `R_mics/00a_extract.R` and then `00`–`11`
   ([MICS pipeline](../../R_mics/README.md)). It needs the MICS downloads in `data/MICS_Datasets/` and the admin-1
   shapefile named in `R_mics/config/paths.R` (steps 04–06); no step uses the network. The primary reads:
   - `data/derived_mics/cbh_build/manifest.rds` and `child_bands/`: MICS shards written by `R_mics/08_build_child_bands.R`
     with the unchanged builder;
   - `data/derived_mics/regional_covariates_wide_mics.csv`: MICS regional overlay from `R_mics/09_regional_covariates.R`;
   - `results/mics_inventory/exclusion_attribution_mics_by_survey.csv` (`R_mics/10`) and `survey_inventory.csv`
     (`R_mics/01`), read by the study flow;
   - `data/derived_mics/mics_region_centroids.csv` (`R_mics/11`), read by the Sahel subgroup refit.

Rerunning `R_cbh/01_make_analysis_data.R` or `R_mics/08_build_child_bands.R` rewrites that build's `manifest.rds`, with
new build times, even when every shard is reused from cache. The manifest's hash is part of the v5 preparation signature,
so the next run re-prepares v5, the fit caches no longer match, and the subgroup, no-nutrition and v6 sensitivities must
be refitted as well (`Rscript run_all.R --primary`, then `--sensitivities`). A changed MICS overlay has the same effect.
`run_all.R --mics-build` therefore skips step 08 when the MICS manifest exists, unless `--force` is given. Editing any
`.R` file directly in `R_cbh/` or `R_cbh/R/` changes every shard signature and also forces the rebuild and refits.

## Stages

1. `00_prepare_regional.R` assembles the DHS child-band shards with the audited `planned17_audit` overlay and the fixed
   HIV incidence panel, and stops unless the DHS part reproduces v3 exactly (5,465,305 records, 75,726 deaths, 916
   survey-regions). It then adds the MICS shards listed in the MICS `manifest.rds`, joined to the MICS overlay. A
   model-ready MICS record without an overlay row is an error, not a silent drop. The same complete-case rule applies to
   both programmes. The stage verifies unique child-band keys, rescales the confounders on the combined sample
   (`covariate_scaling.csv`) and caches the modelling dataset with source hashes (`preparation_provenance.csv`,
   `prepared_sample_dhs_mics.csv`, `prepared_selection_by_survey.csv`).
2. `01_fit.R` fits seven separate binomial/cloglog models with `gamma=2`, PfPR `cr/k=5`, calendar time `cr/k=6`,
   reference MAP knots, survey/country/region random intercepts and a full-band-width offset. Each age band has its own
   covariate coefficients, time spline and random-effect variances.
3. `02_effects.R` calculates supported contrasts, zero-PfPR contrasts and 2005/2015/2024 country estimates using unchanged
   national MAP/IHME inputs.
4. `04_diagnostics.R` verifies stored fitted inputs and aggregate outcomes.
5. `05_compare_regional.R` produces the comparison with v3. `03_report.R` owns the primary curve, diagnostic and burden
   plots and the age-band table.
6. The `reporting/` stages combine the DHS and MICS ledgers for the inclusion flow (`01_study_flow.R`) and survey map, and
   produce the three paper figures:
   - Figure 1, survey map and timing (`03_survey_map.R`);
   - Figure 2, attributable fraction by age (`06_attributable_fraction_by_age.R`);
   - Figure 3, the combined burden comparison: 2024 countries, Nigerian states and the 2000–2024 annual trend as panels
     A–C (`11_burden_comparison_figure.R`; the source file keeps its historical name `burden_comparison/fig4_burden_comparison.png`).

   They also write the national, annual and state tables, the source comparison table, the three-figure paper
   manifest (`04_paper_figures.R`), the results index and the validation. `burden/04_annual_comparison.R` and
   `05_nigeria_state_burden.R` recalculate the associated burden from the saved fits. The PfPR curves by age
   (`pfpr_splines.png`), the inclusion flow and the count version of Figure 3 are supplementary figures.

## Version routing and frozen defaults

`run_regional.R` sets `CBH_PRIMARY_VERSION` for its subprocesses:

| Call | Version | Output |
|---|---|---|
| no flag (default; `--mics` alone is equivalent) | `regional_mics` | `primary_map_regional17_dhsmics_gamma2_v5` |
| `--dhs-only` | `regional` | DHS-only benchmark `primary_map_regional17_gamma2_v3` |
| `--imputed` | `regional_imputed` | DHS-only imputed-covariate history `primary_map_regional17_imputed_gamma2_v4`; model stages only |

The reporting stages, `primary/03_report.R`, `05_compare_regional.R`, `burden/04_annual_comparison.R` and
`05_nigeria_state_burden.R` default to `regional_mics` when run directly. Four entry points keep older defaults:
`00_prepare_regional.R` defaults to `regional` (v3), and `01_fit.R`, `02_effects.R` and `04_diagnostics.R` call
`cbh_primary_settings()`, whose default in `settings.R` is `legacy` (`primary_map_gamma2_v1`). `02_effects.R` and
`04_diagnostics.R` stop when the variable is unset; `00_prepare_regional.R` and `01_fit.R` do not. Run these four
directly only with `CBH_PRIMARY_VERSION=regional_mics` set, or use the runner.

These defaults are frozen deliberately. v5's `code_provenance.csv` hashes `primary/01_fit.R`, `settings.R`,
`reference_knots.csv`, `specification.R`, `00_prepare_regional.R`, `analysis/model.R`, `sensitivity/model.R` and
`covariates/regional.R`, and the fit-cache signatures in `01_fit.R` include the same files. The subgroup and no-nutrition
refits also record `settings.R` and `specification.R` in their `fit_input_provenance.csv`. Editing any of these files
invalidates the v5 fit caches, makes `reporting/10_validate_reporting.R` report stale provenance, and requires fresh v5
and sensitivity fits.

Neither history run executes the results-index stage, so `Key results/README.md` keeps describing v5.

## Sensitivity analyses on the v5 sample

Each runs after the v5 fits exist and is not part of this runner; `Rscript run_all.R --sensitivities` runs all three in
the order below.

- Subgroups (Sahel, Eastern Africa, surveys up to and after the median year 2014): `R_cbh/sensitivity/subgroups/01_fit.R`
  then `02_report.R`, writing `results/cbh/subgroups_dhsmics_map_gamma2_v2/`.
- Adjustment without wasting and stunting: `R_cbh/sensitivity/nutrition/01_fit.R` then `02_report.R`, writing
  `results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v2/`.
- Imputed covariates: `R_cbh/sensitivity/imputation_dhsmics/01_impute_national.R` to `06_report.R`, writing
  `results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/` and `results/cbh/covariate_imputation_dhsmics_v6/`
  ([instructions](../sensitivity/imputation_dhsmics/README.md)). It keeps its own settings so that v5's provenance is not
  touched. Its multiple-imputation check and report are pending. <!-- V6_MI_PENDING -->

`--dhs-only` on both subgroup scripts and both nutrition scripts reproduces the DHS-only refits of v3
(`subgroups_map_gamma2_v1`, `nutrition_adjustment_map_gamma2_v1`). The DHS-only imputed version v4 is built by
`R_cbh/covariates/15_impute_national.R` and `16_impute_regional.R`, fitted by `run_regional.R --imputed` and checked by
`R_cbh/sensitivity/imputation/`. The 18-variable benchmark `primary_map_regional18_gamma2_v2` (settings version
`regional18`) and the eleven-variable `primary_map_gamma2_v1` are preserved history. See [sensitivity
fits](../sensitivity/README.md).

## Outputs and scope

No DHS or MICS extraction, HIV refitting, sensitivity fitting or `.tex` file writes are part of this runner. Annual
national MAP extraction uses existing local rasters. LaTeX table source is provided as `.latex.txt` for author insertion.
The older runner `run.R` explicitly selects the preserved legacy version.

Model objects and the prepared dataset remain under ignored
`data/derived_cbh/models/primary_map_regional17_dhsmics_gamma2_v5/`. Only aggregate results, code and provenance are
exported. Model cache signatures include data, formula, knots, code and software versions; invalid caches are refitted.

The combined sample contains the v3 sample entirely. Differences from v3 reflect the added MICS records and the rescaling
of confounders on the larger sample, and cannot be attributed to individual surveys or covariates. Both versions use the
same fixed HIV imputation and conditional spline uncertainty. Uncertainty in exposure, filled covariates and smoothing
parameters is not propagated into the reported intervals.

Paper export (after validation): `python3 R_cbh/reporting/export_paper.py --destination PATH` writes a local review
manifest. Add `--copy` to copy only the listed primary and supplementary figures, captions and table-source files. The
asset list includes the v6 supplementary figure and caption written by `imputation_dhsmics/06_report.R`, so the export
stops until they exist. The exporter backs up replaced assets and checks that all TeX hashes remain unchanged.

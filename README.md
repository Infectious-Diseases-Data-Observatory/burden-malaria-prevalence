# Malaria prevalence and child mortality

> **Current analysis (23 September 2026):** DHS and UNICEF MICS complete birth histories, seven separate child age-band MAP models with `gamma=2` and 17 regional/annual covariates. Primary Figures 1–3, the inclusion flow and the age-band, national, annual and state tables use `primary_map_regional17_dhsmics_gamma2_v5`: 132 surveys (95 DHS, 37 MICS) in 36 countries, 7,498,459 child-band records and 102,282 deaths; [current outputs](results/cbh/primary_map_regional17_dhsmics_gamma2_v5/RESULTS.md). The DHS-only primary `primary_map_regional17_gamma2_v3` is kept as the benchmark. Start with the [analysis plan](docs/ANALYSIS_PLAN.md), the [primary pipeline](R_cbh/primary/REGIONAL_REFIT.md), the [current result index](<Key results/README.md>), [supplementary results](<Supplementary results/README.md>) and [open issues](docs/OPEN_ANALYSIS_ISSUES.md). The [code audit](docs/CODE_AUDIT.md) is a dated 15 September snapshot. `Rscript run_all.R --help` lists explicit current commands. The survey-region description further below is historical; superseded scripts are in [the archive](archive/2026-09-15-code-audit/).

## Current analysis

Each child contributes one record per completed-month age band (<1, 1–5, 6–11, 12–23, 24–35, 36–47, 48–59) that it entered alive in the 60 months before interview and whose full width ended by interview. Band death is modelled with a binomial/cloglog GAM per age band: a smooth in the regional MAP PfPR₂₋₁₀ at band entry, a calendar-time smooth, 13 survey-region and four national annual confounders, and survey, country and region random intercepts. The fitted curves give attributable fractions and national malaria-attributable under-five deaths for 2000–2024, compared with IHME and UN IGME.

| Folder | Content |
|---|---|
| `R_cbh/` | Child age-band builder (`01_make_analysis_data.R`), regional covariates, HIV incidence imputation, primary fits (`primary/`), burden, reporting and sensitivity analyses |
| `R_mics/` | MICS inventory, region mapping, MAP extraction, conversion to the Births Recode layout and regional covariates ([README](R_mics/README.md)) |
| `results/cbh/` | Aggregate results by versioned folder; the current primary is `primary_map_regional17_dhsmics_gamma2_v5/` |
| `results/mics_inventory/` | MICS inventory, [eligibility assessment](results/mics_inventory/ELIGIBILITY.md) and exclusion attribution |
| `docs/` | Analysis plan, open issues, code audit |
| `Key results/`, `Supplementary results/` | Index of current primary outputs (generated) and of supplementary results (maintained by hand) |

### Reproduction sequence

Run from the project root; `Rscript run_all.R --help` prints the same order. Microdata and model objects stay in the git-ignored `data/` tree. Steps 1–4 use local files only. Step 0 reuses cached public-data snapshots; only its explicit retrieval scripts (`R_cbh/covariates/01_fetch_published.R`, `12_fetch_nutrition.py`) query the DHS API.

```sh
# 0. Upstream DHS inputs (not run by any run_all mode): child-band shards, the child HIV
#    incidence panel and the 17-variable regional overlay (R_cbh/covariates/README.md), in
#    that order: the covariate audits (03, 13, 14) read the HIV panel.
Rscript run_all.R --check-inputs         # local DHS dataset prerequisites
Rscript R_cbh/01_make_analysis_data.R
Rscript R_cbh/hiv/01_fit_incidence.R     # rerun with --extended for the imputed-covariate panel
Rscript R_cbh/covariates/run.R           # then 11, 12 (python3; queries the DHS API, reuses cached pages), 13 and 14

# 1. MICS build (R_mics/README.md): R_mics/00a_extract.R, then 00 to 11 in order.
#    Skips 08_build_child_bands.R when its manifest exists; --force runs it.
Rscript run_all.R --mics-build

# 2. Primary v5: prepare, fresh fits, effects, burden, figures, tables and validation
Rscript run_all.R --primary              # R_cbh/primary/run_regional.R
Rscript R_cbh/primary/run_regional.R --report-only   # from saved fits

# 3. Sensitivity analyses on the v5 sample: subgroups, no nutrition, imputed covariates (v6).
#    The imputed-covariate fits take hours; run detached.
Rscript run_all.R --sensitivities

# 4. Paper assets (validates by default; --copy copies figures, captions and table source, never TeX)
python3 R_cbh/reporting/export_paper.py --destination <Overleaf folder>
```

Rerunning step 0's builder or `R_mics/08_build_child_bands.R` rewrites the build manifest and forces fresh v5 and sensitivity fits (steps 2 and 3); see [REGIONAL_REFIT.md](R_cbh/primary/REGIONAL_REFIT.md).

## Historical survey-region pipeline (R_dhs)

This section and those below describe the earlier `R_dhs` analysis and are kept for reference. Most of its scripts from 04 onwards are in [`archive/2026-09-15-code-audit/R_dhs/`](archive/2026-09-15-code-audit/). `R_dhs/` keeps the data-access and extraction scripts whose outputs (`data/derived_dhs/`) the current pipeline still reads, and two wrappers (`11_study_flow.R`, `33_key_results.R`) that call `R_cbh/reporting/`.

Relationship between *Plasmodium falciparum* parasite prevalence (MAP PfPR₂₋₁₀)
and all-cause child mortality across sub-Saharan Africa, estimated at the
DHS/MIS survey-region level and extrapolated to a national malaria-attributable
burden.

📄 **[SUMMARY.md](SUMMARY.md)** / **SUMMARY.pdf** — methods and results of this
historical pipeline with the key figures embedded.

### What the analysis did

Each DHS or MIS survey region contributes one observation: all-cause
post-neonatal (1 month to 5 years) and neonatal mortality over the **12 months
before interview**, computed from the birth histories with the DHS
synthetic-cohort life table (`DHS.rates::chmort`), paired with the
population-weighted MAP PfPR₂₋₁₀ estimate for that region in the survey year and
18 standardised covariates (regional DHS indicators, national WUENIC vaccine
coverage, child HIV prevalence, World Bank economic series).

Approximate death counts (rate × exposure) are modelled with a ridge-penalised
negative-binomial GAM: a country random intercept, an offset, one ridge block
for the covariates, and one of five prevalence-by-time structures chosen by AIC.
Post-neonatal mortality is the primary outcome and neonatal mortality the
negative control. The prevalence-by-time surface is also fitted in Stan
(`brms`) so that uncertainty in its smoothness propagates into every
attributable fraction, and the additive, linear-in-time and full-surface
structures are compared by PSIS-LOO and by survey-grouped 10-fold
cross-validation. Population-average attributable fractions are then applied to
national MAP prevalence and IGME all-cause mortality to give a burden series for
2000–2024, compared with IHME/GBD and WHO.

### Repository layout (historical pipeline)

```
R_dhs/00_config.R            paths, constants (CHMORT_PERIOD, AF_REFERENCE, ...) and shared helpers
R_dhs/01-03                  data access, MAP extraction, panel assembly
R_dhs/04-13                  main models, plots, sensitivities, burden, triangulation      [mostly archived]
R_dhs/14-29                  manuscript figures and further sensitivity analyses           [archived]
R_dhs/30-37                  Bayesian (brms) refits, model comparison, subgroups, burden trend [mostly archived]
R_dhs/run_all.R              runs the pipeline in order
R_dhs/COVARIATE_SOURCES.md   provenance of every covariate
tests/                       unit tests for the region-name matching
results/dhs_rebuild/         figures (.png) and tables (.csv)                    [tracked]
data/                        raw inputs, derived panel, cached model fits       [git-ignored]
R/                           legacy three-component pipeline, superseded          [kept for reference]
archive/                     superseded scripts                                  [dated audit folders tracked]
```

### Running the historical pipeline

From the repository root:

```bash
Rscript R_dhs/run_all.R                   # full rebuild from local recodes and cached MAP
Rscript R_dhs/run_all.R --analysis-only   # skip data access and MAP extraction
Rscript R_dhs/run_all.R --kfold           # also run the hour-long survey-grouped k-fold (script 35)
```

Every `brms` script caches its sampled fit under `data/` with a fingerprint of
the panel it was fitted to and refits only when the panel changes, so a warm
re-run of the whole pipeline takes about half an hour; a cold one adds roughly
five minutes per Bayesian model. Rendering the summary:

```bash
Rscript -e 'rmarkdown::render("SUMMARY.md", output_format = "pdf_document")'
```

R packages: `rdhs, DHS.rates, terra, sf, malariaAtlas, countrycode, mgcv, brms,
loo, posterior, ggplot2, rmarkdown` (and a working Stan toolchain for `brms`).

#### Data inputs

| Input | Where | How to obtain |
|---|---|---|
| DHS/MIS Births Recodes | `data/dhs/*.rds` | DHS account and registered project; `01_access_dhs_data.R` downloads via `rdhs` using the cached login |
| DHS admin-1 boundary files | `data/dhs_boundaries/` | DHS Spatial Data Repository, fetched by script 01 |
| MAP PfPR₂₋₁₀ rasters, GPW population | `data/map/` | `malariaAtlas` and GPW, fetched by script 02 |
| UNICEF WUENIC coverage | `data/` | UNICEF global workbook, extracted by `02b` |
| DHS StatCompiler indicators (improved water, improved sanitation, wasting) | `data/derived_dhs/statcompiler_covariates.csv` | DHS Program API via `rdhs`, pulled by `02c` at national and region level for every registry survey |
| UNAIDS child HIV numbers, World Bank series | `data/` | fetched in script 03 (World Bank API) and from the UNAIDS 2025 workbook |
| IHME/GBD under-5 malaria deaths | `data/ihme_malaria_u5_deaths_by_country.csv` | GBD Results tool (free login); needed for the burden comparison only |
| WHO WMR 2025 deaths | fetched by script 10 or supplied via `WHO_JSON_PATH` | WHO GHO API |

### Pipeline map

| Script | Does |
|---|---|
| `01_access_dhs_data.R` | Survey registry, recode and boundary downloads |
| `02_extract_map_pfpr.R`, `02b_extract_unicef_immunisation.R`, `02c_fetch_statcompiler_covariates.R` | Population-weighted MAP PfPR₂₋₁₀ per survey region; WUENIC panel; DHS API (StatCompiler) household covariates |
| `03_build_analysis_dataset.R` | Regional mortality (12-month window), region-name reconciliation against boundaries, covariates, inclusion flags; writes `region_merge_quality.csv` and `analysis_inclusion_counts.csv`. Improved water, improved sanitation and wasting come from the StatCompiler indicators matched to the survey regions by name (`statcompiler_region_matching.csv`); wasting is imputed from a country-year model where a survey measured no anthropometry, and those rows count as imputed. `--refresh-covariates` re-attaches them without rebuilding mortality |
| `04_fit_main_models.R`, `05_make_main_plots.R` | Five-specification AIC comparison, REML refits for the three outcomes, attributable fractions, figures 1–5 |
| `06`–`08`, `15`–`20`, `22`–`24` | Sensitivities: samples, model structure, likelihood, specification forest, country slope, timing, ages 5–14, period and sub-region strata, fieldwork season, deprivation proxy, intervention targeting |
| `09_validate_reproduction.R` | Comparison against the legacy aggregate (informational) |
| `10_time_surface_and_burden.R`, `25_national_burden_lagged.R` | Three time structures extrapolated to national burden 2000–2024; 2024 against IHME and WHO |
| `11_study_flow.R` | Flow diagram, every count read from the pipeline outputs |
| `12_triangulation_rct.R`, `13_ihme_share.R`, `14_manuscript_figures.R` | Bednet-trial triangulation, IHME share comparison, main-text manuscript figures |
| `21_map_vs_measured_prevalence.R` | MAP against DHS-measured parasitaemia, by era |
| `26_audit_survey_coverage.R` | Every registry survey accounted for, with the reason for any omission |
| `27_horizon_lag_selection.R` | Mortality window 12–60 months × prevalence lag 0–2 years, with a cluster-level jackknife life table |
| `28_survey_map.R`, `29_subgroup_fits.R` | Survey map and timeline; mgcv subgroups by era and region |
| `30`–`32` | brms tensor surface, brms subgroups, three-way region model (`BRMS_OUTCOME=neonatal` refits 30 and 31 to the negative control) |
| `33_key_results.R` | Copies the curated figures into `Key results/` |
| `34_brms_model_ladder.R`, `35_brms_ladder_kfold.R` | Additive vs linear-in-time vs full surface: PSIS-LOO, stacking; survey-grouped 10-fold CV |
| `36_neonatal_as_covariate.R` | Neonatal mortality (12- and 60-month) as a covariate; neonatal vs post-neonatal scatter |
| `37_burden_trend_brms.R` | National burden 2000–2024 under the Bayesian ladder with credible bands |
| `40_build_person_time.R`, `41_extract_map_window_years.R` | Person-time design: deaths and person-months by survey region, 12-month window (five per survey) and age segment; MAP prevalence for each window's own year |
| `42_person_time_models.R`, `43_person_time_brms.R`, `44_deaths_by_age.R` | One negative-binomial model per age band with a smooth prevalence effect, reported as hazard-ratio curves and attributable fractions against 0% prevalence at 10, 30 and 50% (no per-10-point slopes: the response saturates above about 30%); Stan refits; age-at-death histogram |
| `45`–`47` | Data-and-curve figures for the six-band splits, covariate forest by band, the joint model with band-specific random effects and ridge blocks against the separate fits |
| `48_person_time_burden_nga_cod.R`, `49_person_time_burden_ssa.R` | Malaria-attributable under-5 deaths 2005–2025 from admin-1 prevalence, the Stan attributable fractions by age band (script 43) and UN IGME or IHME all-cause deaths: Nigeria and DRC, then every sub-Saharan country against IHME by age block |
| `50_age_specific_attributable.R` | Attributable fraction of all-cause mortality and attributable deaths per 1,000 child-years by age band as functions of PfPR₂₋₁₀, from the Stan posteriors of script 43 with the mgcv curves alongside in the table (figure 36) |
| `51_snow_polygon_prevalence.R` | **Archived 18 September 2026** (`archive/2026-09-18-snow-comparison/`): Snow et al. (2017) polygon PfPR₂₋₁₀ population-weighted onto survey regions. The Snow comparison was removed from the analysis plan. |

Two constants in `00_config.R` define the primary analysis: `CHMORT_PERIOD`
(12 months) and `AF_REFERENCE` (the 1% counterfactual prevalence). Country
random PfPR slopes are off (`INCLUDE_COUNTRY_PFPR_SLOPE`) and examined only as a
sensitivity.

### Legacy pipeline

The `R/` folder keeps the remaining scripts of the earlier three-component
analysis (data fetch, the RDT-to-microscopy conversion of component 3 and the
component-2 malaria-death time series `11_malaria_deaths_timeseries.R`); components 1
and 2 (country-level malaria share against PfPR, and a DHS mixed model on the log rate
with country random slopes) and the former root runner are in
`archive/2026-09-15-code-audit/`; the root `run_all.R` is now the explicit router for the
current analysis. The legacy pipeline is superseded and kept for reference; its RDT-to-microscopy
conversion (`results/rdt_microscopy_conversion.csv`, microscopy ≈ 0.74 × RDT)
was read by `21_map_vs_measured_prevalence.R` (archived).

### Caveats of the historical pipeline

- The national burden is not a validation against a common estimand: the model
  includes indirect malaria-associated deaths, whereas IHME and WHO assign
  deaths specifically to malaria.
- Regional mortality over a 12-month window rests on few deaths; the
  cluster-level jackknife in script 27 gives the sampling error, and the
  neonatal rate at that window is mostly noise.
- The 2000–2024 trend is not well identified: structures separated by well
  under one AIC point imply changes from −53% to +4%, and the Bayesian bands
  on the year-2000 burden under the time-varying structures are very wide.
- MAP PfPR₂₋₁₀ is available through 2024; 2025 inputs are not.

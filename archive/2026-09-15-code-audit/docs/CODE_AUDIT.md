# Code audit and proposed disposition

7 September 2026 · Source baseline `17828d1` · Companion to the [analysis plan](ANALYSIS_PLAN.md)

## Scope and conclusion

The inventory covers **25 scripts in `R/` and 52 in `R_dhs/`**, approximately 17,900 source lines, plus the root `run_all.R`. The review inspected script responsibilities, inputs/outputs, model specifications, shared helpers, caching and orchestration, with detailed review of the person-time pipeline and its upstream dependencies. Both READMEs, the covariate-source note, the existing covariate-audit README and the region-matching tests supplied context. Script headers and documentation were checked against executable code where they conflicted.

**Make `R_dhs/40–50` the core of the new workflow, after separating their data preparation, model fitting and plotting.** Retain the useful geography and covariate work in `00–03`. Consolidate sensitivity ideas rather than preserving every historical script. The two current directories are not independent pipelines that can simply be reduced to one by deleting `R/`.

All 78 scripts parsed successfully using `Rscript --vanilla` and base R `parse(file=...)`, without evaluating their contents. This checks syntax only. No raw DHS data, model caches or numerical output tables were read; no statistical model, download script or existing test suite was executed. Findings below distinguish directly observed code behaviour from methodological risks requiring validation. No existing code or outputs were modified. The unrelated untracked `lit.txt` was left untouched.

## Findings that determine the implementation order

Line references below refer to the source baseline, not to a refactored implementation.

| ID | Finding and evidence | Consequence and planned action |
|---|---|---|
| A01 | **Person-time data depend on the survey-region model.** `R_dhs/42_person_time_models.R:76` reads `MODEL_BUNDLE_RDS`; lines 114–118 obtain covariates from `read_analysis_data()`. `43:41–45` and `47:40–55` repeat this dependency. | A primary person-time dataset cannot currently be built independently. Move covariates, year centring and preprocessing into stage 1; do not require script `04` merely to fit person-time models. |
| A02 | **Inherited survey/geographic exclusions.** `40:183–185` and `41:30–32` restrict the registry to surveys present in survey-year MAP output. `03:337–339` requires mortality and covariate aggregation to succeed before retaining a region, and `42` later merges those regions into person-time. | A survey can contain eligible retrospective person-time yet be lost for missing survey-year MAP or life-table failure. Build the geography/covariate frame independently and record person-time-specific eligibility. |
| A03 | **National covariates stay at survey year across five windows.** `03:489–657` joins national predictors using `data$year`; `42:115–118` copies their `_analysis` columns into every window/segment. | Malaria exposure changes with historical year, but national HIV, vaccination and other predictors do not. Build annual panels and time-align them explicitly, with the present join retained as a sensitivity. |
| A04 | **Exposure year weights are pooled across ages.** `40:109–129` computes year composition from all under-5 person-time and `40:239–250` collapses it to region/window/year. `42:87–108` shares that exposure across segments. | Age bands with different year composition receive the same mixture. Add segment-specific year weights and compare with the existing approximation. |
| A05 | **The fitted counts are constructed, not observed.** `42:120–122` rescales a design-weighted rate to unweighted person-months and rounds. `40:75` retains cluster IDs but the saved sums do not contain strata; `42/43` do not use those sums for design-based uncertainty. | Preserve for reproducibility but test weighting, rounding and survey dependence. Ordinary NB/Poisson uncertainty is not automatically design-correct. This is a methodological limitation, not a syntax error. |
| A06 | **Person-time caches can miss input changes.** `40:42–43,192–195` caches by segment count and survey ID. `41:46–47` caches by survey and year-span settings. Most decisively, `43:70–76` checks row count, total deaths, formula and selected sampler settings, but not actual prevalence/covariate values, person-time, priors or preprocessing. | Changing HIV/WASH values while preserving rows and total deaths can reuse an old fit. Fingerprint actual inputs, code/configuration and priors; validate producer IDs in downstream consumers. |
| A07 | **The burden runner has a backward dependency.** `49:72` reads `latest_country_burden_comparison.csv`, produced by `10:483`. `R_dhs/run_all.R` runs `49` before `10`. It also runs `02c` before refreshing/building the registry with `01`. | A warm workspace can hide failures that occur on a fresh build, and covariates may use an old registry. Use explicit prerequisites; generate the burden country universe from inputs. |
| A08 | **Optional Bayesian fits feed unconditionally attempted outputs.** The runner only fits `43` with `--person-time-brms`, but attempts `48/49`, whose shared AF helper requires those Stan caches. `50:50–52,148` selects its reporting engine from cache-file existence. | Optional work can fail or change the displayed estimator depending on workspace state. Declare the engine and required outputs in the run specification; report skipped/failed outputs accurately. |
| A09 | **Unavailable years are treated inconsistently.** Script `02` declines survey-year substitution. `42:88` clamps both exposure and lag years to 2000–2024, with a >50% pre-2000 exclusion based on the unlagged calendar shares. `49:83–91,180–187` fills any missing year using the last available/configured year. | Historical gaps as well as future years can receive copied values. Require observed coverage for the primary analysis; preserve labelled borrowing/forecast scenarios and evaluate coverage separately by lag. |
| A10 | **AF truncation differs between figures and burden.** `43:118–125` and `50:61–74` retain signed AFs; `00_config.R:1329` returns `pmax(1 - exp(-est), 0)` for burden. | Near-null and negative associations can contribute positive burden after draw-wise truncation. Use a common signed estimator and a separately labelled non-negative scenario. |
| A11 | **Age-combined intervals are averaged endpoints.** `42:282–284` weights the separate bands' AF estimates and lower/upper limits by observed death shares. | These endpoint averages are not a calibrated interval for the sum. Aggregate draws/replicates, declare the weighting population and address shared age-band sampling uncertainty. |
| A12 | **Missingness differs from the documented general rule.** `03:772–832,857–880` imputes wasting before generic eligibility, then restores an imputation flag. `00_config.R:750–791` assesses missingness after that fill. WASH national fallbacks are values, not counted as missing by this rule. | Pre-imputation missingness and regional/national availability must be reported separately. Specify wasting imputation as an exception rather than saying all model covariates have ≤5% original missingness. |
| A13 | **Geographic/temporal provenance can be lost in external joins.** `00_config.R:723–730` returns a nearest-year value with no maximum gap or source year. `03:611` keeps the first duplicate HIV country-year. `02b` assumes zero for pre-series and absent vaccine series. | Save source year/distance, missingness reason and duplicate-resolution evidence. Test vaccine-zero and censored-HIV assumptions; validate every country-year key. |
| A14 | **The two fitting engines use different regularisation.** `42:153–171` estimates a ridge penalty through fREML. `43:56–61` assigns fixed Normal(0,1) coefficient priors, including other parametric terms. | Describe the Stan fit as a related Bayesian model, not an identical ridge refit with “full uncertainty.” Freeze intended priors and compare on identical datasets. |
| A15 | **Diagnostics are reported but do not gate downstream inference.** `43:110–118` summarises divergences, R-hat and bulk ESS; its sampling settings comments acknowledge prior divergences. There is no stop before saving/reporting problematic fits. `50:108–110` only warns on a row-count mismatch. | Add reviewed convergence criteria and input-hash checks before figures/burden are final. Actual current convergence has not been assessed in this audit. |
| A16 | **Data and analysis are interleaved with plotting.** `42` builds the modelling table, fits multiple groupings/sensitivities and plots. `43/47–50` also combine inference and plotting. `46:24` obtains labels from the survey-region model's coefficient table. | Move modelling-table creation into stage 1, inference into stage 2 and visualisation into stage 3. Labels belong in the covariate dictionary. |
| A17 | **Burden construction repeats calculations and strong assumptions.** `48/49` duplicate regional MAP extraction and allocation logic; `49:122–140` treats absent age columns as zero and clips negative infant remainders; `49:215–260` uses pooled or latest-survey death shares. | One parameterised burden function; explicit source schemas; age-total conservation; uncertainty/sensitivity for sparse or fixed death shares. Do not treat all missing age cells as structural zeros. |
| A18 | **Coverage and descriptive plots are still tied to different samples.** `26:40–41` requires the legacy full panel and classifies omissions relative to it. `44:34–37` counts full-survey deaths without the region matching/MAP-window filters used for modelling. `11/28` describe the survey-region sample. | Generate inclusion and age-at-death tables from the primary person-time build. Keep a separate legacy comparison, not a legacy-dependent primary flow diagram. |
| A19 | **Several documentation claims are stale.** READMEs describe survey-region analysis as primary, and say `13` reads the legacy Method-1 country panel; the actual `13` builds a country-year IHME/IHME share model from external files. The covariate-audit README says the pipeline was unchanged, but `03` now attaches StatCompiler. | Rebuild the README and results index from the final architecture. Treat source code as authoritative for this plan; do not transfer old headline counts or conclusions into new documentation. |
| A20 | **Legacy scripts can write outside the repository.** `R/34_fig2_penalized.R:47–48` and `R/39_penalised_main.R:79–81` overwrite files in a Dropbox/Overleaf directory. | Remove these actions from the active pipeline. Plotting should export locally; external publication/copying should be a separate explicit action. No such script was executed here. |

Additional methodological checks belong in the implementation rather than being asserted as observed numerical errors: full-month death exposure and age heaping; survey-specific calendar conversions; density-to-population raster weighting; heuristic geographic matches; household/mother/birth denominators in regional covariates; and preprocessing leakage in grouped cross-validation. The plan specifies how to examine these.

## Dependencies to transfer before retiring legacy code

| Existing producer | Remaining use | Replacement |
|---|---|---|
| `R/03_component3_rdt_microscopy.R` | `data/dhs_prevalence_by_region.csv` and `results/rdt_microscopy_conversion.csv` feed `R_dhs/21` and `22` | Optional biomarker input builder and one conversion/age-standardisation helper |
| `R/01_fetch_data.R`, `R/11_malaria_deaths_timeseries.R` and related national extraction code | Public national input tables such as `pfpr_by_country_year.csv`, mortality/birth panels and comparator tables are consumed throughout the rebuild | Dedicated external-input and national-exposure preparation; eliminate dependence on legacy model execution |
| `R/30_full_covariates.R`, `R/33_penalized_model.R`, `R/34_fig2_penalized.R` | Legacy aggregate/fit artifacts feed migration mode, validation and the current coverage audit | Freeze a documented migration reference; primary audit derives from current inputs |
| `R_dhs/04_fit_main_models.R` | Main bundle supplies person-time preprocessing; associated coefficients supply figure labels | Independent preprocessing object and covariate dictionary |
| `R_dhs/10_time_surface_and_burden.R` | Older burden result table supplies script `49`'s country universe | Explicit country/input registry |

`R_dhs/13_ihme_share.R` is **not** currently a reason to preserve `R/02_component1_country.R`: despite README claims, its executable code no longer reads that panel.

## Complete script disposition

“Retire” means remove from the future active workflow after replacement/reproduction checks. It does not mean delete now. “Extract” means preserve useful functions/data construction, then retire the current mixed-purpose script. “Optional” means excluded from the default person-time run.

### `R/` — 25 scripts

| Script | Current responsibility | Proposed disposition |
|---|---|---|
| `R/00_utils.R` | Legacy paths, geography, public inputs, mortality and Component 1/2 formulas | Extract only unique external-input/biomarker helpers; use the newer geography implementation; retire legacy model configuration |
| `R/01_fetch_data.R` | MAP/population/boundaries, national series, biomarker-survey universe and recodes | Extract still-needed input preparation; replace the biomarker-dependent survey frame |
| `R/02_component1_country.R` | Country malaria share vs prevalence/GDP/DTP3 | Retire from primary; retain only as historical comparison, superseded by the optional `R_dhs/13` |
| `R/03_component3_rdt_microscopy.R` | Measured regional parasitaemia and RDT/microscopy conversion | Extract optional validation-data builder before retiring |
| `R/04_component2_dhs.R` | Measured-prevalence mortality panel, mixed/count/spline models and negative control | Retire; historical reference only |
| `R/05_prediction_10_to_30.R` | Cross-component mortality predictions | Retire; replace any retained contrast with shared person-time prediction functions |
| `R/06_attributable_fraction.R` | Linear primary-model AF curves | Retire; replaced by common nonlinear person-time contrasts |
| `R/07_ng_drc_malaria_deaths.R` | National/subnational Nigeria and DRC burden under two old methods | Retire after generic person-time burden replacement |
| `R/08_country_malaria_deaths.R` | Country burden under the old methods | Retire after national-input transfer |
| `R/09_method_comparison.R` | Old Component 1/2 comparison figure | Retire; optional external comparison uses new saved results |
| `R/10_triangulation_rct.R` | Trial data, age conversion, old-model predictions and plots | Retire duplicate; preserve trial provenance through optional `R_dhs/12` |
| `R/11_malaria_deaths_timeseries.R` | National time-series preparation, old-model burden and WHO/IHME comparators | Extract national-input/comparator preparation; retire model and plotting branches |
| `R/12_paper_figure1.R` | Legacy scatter and linear NB response figure | Retire |
| `R/14_paper_figure_ihme_share.R` | Legacy cross-sectional IHME share figure | Retire; optional `R_dhs/13` replaces it |
| `R/15_ihme_share_time.R` | Several country-year IHME share models/figures | Retire duplicate exploratory models; preserve input parsing only if still needed |
| `R/18_expanded_panel.R` | Expanded MAP-based DHS panel and initial models | Retire after comparing geography/frame against the new independent builder |
| `R/23_region_health.R` | Regional health covariates, access index and refits | Retire duplicate; retain relevant variable definitions/provenance in the new dictionary |
| `R/25_combined_model.R` | Linear/nonlinear × national/regional-covariate model ladder | Retire; historical model-selection exploration |
| `R/26_neonatal_outcome.R` | Legacy neonatal/post-neonatal model comparisons | Retire; use primary age-specific person-time fits |
| `R/29_paper_figures.R` | Legacy manuscript figures, embedded refits and burden | Retire; one plotting layer |
| `R/30_full_covariates.R` | Expanded predictor table and imputation | Retire after freezing any migration aggregate; definitions superseded by `R_dhs/03` audit |
| `R/32_full_model.R` | Unpenalised full covariate models and coefficient forest | Retire |
| `R/33_penalized_model.R` | Legacy ridge models and saved headline fits | Retire from execution; preserve migration reference metadata |
| `R/34_fig2_penalized.R` | Linear-outcome comparison and external figure copy | Retire; eliminate automatic Overleaf copy |
| `R/39_penalised_main.R` | Penalised headline figures, temporal burden and external copies | Retire; eliminate duplicate figures and automatic copies |

### `R_dhs/` — 52 scripts

| Script | Current responsibility | Proposed disposition |
|---|---|---|
| `R_dhs/00_config.R` | Paths/constants, geography, mortality, external joins, imputation, GAMs, age definitions and AF helpers | Split into configuration and focused shared modules; retain geography tests; no top-level directory creation in pure helper files |
| `R_dhs/01_access_dhs_data.R` | Registry and optional recode/boundary retrieval | Retain as data-input functionality; make local-snapshot processing independent of an online registry refresh |
| `R_dhs/02_extract_map_pfpr.R` | Survey-year MAP extraction | Merge into one boundary/year extractor; remove survey-year eligibility gate for person-time |
| `R_dhs/02b_extract_unicef_immunisation.R` | WUENIC panel extraction and source/zero assumptions | Retain in external-input module; make introduction assumptions explicit |
| `R_dhs/02c_fetch_statcompiler_covariates.R` | National/regional water, sanitation and wasting API extracts | Retain optional retrieval plus validated local panel; run after the intended registry is established |
| `R_dhs/03_build_analysis_dataset.R` | Life-table outcomes, regional/national predictors, imputation and sample flags | Split: independent covariates/joins into primary data stage; life-table branch becomes sensitivity; migration branch leaves default run |
| `R_dhs/04_fit_main_models.R` | Five survey-region specifications, ML/AIC selection and outcome refits | Optional region-mortality sensitivity; stop using it as a person-time prerequisite |
| `R_dhs/05_make_main_plots.R` | Survey-region figures from saved fits | Retire default figures; retain useful plot conventions for optional compatibility results |
| `R_dhs/06_sensitivity_samples.R` | Prevalence thresholds and complete cases | Port relevant sample specifications to person-time sensitivity registry |
| `R_dhs/07_sensitivity_model_structure.R` | Country slopes, covariate levels and spline dimensions | Port selected specifications; avoid separate duplicate fit helpers |
| `R_dhs/08_sensitivity_likelihood.R` | Log-Gaussian survey-region rate model | Optional legacy-estimator sensitivity; use NB/Poisson checks for person-time |
| `R_dhs/09_validate_reproduction.R` | Rebuild vs legacy headline models | Retain as migration-only tool; replace significance-preservation criteria with numerical comparisons and explained changes |
| `R_dhs/10_time_surface_and_burden.R` | Survey-region time structures and national burden/comparators | Extract shared comparator/input handling; retire default model/burden calculations; optional temporal-estimator comparison |
| `R_dhs/11_study_flow.R` | Survey-region flow diagram | Retain plotting concept, driven by new person-time inclusion ledger |
| `R_dhs/12_triangulation_rct.R` | Hardcoded trial data and survey-region model predictions | Optional; move trial data/provenance to an input table and adapt prediction to age-appropriate person-time contrasts |
| `R_dhs/13_ihme_share.R` | Country-year IHME malaria/all-cause share model with national covariates | Optional external comparison; separate its analysis and plots; standardise external input parsing |
| `R_dhs/14_manuscript_figures.R` | Survey-region manuscript panels plus independent burden calculations | Retire repeated computation; manuscript assembly should consume saved canonical panels/results |
| `R_dhs/15_specification_forest.R` | Additional linear refits assembled into a sensitivity forest | Retain forest format; remove embedded refits and per-10-point summaries as the main person-time metric |
| `R_dhs/16_sensitivity_random_effects.R` | Country slopes and conditional/population burden comparisons | Port relevant heterogeneity/standardisation sensitivity; retire duplicate national burden computation |
| `R_dhs/17_sensitivity_timing.R` | Lagged MAP, multiple mortality horizons, refits and plots | Merge annual MAP extraction; retain selected timing comparisons in specification registry; life-table construction optional |
| `R_dhs/18_age_windows_5to14.R` | Extended life-table mortality at ages 5–14 | Optional extension outside the under-5 core |
| `R_dhs/19_period_stratified.R` | Period-specific survey-region fits and descriptive panels | Port a limited, declared era sensitivity; retire multiple overlapping plotting/fitting implementations |
| `R_dhs/20_subregion_stratified.R` | Regional and region/time survey-region comparisons | Port selected regional heterogeneity question; keep sample/specification consistent with person-time |
| `R_dhs/21_map_vs_measured_prevalence.R` | Measured-vs-MAP agreement and model comparison | Optional exposure validation; preserve legacy biomarker producer first and separate data from fitting |
| `R_dhs/22_seasonality_of_fieldwork.R` | Interview-month/latitude panels and seasonal comparisons | Optional measured-prevalence sensitivity; extract reusable fieldwork metadata in data stage |
| `R_dhs/23_deprivation_proxy.R` | Within-survey deprivation/prevalence correlations | Optional confounding diagnostic; one descriptive-analysis function, no implied causal resolution |
| `R_dhs/24_intervention_targeting.R` | ITN/ACT extraction and within-survey targeting correlations | Optional intervention analysis; data extraction reusable, excluded from automatic primary adjustment |
| `R_dhs/25_national_burden_lagged.R` | Lagged survey-region model and national burden | Retire default run; if retained, implement as a burden specification on shared input tables |
| `R_dhs/26_audit_survey_coverage.R` | Registry/legacy/rebuild coverage and omission classification | Replace with primary input/inclusion audit; legacy comparison becomes optional |
| `R_dhs/27_horizon_lag_selection.R` | Multiple life-table horizons, cluster jackknife, lag comparison and plots | Extract validation/resampling ideas; retain only declared life-table sensitivities, not exhaustive default selection |
| `R_dhs/28_survey_map.R` | Survey coverage map and timeline | Retain; consume primary inclusion and boundary tables |
| `R_dhs/29_subgroup_fits.R` | mgcv era/region subgroup fits and interaction checks | Consolidate into a small person-time heterogeneity menu |
| `R_dhs/30_brms_tensor_model.R` | Bayesian survey-region prevalence/time tensor model | Retire default fitting; optional estimator/structure comparison using shared model interface |
| `R_dhs/31_brms_subgroup_fits.R` | Bayesian survey-region subgroup refits | Retire duplicate workflow; selected substantive questions move to person-time sensitivities |
| `R_dhs/32_brms_three_way.R` | Bayesian prevalence × time × region ladder | Optional exploratory extension; outside default workflow |
| `R_dhs/33_key_results.R` | Copies selected figures and writes a second index | Replace with one authoritative results index; avoid stale copied collections |
| `R_dhs/34_brms_model_ladder.R` | Bayesian additive/linear-time/full-surface comparison | Extract limited model-comparison ideas; no obligatory parallel survey-region ladder |
| `R_dhs/35_brms_ladder_kfold.R` | Survey-grouped Bayesian cross-validation | Retain grouping principle; rebuild for selected person-time models with training-fold preprocessing |
| `R_dhs/36_neonatal_as_covariate.R` | Neonatal outcome used to adjust post-neonatal mortality | Optional exploratory sensitivity; not primary adjustment or proof of confounding control |
| `R_dhs/37_burden_trend_brms.R` | National burden under survey-region Bayesian ladder | Retire default run; common burden interface if estimator comparison retained |
| `R_dhs/40_build_person_time.R` | Death/person-time cells, year shares, cluster sums and validation | Retain as core data layer; add calendar metadata, segment-year weights, strata, robust caches and inclusion QC |
| `R_dhs/41_extract_map_window_years.R` | Retrospective-year regional MAP extraction | Merge with other extractors into core data layer |
| `R_dhs/42_person_time_models.R` | Modelling table, four groupings, primary/sensitivity fits, summaries and plots | Split across all three stages; one main grouping plus declared alternative groupings |
| `R_dhs/43_person_time_brms.R` | Six-band Bayesian refits, diagnostics, AFs and comparison plot | Retain as proposed reporting engine; robust input hashes, reviewed diagnostics, separate predictions/plots |
| `R_dhs/44_deaths_by_age.R` | Separate raw-data read for age-at-death distribution | Merge aggregation into `40` replacement; retain plot from saved tables on declared sample |
| `R_dhs/45_person_time_data_and_curves.R` | Raw rates beside age-band HR curves | Retain in plotting stage using canonical tables and common sample labels |
| `R_dhs/46_person_time_covariate_forest.R` | Conditional covariate association plot by age band | Retain; save coefficients in analysis stage and use independent dictionary labels |
| `R_dhs/47_person_time_joint_model.R` | Joint model with age-specific nuisance terms vs separate fits | Retain as targeted age/dependence sensitivity; acknowledge shared dispersion; separate fitting and plotting |
| `R_dhs/48_person_time_burden_nga_cod.R` | Nigeria/DRC burden, IGME retrieval, regional MAP, allocation and figures | Merge with `49` into one burden implementation; country-specific outputs become filters/views |
| `R_dhs/49_person_time_burden_ssa.R` | SSA burden under IGME/IHME denominators and comparisons | Retain secondary analysis; extract external inputs/MAP into data stage and figures into plotting stage |
| `R_dhs/50_age_specific_attributable.R` | mgcv/Stan AF and attributable-rate draws plus principal AF figure | Retain shared effect-calculation and figure functionality; declare reporting engine, eliminate cache-driven selection |
| `R_dhs/run_all.R` | Runs nearly every historical and current branch | Replace with explicit data/analyses/plots stages and named optional specifications |

The **root `run_all.R`** is the 78th script checked. Replace its legacy three-component orchestration with the three-stage entry point only after the new modules pass their reproduction checks.

## Supporting code and documentation

- Retain `tests/test_region_match.R` and `tests/test_rkey.R`, separating optional real-data comparisons from pure fixtures as needed. These cover geography, not person-time accounting or model validity.
- `results/dhs_rebuild/covariate_audit/` contains four additional R scripts outside the two requested folders. Its README explains changes relevant to `03`. Preserve this audit trail; inspect and move any reusable tests into `tests/` during implementation rather than allowing production logic to live under `results/`. Those four scripts were inventoried by filename, not included in the 78-script parse/review scope.
- Rewrite the main and `R_dhs` READMEs around the person-time contract after implementation. Preserve covariate-source provenance but replace fixed historical sample counts with generated reports.
- Keep the current `archive/` contents outside the active pipeline. Do not blanket-delete them or raw/cache directories as part of the first refactor.

## Acceptance checks for removing a script

A script is surplus only once its scientific question has been deliberately omitted or its required outputs have a named replacement. Verify that no retained code reads its outputs, the new pipeline runs without its cached artifacts, and the agreed figures/tables still build. For migrated calculations, compare keyed aggregate inputs, model specifications, effects at fixed prevalence anchors and draw-level burden totals; document intended differences rather than demanding unchanged p-value classifications.

For the immediate next step, the smallest useful implementation is the **independent person-time dataset build**. It addresses A01–A06 and provides a reliable basis for deciding which remaining fits and plots are truly redundant.

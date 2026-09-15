# Analysis plan: malaria prevalence and child mortality

Updated 15 September 2026 (implementation audit and primary-only rerun; primary specification unchanged). This is a working analysis plan revised after inspection of exploratory results, not a prospective preregistration. The [previous plan](ANALYSIS_PLAN_BEFORE_MAP_GAMMA2.md) preserves the original audit, superseded specifications and decision history; the [code audit](CODE_AUDIT.md) compares each component with the current implementation and records the archival changes.

**Current primary analysis:** seven separate child–age-band mortality models using **annual MAP PfPR₂–₁₀ at band entry**, fitted with **`mgcv::bam`, `gamma = 2`, PfPR `bs = "cr", k = 5`, and calendar-year `bs = "cr", k = 6`**. Use the full eligible MAP sample, including eligible observations after 2015. Each age band has its own time spline, confounder coefficients and survey/country/region random effects and variance parameters. Section 3.1 gives the complete specification.

**Supplementary results:** comparisons with the annual Snow prevalence estimates, including matched MAP–Snow fits at gamma=2, all Snow penalisation analyses and Snow-based national mortality scenarios. The Snow availability cutoff at 2015 applies only to those analyses and their matched MAP comparators. Snow is not the primary exposure. The [supplementary results index](<../Supplementary results/README.md>) links to the saved figures and tables.

**Primary output selection:** use the [primary-only gamma=2 manifest](../results/cbh/primary_map_gamma2_v1/fit_manifest.csv) and associated result tables in `results/cbh/primary_map_gamma2_v1/`. All seven entries are `map_full`. The preceding mixed MAP–Snow outputs remain unchanged as historical/supplementary comparisons. Run `Rscript run_all.R --primary` for the primary analysis from the saved prepared sample, without data setup. See [pipeline instructions](../R_cbh/primary/README.md).

## 1. Scientific questions and reporting choices

The primary causal question is: **if PfPR₂–₁₀ changes from X% to Y%, how does all-cause mortality in age group g change?** Estimate this separately for **<1, 1–5, 6–11, 12–23, 24–35, 36–47 and 48–59 completed months**. All seven age-specific effects, including neonatal mortality, are primary quantities of interest. Retain estimates regardless of direction or significance.

The exposure is *P. falciparum* prevalence standardized to ages 2–10 years, expressed on the 0–100 percentage scale. Regional annual prevalence is the same exposure definition for every age group. Values assigned to a child's successive bands can differ because their entry years differ. The outcome is death during a band, conditional on reaching it alive; the cloglog model estimates an annualized mortality hazard.

Report relative and absolute mortality changes for specified X%→Y% contrasts in the same target population. National malaria-attributable deaths are a **secondary extrapolation of the primary MAP model**. They additionally require assumptions about transportability, national exposure and the counterfactual intervention. Agreement with IHME, WHO or trials does not establish causal identification. The plan does not seek a null neonatal association or interpret one as proof of no confounding.

The age-band estimand conditions on reaching that band alive. An intervention beginning at birth can change who survives into later bands; conditional band effects are not automatically the total effect of a lifelong intervention. Keep this distinction explicit in causal interpretation and synthetic-cohort survival figures.

## 2. Stage 1 — make the analysis datasets

### 2.1 Inputs and provenance

The active builder is [R_cbh/01_make_analysis_data.R](../R_cbh/01_make_analysis_data.R), with [configuration](../R_cbh/00_config.R) and [dataset documentation](../R_cbh/README.md). It processes authorized local DHS Births Recodes and joins reviewed geography and external-data snapshots. Record source release, units, coverage, file hashes, survey versions and exclusions. Optional input retrieval is separate from local processing; do not refresh online sources because a cache is absent.

Keep every in-scope survey in the registry with a status for missing recodes, boundaries, dates, exposure or covariates. Survey-year MAP availability alone must not determine whether earlier band entries are usable. Configured entry years are 2000–2024; the primary sample has no Snow-derived 2015 cutoff.

### 2.2 Geography and exposure

Use survey/boundary-version-specific region keys, explicit normalization, curated synonyms and reviewed aggregation rules. Do not infer a match by elimination or accept an unreviewed fuzzy match. Keep unresolved records in the exclusion ledger. Current region identifiers include the survey/boundary version; they are not automatically longitudinal regional identities.

Assign regional annual MAP PfPR₂–₁₀ from the **calendar year containing band entry**. Hold this value fixed over the band, independently of death timing. The current analysis does not average prevalence over a child's realized person-time or the full band, interpolate between years, borrow the nearest year, or carry the last exposure forward. Preserve genuine zero prevalence and distinguish it from missing exposure. `calendar_year` is fractional entry time; `entry_year` is the integer source year.

The region is the mother's residence at interview; historical migration is not reconstructed. Record coverage and geographic-source flags. Consolidation and validation of the inherited DHS regional raster extracts remain necessary; the cell-area correction used in national burden extraction does not by itself validate the older regional extraction.

### 2.3 Child–age-band records and eligibility

Use one row per child and age band, with disjoint intervals `[0,1), [1,6), [6,12), [12,24), [24,36), [36,48), [48,60)` months. Predetermined band widths are 1/12, 5/12, 6/12, 1, 1, 1 and 1 year.

1. The child must reach the band alive; exclude bands after death.
2. Entry must satisfy `interview_cmc - 60 <= band_entry_cmc < interview_cmc`.
3. The full potential band must end by interview, **for deaths and survivors alike**. Exclude a death in a band that would otherwise be incomplete.
4. The entry year must lie within the configured exposure period. Apply geography, exposure and selected-covariate availability rules with recorded exclusion reasons.

`death` is 1 for death during the band and 0 for survival through it. `band_years` is the full predetermined width, **not observed time until death**. The earlier aggregate Poisson/negative-binomial person-time dataset is not the primary analysis dataset.

Retain recorded death-age units/days and B6/B7 disagreement flags. The current DHS month-based neonatal coding includes reported days 0–30; its boundary differs from the exact 28-day definition in the Burstein text and the IHME neonatal export. Review ambiguous intervals, impossible dates, duplicate identifiers, age heaping and calendar conversions. Older Ethiopian recodes use the declared +92 CMC approximation; ET8AFL is already Gregorian. Do not replace documented rules with generic calendar heuristics.

Validate synthetic birth histories, age boundaries, five-year eligibility, complete-band selection, annual joins and exclusion accounting. Preserve survey, mother, cluster and stratum identifiers under the ignored data directory for dependence/design checks. [Build validation](../R_cbh/BUILD_VALIDATION.md) and [open issues](OPEN_ANALYSIS_ISSUES.md) record current evidence and limitations.

### 2.4 Covariates, HIV imputation and missing data

The primary adjustment set is **sex, multiple birth, birth order, maternal age at birth, maternal education, household wealth quintile, urban residence, log child HIV incidence, log GDP per capita, log health expenditure per capita and political stability**. Wealth is categorical; continuous confounders enter linearly after the saved centering/scaling. Retain the existing scaling for the current fitted comparisons. Do not select covariates by significance or introduce a generic missingness threshold to redefine this set.

Annual national covariates are assigned at band-entry year. Wealth, urban residence, education and other survey measurements describe interview-time conditions; do not relabel them as reconstructed historical measurements. Vaccines, maternal death fraction, hypothetical children, wasting, WASH and intervention coverage are not part of the current primary adjustment set. Optional candidate fields can remain in the base dataset without triggering complete-case exclusion.

Use **child HIV incidence at ages 0–14 per 1,000 uninfected population**, not HIV prevalence or infection counts. The [HIV imputation model](../R_cbh/hiv/README.md) uses child and adolescent incidence jointly, calendar-year trends, geographic/country effects and within-country temporal correlation. Preserve reported child rates, handle censored values as censored, and predict missing child series. Adolescent incidence informs imputation; it is not a separate predictor in the mortality model.

The current primary fits use **one fixed posterior-median child-incidence imputation** and complete cases for the other selected model variables. Keep censoring, imputation and missing-source flags; unsupported country-year values remain missing. The historical HIV prevalence field in base shards is provenance only. There is no country/overall-median imputation for the other primary covariates.

Report missingness and exclusions by country, survey, year, age and exposure. Missing HIV/source data, unresolved geography and complete-case selection remain substantive issues. Retain posterior HIV trajectories for future uncertainty propagation through all seven primary models. For validation with held-out surveys/countries, estimate preprocessing and imputation within training folds where appropriate.

### 2.5 Analysis-ready outputs

| Output | Role |
|---|---|
| `data/derived_cbh/manifest.rds` and `child_bands/<survey>.rds` | Authoritative build manifest and eligible child-band shards with audits |
| `data/derived_cbh/hiv_incidence/child_incidence_country_year.csv` | Incidence estimates and source/imputation flags |
| `data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/complete_case_dataset.rds` | Saved incidence-adjusted primary modelling sample and preprocessing |
| `results/cbh/primary_map_gamma2_v1/fit_manifest.csv` | Seven freshly fitted primary MAP gamma=2 objects; all entries are `map_full` |
| Aggregate selection, missingness, geography and model diagnostics | Inspectable accounting without microdata |

The primary fitting sample contains **1,817,912 distinct children contributing 5,885,022 child-band records and 82,415 deaths from 105 surveys in 34 countries**. The [15 September verification](../results/cbh/code_audit_2026_09_15/primary_fit_checks.csv) compared every primary fitted outcome, predictor and offset with the saved prepared sample and verified the selected model hashes and knots. Keep raw/derived microdata and fitted objects under ignored `data/`; export only approved aggregates to results. Validate unique child-band keys and input hashes. The primary-only fitter now reads this prepared dataset directly and checks every outcome, predictor and offset. Its reference MAP knots are preserved in a committed snapshot. A future consolidated data stage must still operate independently of existing outcome models and results tables.

## 3. Stage 2 — primary and sensitivity analyses

### 3.1 Primary working model

**Selected on 10 September 2026: separate age-band MAP models with gamma=2.** For each of the seven bands, fit:

```r
form <- death ~
  s(pfpr_pct, bs = "cr", k = 5) +
  s(calendar_year, bs = "cr", k = 6) +
  sex + multiple_birth + z_birth_order + z_maternal_age_birth +
  z_maternal_education_years + wealth_quintile + urban +
  z_log_hiv_incidence + z_log_gdp_pc +
  z_log_health_expenditure_pc + z_political_stability +
  s(survey, bs = "re") + s(country, bs = "re") +
  s(region, bs = "re") + offset(log(band_years))

# d_g: the full eligible MAP sample for this one age band.
# knots_g: saved PfPR and calendar-year knots from the MAP reference fit.
fit_g <- mgcv::bam(
  form, data = d_g, knots = knots_g,
  family = binomial(link = "cloglog"),
  method = "fREML", discrete = TRUE,
  gamma = 2, select = FALSE,
  na.action = na.fail, control = mgcv::gam.control(maxit = 100)
)
```

For child i and band g, with width Δg in years:

\[
D_{ig}\sim\operatorname{Bernoulli}(q_{ig}),\qquad
q_{ig}=1-\exp(-\Delta_g\lambda_{ig}),
\]

\[
\log\lambda_{ig}=\alpha_g+f_g(P_{r(i),y_{ig}})+h_g(t_{ig})+
\mathbf X_{ig}^{\mathsf T}\boldsymbol\beta_g+
u_{s(i),g}+v_{c(i),g}+b_{r(i),g}.
\]

Here P is regional annual MAP PfPR₂–₁₀, y is band-entry year, and t is fractional calendar year at entry. The exposure surface has no age-band-specific definition; the function f and the entry year depend on the band. Equivalently, `cloglog(q_ig) = log(Δg) + log(λ_ig)`.

Each fit estimates its own intercept, PfPR and time curves, confounder coefficients and survey/country/region random intercepts. Within each model, the random effects are independent zero-mean Gaussians, with a separate estimated variance for each grouping type. No time spline, random intercept or variance parameter is shared across age bands. Separate fitting does not imply independent sampling errors across ages because children and surveys contribute to multiple bands.

`gamma = 2` changes smoothing selection throughout the model, including time and random effects; it is not a PfPR-only multiplier. The primary PfPR basis remains **cr**, not cs. No monotonicity constraint is imposed. Reuse the reference MAP knots for these promoted fits; they are recorded and verified in the fitting pipeline.

The unweighted likelihood, fixed HIV imputation and conditional model covariance are current working choices. Address survey design, residual dependence and imputation/exposure uncertainty before claiming comprehensive uncertainty. The primary decision follows review of results and user instruction; it does not establish superior causal identification or formal model-selection evidence.

**Authoritative implementation:** [primary/01_fit.R](../R_cbh/primary/01_fit.R), configured in [settings.R](../R_cbh/primary/settings.R), with `map_full_age_1` through `map_full_age_7` in the [primary manifest](../results/cbh/primary_map_gamma2_v1/fit_manifest.csv). Fitted objects are in `data/derived_cbh/models/primary_map_gamma2_v1/`. The fitting sample is read from the existing prepared dataset, with saved scaling and [reference knots](../R_cbh/primary/reference_knots.csv). The older combined MAP–Snow workflow remains for supplementary reproduction; its `map_full` results are the numerical comparison for this rerun. The older `sensitivity_single_imputation/age_*.rds` fits used gamma=1.

### 3.2 Effect summaries

For each of the **seven** age groups and any specified X%→Y% contrast:

\[
HR_g(X\to Y)=\exp\{f_g(Y)-f_g(X)\},\qquad
\%\Delta_g=100[HR_g(X\to Y)-1].
\]

Absolute changes require an explicit target population or calibrated baseline: `rate_g(Y) - rate_g(X)`. Express rates with declared person-year units. For band death probabilities use `q = 1 - exp(-Δg × rate)`; a hazard ratio is not a death-probability ratio. Random effects and other additive predictors cancel in a within-target hazard ratio, but not in marginal absolute-rate standardization.

Report full curves, 40%→20% contrasts and exposure support; additional AF anchors at 10%, 30% and 50% remain useful. For reference p0, `AF_g(p,p0) = 1 - exp[f_g(p0)-f_g(p)]`. Keep the zero reference explicit, assess a 1% reference as a sensitivity, and flag extrapolation. Retain signed effects and intervals, including negative values. The current CBH burden calculation does not truncate negative contributions.

Age-band intervals currently condition on smoothing parameters, one HIV imputation and point-estimate exposure. Do not sum marginal confidence limits or assume independent age-band errors to construct national intervals. Future uncertainty propagation should reuse common survey/design replicates and HIV trajectories across bands and countries before aggregating.

### 3.3 Diagnostics and validation

Check convergence, coefficient rank, finite covariance, smoothing gradients/Hessian, fitted outcomes, residual patterns, exposure support and influential surveys/countries. Verify every stored outcome, predictor and offset against intended inputs. Check compact spline contrasts/variances against full prediction matrices, including zero-PfPR extrapolation. Keep exact model hashes, sample definitions, knots, software versions and warning records.

Numerically flagged fits may receive a tighter-tolerance restart under the **same gamma and specification**; record original and selected attempts. Do not use fREML values to rank models across different gamma values or exposure sources. All seven currently selected MAP gamma=2 fits passed the recorded numerical checks.

Use survey-grouped cross-validation for flexible-structure comparisons and country-held-out validation for transportability. Keep all records from each survey together. Examine within-region versus between-region exposure information; the existing coding does not itself identify a within-region causal effect. Assess smoothing-uncertainty-corrected covariance when available and record when it is unavailable. Bayesian sampling checks apply only if a Bayesian extension is undertaken, not as requirements for the selected bam fit.

### 3.4 Supplementary and sensitivity analyses

Each sensitivity needs a declared change, sample, exposure source, gamma, formula, diagnostics and uncertainty method. Report the primary model on the corresponding restricted sample when necessary. Keep specification choices independent of effect direction or significance; do not automatically cross every sensitivity with every other one.

| Priority / role | Comparison | Status and specification |
|---|---|---|
| **Required supplementary: prevalence source** | MAP versus annual Snow estimates | Completed gamma=2 matched fits through 2015; full-sample comparison for context. Section 3.4.2. |
| **Required key sensitivity** | Joint versus separate-age structure; geography; survey period | Future comparisons use MAP and gamma=2. Earlier gamma=1 exploratory results remain historical. Section 3.4.1. |
| Supplementary: smoothing | Primary MAP gamma=2 versus MAP gamma=1; Snow gamma=1, 1.4, 2 and cs | Existing fits/results retained. Snow cs used gamma=1 and changed PfPR only; it is not the primary specification. |
| Required: eligibility and timing | Shorter lookback; entry-year lag; complete-band versus appropriate partial-follow-up alternatives | Adapt to current child-band records; no silent annual borrowing or return to the superseded window-count likelihood. |
| Required: age coding | Recorded days, month/year heaping, B6/B7 disagreement and alternative age boundaries | Primary has seven DHS completed-month bands; the historical 3-versus-4-month split is not the primary definition. |
| Required: missingness and adjustment | Selection patterns, HIV censoring/absent series, imputation uncertainty and reduced adjustment blocks | Other covariates remain complete-case in primary. Adding vaccines, wasting, WASH or intervention measures requires a separate causal rationale. |
| Required: design and dependence | Survey weights; common cluster/stratum bootstrap; residual regional and within-child dependence | Retain unweighted primary point fit as reference; do not substitute rounded effective counts without a separate model definition. |
| Required: dose response / reference | Basis dimensions, linear comparison, 0% versus 1%, exposure restrictions | Match samples; report all supported contrasts and low-PfPR sensitivity. |
| Secondary: heterogeneity | Country omission, calendar interaction, within/between exposure and covariate timing | Separate, declared extensions of the primary structure. |
| Secondary: burden assumptions | Geographic coverage, population weights, national-mean versus subnational integration, age allocation and denominator sources | IHME all-cause is current baseline; IGME/WHO comparisons are supplementary alternatives. |
| Optional independent checks | Measured parasitaemia, fieldwork season, trial triangulation and alternative mortality estimators | Retain only for a stated scientific question; historical survey-region/NB/brms work remains archived. |

### 3.4.1 Key PfPR spline sensitivity analysis — current seven-band CBH model

<a id="key-pfpr-spline-sensitivity-analysis-current-seven-band-cbh-model"></a>

Use the seven full-sample MAP gamma=2 fits as the primary reference. Hold the median child HIV imputation, covariate selection/scaling and entry-year exposure timing fixed. Re-estimate smoothing parameters in each new fit at **gamma=2**, retaining the basis types and dimensions; record subgroup-specific knots and support.

- **Joint versus separate:** compare with a joint model having age-specific PfPR curves and confounder coefficients but one shared calendar-year spline and the historical joint random-effect structure. The completed joint/separate comparison used gamma=1; a matched gamma=2 joint fit remains to be run. Do not relabel the old comparison as gamma=2.
- **Geography:** fit each age separately in UN M49 Western Africa and Eastern plus Middle (Central) Africa. Exclude Southern Africa only for this comparison: Namibia, Eswatini and South Africa. The current groups contain 13 and 18 countries. Each subgroup-age fit has its own time curve and random-effect variances. Separate-age gamma=2 subgroup fits remain planned.
- **Survey period:** fit each age separately before/after the median survey year, counting each survey once and assigning ties to the early group. Current median: 2012; 53 surveys in 2003–2012 and 52 in 2013–2024. Keep whole surveys together and all countries. Exposure years may overlap. Separate-age gamma=2 period fits remain planned.

Overlay `f_g(P)-f_g(20)` using common axes, pointwise intervals and each fit's central exposure range. Report 40%→20% contrasts alongside the full curves. The earlier [key figure](<../Key results/pfpr_spline_sensitivity_by_age_geography_period.png>) and [report](../results/cbh/age_band_hiv_incidence_shared_time_v3/sensitivity_single_imputation/REPORT.md) retain their original gamma=1, joint-reference labels. Their geographic/period findings motivate these checks but are not results of the planned primary-compatible refits. Historical joint-reference and West Africa fits had small negative smoothing-Hessian eigenvalues; those checks remain relevant to their interpretation.

### 3.4.2 Supplementary Snow prevalence comparison

Use the supplied **annual re-fit** of the Snow survey database, not the published five-year estimates: 520 polygons for 2000–2015, microscopy-equivalent PfPR₂–₁₀ posterior means. [Script 51](../R_dhs/51_snow_polygon_prevalence.R) and [the extraction audit](../R_cbh/snow/00_audit_extraction.R) document same-country polygon overlays, exact overlap and GPW 2020 population weighting. Join by reviewed region key and band-entry year. Exclude post-2015 surveys/interviews, bands finishing after 2015 and regional Snow population coverage below 50%; do not carry 2015 exposure forward.

The central supplementary comparison fits MAP and Snow on identical records at gamma=2, with the same covariates, offsets and scaling. Retain source-specific PfPR knots. The matched sample has **3,428,544 records and 55,254 deaths**. Full Snow retains **3,429,310 records and 55,257 deaths**; its 766 additional records lack MAP. Comparing full-period MAP with full Snow changes both source and sample/period, and must be labeled accordingly.

Report curves, EDF, supported 40%→20% and current-to-zero contrasts, diagnostics and national-mortality comparisons. All Snow comparisons and penalisation trials belong in supplementary results. Source uncertainty remains fixed. Retain the supplied MCMC/sparse-year limitations: weak information in 2014–2015, no transmission-limits mask, maximum R-hat 1.0595 and 697 maximum-tree-depth events. Do not average marginal polygon quantiles to create region posterior intervals.

For the completed 2005/2015/2024 mortality comparison, both fitted relationships were evaluated at the **same national MAP prevalence**. The Snow-fitted 2024 result is an explicit transport scenario using MAP exposure, not a Snow prevalence estimate for 2024. This requires interchangeability of the exposure scales as well as temporal transport. These Snow scenarios are supplementary and do not replace primary MAP burden estimates. [Supplementary figures and tables](<../Supplementary results/README.md>).

### 3.5 Secondary national burden from the primary MAP model

For **2005, 2015 and 2024**, use the full-sample MAP gamma=2 curves, annual population-weighted national MAP PfPR and the finer-age **IHME all-cause deaths/rates** export dated 9 September 2026. National population weights use GPW 2020 density × cell area and fractional country overlap. Keep the fixed population-weight year, partial coverage and absent exposure explicit; covered-area prevalence applied nationally remains an assumption. Do not fabricate zero prevalence for missing countries.

For country c and band g:

\[
HR_{cg}^{0}=\exp\{f_g(0)-f_g(P_{c,y})\},\quad
m_{cg}^{0}=m_{cg}^{IHME}HR_{cg}^{0},\quad
m_{cg}^{attr}=m_{cg}^{IHME}(1-HR_{cg}^{0}).
\]

Apply the same fraction to IHME annual death counts, holding annual person-time fixed. These are calibrated reductions in all-cause mortality under zero PfPR. Baseline all-cause mortality is supplied by IHME, not independently predicted by the DHS model.

Combine early/late neonatal inputs for the <1-month model effect, recording the day/month boundary approximation. Use direct matches for 1–5, 6–11 and 12–23 months. For ages 2–4, retain the user-selected **same IHME rate for all three bands and equal division of deaths/person-time**. Check that age inputs sum to under-five totals and that observed = attributable + counterfactual rates/counts. The former IGME-primary, regional death-share and 2025 carry-forward specifications are historical alternatives, not this calculation.

Use the primary-only [country/age estimates](../results/cbh/primary_map_gamma2_v1/burden/country_age_estimates.csv), [country totals](../results/cbh/primary_map_gamma2_v1/burden/country_totals.csv) and [year totals](../results/cbh/primary_map_gamma2_v1/burden/year_summary.csv). These are calculated by [primary/02_effects.R](../R_cbh/primary/02_effects.R) from the newly fitted primary splines and the existing national MAP/IHME baseline inputs. The aggregate primary estimates across 42 estimable countries are **973,148 deaths in 2005, 694,534 in 2015 and 592,353 in 2024**. All 45 country rows remain in outputs, with missing estimates for Cape Verde, Lesotho and São Tomé and Príncipe. The previous mixed-comparison output directory is preserved; [country/year differences](../results/cbh/primary_map_gamma2_v1/burden/comparison_with_previous_primary.csv) document numerical reproduction.

Retain signed contributions, zero/current-prevalence support flags and conditional age-band intervals. National totals remain point estimates because cross-age sampling covariance is not estimated. Source, design, HIV-imputation, smoothing, transport and age-allocation uncertainty are not fully propagated. A nonlinear curve evaluated at national mean prevalence does not equal population-weighted subnational burden. Compare IHME malaria deaths only with explicit age/year/geographic alignment: cause-specific malaria deaths and the modeled reduction in all-cause mortality are different estimands, and export-release alignment remains unverified.

## 4. Stage 3 — plot saved results and separate primary from supplementary reporting

| Reporting location | Content | Authoritative selection |
|---|---|---|
| Main results | Study flow, survey/sample coverage and seven age-band descriptive summaries | Full MAP sample and its exclusion ledger |
| Main results: Figure 1 | Survey map and timing for the 105 included surveys in 34 countries | `survey_map/survey_map_and_timing.png` under current primary results |
| Main results | Seven primary PfPR curves, X%→Y% contrasts and diagnostics | `map_full`, gamma=2, from `results/cbh/primary_map_gamma2_v1/` |
| Main results: secondary burden | Country estimates and age contributions for 2005/2015/2024 using the primary model | Primary-only burden files linked above |
| Supplementary results | Matched MAP–Snow curves and contrasts; available-sample comparisons; all Snow smoothing trials and Snow mortality scenarios | [Supplementary index](<../Supplementary results/README.md>) |
| Supplementary sensitivity results | Structure, geography, survey period, timing, design, missingness and other declared checks | Label old gamma=1 exploratory results separately from planned gamma=2 refits |
| Historical archive | Six-band person-time/count models, old joint fits and earlier burden outputs | [Archived plan](ANALYSIS_PLAN_BEFORE_MAP_GAMMA2.md), original result directories and code audit |

The paper's Figure 2 is the primary PfPR-by-age curve figure; Figure 3 is the country comparison with IHME for 2005/2015/2024, using true base-10 logarithmic scales on both axes with the same limits. All 42 estimable countries per year have positive values and remain visible. The three paper figures omit titles, subtitles and embedded captions; separate caption text and filename mappings are in the [paper figure manifest](../results/cbh/primary_map_gamma2_v1/paper_figures/CAPTIONS.md). Preparing/copying figures does not edit manuscript TeX.

Plotting-only stages must read saved aggregate estimates and must not refit, change exposure definitions or select models according to available cache files. Label source, sample period, gamma, age bands, contrast and uncertainty. A combined MAP–Snow figure is supplementary; a main-results figure must select only the full MAP gamma=2 series. Keep links to authoritative saved artifacts instead of copying figures into competing locations.

## 5. Implementation status and next steps

**Primary-only rerun:** [R_cbh/primary/](../R_cbh/primary/README.md) now separates fresh fitting, effect/burden calculation, aggregate diagnostic checks and MAP-only plotting. It regenerates the study inclusion flow and key-results index. The prepared dataset, HIV imputation, national exposure inputs and IHME inputs remain fixed; supplementary fits are not run. New outputs are saved separately from the preceding comparison workflow, with numerical differences documented in the [rerun report](../results/cbh/primary_map_gamma2_v1/REPORT.md).

**15 September audit and cleanup:** the [component-by-component audit](CODE_AUDIT.md) verifies that the saved primary specification and data match this plan, while distinguishing unimplemented sensitivities and unresolved upstream issues. Sixty-eight superseded scripts were moved to [the archive](../archive/2026-09-15-code-audit/), along with preserved copies of three replaced entry points. Existing data, fitted objects and numerical results were not changed.

`Rscript run_all.R --primary` forces fresh primary fitting and regenerates all implemented primary outputs without data setup. `Rscript R_cbh/primary/run.R --resume` reuses only valid completed fit caches; `--report-only` uses saved fits for effects, diagnostics and plotting. `Rscript run_all.R --audit` retains verification of the preserved pre-rerun primary snapshot; `--check-inputs` checks local dataset prerequisites. The [results index](<../Key results/README.md>) links the new authoritative results.

**Remaining command-routing limitation:** older [sensitivity fitting](../R_cbh/sensitivity/01_fit.R) and [burden defaults](../R_cbh/burden/settings.R) still refer to gamma=1-era fits. `analysis/01_fit_complete_case.R` should be used with `--prepare-only` for dataset preparation; its default fit is historical. These scripts remain because they provide reference frames/knots or input tables used by the selected workflow. Do not use those defaults as a gamma=2 primary run. The new primary-only pipeline bypasses these defaults. National burden still reads explicit baseline columns from the earlier saved national input tables; consolidating their upstream producers remains a data-stage task.

Retain the three-stage workflow: **make datasets → fit declared models and calculate effects → plot saved results**. Priorities are:

1. Preserve the current primary/supplementary separation and input provenance. Primary-only routing and plotting are now implemented; in-sample fitted-outcome checks supplement the numerical diagnostics, but are not cross-validation or influence analysis.
2. Complete primary-compatible gamma=2 joint, geographic and period sensitivity fits; preserve historical comparisons with their actual settings.
3. Resolve the [open data and inference issues](OPEN_ANALYSIS_ISSUES.md), especially regional exposure coverage, missing-data selection, calendar/age coding, imputation uncertainty and survey/within-child dependence.
4. Consolidate the independent data build and prediction interfaces. Rebuild from a declared input snapshot and verify cache hashes. Fresh primary fitting now uses prepared data and committed reference knots; baseline national input tables still need an independent producer interface.
5. Finish separating input-only functions from the retained legacy bridges, then archive their remaining surplus branches. The 15 September cleanup removed the superseded independent modelling/plotting scripts; still-needed producers were retained with explicit reasons. Keep historical code/results traceable and do not execute legacy external-copy actions during migration.

Completion requires auditable survey accounting, validated child-band eligibility and joins, consistent primary specification and artifact selection, reviewed diagnostics, declared uncertainty and reproducible main/supplementary figures. Historical six-band/count-model choices are preserved in the archive and must not override this plan.

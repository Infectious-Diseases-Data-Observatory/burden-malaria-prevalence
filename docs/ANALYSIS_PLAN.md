# Analysis plan: malaria prevalence and child mortality

Draft for review · 7 September 2026 · Based on repository commit `17828d1`

**Subsequent model trial:** the user has requested Burstein's seven age bands and exclusion of maternal death fraction. The [proposed joint PfPR age-band model](PFPR_AGE_BAND_MODEL.md) specifies this trial, including a complementary log–log rate parameterization and age-specific confounder effects for review. The six-band separate-fit specification below records the earlier primary plan and is not the specification for this new trial.

**Use the newer person-time analysis as the primary analysis**, as confirmed during this review. Build one pipeline with three stages: **make analysis datasets → fit primary and sensitivity analyses → plot saved results**. The survey-region mortality approach becomes a limited sensitivity and migration reference. The earlier three-component pipeline in `R/` should leave the active workflow once its remaining data-building functions have been transferred.

This is a plan, not a refactor or a new results report. The review covered the structure, dependencies and analytical specifications of all 77 R scripts in `R/` and `R_dhs/`, plus the root runner, with detailed inspection of the person-time and shared data-processing code. All 78 scripts parsed under R 4.6.0. No raw DHS records or fitted-model objects were opened, and no analysis was rerun. Current sample sizes, model diagnostics and numerical results therefore remain unverified. The [code audit and complete script disposition](CODE_AUDIT.md) provide the supporting evidence.

## 1. Scientific questions and reporting choices

The primary causal question is: **if PfPR₂–₁₀ changes from X% to Y%, how does all-cause mortality in age group gᵢ change?** Estimate this effect separately for each of six age groups. The exposure is MAP *P. falciparum* prevalence standardised to ages 2–10 years (PfPR₂–₁₀), expressed as a percentage. With the person-time design, the outcome is the all-cause mortality rate within the age group.

Fit the data for each age group separately, following the most recent person-time models. Each group has its own prevalence–mortality curve, covariate effects and other model parameters. All six age-specific effects, including the effect below one month of age, are primary quantities of interest. Report the relative and absolute changes in mortality for a specified X% → Y% prevalence contrast, evaluated in the same target population under both scenarios.

The aim is causal estimation from observational data. Interpreting the fitted X% → Y% contrasts causally requires a defined exposure-change scenario, adequate control of confounding and appropriate exposure measurement. National attributable deaths are a **secondary extrapolation** of these age-specific effects and additionally require transportability to the populations receiving the estimates. Agreement with IHME, WHO or trial estimates can inform assessment of the results but does not establish these assumptions.

| Decision | Proposed specification | Relation to existing code |
|---|---|---|
| Primary design | Deaths and person-months by survey × region × retrospective window × age segment | Retain `40–43` |
| Survey frame | SSA DHS/MIS surveys with usable full birth histories, registry years 2000–2025; include supported retrospective windows independently of survey-year MAP availability | Broaden the current dependency on script `02` |
| Follow-up | Five non-overlapping 12-month windows before each child's interview; exclude the incomplete interview month | Retain `40` |
| Primary age groups | <1, 1–3, 4–11, 12–23, 24–35, 36–59 completed months | Confirmed 4-month boundary; retain `AGE6B`, used by the latest separate age-group models |
| Exposure | Person-time-weighted annual MAP prevalence for the years within each window, using weights specific to the age segment | Refine `40–42`, which currently use pooled under-5 year weights |
| Prevalence eligibility | Retain all finite prevalence values, including zero; no country or region prevalence floor | Retain person-time intent; remove inherited survey-region restrictions |
| Primary contrast | Change PfPR from any specified X% to Y%, estimating the mortality change separately in each age group | Generalise the existing prevalence-response predictions; 0% is a special reference for attributable-fraction summaries |
| Main model structure | Separate negative-binomial additive models by age band; smooth prevalence and year, segment and window effects, country and survey intercepts, regularised covariates | Retain `42–43` |
| Reporting engine | Recommend `brms` for final primary curves and burden draws, with `mgcv` for replication and efficient sensitivity screening | Makes the latest `43/48–50` reporting choice explicit; currently engine selection partly depends on cache availability |
| Historical burden | 2005 through the last year with the required observed inputs; currently configured as 2024 | Remove unlabelled carry-forward from the main series; later-year scenarios separate |

The person-time design, the age-specific X% → Y% causal question, the six age groups with a 4-month boundary, and separate fitting by age group have been confirmed by the user. The remaining choices are proposed defaults for the streamlined implementation. In particular, the final engine, count-weighting approximation, covariate block and missing-data rules should be settled in the specification before refitting. This is a revised analysis plan informed by existing exploratory work, not a claim of prospective preregistration.

## 2. Stage 1 — make the analysis datasets

### 2.1 Input manifest and survey registry

Separate optional input retrieval from processing local inputs. Record source, release, file checksum, available years, units and geographic coverage for every input. The normal data-build command should operate on a declared input snapshot; it should not refresh online sources because a cache happens to be absent.

The registry should retain every in-scope survey and distinguish unavailable recodes, unavailable boundaries, unresolved geography, unusable dates and missing exposure. Choose recode versions explicitly rather than taking the first duplicate survey key. Retain the recode variable dictionary and label-decoding information needed to handle different DHS phases.

Core inputs are authorised Births Recodes, survey boundaries, annual MAP rasters, the GPW population surface, national UNICEF/WUENIC and UNAIDS files, World Bank series and DHS StatCompiler aggregates. PR/HR/KR/IR recodes are optional sources for biomarker validation or indicators that the Births Recode cannot represent adequately. Their absence should not block the core mortality analysis.

A 2025 survey may contain eligible exposure and mortality during 2020–2024 even if there is no 2025 MAP raster. Conversely, older retrospective windows may predate MAP coverage. Determine eligibility at the region-window/segment level after constructing follow-up, not by requiring MAP for the registry's survey year first.

### 2.2 Geographic harmonisation

Use one geography module for mortality, covariates and MAP extraction. Retain the current encoding repairs, synonyms, language normalisation and documented aggregation from finer recode units to coarser survey boundaries. Save a reusable crosswalk with the survey ID, original region variable and label, boundary key, match method and any donor survey used for aggregation.

Require one-to-one matches unless an explicit aggregation rule applies. Review heuristic and elimination matches; matching the last remaining pair is not itself evidence that the geographies are equivalent. Record unresolved regions rather than silently dropping them. Qualify keys by survey/boundary version: a normalised region name alone does not establish a stable longitudinal region.

### 2.3 Deaths and time at risk

Retain the nine primitive age segments already implemented, in months:

`[0,1), [1,3), [3,4), [4,6), [6,12), [12,24), [24,36), [36,48), [48,60)`.

Assign these segments to the six primary groups `g₁` through `g₆`: `[0,1), [1,4), [4,12), [12,24), [24,36), [36,60)` months, retaining segment rows within each group for fitting. These intervals are disjoint and correspond to the completed-month labels in Section 1.

Construct follow-up from birth date, interview date, survival status and recorded/imputed age at death. Split each child's contribution across the five retrospective windows and the age segments. Sum both weighted and unweighted deaths and person-months. Keep zero-death cells with positive time at risk; distinguish these from cells with no observed follow-up.

Before accepting the build, check survival-status/age-at-death consistency, missing and special codes, deaths before birth or after interview, age limits, survey calendar conventions and boundary-month assignment. Replace the generic “7–9 year gap implies Ethiopian calendar” heuristic with survey-specific documented calendar metadata. The current convention assigns a full month at risk in the death month; retain it for reproduction, then assess fractional-month exposure and available day-of-death information, particularly for infancy.

Save calendar-year person-time contributions **by segment as well as region and window**. Retain survey, cluster and stratum identifiers in an ignored internal aggregate for design-based checks. The current cluster sums omit strata and are not used to adjust the model uncertainty.

Create the age-at-death descriptive table during this same pass through the recodes. Produce both the full eligible-history summary and the exact fitted-sample summary so they can be reconciled. Script `44` currently reads the recodes separately and can count deaths outside the model's matched regions or retained windows.

Validate aggregation with synthetic birth histories, conservation of deaths and exposure across cells, and comparisons against DHS life-table estimates for a deliberately varied set of surveys. A piecewise-exponential rate model and a synthetic-cohort probability estimator need not agree exactly. The DHS documentation identifies interview date, birth date, survival, age at death and sampling weight as mortality inputs; its probability estimator should remain a reference rather than be relabelled as a person-time model. [DHS Guide to Statistics](https://www.dhsprogram.com/Data/Guide-to-DHS-Statistics/Early_Childhood_Mortality.htm).

### 2.4 MAP exposure

Consolidate survey-year, lagged, window-year and burden extraction into one function accepting a boundary and a vector of years. Record the MAP release and layer, polygon coverage, missing-pixel coverage and population denominator. Verify that population-density weighting includes the appropriate cell-area contribution where the raster grid requires it; reproducing a second implementation of the same formula is not an independent check of the weighting.

For each age segment in each region-window, construct `PfPR = sum(year_share × annual_PfPR)` over supported years. Retain the contributing years, weights and exposure lag. Use lag 0 for the primary model and lag 1 for the planned sensitivity.

The primary dataset should require complete MAP coverage for the included cell's exposure years. Preserve the current boundary-year borrowing rule as an explicitly labelled sensitivity, including the fraction of person-time borrowed and source year. Never silently clamp 1999 to 2000 or 2025 to 2024. Do not remove available years simply because a newer survey-year surface is missing.

### 2.5 Covariates and external joins

Build a survey-region covariate table independently of whether a life-table mortality estimate can be produced. Build national country-year panels independently of the DHS outcomes. Join national annual covariates to the **calendar years represented by person-time**, then aggregate on the declared natural scale and transform for modelling. Preserve the current survey-year joins as a sensitivity.

Survey-derived covariates generally describe interview-time conditions or their own retrospective eligibility windows. Keep their measurement period explicit when reusing them across the five mortality windows; they should not be presented as reconstructed historical regional measurements.

| Candidate block | Existing variables | Processing and proposed treatment |
|---|---|---|
| Settlement and resources | `pct_urban`, `educ_yrs`, `wealth_q`, `elec_dhs` | Audit denominators, deduplication and code/label decoding; distinguish household from mothers-with-birth-history populations |
| Routine vaccination | `dtp3_reg`, `measles`, `pentavalent3_reg`, `pcv3_reg`, `rotavirus_complete_reg` | Preserve eligible living-child denominators and survey-specific schedules; sparse fields remain available for restricted sensitivities |
| Maternal/birth characteristics | `facility`, `birth_int`, `mage1` | State whether the denominator is recent births, all recorded births or unique mothers; harmonise intentionally |
| Feeding and nutrition | `excl_bf`, `stunting`, `underweight`, `wasting` | Audit phase-specific availability, anthropometric standards and valid codes; use the documented StatCompiler wasting measure where available |
| Water and sanitation | `imp_water`, `imp_sanit` | Prefer the corrected StatCompiler regional definitions; retain national fallback and recode estimates under separate source flags |
| National immunisation | `hib3_wuenic`, `pcv3_wuenic`, `rotac_wuenic` | Join country-year WUENIC; distinguish reported zero, documented pre-introduction, absent series and internal gaps |
| Paediatric HIV | `log_hiv_prev` | Derive 0–14 prevalence from UNAIDS PLHIV counts / World Bank 0–14 population; preserve numerator, denominator, release, year and censoring flags |
| National context | `log_gdp`, `log_hexp_pc`, `polstab` | Retain the existing candidate measures; use exact-year values where possible and an explicit, bounded fallback with source year recorded |

This accounts for all 25 current candidates. Eligibility by missingness is not a sufficient causal rationale for adjustment. Freeze the intended confounder block before examining revised malaria effects; identify nutrition and intervention measures whose adjustment may change the interpretation. Do not select covariates by significance. ITN/ACT targeting and any future SMC series belong in explicitly labelled secondary analyses rather than being added automatically.

For HIV, keep the existing log-scale construction and record censored values such as `<500` separately from exact counts. The existing midpoint substitution is an assumption to test. An unavailable child series must remain missing rather than be replaced by an adolescent series. Retain lower/upper estimates when available for sensitivity work.

For vaccines, lack of a WUENIC series is not by itself proof of non-introduction. Use zero where reported or supported by introduction metadata; otherwise classify it as unknown. Quantify the effect of reproducing the existing assumption that every pre-series or absent-series year is zero.

### 2.6 Missingness and preprocessing

Preserve the current ≤5% missingness rule as the starting general rule, evaluated before imputation on a clearly defined, unique survey-region base rather than after expansion into age/window rows. Also report missingness by country, survey, cell and person-time so this threshold cannot conceal systematic gaps.

Variables passing the rule can use the existing country-median, then overall-median imputation for a reproducible baseline. Keep original values, analysis values, imputation method and geographic fallback flags separately. Estimate transformations and standardisation on the primary training sample and save them as a preprocessing object independent of any fitted mortality model.

Wasting is currently a substantive exception: a country/year/region model fills missing values **before** the generic missingness rule is applied. Preserve and report its original missingness. Treat this as a declared model-based imputation choice, compare with excluding wasting and with observed-only data, and assess uncertainty rather than presenting it as ordinary ≤5% median imputation.

For grouped cross-validation, refit imputation models, medians, scaling and any data-driven covariate eligibility within training folds. Do not borrow held-out survey information through the existing complete-panel preprocessing.

### 2.7 Analysis-ready outputs

| Output | Unit/key | Purpose |
|---|---|---|
| `survey_registry` and `region_crosswalk` | Survey; survey × source region | Input accounting and common geography |
| `person_time_cells` | Survey × region × window × primitive segment | Weighted/unweighted deaths and exposure; demographic validity flags |
| `person_time_year_weights` | Same cell × calendar year | Temporal matching for prevalence and national covariates |
| `region_covariates` and `national_covariates` | Survey × region; country × year | Reusable covariates with denominators and source flags |
| `analysis_person_time` | Survey × region × window × segment | Main model inputs, lag alternatives, age-band mappings and inclusion flags |
| `preprocessing` and `analysis_specifications` | Run × specification | Covariate transformations, scaling, model settings and sample definitions |
| `analysis_region_mortality` | Survey × region × horizon | Optional 12-/60-month life-table comparison, generated only when requested |
| `analysis_burden` | Country × region × year × age band × denominator source | National deaths, allocation shares, MAP and provenance for extrapolation |
| `analysis_measured_prevalence` | Survey × region | Optional biomarker/seasonality comparison |

Use RDS for typed analysis objects and CSV for inspectable aggregate tables where appropriate. Keep raw records, cluster-level material and fitted models containing analysis data under ignored `data/`. Every analysis dataset needs a dictionary, unique-key check, inclusion ledger and input manifest. None should require a previously fitted outcome model or a table from `results/` to be built.

## 3. Stage 2 — primary and sensitivity analyses

### 3.1 Primary working model

Fit a separate model to the data in each of the six age groups, using the following structure:

```text
deaths ~ smooth(PfPR, k=5) + segment + window
       + smooth(calendar_year, k=8) + regularised_covariates
       + country_intercept + survey_intercept
       + offset(log(person_months))
```

The segment term is omitted for a one-segment band. Each age band has its own nuisance effects, smooths and negative-binomial dispersion. The primary prevalence association is constant across calendar time, with the smooth year term adjusting the baseline; time-varying prevalence effects are a planned sensitivity.

**The current outcome is a constructed count**, `round(weighted_deaths / weighted_person_months × unweighted_person_months)`, with unweighted person-months as the offset. It is not the observed death count, and unweighted exposure is not a survey-design effective sample size. Preserve this approximation for reproduction, quantify its rounding effects and compare it with unweighted observed-count fits and a survey-design-aware analysis before treating its intervals as final. The minimum design check should resample clusters within survey strata, preserve all contributions from each sampled cluster across ages/windows, and refit the complete age-band set using a common replicate. Model-based country/survey intercepts do not replace this check.

Recommend one named `brms` primary fit per band for the final figures and burden predictions. Retain the current priors as a reproducibility baseline, document them explicitly and conduct prior predictive/sensitivity checks. The existing Normal(0,1) coefficient priors are **not identical** to `mgcv`'s empirically estimated ridge penalty; engine comparisons test regularisation as well as computation. Use `mgcv` as a quick reproducibility and sensitivity engine, reporting material discrepancies and refitting consequential alternatives in the reporting engine.

### 3.2 Effect summaries

Let `rate_i(p)` denote the predicted all-cause mortality rate in age group `gᵢ` under prevalence `p`, evaluated for the same target population and adjustment settings. For any specified change from `X%` to `Y%`, use a single prediction routine to calculate:

```text
Rate_ratio_i(X → Y) = rate_i(Y) / rate_i(X)
Percentage_change_i(X → Y) = 100 × [Rate_ratio_i(X → Y) - 1]
Absolute_change_i(X → Y) = rate_i(Y) - rate_i(X)
```

Report each contrast for all six groups, with uncertainty. Express absolute changes in deaths per 1,000 child-years. A negative percentage or absolute change denotes a reduction in mortality. Because the prevalence-response curve is nonlinear, both X and Y must be specified; a single per-10-percentage-point coefficient does not describe all possible contrasts.

Attributable fractions and attributable rates are additional summaries of the same fitted curves. For prevalence `p` and reference `p0 = 0`, calculate:

```text
HR_age(p, p0) = rate_age(p) / rate_age(p0)
AF_age(p, p0) = 1 - rate_age(p0) / rate_age(p)
Attributable_rate_age(p, p0) = rate_age(p) - rate_age(p0)
```

Report the age-specific curves and X% → Y% contrasts with their exposure support, flagging unsupported extrapolation. Retain attributable-fraction anchors at 10%, 30% and 50% as additional summaries. For those summaries, keep the 0% reference visible and compare with 1% on the same sample.

Retain signed contrasts and intervals, including negative values. A non-negative burden scenario, if wanted, must have an explicit label and be compared with the untruncated calculation. Current effect tables and burden routines apply different truncation rules.

Specify the prediction population. Setting random effects to zero and covariates to their means defines a reference cell; it is not automatically a marginal population-average mortality rate. For descriptive absolute rates, standardise over an explicit covariate/segment distribution. The random intercept and covariate contributions cancel from HRs within the present additive, no-interaction specification, but not necessarily from future interactions or marginal absolute rates.

Combine age-band and geographic estimates at the **draw or bootstrap-replicate level**, then summarise the total. Do not average separate lower and upper confidence limits. Separate Bayesian band models omit shared sampling dependence between ages; use the common design replicates or a suitable joint model to assess the effect on aggregate uncertainty.

### 3.3 Diagnostics and model comparison

Check convergence, dispersion, residual patterns, zero counts, fitted rates by age/window, exposure-response support, influential surveys/countries and sensitivity to spline dimensions. Examine whether repeated observations from the same survey-region need an additional intercept or correlated structure. Add a within-region versus between-region exposure analysis: the current within/between variables are descriptive and do not identify a within-region effect in the fitted model.

For Bayesian reporting, require reviewed sampling diagnostics: R-hat below 1.01, adequate bulk/tail effective sample sizes and Monte Carlo precision for the reported contrasts, and investigation/resolution of divergent transitions. Record tree-depth and energy diagnostics. Merely writing these values into a CSV should not make a fit eligible for final reporting. [Stan diagnostic guidance](https://mc-stan.org/learn-stan/diagnostics-warnings.html).

Use survey-grouped cross-validation when comparing flexible structures; all regions, age segments and windows from a survey stay in one fold. Include country-held-out checks where transportability is the question. Pointwise leave-one-cell-out prediction is insufficient for either task. Use a fixed, small model menu and common comparison samples; keep the additive primary specification unless the plan is explicitly revised. Choose lookback periods and covariate sets independently of the estimated effect direction or significance in any age group.

For `mgcv` uncertainty, evaluate the smoothing-uncertainty-corrected covariance when available and record when it is unavailable. The package supports this option; moving to Stan is not the only possible improvement over conditional covariance estimates. [mgcv prediction documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/predict.gam.html).

### 3.4 Planned sensitivities

Each sensitivity should have a specification ID, one stated change, declared sample, saved estimates and diagnostics. Report the primary fit on the corresponding restricted sample where needed to distinguish a specification effect from a sample-composition effect. Avoid automatically crossing every choice with every other choice.

| Priority | Question | Prespecified comparisons | Existing code to reuse |
|---|---|---|---|
| Required | Recall window | All five windows vs windows 1–4 vs window 1 only | `42` |
| Required | Exposure timing | Lag 0 vs lag 1; complete MAP coverage vs current boundary-year borrowing; segment-specific vs pooled under-5 year weights | `40–42`; ideas from `17/27` |
| Required | Covariate timing | Person-time-year national values vs survey-year values | External joins in `03`, currently inherited by `42` |
| Required | Count/weighting assumptions | Current effective counts vs observed unweighted counts; rounding impact; cluster/stratum resampling | `40/42`, validation ideas from `27` |
| Required | Age definition and recorded timing | Main 4-month split vs 3-month split; broader four/five groups; death-month exposure and age-heaping checks | `40/42/44/47` |
| Required | Likelihood and dependence | NB vs Poisson; add survey-region intercept; no survey intercept; country prevalence heterogeneity | `42/47`, selected ideas from `07/16` |
| Required | Missing data and source substitution | Complete cases; exclude model-imputed wasting; regional vs national WASH; vaccine zero assumptions; HIV censoring and absent-series handling | `03/06`, covariate audit |
| Required | Dose-response and reference | Spline basis alternatives, linear comparison, 0% vs 1% reference; separate ≥1% and 5–40% restrictions | `06/07/42` |
| Required | Uncertainty/regularisation | `brms` vs `mgcv`; prior sensitivity; compare conditional and corrected/design-based intervals | `43/50` |
| Secondary | Heterogeneity and confounding | Calendar interaction; region/era strata; leave-one-country-out; within/between exposure decomposition; reduced adjustment blocks | Ideas from `19/20/23/29–35` |
| Secondary | Alternative mortality estimator | Synthetic-cohort post-neonatal/neonatal estimates over 12 and 60 months on comparable geography and supported exposure | `03/04/17/27` |
| Optional extension | Independent comparisons | MAP vs measured parasitaemia/fieldwork season; trial triangulation; IHME share analysis; ITN/ACT targeting | `12/13/21/22/24` |

Older-child mortality (5–14 years), neonatal mortality as a covariate and the extensive three-way region/time model ladder leave the default workflow. They can be revived for a specific scientific question; the person-time primary does not require them.

### 3.5 Secondary national burden

Consolidate Nigeria/DRC and all-SSA estimates into one country-parameterised calculation. Use the latest eligible survey's boundaries where available; report countries using national prevalence and pooled age shares separately. Construct the country universe from an explicit geographic/input manifest, not from an older model's results table.

Use UN IGME all-cause age-block death counts as the primary denominator and IHME all-cause counts as a denominator sensitivity. Harmonise country, year, sex and age definitions before comparison. Validate that age blocks sum to under-5 totals and distinguish missing components from true zero counts. The current “missing column becomes zero” and “negative remainder becomes zero” rules should not silently repair an inconsistent export.

For each draw, calculate `sum_region,age[AF(age, PfPR_region,year) × allocated_allcause_deaths]`. Retain the latest-survey regional and within-block age death shares as the reproduction baseline, but report their sampling uncertainty and the assumption that they are fixed through time. Test alternative geographic allocations and wider age blocks. A zero-death age band needs an explicit allocation rule that still conserves the national total.

Reuse each band's model draw across every country, region, year and denominator. Preserve the dependence induced by shared coefficients. State separately which uncertainty sources are included: model parameters, survey sampling, imputation, prevalence, death shares and external mortality inputs. The existing posterior bands include model uncertainty only.

Show neonatal, 1–11-month and 1–4-year contributions, plus post-neonatal and all-under-5 totals. Compare with IHME/WHO on matched age, year and geographic sets. A WHO under-5 proxy from an all-age total must remain labelled as a proxy. Report later-year carry-forward estimates only as scenarios with an input-year flag; do not merge them into the observed-input historical series.

## 4. Stage 3 — plot saved results

Plotting scripts should read saved descriptive tables, model contrasts and burden estimates. They should not open raw recodes, fetch data, fit models, choose an engine based on available cache files or recompute burden under a separate formula.

| Figure/table | Content | Source to consolidate |
|---|---|---|
| Study flow and survey coverage | Surveys, countries, geography, retained windows, deaths/person-time and reasons for exclusion; map/timeline | `11/26/28`, adapted to person-time |
| Age-specific data and dose-response | Descriptive rates alongside adjusted HR/AF curves for the six main bands, with exposure support | `42/45/50` |
| Primary results table | Rate ratios, percentage changes and absolute mortality-rate changes for specified X% → Y% contrasts in all six age groups; uncertainty method and model diagnostics; additional AF anchors at 10/30/50% | `43/50`, extending their prediction summaries |
| Sensitivity forest | The same X% → Y% contrasts by age group and specification, with sample sizes and required diagnostics; AF summaries where relevant | `15/42/43/47` |
| Recall and age-at-death diagnostics | Window effects, death-age distribution/heaping, full-frame vs fitted-sample accounting | `42/44` |
| Conditional covariate associations | Age-specific coefficients per SD on the declared transformed scale | `46`; labels from a dictionary, not another fitted model |
| Secondary burden | Country/SSA trends, age contributions, denominator comparison and Nigeria/DRC views | `48/49`; replace competing temporal figure calculations |

Each exported figure should have a saved plotting-data table and metadata identifying dataset version, specification, age definition, counterfactual, engine and interval type. Generate a results index with links to the authoritative files; avoid maintaining a second collection of copied figures with potentially stale filenames. Manuscript layouts can assemble these saved panels without repeating analysis.

## 5. Implementation sequence and completion criteria

Use three public entry points, backed by small shared function modules and one configuration/specification file. A proposed layout is:

```text
run_all.R                       explicit data / analyses / plots stages
R/01_make_analysis_data.R        local inputs → analysis-ready datasets and QC
R/02_fit_analyses.R              declared specifications → fits, contrasts, burden
R/03_plot_results.R              saved tables → figures and results index
R/lib/                          geography, person-time, external inputs,
                                covariates, models, predictions, plotting helpers
config/                         settings, covariate dictionary, sensitivity menu
data/                           ignored raw, derived and model artifacts
results/                        approved aggregate tables and figures
docs/                           plan, provenance, code audit and migration record
```

This is the intended future layout; neither current code folder has been moved. No new workflow package is necessary to establish these boundaries. Optional fetch, sensitivity, Bayesian and burden steps must declare their prerequisites and report skipped/failed status accurately.

1. **Freeze a reproduction baseline.** Record source versions, current specification and a small set of aggregate comparison outputs before numerical changes. Replace weak cache checks with hashes covering actual input values, crosswalks, configuration, formulas, priors, preprocessing and relevant code/package versions.
2. **Extract the data layer.** Decouple covariates/geography from life-table fitting; combine the overlapping raster extractors; move construction currently embedded in `42` into the data stage. Prove that the data stage runs without fitted models or legacy result tables.
3. **Consolidate fitting and prediction.** Introduce one specification registry, one model-result format and one HR/AF/burden prediction interface. Reproduce the existing person-time specification first, then implement the analytical changes separately with a numerical change log.
4. **Consolidate figures.** Rebuild the declared figure/table set solely from saved results; verify labels, sample counts, age units, counterfactuals and uncertainty provenance. A plotting-only run must leave fits and data unchanged.
5. **Retire surplus code.** Use the disposition table in the audit. Remove obsolete scripts from execution first; preserve their Git history and any documented reference snapshot. Delete only after all required outputs and dependencies have replacements. Do not execute legacy external-copy actions during migration.

The streamlined pipeline is ready when it can rebuild from the declared local input snapshot in a fresh derived-output directory; account for every survey and excluded cell; conserve deaths/exposure through aggregation; pass key/merge/calendar checks; reproduce or explain changes in the baseline results; enforce reviewed model diagnostics; and regenerate all declared plots without recomputation or stale artifacts. Current region-name tests should be retained, with targeted new tests for person-time accounting, temporal joins, cache invalidation and draw-level aggregation.

**Next implementation priority:** establish the independent person-time data layer and its validation before consolidating models or deleting scripts. That removes the central dependency tangle while making subsequent changes in estimates interpretable.

# Analysis plan: malaria prevalence and child mortality

Draft for review · 7 September 2026 · Based on repository commit `17828d1`

**Primary-model decision, 9 September 2026:** **fit seven separate age-band models as the main mortality analysis**, following review of the spline sensitivity results. Each band has its own PfPR curve, calendar-year curve, confounder coefficients and survey/country/region random effects and variance parameters. Section 3.1 specifies the current primary model. The joint model with one shared calendar-year spline is now a key sensitivity comparator. This decision supersedes the 8 September shared-time primary specification and the historical six-band negative-binomial plan.

**Key sensitivity update, 9 September 2026:** the joint-versus-separate, geographic and survey-period PfPR spline comparisons remain **required key sensitivity analyses**. Future geographic and period comparisons should use the separate-age primary structure within each subgroup. Section [3.4.1](#key-pfpr-spline-sensitivity-analysis-current-seven-band-cbh-model) distinguishes that plan from the completed exploratory joint-model subgroup fits; the original [combined figure](<../Key results/pfpr_spline_sensitivity_by_age_geography_period.png>) remains in Key results.

**Earlier model trial:** the user requested Burstein's seven age bands and exclusion of maternal death fraction. The [joint PfPR age-band specification](PFPR_AGE_BAND_MODEL.md) records the earlier trial and now serves as the joint sensitivity specification. The seven bands and complementary log–log hazard parameterization are retained in the separate-age primary model.

**Historical mortality model revision, 8 September 2026:** the joint seven-band trial used one shared calendar-year spline across all age bands, alongside age-specific PfPR splines and confounder coefficients. This remains the joint sensitivity specification. The new separate-age primary fits each estimate their own calendar-year spline. The HIV imputation model is unchanged.

**Country burden revision, 9 September 2026:** the [national burden stage](../R_cbh/burden/README.md) uses the newly supplied, finer-age IHME all-cause export for **2024**. It applies each country's population-weighted 2024 PfPR-to-zero hazard ratio to the corresponding age-band mortality rate and annual death count. The user selected the same IHME 2-4-year rate for ages two, three and four, with equal death/person-time shares, where the export remains aggregated. This stage reports missing/partial exposure coverage, signed attributable effects and model/HIV-imputation uncertainty; it does not use the old 2025 carry-forward scenario.

**Country burden implementation after the primary-model switch:** updated **2005, 2015 and 2024** estimates now use the seven separate age-band fits, holding one posterior-median HIV-incidence imputation fixed. [Updated country tables and comparisons](../results/cbh/age_band_separate_v1/primary_update_2005_2015_2024/REPORT.md) are stored under `results/cbh/age_band_separate_v1/`. The IHME inputs, national PfPR, fixed GPW 2020 population weights and age allocations are unchanged. Age-band intervals are conditional on the fixed imputation; national totals are point estimates. Previous joint-model results, which pooled ten HIV-imputation fits, remain in their original directory.

**Implementation update:** the new [R_cbh pipeline](../R_cbh/README.md) builds child–age-band data for that trial. It uses annual exposures/covariates at band entry, restricts entry to the five years before interview, and requires the full potential band to end by interview for deaths and survivors alike. The historical audit and six-band plan below describe the earlier review; they do not describe this new dataset build.

**Historical pipeline decision, 7 September 2026:** the newer person-time analysis was chosen over the survey-region mortality approach. The current primary design is the child–age-band CBH design described above. Retain the three-stage workflow: **make analysis datasets → fit primary and sensitivity analyses → plot saved results**. The survey-region mortality approach remains a limited sensitivity and migration reference.

The original 7 September review covered the structure, dependencies and analytical specifications of all 77 R scripts in `R/` and `R_dhs/`, plus the root runner, with detailed inspection of the person-time and shared data-processing code. All 78 scripts parsed under R 4.6.0. At that review stage, no raw DHS records or fitted-model objects were opened and no analysis was rerun, so numerical results were unverified. The dated implementation updates and the current CBH sensitivity results below record subsequent work. The [code audit and complete script disposition](CODE_AUDIT.md) provide the historical review's supporting evidence.

## 1. Scientific questions and reporting choices

The primary causal question is: **if PfPR₂–₁₀ changes from X% to Y%, how does all-cause mortality in age group gᵢ change?** Estimate this effect separately for seven age groups: <1, 1–5, 6–11, 12–23, 24–35, 36–47 and 48–59 completed months. The exposure is MAP *P. falciparum* prevalence standardised to ages 2–10 years (PfPR₂–₁₀), expressed as a percentage and assigned at band entry. The observed outcome is death during the band, conditional on reaching it alive; the complementary log–log model estimates an annualized mortality hazard.

Fit each age group separately using the child–age-band binomial model in section 3.1. Each group has its own prevalence–mortality curve, calendar-time curve, covariate effects and random-effect variance parameters. All seven age-specific effects, including the effect below one month of age, are primary quantities of interest. Report the relative and absolute changes in mortality for a specified X% → Y% prevalence contrast, evaluated in the same target population under both scenarios.

The aim is causal estimation from observational data. Interpreting the fitted X% → Y% contrasts causally requires a defined exposure-change scenario, adequate control of confounding and appropriate exposure measurement. National attributable deaths are a **secondary extrapolation** of these age-specific effects and additionally require transportability to the populations receiving the estimates. Agreement with IHME, WHO or trial estimates can inform assessment of the results but does not establish these assumptions.

The following decision table records the **historical 7 September person-time plan**. Its six-band age definitions, count likelihood and reporting-engine recommendations are superseded for the current CBH primary analysis by section 3.1; it is retained for migration and provenance.

| Decision | Historical proposed specification | Relation to existing code |
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
| Paediatric HIV | `log_hiv_incidence` (active CBH model); `log_hiv_prev` (legacy) | Use reported ages 0–14 incidence per 1,000 uninfected population, with a country/year hierarchical imputation model informed by adolescent incidence; preserve release, year, censoring and imputation flags |
| National context | `log_gdp`, `log_hexp_pc`, `polstab` | Retain the existing candidate measures; use exact-year values where possible and an explicit, bounded fallback with source year recorded |

This accounts for all 25 current candidates. Eligibility by missingness is not a sufficient causal rationale for adjustment. Freeze the intended confounder block before examining revised malaria effects; identify nutrition and intervention measures whose adjustment may change the interpretation. Do not select covariates by significance. ITN/ACT targeting and any future SMC series belong in explicitly labelled secondary analyses rather than being added automatically.

**HIV update, 8 September 2026:** the active CBH analysis uses log child HIV incidence, replacing prevalence at the user's request. Fit the child and adolescent incidence series jointly with separate calendar-year smooths, regional/country effects and within-country temporal correlation. Treat `<0.01` rates as censored, preserve numeric child estimates and predict missing child series. The resulting child rate—not the adolescent predictor—enters the mortality model at band-entry year. Retain posterior trajectories and imputation status; propagate external-covariate uncertainty through mortality fits. Countries without usable adolescent information remain missing. See [implementation and limitations](../R_cbh/hiv/README.md). The previous count-derived prevalence and midpoint assumptions remain part of the archived analysis only.

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

**Current primary specification — seven separate CBH age-band fits (9 September 2026).** Fit one unweighted binomial complementary log–log model per age band with `mgcv::bam`, discrete fitting and fREML smoothing selection. For each band g, use:

```r
death ~
  s(pfpr_pct, bs = "cr", k = 5) +
  s(calendar_year, bs = "cr", k = 6) +
  sex + multiple_birth + z_birth_order + z_maternal_age_birth +
  z_maternal_education_years + wealth_quintile + urban +
  z_log_hiv_incidence + z_log_gdp_pc +
  z_log_health_expenditure_pc + z_political_stability +
  s(survey, bs = "re") + s(country, bs = "re") +
  s(region, bs = "re") + offset(log(band_years))
```

The intercept is estimated separately in each fit. On the hazard scale, `log(lambda_ig) = alpha_g + f_g(P[r(i), entry_year_ig]) + h_g(t_ig) + X_ig beta_g + u_s(i),g + v_c(i),g + b_r(i),g`, and `q_ig = 1 - exp(-band_years_g * lambda_ig)`. Regional annual PfPR₂–₁₀ is the same exposure definition for every age group; values can differ across a child's bands because entry years differ.

Each fit independently estimates its PfPR and calendar-year spline, confounder coefficients and survey/country/region random intercepts. Within each fit, the random effects are independent zero-mean Gaussians with a separate variance for each grouping type. The three random-effect variances and all smoothing parameters are estimated separately by age. No calendar-year function, random intercept or variance parameter is shared across age bands. Region identifiers retain the current survey-specific coding. Country replaces the joint model's country-by-age term because each fit contains one age band. The unweighted likelihood is the current working choice; survey-design uncertainty remains an outstanding issue. Separate fitting does not imply independent sampling errors across age bands: uncertainty for totals across bands must account for shared children, surveys and HIV-imputation trajectories.

Retain the current seven age bands and their full predetermined widths, entry within five years before interview, complete-potential-band eligibility for deaths and survivors alike, annual PfPR and external covariates at band entry, the existing adjustment set and scaling, and exclusion of vaccines and maternal death fraction. HIV adjustment remains log child incidence from the existing imputation model; adolescent incidence is informative in imputation and is not a separate mortality predictor. The initial comparison uses one fixed posterior-median HIV imputation, with other covariates requiring complete cases. Final uncertainty reporting must distinguish these conditional fits from results propagating HIV-imputation and other relevant uncertainty.

**Rationale and model status.** The separate fits remove cross-age restrictions on baseline time trends and residual variation while preserving the broad PfPR curve shapes observed in the sensitivity comparison. Their promotion follows inspection of results and the user's explicit decision; it does not establish superior causal identification or a formal model-comparison result. Retain all seven bands and report their estimates regardless of direction or significance. The joint model with one shared calendar-year spline remains a required sensitivity analysis.

**Implementation transition.** The seven existing single-imputation fits are available as `age_1.rds` through `age_7.rds` under `data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/sensitivity_single_imputation/`, fitted by `R_cbh/sensitivity/01_fit.R`. They implement the selected primary structure. The country burden stage now uses these fits by default and has updated 2005, 2015 and 2024 estimates under the primary result ID `age_band_separate_v1`. Named primary fitting/curve outputs, propagation of HIV uncertainty through the separate fits and separate-age geographic/period refits remain subsequent implementation steps. Existing joint-model results, the original key comparison figure and earlier country burden estimates retain their original provenance.

#### Historical six-band person-time specification (superseded for the CBH primary)

The following model and reporting-engine recommendations record the earlier count-data analysis. Fit a separate model to each of the six historical age groups using:

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

| Priority | Question | Planned comparisons | Existing code to reuse |
|---|---|---|---|
| **Key required — current CBH model** | **Stability of age-specific PfPR curves** | **Separate-age primary vs joint sensitivity; West vs East/Central Africa and early vs late surveys using separate-age fits within each subgroup. Section 3.4.1 distinguishes planned fits from completed exploratory results.** | **`R_cbh/sensitivity/01_fit.R`, `02_report.R`; subgroup structure to be updated** |
| Required — current CBH model | Alternative prevalence source | Annual Snow re-fit posterior mean versus MAP, with seven separate age-band fits and MAP refitted on the same records through 2015. Hold the adjustment set, scaling and median child HIV imputation fixed; compare full spline curves and supported contrasts. | `R_dhs/51_snow_polygon_prevalence.R`; `R_cbh/snow/00_audit_extraction.R`, `01_prepare.R`, `02_fit.R`, `03_report.R` |
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

**Annual Snow exposure sensitivity.** Use the supplied annual re-fit of the Snow parasite-survey database (`data/snow_prevalence_model/`, script 51's `annual_csv` option), which provides microscopy-equivalent PfPR₂–₁₀ posterior means for 520 polygons in 2000–2015. This is distinct from the published five-year Snow estimates. Overlay same-country Snow polygons onto each survey's DHS regional boundaries using exact overlap and GPW 2020 population counts. Join the resulting annual means to the existing reviewed region keys and band-entry calendar year; never carry 2015 values into later years. The first implementation excludes surveys and individual interviews after 2015, requires full potential bands to finish by the end of 2015, and excludes Snow population coverage below 50%, following script 51's flag. Preserve all remaining eligibility rules. Refit MAP on the identical comparison records; retain the original full-period MAP curves as context. Snow fits need not require MAP availability; if the samples differ, add matched Snow fits. Report source uncertainty as fixed, distinguish model-based intervals from exposure credible intervals, and carry forward the supplied source MCMC and sparse-year limitations. The [implementation and results](../results/cbh/age_band_snow_2000_2015_v1/REPORT.md) include extraction checks, selection, diagnostics and curve comparisons.

**Snow penalisation sensitivity, 10 September 2026.** Compare two separate changes against the seven full-sample Snow fits: (1) retain PfPR `cr`, k=5, and increase `bam` gamma from 1 to 1.4; (2) change only PfPR to `cs`, k=5, with gamma=1. Retain calendar-year `cr`, k=6, original knots, identical records and adjustment set, and fixed median HIV imputation. Gamma affects smoothing selection throughout the model; cs adds shrinkage of the PfPR linear component as well as curvature. Neither imposes monotonicity. Compare full curves, EDF and supported 40%→20% contrasts with conditional intervals, without selecting the preferred result by significance or visual agreement. [Code](../R_cbh/snow/README.md) and [results](../results/cbh/age_band_snow_2000_2015_v1/penalty_sensitivity/REPORT.md) are separate from the original Snow and MAP results.

**Additional Snow gamma=2 sensitivity.** Extend the penalisation comparison with seven separate age-band fits using PfPR `cr`, k=5, calendar-year `cr`, k=6, and gamma=2. Reuse the exact full Snow sample, reference knots, adjustment set and fixed HIV imputation. Compare with gamma=1 and 1.4, retaining the separate cs analysis for context. The [extended comparison](../results/cbh/age_band_snow_2000_2015_v1/penalty_sensitivity/gamma2/comparison/REPORT.md) reports curves, EDF, supported contrasts and numerical diagnostics.

**MAP versus Snow with gamma=2.** Fit each exposure on identical eligible records through 2015 to isolate exposure-source differences, retaining each source's original PfPR knots and the same confounders, scaling, offsets, HIV imputation and cr basis dimensions. Also fit full-period MAP models at gamma=2 and retain full-sample Snow gamma=2 as context; that comparison additionally changes period and sample. Compare national attributable mortality for 2005, 2015 and 2024 by applying both fitted relationships to the same annual national MAP exposure and the existing IHME all-cause inputs. Snow has no 2024 exposure: this is an explicit cross-source exposure and temporal transport scenario, not a source-specific Snow burden estimate. Retain the earlier MAP gamma=1 estimates, signed effects, missing-country rows, age-allocation assumptions and support flags. Report conditional age-band intervals without inventing cross-age or cross-model sampling independence. [Fits and curves](../results/cbh/map_snow_gamma2_v1/REPORT.md); [country comparisons](../results/cbh/map_snow_gamma2_v1/burden/REPORT.md).

### 3.4.1 Key PfPR spline sensitivity analysis — current seven-band CBH model

<a id="key-pfpr-spline-sensitivity-analysis-current-seven-band-cbh-model"></a>

**Status and purpose.** Promote this suite to a required part of reporting the current age-specific mortality analysis and assessing subsequent national burden extrapolation. The initial results were inspected on 9 September 2026 before this priority was agreed; this is not a claim of prospective prespecification. Retain and report all comparisons regardless of the direction or statistical significance of the estimated effects. These specifications apply to the seven-band child–age-band model, rather than the historical six-band negative-binomial plan above.

**Primary reference and fixed inputs.** Use the seven separate-age fits in section 3.1 as the reference for future sensitivity reporting. The completed exploratory figure used the previous joint model as its reference and is retained with that labeling. Hold the saved posterior-median child HIV-incidence imputation fixed for a directly comparable initial screen, retaining complete-case selection for other confounders, their scaling and exposure timing at band entry. The full sample contains 5,885,022 child-band records and 82,415 deaths from 105 surveys in 34 countries.

Run the following as **three separate comparisons**, without automatically crossing geography and period:

| Comparison | Required implementation |
|---|---|
| Joint versus separate-age primary | Compare the seven primary fits with the previous joint seven-band model, which has separate PfPR curves and confounder coefficients but one shared calendar-year spline, shared survey/region intercepts and common variance parameters for each random-effect type. Keep the same full sample and HIV imputation. This comparison has already been fitted; only its primary/comparator designation changes. |
| Geography | Fit each age band separately within UN M49 Western Africa and within Eastern plus Middle (Central) Africa: seven fits per geographic group. Exclude Southern Africa from this comparison only: Namibia, Eswatini and South Africa in the current sample. The groups contain 13 and 18 countries. Each subgroup-age fit estimates its own time curve and random-effect variances. The completed exploratory geographic results used joint models and remain preliminary evidence; these separate-age geographic fits are planned, not yet run. |
| Survey period | Fit each age band separately in early and late surveys: seven fits per period. Take the median metadata survey year, counting each included survey once; assign surveys at or below the median to early and those above it to late. The present median is 2012: 53 surveys in 2003–2012 and 52 in 2013–2024. Keep surveys intact, assign ties to early, and retain all countries. Retrospective exposure years can overlap. Each subgroup-age fit has its own time curve and random-effect variances. The completed exploratory period results used joint models; the separate-age period fits are planned, not yet run. |

Re-estimate smoothing parameters in every new fit, retaining the same cubic regression spline basis dimensions; knots adapt to each fitting sample. Record the exact sample, country/survey assignments, formula, exposure support, convergence and smoothing diagnostics for every fit.

**Required comparison figure.** Overlay each sensitivity curve with the separate-age primary curve as `f_g(P) - f_g(20)`, using common axes and the same 20% PfPR reference. Label the model structure and HIV-imputation treatment explicitly. Show pointwise 95% intervals and each fit's central 95% exposure range. Report 40%→20% hazard ratios as a numerical anchor, alongside full curves and support flags; this one contrast does not capture all important shape differences. Compare curves on overlapping exposure support and distinguish changes in sample composition from changes in specification. Preserve the existing joint-reference figure as the initial comparison rather than presenting its joint subgroup results as separate-age results.

**Initial findings to retain (completed joint-reference comparison).** Separate age-band fits leave the overall spline shapes broadly similar to the previous joint reference; the largest absolute difference over common central exposure ranges is approximately 0.10 on the log-hazard-ratio scale. In the completed joint subgroup fits, geographic differences are clearest at low PfPR in older children, particularly 48–59 months, where West Africa has a steeper rise toward 20% prevalence. Period differences are clearest at 24–35 and 48–59 months. For 48–59 months, the joint subgroup 40%→20% hazard ratio is 0.93 (95% interval 0.81–1.08) early versus 1.16 (1.01–1.33) late. These findings motivate the planned separate-age subgroup checks and examination of low-prevalence curve stability before interpreting PfPR-to-zero burden estimates; they are not results from those planned refits.

**Interpretation and remaining checks.** The eleven fits in the initial sensitivity suite converged and had finite coefficients and covariance matrices. The previous joint reference and joint West Africa fits have small negative smoothing-Hessian eigenvalues, so smoothing stability remains an open check for those comparators. The initial intervals condition on one HIV imputation and the fitted smoothing parameters; they omit survey-design, residual-clustering and MAP uncertainty. The comparisons are descriptive, not formal tests of curve differences or proof of causal effect heterogeneity. For final inference, address the outstanding uncertainty and numerical checks, and assess country composition and exposure support when interpreting geographic and period differences.

![Key sensitivity: PfPR splines by age, geography and survey period](<../Key results/pfpr_spline_sensitivity_by_age_geography_period.png>)

Implementation: [sensitivity code and instructions](../R_cbh/sensitivity/README.md). Results: [full report, contrasts and diagnostics](../results/cbh/age_band_hiv_incidence_shared_time_v3/sensitivity_single_imputation/REPORT.md). Figure source: `R_cbh/sensitivity/02_report.R`; curated copy: `Key results/pfpr_spline_sensitivity_by_age_geography_period.png`.

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

# Analysis plan: malaria prevalence and child mortality

Updated 17 September 2026 (18-variable regional adjustment, urban residence retained, and vaccination/regional-mean fallback). This is a working analysis plan revised after inspection of exploratory results, not a prospective preregistration. The [previous plan](ANALYSIS_PLAN_BEFORE_MAP_GAMMA2.md) preserves the original audit, superseded specifications and decision history; the [code audit](CODE_AUDIT.md) compares each component with the current implementation and records the archival changes.

**Current primary analysis:** seven separate child–age-band mortality models using **annual MAP PfPR₂–₁₀ at band entry**, fitted with **`mgcv::bam`, `gamma = 2`, PfPR `bs = "cr", k = 5`, and calendar-year `bs = "cr", k = 6`**. Use the full eligible MAP sample, including eligible observations after 2015. Each age band has its own time spline, confounder coefficients and survey/country/region random effects and variance parameters. Section 3.1 gives the complete specification. **The latest 17 September revision uses 18 survey-region/national annual confounders. Retain urban residence (`urban_pct`), DTP3/measles, facility delivery, short birth interval, WASH and electricity alongside the other retained adjustment concepts. Drop Hib3, PCV, rotavirus and exclusive breastfeeding from both the formula and complete-case requirements. Apply exact country/survey-year UNICEF substitution for missing DTP3/measles and available-region means within the same survey for remaining regional gaps, preserving provenance. This revised specification has not yet been fitted.**

**Supplementary results:** comparisons with the annual Snow prevalence estimates, including matched MAP–Snow fits at gamma=2, all Snow penalisation analyses and Snow-based national mortality scenarios. The Snow availability cutoff at 2015 applies only to those analyses and their matched MAP comparators. Snow is not the primary exposure. The [supplementary results index](<../Supplementary results/README.md>) links to the saved figures and tables.

**Previously fitted benchmark (superseded adjustment set):** the [primary-only gamma=2 manifest](../results/cbh/primary_map_gamma2_v1/fit_manifest.csv) and associated result tables in `results/cbh/primary_map_gamma2_v1/`. All seven entries are `map_full`, but these fits, their sample counts, sensitivity results, burden estimates and paper figures still use the previous eleven-variable individual/country adjustment set. Preserve them as a benchmark. **Do not present them as results of the revised regional adjustment model.** `Rscript run_all.R --primary` currently reproduces that benchmark; a new prepared dataset and versioned fitting/reporting configuration are required before promoting revised results. See [pipeline instructions](../R_cbh/primary/README.md) and the [regional covariate audit](../R_cbh/covariates/README.md).

## 1. Scientific questions and reporting choices

**Manuscript naming rule (code pointer):** call the proposed method **“PfPR-ACM model”** in every main-paper and supplementary figure that identifies the method, including legends, axis labels and captions. Do not use “New method”, “Our model” or “Primary MAP model” as its display name. Plotting code must source [R_cbh/reporting/labels.R](../R_cbh/reporting/labels.R) and use `cbh_paper_model_label()` rather than defining its own label. Where needed, append an exposure or sensitivity qualifier to this shared name. Preserve descriptive outcome labels and the no-title/no-subtitle convention; figures without a model-name label do not need an extra label. This is a presentation change, not a change to the model, estimates or internal data identifiers.

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

**Decision, 17 September 2026:** all confounding predictors are **survey-region summaries**, with no individual child, mother or household predictor in the primary model. Retain urban/rural residence and the other previous concepts, plus DTP3/measles vaccination, facility delivery, birth interval, WASH and electricity. Exclude Hib3, PCV, rotavirus and exclusive breastfeeding from the primary model and its complete-case selection. The outcome remains a child–age-band observation. Regional averaging changes the adjustment estimand and does not substitute for individual-level confounding control.

| Confounder | Regional definition and extraction denominator |
|---|---|
| Sex; multiple birth | Weighted percentages male and from multiple births among live births in the 60 months before interview; each birth counted once, irrespective of survival or eventual model inclusion |
| Birth order; maternal age at birth | Weighted mean birth order and maternal age at birth among those recent births |
| Maternal education; wealth; urban residence | Weighted mean education years, mean wealth-quintile score (1–5), and percentage urban among distinct interviewed mothers with a birth in those 60 months; each mother counted once. These are summaries of mothers represented in the birth history, not all households or all women |
| DTP3; measles vaccination | Published DHS survey-region coverage (normally ages 12–23 months); if missing, substitute the UNICEF WUENIC national estimate for the survey year (`IM_DTP3`, `IM_MCV1`). Preserve regional estimates and mark national substitutions |
| Facility delivery | Published regional percentage delivered in a health facility; prefer the five-year recall window, otherwise retain the reported three-/two-year window with its flag |
| Birth interval | Regional percentage of non-first births in the previous five years with a preceding interval <24 months: sum published 7–17- and 18–23-month categories. First births are outside this coverage denominator; they are **not excluded from mortality analysis** |
| WASH | Two separate published regional percentages: households with an improved drinking-water source and households with an improved sanitation facility |
| Electricity | Published regional percentage of households with electricity |
| HIV incidence, GDP, health expenditure, political stability | Existing national annual measures assigned to each survey-region at band-entry year; log-transform the first three as previously specified |

For recode-derived means use DHS women's sampling weights (`v005/1e6`) and the denominator above, calculated **before exposure/outcome complete-case selection**. Never average child-band rows, which would weight surviving children repeatedly. Published household indicators use their source's household denominator and weighting; do not approximate these using women’s weights on deduplicated birth histories. Preserve the region key's survey/boundary version. National annual variables remain constant across regions within a country-year because no regional source is available; assigning them to a region does not create regional variation.

All 18 scalar confounders enter linearly after centering/scaling on the revised complete-case modelling sample. Mean wealth-quintile score is a numerical regional summary, not the former five-level individual factor. Save new scaling in a new analysis version; do not reuse the old individual-covariate scaling. Within each separately fitted age band, regional covariates have their own coefficients. PfPR and fractional calendar time remain band-entry exposure/time variables rather than being averaged over children.

Survey-region DHS coverage, wealth, urban residence and education describe survey-time conditions or the indicator's stated retrospective window. Do not relabel them as reconstructed historical values at band entry. **Urban/rural residence is required as `urban_pct`, the weighted percentage urban among distinct interviewed mothers with a recent birth in the region.** DTP3 and measles are the only selected vaccine coverage measures. Preserve the excluded vaccine/breastfeeding extraction history, but neither their missingness nor pre-introduction assumptions affect primary eligibility. Maternal death fraction, hypothetical children and wasting remain outside the adjustment set.

**National UNICEF vaccination imputation:** retain every available regional DTP3/measles value; fill only unavailable values with a valid national WUENIC percentage for the **same survey year**. This timing matches those survey-time regional measures. Preserve genuine reported zeros. This is deterministic national substitution, not a fitted imputation model. The historical national Hib3/PCV/rotavirus series and pre-rollout zero scenarios are retained for provenance only and are not used by the current primary adjustment set.

Save the value before substitution, an imputation indicator, national source year and source category for every vaccine. Substitutions remain national values assigned to survey-regions; they do not reconstruct within-country differences. Imputation flags describe provenance and do not automatically add predictor terms. Uncertainty in national estimates or regional departures from them is not propagated by this single substitution. Include an observed-regional-only/reduced-adjustment sensitivity to assess dependence on the fallback. Remaining regional gaps follow the available-region rule below; retained national annual gaps remain missing after the existing HIV imputation.

**Regional missingness fallback:** after the DTP3/measles national fallback, substitute each missing regional covariate with the arithmetic mean of finite available regions in the **same survey**, separately for each variable. Count each donor region once, not once per child-band record. Preserve all observed values and store original values, imputation flags, donor counts and source labels. Missing individual responses are already excluded from within-region weighted means and must not propagate into an otherwise estimable aggregate. If no usable donor region exists in that survey, leave the value missing; do not borrow across surveys or apply regional means to national country-year variables. This carries forward the user-specified within-survey imputation approach; its uncertainty is not propagated. Retain a no-regional-substitution sensitivity.

The [regional pipeline](../R_cbh/covariates/README.md) separately retrieves public indicator snapshots, extracts weighted recode summaries, resolves published regional joins and audits availability. Exposure/region identity must still be resolved; a missing regional vaccination covariate may use the declared UNICEF fallback, without imputing location or PfPR. The [missingness audit](../results/cbh/regional_adjustment_reduced_v3/README.md) distinguishes missing records, entirely lost survey-regions and partially reduced survey-regions. A regional mean can be observed when some individual responses are missing; respondent-level missingness and small denominators are reported separately. Do not drop a mortality record solely because that child's own confounder response is missing.

Use **child HIV incidence at ages 0–14 per 1,000 uninfected population**, not HIV prevalence or infection counts. The [HIV imputation model](../R_cbh/hiv/README.md) uses child and adolescent incidence jointly, calendar-year trends, geographic/country effects and within-country temporal correlation. Preserve reported child rates, handle censored values as censored, and predict missing child series. Adolescent incidence informs imputation; it is not a separate predictor in the mortality model.

The revised primary model retains **one fixed posterior-median child-incidence imputation** and requires complete cases for the selected **regional/annual** model variables **after UNICEF vaccination and available-region substitution**. The currently saved fits implement the earlier adjustment set. Keep censoring, imputation and missing-source flags; unsupported country-year values remain missing. The historical HIV prevalence field in base shards is provenance only. There is no country/overall-median imputation. The only new national substitution is the explicitly specified UNICEF vaccination fallback.

Report missingness and exclusions by country, survey, year, age and exposure. Missing HIV/source data, unresolved geography and complete-case selection remain substantive issues. Retain posterior HIV trajectories for future uncertainty propagation through all seven primary models. For validation with held-out surveys/countries, estimate preprocessing and imputation within training folds where appropriate.

**Availability under the reduced 18-variable policy, 17 September:** start with 6,357,802 MAP-eligible child-band records across 120 surveys. Complete cases after the declared substitutions retain **5,680,117 records, 1,755,838 children and 78,634 deaths in 973 survey-regions, 100 surveys and 34 countries**; 10.66% of the starting records remain excluded. Relative to the previous fitted sample, 973/1,015 survey-regions remain and **42 (4.14%)** are lost; 5,617,520 of its 5,885,022 records remain (4.55% excluded). Urban proportion has 0% missingness in every audited denominator. See the [current audit](../results/cbh/regional_adjustment_reduced_v3/README.md). The [22-variable UNICEF audit](../results/cbh/regional_adjustment_unicef_v2/README.md) and its pre-rollout/region-mean scenarios remain historical and must not supply current sample counts. These are data-availability counts, not fitted model results.

### 2.5 Analysis-ready outputs

| Output | Role |
|---|---|
| `data/derived_cbh/manifest.rds` and `child_bands/<survey>.rds` | Authoritative build manifest and eligible child-band shards with audits |
| `data/derived_cbh/hiv_incidence/child_incidence_country_year.csv` | Incidence estimates and source/imputation flags |
| `data/derived_cbh/regional_adjustment/reduced_v3/regional_covariates_wide.csv` | Current 14-regional-variable overlay after UNICEF and available-region fallback, including original values and imputation flags; urban_pct retained |
| `data/derived_cbh/regional_adjustment/unicef/country_year_estimates.csv` | WUENIC source snapshot: only DTP3/measles are selected for current primary adjustment; the other three vaccine series are historical |
| `results/cbh/regional_adjustment_reduced_v3/complete_case_summary.csv` | Current availability after substitution, including losses relative to the previous fitted sample |
| `data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/complete_case_dataset.rds` | Previous individual-adjustment modelling sample and preprocessing; benchmark only |
| `results/cbh/primary_map_gamma2_v1/fit_manifest.csv` | Seven previous-adjustment MAP gamma=2 objects; all entries are `map_full` |
| Aggregate selection, missingness, geography and model diagnostics | Inspectable accounting without microdata |

The **previously fitted** sample contains **1,817,912 distinct children contributing 5,885,022 child-band records and 82,415 deaths from 105 surveys in 34 countries**. The [15 September verification](../results/cbh/code_audit_2026_09_15/primary_fit_checks.csv) compared every primary fitted outcome, predictor and offset with the saved prepared sample and verified the selected model hashes and knots. Keep raw/derived microdata and fitted objects under ignored `data/`; export only approved aggregates to results. Validate unique child-band keys and input hashes. The primary-only fitter now reads this prepared dataset directly and checks every outcome, predictor and offset. Its reference MAP knots are preserved in a committed snapshot. A future consolidated data stage must still operate independently of existing outcome models and results tables.

## 3. Stage 2 — primary and sensitivity analyses

### 3.1 Primary working model

**Separate age-band MAP models with gamma=2, with the regional adjustment revision of 17 September 2026.** For each of the seven bands, fit:

```r
form <- death ~
  s(pfpr_pct, bs = "cr", k = 5) +
  s(calendar_year, bs = "cr", k = 6) +
  z_male_pct + z_multiple_birth_pct + z_mean_birth_order +
  z_mean_maternal_age_birth + z_mean_maternal_education_years +
  z_mean_wealth_quintile + z_urban_pct +
  z_dtp3_pct + z_measles_pct + z_facility_delivery_pct +
  z_short_birth_interval_pct + z_improved_water_pct +
  z_improved_sanitation_pct + z_electricity_pct +
  z_log_hiv_incidence + z_log_gdp_pc +
  z_log_health_expenditure_pc + z_political_stability +
  s(survey, bs = "re") + s(country, bs = "re") +
  s(region, bs = "re") + offset(log(band_years))

# d_g: revised regional-adjustment complete cases for this one age band.
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
\mathbf X_{r(i),s(i),y_{ig}}^{\mathsf T}\boldsymbol\beta_g+
u_{s(i),g}+v_{c(i),g}+b_{r(i),g}.
\]

Here X contains survey-region summaries from Section 2.4 and national annual covariates assigned to those regions; it contains no individual confounding predictor. P is regional annual MAP PfPR₂–₁₀, y is band-entry year, and t is fractional calendar year at entry. The exposure surface has no age-band-specific definition; the function f and the entry year depend on the band. Equivalently, `cloglog(q_ig) = log(Δg) + log(λ_ig)`.

Each fit estimates its own intercept, PfPR and time curves, confounder coefficients and survey/country/region random intercepts. Within each model, the random effects are independent zero-mean Gaussians, with a separate estimated variance for each grouping type. No time spline, random intercept or variance parameter is shared across age bands. Separate fitting does not imply independent sampling errors across ages because children and surveys contribute to multiple bands.

`gamma = 2` changes smoothing selection throughout the model, including time and random effects; it is not a PfPR-only multiplier. The primary PfPR basis remains **cr**, not cs. No monotonicity constraint is imposed. Reuse the reference MAP knots for these promoted fits; they are recorded and verified in the fitting pipeline.

The unweighted likelihood, fixed HIV imputation and conditional model covariance are current working choices. Address survey design, residual dependence and imputation/exposure uncertainty before claiming comprehensive uncertainty. The primary decision follows review of results and user instruction; it does not establish superior causal identification or formal model-selection evidence.

**Revised formula/data specification:** [covariates/regional.R](../R_cbh/covariates/regional.R) defines `cbh_regional_spec()` and `cbh_regional_formula()`. [covariates/settings.R](../R_cbh/covariates/settings.R) selects the current `regional_adjustment_reduced_v3` 18-variable data policy; [covariates/unicef.R](../R_cbh/covariates/unicef.R) implements national vaccination substitution, and `cbh_regional_mean_fill()` implements the available-region fallback. `cbh_regional_spec()` and `cbh_regional_formula()` default to 18 predictors with `urban_pct` retained; `expanded=TRUE` is explicitly historical. Refit and reporting integration are pending. **Previous benchmark implementation:** [primary/01_fit.R](../R_cbh/primary/01_fit.R), configured in [settings.R](../R_cbh/primary/settings.R), with `map_full_age_1` through `map_full_age_7` in the [primary manifest](../results/cbh/primary_map_gamma2_v1/fit_manifest.csv). Fitted objects are in `data/derived_cbh/models/primary_map_gamma2_v1/`. The fitting sample is read from the existing prepared dataset, with saved scaling and [reference knots](../R_cbh/primary/reference_knots.csv). The older combined MAP–Snow workflow remains for supplementary reproduction; its `map_full` results are the numerical comparison for this rerun. The older `sensitivity_single_imputation/age_*.rds` fits used gamma=1.

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
| Required: missingness and adjustment | Selection patterns, HIV censoring/absent series, imputation uncertainty and reduced adjustment blocks | The expanded regional adjustment set in Section 2.4 is primary. Compare its selection and effects with reduced adjustment blocks and the previous individual-adjustment benchmark; wasting remains an optional separate extension. |
| Required: design and dependence | Survey weights; common cluster/stratum bootstrap; residual regional and within-child dependence | Retain unweighted primary point fit as reference; do not substitute rounded effective counts without a separate model definition. |
| Required: dose response / reference | Basis dimensions, linear comparison, 0% versus 1%, exposure restrictions | Match samples; report all supported contrasts and low-PfPR sensitivity. |
| Secondary: heterogeneity | Country omission, calendar interaction, within/between exposure and covariate timing | Separate, declared extensions of the primary structure. |
| Secondary: burden assumptions | Geographic coverage, population weights, national-mean versus subnational integration, age allocation and denominator sources | IHME all-cause is current baseline; IGME/WHO comparisons are supplementary alternatives. |
| Optional independent checks | Measured parasitaemia, fieldwork season, trial triangulation and alternative mortality estimators | Retain only for a stated scientific question; historical survey-region/NB/brms work remains archived. |

### 3.4.1 Key PfPR spline sensitivity analysis — current seven-band CBH model

<a id="key-pfpr-spline-sensitivity-analysis-current-seven-band-cbh-model"></a>

Once refitted, use the seven revised regional-adjustment MAP gamma=2 fits as the primary reference. Existing subgroup and source comparisons below describe the earlier adjustment set and require matching refits. Hold the median child HIV imputation, covariate selection/scaling and entry-year exposure timing fixed. Re-estimate smoothing parameters in each new fit at **gamma=2**, retaining the basis types and dimensions; record subgroup-specific knots and support.

- **Joint versus separate:** compare with a joint model having age-specific PfPR curves and confounder coefficients but one shared calendar-year spline and the historical joint random-effect structure. The completed joint/separate comparison used gamma=1; a matched gamma=2 joint fit remains to be run. Do not relabel the old comparison as gamma=2.
- **Geography:** fit each age separately in UN M49 Western Africa and Eastern plus Middle (Central) Africa. Exclude Southern Africa only for this comparison: Namibia, Eswatini and South Africa. The current groups contain 13 and 18 countries. Each subgroup-age fit has its own time curve and random-effect variances. Separate-age gamma=2 subgroup fits remain planned.
- **Seasonal-area restriction (Sahel, 12°N):** refit all seven primary MAP gamma=2 age-band models using only survey regions with a cached boundary centroid at or north of 12°N, west of 36°E and outside Ethiopia, Eritrea, Somalia and Djibouti. This operational geographic proxy selects whole regions, not individual child locations or measured transmission seasonality. Retain the same primary eligibility, confounders/scaling, fixed child HIV imputation, basis types/dimensions and duration offsets. Re-estimate smoothing parameters with subgroup-specific knots at quantiles of distinct predictor values. The selected sample contains **378,680 children, 1,230,499 records and 15,999 deaths across 196 survey regions, 25 surveys and seven countries** (Burkina Faso, Chad, The Gambia, Mali, Mauritania, Niger and Senegal). Nigeria's broad northern-zone centroids fall below 12°N and are excluded. Overlay conditional HR curves `exp{f_g(P)-f_g(20)}` and pointwise intervals against the saved full-sample primary curves, drawing each fit over its central exposure range. Compare supported 40%→20% contrasts descriptively: subgroup and full-sample fits overlap, so do not assume independent estimates or treat interval overlap as an interaction test. [Code](../R_cbh/sensitivity/sahel/01_fit.R), [results and geographic audit](../results/cbh/sahel_map_gamma2_v1/REPORT.md), [overlay](../results/cbh/sahel_map_gamma2_v1/sahel_vs_full_pfpr_hazard_ratios.png). This is distinct from the earlier descriptive monthly-mortality analysis using an 11°N cutoff.
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

**Revision status:** all existing figures and numerical burden results below use the previous adjustment set. Regenerate them from the new version after the regional-adjustment fits are completed; retain the figure definitions and naming conventions.

| Reporting location | Content | Authoritative selection |
|---|---|---|
| Main results | Study flow, survey/sample coverage and seven age-band descriptive summaries | Full MAP sample and its exclusion ledger |
| Main results: Figure 1 | Survey map and timing for the 105 included surveys in 34 countries | `survey_map/survey_map_and_timing.png` under current primary results |
| Main results | Seven primary PfPR curves, X%→Y% contrasts and diagnostics | `map_full`, gamma=2, from `results/cbh/primary_map_gamma2_v1/` |
| Main results: secondary burden | Country estimates and age contributions for 2005/2015/2024 using the primary model | Primary-only burden files linked above |
| Supplementary results | Matched MAP–Snow curves and contrasts; available-sample comparisons; all Snow smoothing trials and Snow mortality scenarios | [Supplementary index](<../Supplementary results/README.md>) |
| Supplementary sensitivity results | Structure, geography, survey period, timing, design, missingness and other declared checks | Label old gamma=1 exploratory results separately from planned gamma=2 refits |
| Historical archive | Six-band person-time/count models, old joint fits and earlier burden outputs | [Archived plan](ANALYSIS_PLAN_BEFORE_MAP_GAMMA2.md), original result directories and code audit |

The paper's Figure 2 is the primary PfPR-by-age curve figure. **Figure 3** shows age bands on the x-axis and the malaria-attributable share of all-cause deaths **within each band** on the y-axis (0–100%), with four lines for PfPR of 10%, 20%, 30% and 40%. Calculate `100 × {1 − exp[f_g(0) − f_g(P)]}` from the saved primary MAP gamma=2 models, using zero prevalence as the counterfactual and the same hazard-based, fixed-person-time convention as the national burden analysis. This is not a distribution of total under-five deaths across ages. Connect the seven point estimates without implying a continuous-age fit; retain conditional intervals and zero-exposure extrapolation flags in the source table. [Figure](../results/cbh/primary_map_gamma2_v1/attributable_fraction_by_age.png), [estimates](../results/cbh/primary_map_gamma2_v1/attributable_fraction_by_age.csv), [code](../R_cbh/reporting/06_attributable_fraction_by_age.R).

**Figure 4** (previously Figure 3) is the country comparison with IHME for 2005/2015/2024, using true base-10 logarithmic scales on both axes with the same limits. Both axes start at 1,000 deaths. All 42 estimable countries per year have positive estimates; points below the display minimum on either axis are clipped from the figure, with all underlying estimates retained. All five paper figures omit titles, subtitles and embedded captions; separate caption text and filename mappings are in the [paper figure manifest](../results/cbh/primary_map_gamma2_v1/paper_figures/CAPTIONS.md). Preparing/copying figures does not edit manuscript TeX.

**Figure 5** compares annual malaria-attributable under-five mortality rates from the PfPR-ACM model with IHME cause-specific malaria and UN IGME, for 2004–2024. The figure legend must use “UN IGME”; the caption identifies the underlying CA-CODE 2026 release, which is distinct from the WHO World Malaria Report. Use the same 42 national-burden countries in every year. Recalculate national MAP exposure with population counts (density × cell area) and fixed GPW 2020 spatial weights. Apply the saved primary age-band effects to annual IHME all-cause deaths; preserve the equal-person-time assumption for ages 2–4 and signed age contributions. Sum deaths across countries and divide all three sources by the same annual under-five person-years inferred from the IHME all-cause under-five count/rate pair. Report deaths and rates per 100,000 child-years, not per live birth or an average of country rates. Use the direct under-five CA-CODE 2026 series from the official UN IGME/UNICEF portal, as agreed; do not multiply all-age WHO totals by a constant child fraction. Input coverage is complete for all 882 country-years. Small differences between IHME source-implied population denominators are retained in the audit; exact release alignment remains unresolved. Plot point estimates without summing marginal uncertainty bounds. The method shares IHME all-cause inputs, and attributable all-cause reductions and cause-specific comparator estimates are different estimands. [Figure and audit](../results/cbh/primary_map_gamma2_v1/annual_comparison/README.md), [annual totals](../results/cbh/primary_map_gamma2_v1/annual_comparison/annual_totals_2004_2024.csv), [calculation code](../R_cbh/burden/04_annual_comparison.R).

Plotting-only stages must read saved aggregate estimates and must not refit, change exposure definitions or select models according to available cache files. Label source, sample period, gamma, age bands, contrast and uncertainty. A combined MAP–Snow figure is supplementary; a main-results figure must select only the full MAP gamma=2 series. Keep links to authoritative saved artifacts instead of copying figures into competing locations.

## 5. Implementation status and next steps

**Previous-adjustment primary-only rerun:** [R_cbh/primary/](../R_cbh/primary/README.md) now separates fresh fitting, effect/burden calculation, aggregate diagnostic checks and MAP-only plotting. It regenerates the study inclusion flow and key-results index. The prepared dataset, HIV imputation, national exposure inputs and IHME inputs remain fixed; supplementary fits are not run. New outputs are saved separately from the preceding comparison workflow, with numerical differences documented in the [rerun report](../results/cbh/primary_map_gamma2_v1/REPORT.md).

**15 September audit and cleanup:** the [component-by-component audit](CODE_AUDIT.md) verified that the saved primary specification and data matched the plan at that date, before the 17 September regional-adjustment revision, while distinguishing unimplemented sensitivities and unresolved upstream issues. Sixty-eight superseded scripts were moved to [the archive](../archive/2026-09-15-code-audit/), along with preserved copies of three replaced entry points. Existing data, fitted objects and numerical results were not changed.

`Rscript run_all.R --primary` currently forces fresh fitting of the previous adjustment specification and regenerates its outputs without data setup; it does not yet implement the regional revision. `Rscript R_cbh/primary/run.R --resume` reuses only valid completed fit caches; `--report-only` uses saved fits for effects, diagnostics and plotting. `Rscript run_all.R --audit` retains verification of the preserved pre-rerun primary snapshot; `--check-inputs` checks local dataset prerequisites. The [results index](<../Key results/README.md>) links the new authoritative results.

**Remaining command-routing limitation:** older [sensitivity fitting](../R_cbh/sensitivity/01_fit.R) and [burden defaults](../R_cbh/burden/settings.R) still refer to gamma=1-era fits. `analysis/01_fit_complete_case.R` should be used with `--prepare-only` for dataset preparation; its default fit is historical. These scripts remain because they provide reference frames/knots or input tables used by the selected workflow. Do not use those defaults as a gamma=2 primary run. The new primary-only pipeline bypasses these defaults. National burden still reads explicit baseline columns from the earlier saved national input tables; consolidating their upstream producers remains a data-stage task.

Retain the three-stage workflow: **make datasets → fit declared models and calculate effects → plot saved results**. Priorities are:

1. Resolve flagged regional-covariate extraction issues, prepare the expanded complete-case dataset and refit all seven models under a new version; regenerate sensitivities, burden, sample flow and paper results before promoting them. Preserve the primary/supplementary separation and input provenance. Primary-only routing and plotting are now implemented; in-sample fitted-outcome checks supplement the numerical diagnostics, but are not cross-validation or influence analysis.
2. Complete primary-compatible gamma=2 joint, geographic and period sensitivity fits; preserve historical comparisons with their actual settings.
3. Resolve the [open data and inference issues](OPEN_ANALYSIS_ISSUES.md), especially regional exposure coverage, missing-data selection, calendar/age coding, imputation uncertainty and survey/within-child dependence.
4. Consolidate the independent data build and prediction interfaces. Rebuild from a declared input snapshot and verify cache hashes. Fresh primary fitting now uses prepared data and committed reference knots; baseline national input tables still need an independent producer interface.
5. Finish separating input-only functions from the retained legacy bridges, then archive their remaining surplus branches. The 15 September cleanup removed the superseded independent modelling/plotting scripts; still-needed producers were retained with explicit reasons. Keep historical code/results traceable and do not execute legacy external-copy actions during migration.

Completion requires auditable survey accounting, validated child-band eligibility and joins, consistent primary specification and artifact selection, reviewed diagnostics, declared uncertainty and reproducible main/supplementary figures. Historical six-band/count-model choices are preserved in the archive and must not override this plan.

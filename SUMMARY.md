---
title: "DHS/MIS malaria prevalence and child mortality analysis"
output:
  word_document: default
  pdf_document:
    latex_engine: xelatex
  html_document:
    self_contained: true
header-includes:
  - \usepackage{float}
  - \floatplacement{figure}{H}
geometry: margin=0.9in
---

## Main results

- **Higher malaria prevalence is associated with higher post-neonatal
  mortality, and this holds regardless of covariate adjustment.** Each 10
  percentage-point increase in PfPR2-10 is associated with a 6.9% increase in
  post-neonatal mortality in the linear summary model (95% CI 4.5%–9.4%;
  `p<0.001`), and with a 9.0% increase (6.4%–11.8%) at the reference year 2013
  in the selected model, which lets the slope change over calendar time. The
  estimate is essentially unchanged whether we use the regional covariate block
  alone or the full regional-plus-national block that includes national vaccine
  coverage and child HIV prevalence.

- **The primary analysis now uses a 12-month mortality window rather than the
  DHS default of 60 months.** The exposure is MAP prevalence in the survey year,
  so a 60-month window centres the outcome about 2.5 years before the exposure it
  is regressed on; a 12-month window nearly removes that mismatch. The cost is
  precision (the interval on a regional rate roughly doubles) and 56 region-years
  that record no post-neonatal or no neonatal death within the year, leaving
  **1,055 region-years (118 surveys, 35 countries)**. The dose-response barely
  moves across windows once the model structure is held fixed.

- **Neonatal mortality, the negative control, is null.** The linear neonatal
  estimate is +0.1% per 10 PfPR points (95% CI −2.2% to +2.5%; `p=0.92`), the
  prespecified structure refitted to neonatal mortality gives +1.2% (−1.4% to
  +3.8%; `p=0.38`), and the Bayesian prevalence-by-time surface fitted to
  neonatal deaths sits on zero at every year and in all nine subgroups. At 60
  months on the same panel the control had drifted to +1.7% (`p=0.045`); the
  shorter window removes that. This, with the stability under HIV and vaccine
  adjustment, reduces the concern that the main result reflects broad
  differences in child survival between higher- and lower-prevalence settings.

- **The data support letting the prevalence effect change over time, but do not
  say how.** By AIC, three structures tie within 0.44 points: a linear PfPR term
  with a time interaction, a spline with a `ti` interaction, and the additive
  spline `s(PfPR)+s(year)`. Bayesian refits with shared priors cannot be
  separated by pointwise leave-one-out (differences of 0.2–0.4 against standard
  errors of 1.2–2.3). Holding out whole surveys in a grouped 10-fold
  cross-validation does separate them: the additive model trails the full
  tensor surface by 7.1 in held-out log predictive density (pointwise SE 3.4,
  survey-level SE 4.4) and the linear-in-time model trails it by only 1.1 (SE
  1.7). The gain sits in survey regions with prevalence of 30% or more. The
  attributable fractions at 2013 are nearly identical under all three
  structures; what differs is the historical trend.

- **Population-average malaria-attributable fractions of post-neonatal
  mortality are 7.5%, 22.2% and 34.6% at 10%, 30% and 50% PfPR2-10** (versus a
  1% counterfactual, reference year 2013, selected model). The Bayesian surface,
  which propagates uncertainty in how smooth the curve is, gives 11.6%, 27.8%
  and 36.1% with wider intervals (2.6%–21.4%, 15.1%–39.2%, 22.9%–47.4%).

- **The best-fitting national extrapolation estimates about 527,000
  malaria-attributable post-neonatal deaths in 2024**, 95% interval
  373,000–658,000 (spline with `ti` interaction); the full `te` surface gives
  503,000 and a time-stable spline 435,000. IHME/GBD's 2024 under-five figure is
  about 428,000 and the WHO under-five proxy about 434,000, so the point estimate
  is 23% above IHME with overlapping intervals. This is not a validation against
  a common estimand: the model captures direct and indirect malaria-associated
  mortality, whereas IHME and WHO assign deaths specifically to malaria. The
  2000–2024 trend remains unidentified: the two best structures, 0.3 AIC apart,
  imply changes of −5% and −53%.

![Malaria and post-neonatal mortality. **A:** all-cause post-neonatal mortality
versus PfPR2-10, with point size proportional to the survey region's weighted
birth exposure and the fitted ridge-GAM curve. **B:** the model-estimated
malaria-attributable fraction of post-neonatal deaths as a function of
PfPR2-10.](results/dhs_rebuild/figure1_pfpr_mortality_and_af.png){width=100%}

The sections below give the data, methods and supporting detail.

## Data and methods

**Sample.** The unit of analysis is a DHS/MIS survey-region-year (Births
Recodes, sub-Saharan Africa, 2000–2024; AIS excluded). Each region receives the
population-weighted MAP PfPR2-10 estimate for the survey year; regional
all-under-five and neonatal mortality come from birth histories via the DHS
synthetic-cohort life-table method over the 12 months before interview, and
post-neonatal mortality is their difference. The assembled panel has 1,128
region-years (120 surveys, 36 countries); the primary sample requires
country-mean PfPR2-10 above 1% and regional PfPR2-10 at least 1% (1,111
region-years), plus at least one post-neonatal and one neonatal death in the
window, giving **1,055 region-years (118 surveys, 35 countries)**. The full data
flow is shown in the [analysis-flow diagram](#analysis-flow).

![Where and when the surveys were fielded. Left: countries contributing to the
panel, shaded by number of surveys. Right: fieldwork year of each survey, with
point size giving the number of survey
regions.](results/dhs_rebuild/figure10_survey_map.png){width=100%}

**Mortality window.** Regional mortality estimated over 12 months tracks the
60-month estimate closely for post-neonatal mortality (correlation 0.90 across
1,128 survey regions, median ratio 1.06) and less closely for neonatal mortality
(correlation 0.65), with intervals about 1.9 times wider in both cases. The
model comparison was repeated at 12, 24, 36, 48 and 60 months with MAP
prevalence lagged 0, 1 or 2 years. The selected structure varies across the 15
cells (a time-interaction form in 12 of them). The attributable fraction at 30%
prevalence ranges from 14% to 32% for the selected structure across the 15
cells and from 22% to 33% for the additive spline; at zero lag, the primary
choice, the selected-structure range is 21% to 32%. The window changes precision
far more than it changes the estimate; the prevalence lag matters more than the
window.

![Post-neonatal and neonatal mortality over the 12-month window against the
60-month window, one point per survey region with 95% delete-one-cluster
jackknife intervals.](results/dhs_rebuild/figure8b_horizon_12_vs_60.png){width=100%}

**Panel coverage.** The panel is assembled from the DHS survey registry rather
than inherited from an earlier aggregate. Three faults in that path were
repaired — a Windows-1250 mis-decoding of accented region names, a `v022` that
made `chmort` abort for five surveys (empty in four, and populated but with
clusters not nested within strata in Malawi 2004), and a MAP regional extraction
that had never completed — and recode region names are reconciled against
boundary names across language, word order, spelling, qualifier clauses and
granularity. Together these took the panel from 936 region-years (106 surveys)
to 1,128 (120 surveys).

Of the 1,134 survey regions for which MAP publishes a prevalence estimate,
**1,128 (99.5%) merge into the panel**. Only 995 of those joined on an exact key;
the rest were reconciled by token signature (78), prefix (22), edit distance
(18), a curated synonym table (7) or elimination (5).

This matters for interpretation, not only for coverage: on the 60-month window
the linear no-interaction estimate fell at every step as coverage improved
(+8.5%, +8.30%, +7.61%, +7.33% per 10 PfPR points across the four successive
panels), which suggests the surveys previously dropped were not missing at
random with respect to the association. Six regions across two surveys still
cannot be joined — Mali 2012 did not survey four of its regions, and Uganda
2016's Buganda/Central labelling is genuinely ambiguous — and per-survey merge
detail is recorded in `results/dhs_rebuild/region_merge_quality.csv`.

**Covariates.** Twenty-five candidate covariates are screened before imputation;
those with at most 5% missingness are retained and singly imputed (country
median, overall-median fallback). The **18 retained covariates** enter a single
standardised, ridge-penalised block (proportions logit-transformed):

| Level | Covariates |
|---|---|
| Survey-region (DHS) | urban residence, DTP3, measles, facility delivery, maternal education, exclusive breastfeeding, short birth interval, maternal age at first birth, improved water, improved sanitation, household electricity |
| National, exact year | WUENIC Hib3, PCV-completion, rotavirus-completion; and child (0–14) HIV prevalence (log) |
| National, nearest year | log GDP per capita, log health expenditure per capita, political stability |

Seven candidates are excluded for excess missingness (the three direct DHS
vaccine completion recodes, household wealth, and the three anthropometry
measures). Political stability is read from the current World Bank
`GOV_WGI_PV_EST` series after the `PV.EST` indicator was archived. National
DTP3, electricity access and health-expenditure-share are dropped as redundant
with their regional or per-capita counterparts.

Child HIV prevalence is derived from the UNAIDS 2025 estimates workbook (number
of children 0–14 living with HIV) divided by the World Bank 0–14 population
(`SP.POP.0014.TO`), joined by ISO3 and exact survey year, and entered on the log
scale because it spans a very wide range across countries (country means from
about 0.02% in Madagascar to about 3.8% in Eswatini). Nigeria and Comoros
publish only a 15–19 series (4% of the sample) and are imputed with a provenance
flag. Source choices are in
[`R_dhs/COVARIATE_SOURCES.md`](R_dhs/COVARIATE_SOURCES.md).

**Model.** Approximate death counts (rate × birth exposure) are modelled with a
negative-binomial GAM: `offset(log(exposure))`, the ridge covariate block and
country random intercepts. Country-specific random PfPR slopes, present in
earlier versions, are no longer in the primary model because they have no clear
interpretation; re-admitting them is reported as a sensitivity. Five PfPR/time
specifications are compared by ML/AIC on the same 1,055 observations, and the
selected form is refitted by REML for post-neonatal, all-under-five and neonatal
mortality.

**Bayesian refits.** The prevalence-by-time surface is also fitted in Stan via
`brms`, with `t2(PfPR, year)` in place of `te()`, the ridge block as a
`normal(0, 1)` prior on the standardised covariates, and half-t priors on the
smoothing standard deviations, so that uncertainty in how smooth the surface is
propagates into every attributable fraction rather than being fixed by an
optimiser. The same programme is refitted to the neonatal outcome and to
subgroups by era and region.

## Association and model selection

From the linear summary models, each 10 percentage-point increase in PfPR2-10
is associated with:

| Outcome | Change in mortality | 95% CI | p-value |
|---|---:|---:|---:|
| Post-neonatal | +6.92% | +4.52% to +9.37% | <0.001 |
| All under-five | +4.67% | +2.77% to +6.61% | <0.001 |
| Neonatal (negative control) | +0.12% | −2.18% to +2.48% | 0.92 |

![Modelled change in mortality versus 1% PfPR for the three outcomes. The
post-neonatal and all-under-five associations rise with prevalence; the
neonatal negative control stays flat and its interval spans
zero.](results/dhs_rebuild/figure3_negative_control_curves.png){width=88%}

Five specifications are compared by AIC on the 1,055-observation post-neonatal
sample:

| Model | Specification | AIC | ΔAIC |
|---|---|---:|---:|
| Linear, interaction | `PfPR + PfPR × year + s(year)` | 8144.70 | 0.00 |
| Spline, `ti` interaction | `s(PfPR) + s(year) + ti(PfPR, year)` | 8144.80 | 0.11 |
| Spline, no interaction | `s(PfPR) + s(year)` | 8145.14 | 0.44 |
| Full tensor surface | `te(PfPR, year)` | 8150.71 | 6.01 |
| Linear, no interaction | `PfPR + s(year)` | 8151.06 | 6.37 |

The top three are indistinguishable by AIC. In the selected model the
interaction is +0.0048 per year (`p=0.0003`): the proportional effect of
prevalence is estimated to have strengthened over 2000–2024. An independent AIC
comparison for neonatal mortality also lands on the linear time-interaction form
(by 0.96 AIC over the no-interaction form), but with a null prevalence
coefficient (+1.1% per 10 points, 95% CI −1.4% to +3.7%; `p=0.39`).

The selected specification gives the following malaria-attributable fractions of
post-neonatal mortality (country random effects set to zero; reference year
2013):

| PfPR2-10 | Attributable fraction | 95% interval |
|---:|---:|---:|
| 10% | 7.5% | 5.4%–9.5% |
| 30% | 22.2% | 16.4%–27.6% |
| 50% | 34.6% | 26.2%–42.0% |

### Does the prevalence effect change over time? A Bayesian comparison

Three models were fitted in `brms` with identical priors and differ only in how
prevalence and calendar time enter: **A**, additive `s(PfPR) + s(year)`; **C**,
A plus `s(PfPR, by = year)`, so the prevalence curve shifts linearly with time
(pinned to zero at the 1% reference so it does not duplicate the linear year
term); and **D**, the full `t2(PfPR, year)` surface. `brms` has no `ti()`, so D
is not nested above A or C and the comparison is predictive rather than a test
of one term.

| Model | Pointwise LOO elpd | Survey-grouped 10-fold elpd | Difference from best (SE) | Stacking weight (10-fold) |
|---|---:|---:|---:|---:|
| D: full surface | −4072.3 | −4129.6 | 0.0 | 0.86 |
| C: curve shifts linearly in time | −4072.1 | −4130.6 | −1.1 (1.7) | 0.14 |
| A: additive | −4072.4 | −4136.7 | −7.1 (3.4) | 0.00 |

Pointwise leave-one-out cannot separate the three: the other regions of the same
survey stay in the training set, so the country intercept and the year smooth
already pin down that survey's level. Holding out whole surveys (10 folds of
11–12 surveys, identical across models) asks whether the time structure
generalises to a survey the model has never seen, and there the additive model
falls behind by 7.1 (about two pointwise standard errors; 1.6 when the 118
surveys rather than the 1,055 region-years are treated as the exchangeable
units). The gain is consistent (70 of 118 surveys favour D; removing the three
most favourable still leaves 3.3) and sits almost entirely in survey regions with
prevalence of 30% or more (+6.7, SE 2.4) and in West Africa (+5.5, SE 2.0). C and
D are indistinguishable.

Under C and D the attributable fraction at 50% prevalence is higher in 2020 than
in 2005 by 17 and 23 percentage points respectively (95% intervals 2–32 and
6–40); at 30% the intervals include zero (11, −3 to 24; 15, −1 to 31). The
interaction is a real feature of the posterior, but it lives at high prevalence
in recent years, where the data are thinnest.

![The three models with shared priors, curves at 2005, 2013 and 2020 with 95%
credible intervals. In A the three years coincide by
construction.](results/dhs_rebuild/figure17_brms_ladder.png){width=100%}

![Survey-grouped 10-fold cross-validation: held-out log predictive density per
survey, relative to the additive model. Points above zero favour the
time-structured model.](results/dhs_rebuild/figure18_brms_ladder_kfold.png){width=100%}

### Subgroups by era and region

The same Bayesian structure fitted separately by era (before/from 2013, the
median survey year) and region (West Africa versus Central and East) gives the
following attributable fractions at 30% prevalence (95% credible intervals; n is
surveys):

| Subgroup | Surveys | AF at 30% PfPR |
|---|---:|---:|
| All data | 118 | 28.3% (17.7%–37.6%) |
| Before 2013 | 66 | 19.9% (6.4%–34.5%) |
| 2013 onwards | 52 | 35.3% (21.5%–47.8%) |
| West Africa | 51 | 31.8% (11.7%–49.8%) |
| Central and East | 67 | 21.3% (7.0%–34.0%) |
| West, before 2013 | 25 | 47.2% (16.1%–66.9%) |
| Central and East, before 2013 | 41 | 8.0% (0.0%–24.8%) |
| West, 2013 onwards | 26 | 20.3% (0.0%–45.9%) |
| Central and East, 2013 onwards | 26 | 38.3% (22.5%–52.5%) |

The era contrast is the clear signal and is consistent with the time
interaction above. The regional contrast is not: the credible intervals overlap
along the whole prevalence range, a three-way model with region-specific
surfaces buys nothing in leave-one-out (elpd difference 0.3, SE 2.7), and the
2×2 cells occupy different prevalence ranges, so the apparent crossover between
West and Central–East across eras should not be over-read. The neonatal control
fitted the same way straddles zero in every subgroup.

![Subgroup dose-response curves, fitted in Stan. Dashed lines are West Africa,
solid Central and East; light colours are before 2013, dark from 2013. Each
curve spans only the prevalence range its subgroup
observes.](results/dhs_rebuild/figure15_brms_subgroup_curves.png){width=100%}

## Robustness

Robustness is the central message: the association barely moves across covariate
blocks, samples, model structures and likelihood.

- **Covariate adjustment.** Child HIV prevalence is associated with +8.2% higher
  post-neonatal mortality per standard deviation (`p=0.012`), behind regional
  DTP3 coverage (−9.5%) and maternal age at first birth (−7.6%), yet adding the
  national block leaves the malaria estimate essentially unchanged: the regional
  block alone and the full block give attributable fractions of 24.0% and 22.3%
  at 30% PfPR.

![Ridge-standardised conditional associations for the 18 covariates in the
selected post-neonatal model (percentage change per 1 SD after transformation;
ridge-shrunk, not independent causal
effects).](results/dhs_rebuild/figure5_covariate_forest.png){width=90%}

- **Samples, structure and likelihood.** The attributable fraction at 30% PfPR
  stays near 22% across alternative samples (22.2% including sub-1% regions,
  23.6% on complete cases) and model-structure choices (23.3% with
  country-specific PfPR slopes re-admitted, 24.0% with regional covariates
  only); it rises on the 5–40% prevalence range (30.0%, a smaller sample of 637)
  and in a national-covariates-only model that fits far worse (31.0%, AIC +176).
  A log-Gaussian rate model gives +5.6% per 10 PfPR points (95% CI +0.4% to
  +10.9%; `p=0.034`), less precise than the count model but now excluding zero.

- **Country random effects.** Re-admitting the country-specific PfPR slope moves
  the national burden by at most 7.7% in any year, and including each in-sample
  country's own random effects rather than the population average moves it by at
  most 7.1% in the other direction.

## External validation: bednet trials

An independent triangulation against eight insecticide-treated-net/curtain
cluster-RCTs — which measured both parasite prevalence and child mortality —
gives a slope of the same order as the observational model: a weighted
meta-regression gives about an 8.6% mortality reduction per 10-point drop in
age-standardised PfPR2-10, against the model's linear 6.9% (the trial estimate
is imprecise — `p=0.12`, eight heterogeneous trials — and is wholly external to
the panel). Predicted and observed trial effects broadly track each other.

![RCT triangulation. **A:** observed under-five mortality reduction versus the
age-standardised prevalence reduction across trial arms (weighted fit). **B:**
mortality reductions predicted from the model's prevalence slope versus observed
(dotted line = 1:1).](results/dhs_rebuild/figure6_triangulation_rct.png){width=100%}

## National burden

### Latest-year comparison {#aggregate-burden}

The comparison is not perfectly contemporaneous: the model surface is evaluated
at 2024 using the latest available 2024 MAP PfPR, IGME all-cause mortality and
live births for 43 matched sub-Saharan African countries; IHME is 2024
under-five malaria deaths; WHO WMR 2025 reports all-age country deaths for 2024,
so the under-five comparator applies WHO's ~75% age-share statement. Three time
structures are extrapolated; the spline with a `ti` interaction fits best of the
three (AIC 8144.8, against 8145.1 for the time-stable spline and 8150.7 for the
`te` surface).

| Source | Deaths | 95% interval |
|---|---:|---:|
| Spline, `ti` interaction (best of three) | **527,000** | 373,000–658,000 |
| Full `te` surface | 503,000 | 304,000–666,000 |
| Spline, no time interaction | 435,000 | 291,000–565,000 |
| IHME/GBD SSA under-five, 2024 | 428,000 | 392,000–469,000 |
| WHO African-region all-age, 2024 | 579,000 | 531,000–706,000 |
| WHO African-region under-five proxy, 2024 | 434,000 | 398,000–530,000 |

The model intervals propagate both the fitted attributable-fraction uncertainty
and the all-cause post-neonatal mortality uncertainty (from the IHME/GBD
relative interval).

### Temporal trend

The historical trend is the quantity most sensitive to whether a time
interaction is admitted. Applied to national PfPR and all-cause mortality
series, estimated post-neonatal malaria deaths from 2000 to 2024 change by:

| Model | 2000 | 2024 | Change |
|---|---:|---:|---:|
| Spline `ti` interaction (best of three) | 553,000 | 527,000 | −5% |
| Full `te` surface | 486,000 | 503,000 | +4% |
| Spline, no time interaction | 925,000 | 435,000 | −53% |

The spread is the headline uncertainty in the trend, and it is wide: a
time-stable attributable fraction implies the burden more than halved over the
period, whereas the time-varying structures, which let the attributable fraction
rise as prevalence fell, imply little change. The two best structures are 0.3
AIC apart. The survey-grouped cross-validation above favours the time-varying
structures over the additive one, but modestly, and mainly through
high-prevalence regions. The trend should be read as a range, not a point.

![National malaria-attributable post-neonatal deaths over time under the three
time structures.](results/dhs_rebuild/time_surface_burden_timeseries.png){width=90%}

### Country differences from IHME {#country-differences}

The aggregate is about 23% above IHME (model ≈527,000 vs IHME ≈428,000 in 2024,
a ratio of 1.23), and the country-level differences are larger in both
directions (spline `ti` model; ten countries with the largest absolute
difference from IHME/GBD, all for 2024):

| Country | Model | IHME 2024 | Model/IHME |
|---|---:|---:|---:|
| DR Congo | 125,600 | 61,800 | 2.03 |
| Nigeria | 192,100 | 138,100 | 1.39 |
| Uganda | 14,300 | 26,100 | 0.55 |
| Burundi | 4,200 | 14,400 | 0.29 |
| Chad | 13,000 | 4,900 | 2.67 |
| Ethiopia | 4,300 | 8,900 | 0.48 |
| Tanzania | 5,400 | 9,900 | 0.55 |
| Mozambique | 14,200 | 10,300 | 1.38 |
| Cameroon | 13,100 | 16,300 | 0.80 |
| Ghana | 3,800 | 6,700 | 0.56 |

The model sits well above IHME in DR Congo, Nigeria and Chad and well below it
in Uganda, Burundi, Ethiopia and Tanzania. DR Congo alone accounts for about
64,000 of the 99,000 aggregate excess. (IHME/GBD malaria estimates are
unavailable for Sudan and South Africa, so the IHME total covers 41 of the 43
countries; both contribute negligibly.)

## Interpretation and limitations

The data provide fairly robust evidence of a positive cross-sectional
association between malaria prevalence and post-neonatal mortality that survives
adjustment for national vaccine rollout and child HIV prevalence, the neonatal
negative control, alternative samples, model-structure choices and likelihood.
The main remaining limitations are:

- the aggregate burden is not a validation against a common estimand — the model
  includes indirect malaria-associated deaths, whereas IHME and WHO assign
  deaths specifically to malaria, and national extrapolation uses
  population-average effects that ignore country random effects;
- the 12-month window buys alignment between outcome and exposure at the price
  of precision, and 56 region-years with no death of one type in the window are
  excluded from the shared sample;
- the historical trend is not well identified. Three specifications are
  separated by 0.44 AIC yet imply changes from −53% to +4% over 2000–2024, and
  the evidence for a time-varying effect, while consistent across held-out
  surveys, rests on high-prevalence regions in recent years where data are
  thin. The trend should be read as a range, not a point;
- 2025 exposure and all-cause mortality inputs are unavailable, and the WHO
  under-five comparison applies a regional 75% age-share to every country;
- child HIV prevalence is derived from UNAIDS 0–14 counts and a World Bank
  denominator rather than a published rate, and Nigeria and Comoros (about 4% of
  the sample, including the largest-burden country) lack an under-15 series and
  are imputed;
- the Bayesian fits fix the ridge penalty on the covariate block at its
  standardised scale rather than estimating it, and their smoothing parameters
  are weakly identified (posterior close to prior), which is why their intervals
  are wider than the `mgcv` ones;
- on the 60-month window the dose-response estimate fell monotonically as survey
  coverage was repaired (+8.5% to +7.3% per 10 PfPR points across four
  successive panels). That pattern is not noise-shaped and suggests the
  previously dropped surveys differed systematically; further coverage work
  could move it again.

**Changes from the 28 August version.** The mortality window moved from 60 to
12 months (sample 1,111 to 1,055); country-specific random PfPR slopes were
removed from the primary model; the neonatal control moved from marginally
positive to null; the AIC selection became a three-way tie rather than a clear
preference for the time interaction; Bayesian refits, a survey-grouped
cross-validation and era-by-region subgroups were added; the 2024 burden moved
from 546,000 to 527,000 and the country pattern against IHME is again mixed
rather than one-directional.

## Analysis flow {#analysis-flow}

![Data and analysis flow. External data (MAP PfPR2-10 rasters, UNICEF WUENIC
vaccine coverage, UNAIDS/World Bank child HIV prevalence, and World Bank
economic series) feed the assembled survey-region panel, which is filtered to
the primary sample, analysed with the ridge-penalised model, and extrapolated to
a national burden.](results/dhs_rebuild/study_flow_diagram.png){width=100%}

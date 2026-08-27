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
  percentage-point increase in PfPR2-10 is associated with about a 7.5% increase
  in post-neonatal mortality (ridge-linear summary; 95% CI 4.2%–10.9%;
  `p<0.001`). The estimate is essentially unchanged whether we use the regional
  covariate block alone or the full regional-plus-national block that now
  includes national vaccine coverage and child HIV prevalence.

- **Neonatal mortality, the negative control, shows no such association.** The
  ridge-linear neonatal estimate is +1.3% per 10 PfPR points (95% CI −0.7% to
  +3.4%; `p=0.20`), and the selected prevalence term for neonatal mortality is
  effectively flat. This, together with the stability under HIV and vaccine
  adjustment, reduces the concern that the main result simply reflects broad
  differences in child survival between higher- and lower-prevalence settings.

- **A prevalence-by-time interaction now fits best — a change from the earlier,
  smaller panel.** The linear PfPR term with a time interaction is selected by
  AIC ahead of the full `te(PfPR, year)` surface (ΔAIC 3.5), a spline
  time-interaction (ΔAIC 4.0), the additive spline `s(PfPR)+s(year)` (ΔAIC 6.1)
  and the linear no-interaction form (ΔAIC 12.9). The interaction is +0.0042 per
  year (`p=2.4e-05`), where on the 936-region-year panel it had been borderline
  and the time-stable spline was selected. Population-average
  malaria-attributable fractions of post-neonatal mortality are 7%, 22% and 34%
  at 10%, 30% and 50% PfPR2-10 (versus a 1% counterfactual), down from 14%, 34%
  and 40%. Because the historical trend depends strongly on whether a time
  interaction is admitted, this selection matters more than the modest change in
  the dose-response itself.

- **The preferred model estimates about 543,000 malaria-attributable
  post-neonatal deaths in 2024 — the same order as IHME and WHO.** The 95%
  uncertainty interval is 336,000–716,000, propagating both the fitted
  attributable-fraction (model-parameter) uncertainty and the all-cause
  post-neonatal mortality uncertainty (from the IHME/GBD relative interval); the
  2000 estimate is 640,000. Across the three time structures the 2024 total
  ranges from 477,000 to 547,000. The point estimate is about 27% higher than
  the IHME/GBD 2024 estimate (about 428,000) and about 22% higher than an
  approximate WHO under-five figure (about 446,000); the intervals
  overlap. This is not a
  validation against a common estimand: the model captures direct and indirect
  malaria-associated post-neonatal mortality, whereas IHME and WHO assign deaths
  specifically to malaria.

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
synthetic-cohort life-table method, and post-neonatal mortality is their
difference. The assembled panel has 1,125 region-years (119 surveys, 36
countries); the primary sample requires country-mean PfPR2-10 above 1% and
regional PfPR2-10 at least 1%, plus positive outcomes/exposure, giving **1,108
region-years (117 surveys, 35 countries)**. The full data flow is shown in the
[analysis-flow diagram](#analysis-flow).

**Panel coverage.** The panel is now assembled from the DHS survey registry
rather than inherited from an earlier aggregate. Three faults in that path were
repaired — a Windows-1250 mis-decoding of accented region names, an empty `v022`
that made `chmort` abort for four surveys, and a MAP regional extraction that
had never completed — and recode region names are now reconciled against
boundary names across language, word order, spelling, qualifier clauses and
granularity. Together these took the panel from 936 region-years (106 surveys)
to 1,125 (119 surveys).

This matters for interpretation, not only for coverage: the linear
no-interaction estimate fell at every step as coverage improved (+8.5%, +8.30%,
+7.61%, +7.33% per 10 PfPR points across the four successive panels), which
suggests the surveys previously dropped were not missing at random with respect
to the association. Six regions across two surveys still cannot be joined — Mali
2012 did not survey four of its regions, and Uganda 2016's Buganda/Central
labelling is genuinely ambiguous — and per-survey merge detail is recorded in
`results/dhs_rebuild/region_merge_quality.csv`.

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
measures). Political stability is now retained: the World Bank had archived the
`PV.EST` indicator, so the series was silently empty and the covariate was
dropped for 100% missingness; it is read from the current `GOV_WGI_PV_EST`
series. National DTP3, electricity access and health-expenditure-share are
dropped as redundant with their regional or per-capita counterparts.

Child HIV prevalence is derived from the UNAIDS 2025 estimates workbook (number
of children 0–14 living with HIV) divided by the World Bank 0–14 population
(`SP.POP.0014.TO`), joined by ISO3 and exact survey year, and entered on the log
scale because it spans a very wide range across countries (country means from
about 0.02% in Madagascar to about 3.8% in Eswatini). Nigeria and Comoros
publish only a 15–19 series (4% of the
sample) and are imputed with a provenance flag. Source choices are in
[`R_dhs/COVARIATE_SOURCES.md`](R_dhs/COVARIATE_SOURCES.md).

**Model.** Approximate death counts (rate × birth exposure) are modelled with a
negative-binomial GAM: `offset(log(exposure))`, the ridge covariate block,
country random intercepts and country-specific random linear PfPR slopes. Five
PfPR/time specifications are compared by ML/AIC on the same 1,108 observations,
and the selected form is refitted by REML for post-neonatal, all-under-five and
neonatal mortality.

## Association and model selection

From the ridge-linear summary models, each 10 percentage-point increase in
PfPR2-10 is associated with:

| Outcome | Change in mortality | 95% CI | p-value |
|---|---:|---:|---:|
| Post-neonatal | +7.49% | +4.17% to +10.90% | <0.001 |
| All under-five | +5.26% | +2.73% to +7.86% | <0.001 |
| Neonatal (negative control) | +1.34% | −0.72% to +3.45% | 0.20 |

![Modelled change in mortality versus 1% PfPR for the three outcomes. The
post-neonatal and all-under-five associations rise with prevalence; the
neonatal negative control stays flat and its interval spans
zero.](results/dhs_rebuild/figure3_negative_control_curves.png){width=88%}

Five specifications are compared by AIC on the 1,108-observation post-neonatal
sample:

| Model | Specification | AIC | ΔAIC |
|---|---|---:|---:|
| Linear, interaction | `PfPR + PfPR × year + s(year)` | 8556.01 | 0.00 |
| Full tensor surface | `te(PfPR, year)` | 8559.47 | 3.46 |
| Spline, `ti` interaction | `s(PfPR) + s(year) + ti(PfPR, year)` | 8560.01 | 4.00 |
| Spline, no interaction | `s(PfPR) + s(year)` | 8562.08 | 6.07 |
| Linear, no interaction | `PfPR + s(year)` | 8568.90 | 12.89 |

An independent AIC comparison for neonatal mortality selects an effectively
linear, flat prevalence term; it selects the linear no-interaction form, whose
prevalence coefficient is null (`p=0.21`), consistent with the negative-control
table above.

The selected specification gives the following malaria-attributable fractions of
post-neonatal mortality (country random effects set to zero; reference year
2013):

| PfPR2-10 | Attributable fraction | 95% interval |
|---:|---:|---:|
| 10% | 7.3% | 4.6%–9.9% |
| 30% | 21.7% | 14.2%–28.6% |
| 50% | 33.9% | 22.8%–43.4% |

## Robustness

Robustness is the central message: the association barely moves across covariate
blocks, samples, model structures and likelihood.

- **Covariate adjustment.** Child HIV prevalence is itself the strongest single
  predictor of post-neonatal mortality in the block (+14.1% per standard
  deviation), yet adding it leaves the malaria estimate essentially unchanged.
  The regional block alone and the full regional-plus-national block give a
  similar attributable fraction (22.5% and 21.5% at 30% PfPR).

![Ridge-standardized conditional associations for the 17 covariates in the
best-fitting post-neonatal model (percentage change per 1 SD after
transformation; ridge-shrunk, not independent causal
effects).](results/dhs_rebuild/figure5_covariate_forest.png){width=90%}

- **Samples, structure and likelihood.** The attributable fraction at 30% PfPR
  stays near 22% across alternative samples (21.4% including sub-1% regions,
  25.8% on the 5–40% range, 23.1% on complete cases) and across model-structure
  choices (18.7% without country-specific PfPR slopes, 22.5% with regional
  covariates only); it rises materially only in a national-covariates-only model
  that fits far worse (33.6%). A log-Gaussian rate model is less precise and its
  interval spans zero (+5.3% per 10 PfPR points, 95% CI −1.6% to +12.8%;
  `p=0.14`).

## External validation: bednet trials

An independent triangulation against eight insecticide-treated-net/curtain
cluster-RCTs — which measured both parasite prevalence and child mortality —
gives a slope of the same order as the observational model: a weighted
meta-regression gives about an 8.6% mortality reduction per 10-point drop in
age-standardised PfPR2-10, against the model's ridge-linear 7.5% (the trial
estimate is imprecise — `p=0.12`, eight heterogeneous trials — and it is
unchanged by the panel rebuild, being wholly external to it). Predicted and observed trial
effects broadly track each other.

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
so the under-five comparator applies WHO's ~75% age-share statement.

| Source | Deaths |
|---|---:|
| Selected model (full `te` surface) | **543,000** |
| Spline, `ti` interaction | 547,000 |
| Spline, no time interaction | 477,000 |
| IHME/GBD SSA under-five, 2024 | 428,000 |
| WHO African-region all-age, 2024 | 579,000 |
| WHO African-region under-five proxy, 2024 | 434,000 |

### Temporal trend

The historical trend is the quantity most sensitive to whether a time
interaction is admitted, and on the rebuilt panel the interaction is now
favoured. Applied to national PfPR and all-cause mortality series, estimated
post-neonatal malaria deaths from 2000 to 2024 fall by:

| Model | 2000 | 2024 | Change |
|---|---:|---:|---:|
| Full `te` surface (selected) | 640,000 | 543,000 | −15% |
| Spline `ti` interaction | 693,000 | 546,000 | −21% |
| Spline, no time interaction | 961,000 | 476,000 | −50% |

The spread across these three structures is the headline uncertainty in the
trend, and it is wide: a time-stable attributable fraction implies the burden
halved over the period, whereas the time-varying structures — which let the
attributable fraction rise gradually — imply a decline of only 15–21%. On the
936-region-year panel the time-stable form was selected and the −52% figure was
the headline; on the fuller panel it is the worst-fitting of the five
specifications, so the smaller declines should now carry more weight.

![National malaria-attributable post-neonatal deaths over time under the three
time structures.](results/dhs_rebuild/time_surface_burden_timeseries.png){width=90%}

### Country differences from IHME {#country-differences}

The aggregate is now about 27% above IHME (model ≈543,000 vs IHME ≈428,000 in
2024, a ratio of 1.27), and the country-level differences are larger still
(AIC-selected full `te` surface; ten countries with the largest absolute
difference from IHME/GBD, all for 2024):

| Country | Model | IHME 2024 | Model/IHME |
|---|---:|---:|---:|
| DR Congo | 126,000 | 61,800 | 2.05 |
| Nigeria | 200,000 | 138,100 | 1.45 |
| Chad | 13,400 | 4,900 | 2.73 |
| Mozambique | 14,900 | 10,300 | 1.44 |
| Benin | 9,900 | 7,400 | 1.34 |
| Guinea | 9,500 | 7,100 | 1.35 |
| Niger | 22,700 | 20,300 | 1.11 |
| Central African Republic | 6,100 | 3,900 | 1.58 |
| South Sudan | 6,600 | 4,800 | 1.38 |
| Togo | 2,800 | 1,900 | 1.48 |

The pattern differs from the earlier panel in a way worth noting: previously the
model's excess over IHME in DR Congo and Nigeria was offset by lower estimates
elsewhere, which kept the aggregates close. On the rebuilt panel the ten largest
differences all run in the same direction, so the aggregate sits well above IHME
rather than near it. DR Congo alone accounts for about 65,000 of the 115,000
aggregate excess. (IHME/GBD malaria estimates are unavailable for Sudan and
South Africa, so the IHME total covers 41 of the 43 countries; both contribute
negligibly.)

## Interpretation and limitations

The data provide fairly robust evidence of a positive cross-sectional
association between malaria prevalence and post-neonatal mortality that survives
adjustment for national vaccine rollout and child HIV prevalence, the neonatal
negative control, alternative samples, and model-structure choices. The main
remaining limitations are:

- the aggregate burden is not a validation against a common estimand — the model
  includes indirect malaria-associated deaths, whereas IHME and WHO assign
  deaths specifically to malaria, and national extrapolation uses
  population-average effects that ignore country random effects;
- 2025 exposure and all-cause mortality inputs are unavailable, and the WHO
  under-five comparison applies a regional 75% age-share to every country;
- child HIV prevalence is derived from UNAIDS 0–14 counts and a World Bank
  denominator rather than a published rate, and Nigeria and Comoros (about 4% of
  the sample, including the largest-burden country) lack an under-15 series and
  are imputed;
- the log-Gaussian likelihood sensitivity is positive but its interval spans
  zero, so the exact effect size remains uncertain;
- the historical trend is not well identified. The five specifications are
  separated by only ~13 AIC points yet imply declines from 15% to 50% over
  2000–2024, and the selection flipped from the time-stable spline to a
  time-interaction form when the panel was completed. The trend should be read
  as a range, not a point;
- the dose-response estimate fell monotonically as survey coverage was repaired
  (+8.5% to +7.3% per 10 PfPR points across four successive panels). That
  pattern is not noise-shaped and suggests the previously dropped surveys
  differed systematically; it is possible that further coverage work would move
  it again.

## Analysis flow {#analysis-flow}

![Data and analysis flow. External data (MAP PfPR2-10 rasters, UNICEF WUENIC
vaccine coverage, UNAIDS/World Bank child HIV prevalence, and World Bank
economic series) feed the assembled survey-region panel, which is filtered to
the primary sample, analysed with the ridge-penalised model, and extrapolated to
a national burden.](results/dhs_rebuild/study_flow_diagram.png){width=100%}

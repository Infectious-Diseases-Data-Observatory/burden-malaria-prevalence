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
  percentage-point increase in PfPR2-10 is associated with about a 8.6% increase
  in post-neonatal mortality (ridge-linear summary; 95% CI 5.7%–11.7%;
  `p<0.001`). The estimate is essentially unchanged whether we use the regional
  covariate block alone or the full regional-plus-national block that now
  includes national vaccine coverage and child HIV prevalence.

- **Neonatal mortality, the negative control, shows no such association.** The
  ridge-linear neonatal estimate is +1.0% per 10 PfPR points (95% CI −1.2% to
  +3.2%; `p=0.37`), and the selected prevalence term for neonatal mortality is
  effectively flat. This, together with the stability under HIV and vaccine
  adjustment, reduces the concern that the main result simply reflects broad
  differences in child survival between higher- and lower-prevalence settings.

- **A time-stable nonlinear prevalence-mortality relationship fits best.** The
  additive spline `s(PfPR)+s(year)` is selected by AIC over a spline
  time-interaction (ΔAIC 4.9; interaction `p=0.16`), a linear form (ΔAIC 10.2),
  and a full `te(PfPR, year)` surface (ΔAIC 11.2). It implies population-average
  malaria-attributable fractions of post-neonatal mortality of 14%, 34% and 40%
  at 10%, 30% and 50% PfPR2-10 (versus a 1% counterfactual).

- **The preferred model estimates about 482,000 malaria-attributable
  post-neonatal deaths for the latest period — close to IHME.** The
  model-parameter 95% interval is 348,000–606,000. The point estimate is about
  3% higher than the IHME/GBD 2025 estimate and about 11% higher than an
  approximate WHO under-five figure; the intervals overlap. This is not a
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
difference. The assembled panel has 936 region-years (106 surveys, 35
countries); the primary sample requires country-mean PfPR2-10 above 1% and
regional PfPR2-10 at least 1%, plus positive outcomes/exposure, giving **921
region-years (105 surveys, 34 countries)**. The full data flow is shown in the
[analysis-flow diagram](#analysis-flow).

**Covariates.** Twenty-five candidate covariates are screened before imputation;
those with at most 5% missingness are retained and singly imputed (country
median, overall-median fallback). The **17 retained covariates** enter a single
standardised, ridge-penalised block (proportions logit-transformed):

| Level | Covariates |
|---|---|
| Survey-region (DHS) | urban residence, DTP3, measles, facility delivery, maternal education, exclusive breastfeeding, short birth interval, maternal age at first birth, improved water, improved sanitation, household electricity |
| National, exact year | WUENIC Hib3, PCV-completion, rotavirus-completion; and child (0–14) HIV prevalence (log) |
| National, nearest year | log GDP per capita, log health expenditure per capita |

Eight candidates are excluded for excess missingness (the three direct DHS
vaccine completion recodes, household wealth, the three anthropometry measures,
and political stability). National DTP3, electricity access and
health-expenditure-share are dropped as redundant with their regional or
per-capita counterparts.

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
PfPR/time specifications are compared by ML/AIC on the same 921 observations,
and the selected form is refitted by REML for post-neonatal, all-under-five and
neonatal mortality.

## Association and model selection

From the ridge-linear summary models, each 10 percentage-point increase in
PfPR2-10 is associated with:

| Outcome | Change in mortality | 95% CI | p-value |
|---|---:|---:|---:|
| Post-neonatal | +8.65% | +5.69% to +11.68% | <0.001 |
| All under-five | +6.37% | +4.01% to +8.79% | <0.001 |
| Neonatal (negative control) | +0.99% | −1.18% to +3.21% | 0.37 |

![Modelled change in mortality versus 1% PfPR for the three outcomes. The
post-neonatal and all-under-five associations rise with prevalence; the
neonatal negative control stays flat and its interval spans
zero.](results/dhs_rebuild/figure3_negative_control_curves.png){width=88%}

Five specifications are compared by AIC on the 921-observation post-neonatal
sample:

| Model | Specification | AIC | ΔAIC |
|---|---|---:|---:|
| Spline, no interaction | `s(PfPR) + s(year)` | 7100.01 | 0.00 |
| Spline, `ti` interaction | `s(PfPR) + s(year) + ti(PfPR, year)` | 7104.92 | 4.91 |
| Linear, interaction | `PfPR + PfPR × year + s(year)` | 7109.24 | 9.23 |
| Linear, no interaction | `PfPR + s(year)` | 7110.21 | 10.19 |
| Full tensor surface | `te(PfPR, year)` | 7111.23 | 11.21 |

An independent AIC comparison for neonatal mortality selects an effectively
linear, flat prevalence term; in the prespecified REML refit the neonatal PfPR
smooth is linear (edf 1.00) and null (`p=0.37`), consistent with the
negative-control table above.

The selected time-stable spline gives the following malaria-attributable
fractions of post-neonatal mortality (country random effects set to zero;
reference year 2014):

| PfPR2-10 | Attributable fraction | 95% interval |
|---:|---:|---:|
| 10% | 14.1% | 7.9%–20.0% |
| 30% | 34.4% | 26.0%–41.9% |
| 50% | 40.2% | 30.8%–48.4% |

## Robustness

Robustness is the central message: the association barely moves across covariate
blocks, samples, model structures and likelihood.

- **Covariate adjustment.** Child HIV prevalence is itself a strong predictor of
  post-neonatal mortality (+15.0% per standard deviation, `p<0.001`), yet adding
  it leaves the malaria estimate essentially unchanged and only weakens the
  PfPR-by-time interaction (`p` 0.08→0.16). The regional block alone and the
  full regional-plus-national block give a similar attributable fraction (about
  34–36% at 30% PfPR).

![Ridge-standardized conditional associations for the 17 covariates in the
best-fitting post-neonatal model (percentage change per 1 SD after
transformation; ridge-shrunk, not independent causal
effects).](results/dhs_rebuild/figure5_covariate_forest.png){width=90%}

- **Samples, structure and likelihood.** The attributable fraction at 30% PfPR
  stays near 34% across alternative samples (including sub-1% regions, the 5–40%
  range, and complete cases) and across model-structure choices; it rises
  materially only in a national-covariates-only model that fits far worse. A
  log-Gaussian rate model is less precise but still positive (+5.7% per 10 PfPR
  points; `p=0.14`).

## External validation: bednet trials

An independent triangulation against eight insecticide-treated-net/curtain
cluster-RCTs — which measured both parasite prevalence and child mortality —
recovers essentially the same slope as the observational model: a weighted
meta-regression gives about an 8.6% mortality reduction per 10-point drop in
age-standardised PfPR2-10, matching the model's ridge-linear slope (though
imprecise — `p=0.12`, eight heterogeneous trials). Predicted and observed trial
effects broadly track each other.

![RCT triangulation. **A:** observed under-five mortality reduction versus the
age-standardised prevalence reduction across trial arms (weighted fit). **B:**
mortality reductions predicted from the model's prevalence slope versus observed
(dotted line = 1:1).](results/dhs_rebuild/figure6_triangulation_rct.png){width=100%}

## National burden

### Latest-year comparison {#aggregate-burden}

The comparison is not perfectly contemporaneous: the model surface is evaluated
at 2025 using the latest available 2024 MAP PfPR, IGME all-cause mortality and
live births for 43 matched sub-Saharan African countries; IHME is 2025
under-five malaria deaths; WHO WMR 2025 reports all-age country deaths for 2024,
so the under-five comparator applies WHO's ~75% age-share statement.

| Source | Deaths |
|---|---:|
| Selected time-stable model | **482,000** |
| IHME/GBD SSA under-five, 2025 | 470,000 |
| WHO African-region all-age, 2024 | 579,000 |
| WHO African-region under-five proxy, 2024 | 434,000 |

### Temporal trend

The additive spline is time-stable, but the time-varying alternatives — although
not selected — remain the main sensitivity for the historical trend. Applied to
national PfPR and all-cause mortality series, estimated post-neonatal malaria
deaths from 2000 to 2024 fall by:

| Model | 2000 | 2024 | Change |
|---|---:|---:|---:|
| No interaction (selected) | 1,012,000 | 482,000 | −52% |
| Spline `ti` interaction | 878,000 | 527,000 | −40% |
| Full `te` surface | 802,000 | 551,000 | −31% |

Under the selected model the attributable fraction is time-stable; the
time-varying alternatives instead let it rise gradually over the period.
Adjusting for child HIV prevalence — itself a large time-varying driver of child
mortality through ART rollout — brings the time-varying declines closer to the
time-stable one than in the unadjusted analysis.

![National malaria-attributable post-neonatal deaths over time under the three
time structures.](results/dhs_rebuild/time_surface_burden_timeseries.png){width=90%}

### Country differences from IHME {#country-differences}

The aggregate agreement with IHME masks large country-level differences (AIC-selected
time-stable spline; "very different" = at least two-fold and ≥1,000 deaths):

| Country | Model | IHME 2025 | Model/IHME |
|---|---:|---:|---:|
| DR Congo | 110,000 | 68,300 | 1.61 |
| Nigeria | 178,000 | 150,500 | 1.18 |
| Chad | 12,300 | 5,500 | 2.23 |
| Sudan | 3,900 | 1,400 | 2.89 |
| Uganda | 13,400 | 28,500 | 0.47 |
| Burundi | 3,900 | 15,700 | 0.25 |
| Tanzania | 5,100 | 10,700 | 0.47 |
| Ethiopia | 4,000 | 9,300 | 0.42 |
| Ghana | 3,500 | 7,300 | 0.49 |
| Rwanda | 300 | 2,500 | 0.11 |

DR Congo and Nigeria account for most of the model's upward difference from
IHME, offset by lower estimates elsewhere; this is why the aggregate totals are
close despite the country-level spread.

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
- the log-Gaussian likelihood sensitivity is positive but imprecise, so the
  exact effect size and the temporal trend remain uncertain.

## Analysis flow {#analysis-flow}

![Data and analysis flow. External data (MAP PfPR2-10 rasters, UNICEF WUENIC
vaccine coverage, UNAIDS/World Bank child HIV prevalence, and World Bank
economic series) feed the assembled survey-region panel, which is filtered to
the primary sample, analysed with the ridge-penalised model, and extrapolated to
a national burden.](results/dhs_rebuild/study_flow_diagram.png){width=100%}

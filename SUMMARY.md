---
title: "DHS/MIS malaria prevalence and child mortality analysis"
output:
  html_document:
    self_contained: true
  pdf_document:
    latex_engine: xelatex
  word_document: default
---

# DHS/MIS malaria prevalence and child mortality analysis

## Executive summary

The rebuilt analysis uses 921 DHS/MIS survey-region-years from 105 surveys in
34 countries. The primary negative-binomial ridge-GAM finds a strong positive
association between MAP PfPR2-10 and all-cause post-neonatal mortality, with no
corresponding association in the neonatal negative control.

The best-fitting specification is:

```text
deaths ~ s(PfPR) + ti(PfPR, calendar year) + s(calendar year)
       + ridge-penalised covariates
       + country random intercept
       + country random PfPR slope
       + offset(log birth exposure))
```

This ANOVA-style `ti` model has lower AIC than the full `te(PfPR, year)` surface
(ΔAIC 5.91) and the spline model without a time interaction (ΔAIC 3.30).
However, the evidence for the nonlinear `ti` interaction itself is borderline
(`p=0.057`). The model-selection result therefore supports allowing the
prevalence-mortality relationship to vary over time, but it does not establish
the exact form of that change with high confidence.

The downstream 2025-surface estimate is approximately **611,000
malaria-attributable post-neonatal deaths** across 43 matched sub-Saharan
African countries (model-parameter 95% interval 441,000–749,000). This estimate
uses 2024 MAP prevalence, IGME mortality and live-birth inputs because 2025
values are not yet available for those inputs. It is approximately:

- 30% higher than the IHME/GBD 2025 SSA estimate of 470,000 under-five malaria
  deaths;
- 41% higher than the approximate WHO under-five total of 434,000, obtained by
  multiplying the WHO African-region 2024 all-age total of 579,000 by 75%.

The model interval overlaps both comparators. It covers uncertainty in the
fitted DHS/MIS relationship, but not uncertainty in MAP prevalence, IGME
all-cause mortality or live births.

## Analysis structure

### Data selection

The unit of analysis is a DHS/MIS survey-region-year.

1. Eligible DHS and MIS Births Recodes are surveys conducted from 2000 through
   2025 in sub-Saharan Africa. AIS surveys are not included.
2. Survey regions are matched to DHS first-level survey boundaries.
3. Each region receives the population-weighted MAP PfPR2-10 estimate from the
   survey's calendar year. A missing 2025 MAP surface is flagged; 2024 is not
   silently substituted.
4. Regional all-under-five mortality and neonatal mortality are estimated from
   birth histories using the DHS synthetic-cohort life-table method.
   Post-neonatal mortality is calculated as U5MR minus NNMR. The birth-exposure
   denominator is retained for the count model offset.
5. The assembled validation panel contains 936 survey-region-years from 106
   surveys and 35 countries. The shared primary sample requires:
   - country mean PfPR2-10 greater than 1%;
   - regional PfPR2-10 at least 1%;
   - positive finite birth exposure, U5MR, NNMR and post-neonatal mortality.
6. The resulting primary sample contains 921 region-years, 105 surveys and 34
   countries, spanning 2000–2024.

Both the country-level and region-level prevalence flags remain in the analysis
dataset. Sensitivity analyses separately include sub-1% regions and restrict
the prevalence range to 5–40%.

### Covariate processing

Twenty-one candidate covariates are evaluated before imputation. Variables
with at most 5% missingness are retained and singly imputed using the
within-country median, with the overall median as a fallback. Variables with
more than 5% missingness are flagged and excluded.

The 16 included covariates are:

| Level | Included covariates |
|---|---|
| Survey-region | Urban residence, regional DTP3 coverage, measles vaccination, facility delivery, maternal education, exclusive breastfeeding, birth interval under 24 months, maternal age at first birth, improved water, improved sanitation and household electricity |
| National, nearest year | DTP3 coverage, log GDP per capita, health expenditure as a share of GDP, log health expenditure per capita and electricity access |

The five excluded variables are household wealth (5.56% missing), stunting
(25.2%), underweight (27.5%), wasting (27.5%) and political stability (100%).

Proportions are logit transformed, continuous variables remain on their
declared scale, and all included variables are standardised. They enter the
model as one ridge-penalised matrix block; the amount of shrinkage is estimated
from the data.

### Main outcome model

For each survey region, the mortality rate is converted to an approximate death
count using the DHS birth exposure. Models use:

- a negative-binomial likelihood with log link;
- `offset(log(exposure))`;
- a ridge penalty on the full eligible covariate block;
- a smooth calendar-year term;
- country random intercepts;
- country-specific random linear PfPR slopes.

Four prespecified models are compared using maximum likelihood and AIC on the
same 921 observations:

| Model | Prevalence and time specification | AIC | ΔAIC |
|---|---|---:|---:|
| Spline with `ti` interaction | `s(PfPR) + s(year) + ti(PfPR, year)` | 7137.68 | 0.00 |
| Linear with interaction | `PfPR + PfPR × year + s(year)` | 7139.85 | 2.17 |
| Spline without interaction | `s(PfPR) + s(year)` | 7140.98 | 3.30 |
| Linear without interaction | `PfPR + s(year)` | 7149.66 | 11.99 |

The selected specification is refitted with REML for post-neonatal mortality,
all-under-five mortality and neonatal mortality. The linear no-interaction
model is also refitted for all three outcomes to give an interpretable
percentage change per 10 PfPR percentage points.

### Full `te(PfPR, year)` check

A fifth model was fitted as a sensitivity analysis:

```text
deaths ~ te(PfPR, year)
       + ridge covariates
       + country random intercept and PfPR slope
       + offset(log exposure))
```

The full tensor includes the PfPR main effect, year main effect and their
interaction in one surface. It has AIC 7143.59 (ΔAIC 5.91) with 11.46 effective
degrees of freedom. It is more flexible but is not preferred to the more
parsimonious `s(PfPR) + s(year) + ti(PfPR, year)` decomposition.

## Main results

### DHS/MIS association

From the ridge-linear summary models, each 10 percentage-point increase in
PfPR2-10 is associated with:

| Outcome | Change in mortality | 95% CI | p-value |
|---|---:|---:|---:|
| Post-neonatal | +9.04% | +5.98% to +12.19% | <0.001 |
| All under-five | +6.63% | +4.27% to +9.04% | <0.001 |
| Neonatal negative control | +0.95% | −1.19% to +3.14% | 0.388 |

The neonatal result supports the interpretation that the main association is
not simply a marker of uniformly worse child survival in higher-prevalence
settings.

### Attributable fraction

The selected `ti` model estimates the following population-average fractions
of post-neonatal mortality attributable to malaria relative to a 1% PfPR
counterfactual:

| PfPR2-10 | 2014 reference year | 2025 surface |
|---:|---:|---:|
| 10% | 14.6% | 18.6% |
| 30% | 35.6% | 43.0% |
| 50% | 43.0% | 51.1% |

Country random effects are set to zero for these population-average
predictions and for the national burden estimates.

![Attributable-fraction curves under the alternative time
structures](results/dhs_rebuild/time_surface_af_comparison.png)

## Does the time interaction materially affect predictions?

Yes, particularly for historical trends and the latest aggregate burden, but
the alternatives remain within overlapping uncertainty intervals.

At 30% PfPR, the attributable fraction is:

| Model | 2000 | 2014 | 2024 | 2025 |
|---|---:|---:|---:|---:|
| No time interaction | 35.7% | 35.7% | 35.7% | 35.7% |
| Selected `ti` interaction | 24.8% | 35.6% | 42.4% | 43.0% |
| Full `te` surface | 28.4% | 32.2% | 40.3% | 39.1% |

Applied to the same national prevalence and all-cause mortality series, the
estimated post-neonatal malaria burden changes from 2000 to 2024 as follows:

| Model | 2000 deaths | 2024 deaths | Change |
|---|---:|---:|---:|
| No interaction | 1,049,000 | 503,000 | −52.1% |
| Selected `ti` interaction | 753,000 | 601,000 | −20.2% |
| Full `te` surface | 867,000 | 570,000 | −34.2% |

This is the most important sensitivity for the manuscript's temporal bottom
line. The claim of an approximately 50–57% decline follows the time-stable
relationship. The AIC-selected time-varying relationship implies a much
smaller decline of about 20%, while the full `te` surface lies between them.

![National burden trajectories under the alternative time
structures](results/dhs_rebuild/time_surface_burden_timeseries.png)

For the latest 2025-surface comparison, the estimates are:

| Model | Estimated deaths | Model-parameter 95% interval |
|---|---:|---:|
| No interaction | 503,000 | 370,000–621,000 |
| Selected `ti` interaction | 611,000 | 441,000–749,000 |
| Full `te` surface | 551,000 | 218,000–804,000 |

Thus the selected interaction raises the latest point estimate by about
108,000 deaths, or 21%, compared with no interaction. The full `te` estimate is
about 60,000 lower than the selected `ti` estimate.

## Other sensitivity analyses

### Alternative samples

All sample sensitivities retain a positive nonlinear prevalence relationship:

- Primary sample: AF at 30% PfPR = 35.6%.
- Include regions below 1% PfPR: 33.8%.
- Restrict PfPR to 5–40%: 32.5%.
- Complete cases only: 35.6%.

These choices do not materially change the cross-sectional attributable
fraction.

### Model structure and covariate blocks

- Removing the country random PfPR slope worsens fit by about 16.7 AIC units.
- Regional covariates alone have the lowest AIC in the structural sensitivity
  set, 3.1 units below the full regional-plus-national block.
- National covariates alone fit very poorly (ΔAIC 219).
- Changing the PfPR spline basis from `k=6` to `k=4` or `k=8` changes AF at 30%
  only from about 34.3% to 35.1%.
- Across structural sensitivities, AF at 30% ranges from 33.2% to 42.8%; the
  upper value occurs in the national-covariate-only model, which fits poorly.

### Alternative likelihood

The log-Gaussian rate sensitivity is less precise than the primary
negative-binomial count model:

- ridge-linear post-neonatal estimate: +7.13% per 10 PfPR points
  (95% CI −0.60% to +15.47%; `p=0.072`);
- nonlinear prevalence term: `p=0.049`;
- its post-neonatal time interaction is not supported (`p=0.599`).

This weakens certainty about the exact effect size and the time interaction,
but the point estimate remains positive.

## Latest national burden comparison

### Aggregate bottom line

The closest available comparison is not perfectly contemporaneous:

- **Model:** surface evaluated at 2025 using 2024 MAP PfPR, IGME mortality and
  live births for 43 matched sub-Saharan African countries.
- **IHME:** 2025 under-five malaria deaths.
- **WHO:** WMR 2025 estimates for 2024. WHO reports all-age country deaths, so
  the under-five comparator is an approximate 75% share.

| Source | Deaths |
|---|---:|
| Selected `ti` prevalence/all-cause model | **611,000** |
| IHME/GBD SSA under-five, 2025 | 470,000 |
| WHO African-region all-age, 2024 | 579,000 |
| WHO African-region under-five proxy, 2024 | 434,000 |

WHO's 2025 annex explicitly supplies country estimates through 2024, not 2025.
WHO also states that children under five account for approximately 75% of
malaria deaths in the African Region. Sources: [WMR 2025 Annex
4H](https://www.who.int/publications/m/item/annexes-world-malaria-report-2025),
[WHO malaria fact
sheet](https://www.who.int/news-room/fact-sheets/detail/malaria%E2%9C%85), and
the WHO GHO `MALARIA_EST_DEATHS` country indicator.

### Countries with the largest differences from IHME

The table below uses the AIC-selected `ti` model. “Very different” is defined in
the generated results as at least a two-fold difference and an absolute
difference of at least 1,000 deaths.

| Country | Model | IHME 2025 | Model/IHME | Interpretation |
|---|---:|---:|---:|---|
| Nigeria | 224,600 | 150,500 | 1.49 | Largest absolute excess: +74,100 |
| DR Congo | 137,300 | 68,300 | 2.01 | Model approximately twice IHME |
| Chad | 15,900 | 5,500 | 2.89 | Model substantially higher |
| Sudan | 5,200 | 1,400 | 3.88 | Model substantially higher |
| South Sudan | 7,400 | 4,300 | 1.73 | Model higher, below two-fold threshold |
| Uganda | 17,000 | 28,500 | 0.60 | Model lower by about 11,500 |
| Burundi | 4,900 | 15,700 | 0.31 | Model substantially lower |
| Ethiopia | 5,300 | 9,300 | 0.57 | Model lower |
| Rwanda | 400 | 2,500 | 0.15 | Model substantially lower |

Against the approximate WHO country-level under-five proxy, the clearest
two-fold discrepancies are:

- higher under the model: DR Congo;
- lower under the model: Tanzania, Ethiopia, Rwanda, Madagascar and Kenya.

Nigeria and DR Congo together account for most of the model's aggregate excess
over both IHME and the WHO proxy.

### Interpretation

The comparison is not a validation against a common estimand:

- the prevalence/all-cause model is intended to include direct and indirect
  malaria-associated post-neonatal deaths;
- IHME and WHO estimate deaths assigned specifically to malaria;
- the model uses population-average effects and ignores country random effects
  for national extrapolation;
- 2025 exposure and all-cause mortality inputs are unavailable, and WHO's
  latest estimates are for 2024;
- the WHO country under-five comparison applies a regional 75% age-share
  assumption to every country.

The robust conclusion is that malaria prevalence predicts post-neonatal
mortality and implies a large contemporary burden. The exact total and,
especially, the amount of decline since 2000 depend materially on whether the
prevalence-mortality relationship is allowed to change over calendar time.

## Reproducible outputs

- `R_dhs/10_time_surface_and_burden.R` fits the `te` sensitivity and produces
  all time-surface and national comparison tables.
- `results/dhs_rebuild/time_surface_model_comparison.csv` contains the AIC
  comparison.
- `results/dhs_rebuild/time_surface_af_anchors.csv` and
  `time_surface_burden_timeseries.csv` contain the time-interaction comparisons.
- `results/dhs_rebuild/latest_burden_totals.csv` contains aggregate results.
- `results/dhs_rebuild/latest_country_burden_comparison.csv` contains the
  country-level model, IHME and WHO comparisons.

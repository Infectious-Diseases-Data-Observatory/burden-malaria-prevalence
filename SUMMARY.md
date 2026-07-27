---
title: "DHS/MIS malaria prevalence and child mortality analysis"
output:
  html_document:
    self_contained: true
  pdf_document:
    latex_engine: xelatex
  word_document: default
header-includes:
  - \usepackage{float}
  - \floatplacement{figure}{H}
geometry: margin=0.9in
---

## Summary of main results

- **Higher malaria prevalence is associated with higher post-neonatal
  mortality.** In the ridge-linear model, a 10 percentage-point increase in
  PfPR2-10 is associated with 8.8% higher post-neonatal mortality (95% CI:
  5.8%–11.8%). The estimate for neonatal mortality, our negative-control
  outcome, is small and statistically consistent with no association.
  [See association results](#association-results) and
  [neonatal model selection.](#neonatal-selection)

- **After adding national vaccine coverage, a time-stable nonlinear
  prevalence-mortality relationship fits best.** The additive spline has lower
  AIC than the spline interaction (ΔAIC 1.81) and the full
  `te(PfPR, year)` surface (ΔAIC 11.66). The spline interaction has `p=0.078`.
  [See model comparison.](#model-comparison)

- **Time interactions remain an important sensitivity for the historical
  trend, even though they are not selected.** From 2000 to 2024, estimated
  deaths fall by 53% with no interaction, 25% with the spline `ti` interaction
  and 32% with the full `te` surface. [See time-interaction
  results.](#time-interaction-results)

- **Our preferred model estimates about 477,000 malaria-attributable
  post-neonatal deaths for the latest period, but this estimate is uncertain
  and is not directly comparable to IHME or WHO.** The model-parameter 95%
  interval is 342,000–593,000. The point estimate is about 1% higher than the
  IHME/GBD 2025 estimate and 10% higher than an approximate WHO under-five
  estimate; the interval overlaps both. The calculation evaluates the
  time-stable relationship using 2024 prevalence, mortality and birth inputs.
  [See aggregate burden comparison.](#aggregate-burden)

- **The aggregate agreement with IHME masks large country differences.**
  The model is higher for DR Congo, Chad and Sudan, but lower for Uganda,
  Burundi, Tanzania, Ethiopia, Ghana and Rwanda. Nigeria is about 17% higher
  than IHME. [See country-level comparison.](#country-differences)

- **Adding WUENIC vaccine coverage changes the temporal conclusion more than
  the cross-sectional association.** Exact-year national Hib3, PCV-completion
  and rotavirus-completion coverage have no missing values in the analysis
  panel and enter the ridge block. The direct DHS vaccine fields remain too
  sparse. [See covariate processing.](#covariates)

- **Our current bottom line is that the positive cross-sectional association
  remains fairly robust, while likelihood choice remains an important
  limitation.** Alternative samples produce similar attributable fractions,
  but the log-Gaussian sensitivity is positive and statistically imprecise.
  [See sensitivity analyses](#sensitivity-results) and [interpretation and
  limitations.](#interpretation-limitations)

- **The vaccine-adjusted result brings the preferred aggregate estimate close
  to IHME, but does not validate either estimate.** The model measures direct
  and indirect malaria-associated post-neonatal mortality, whereas IHME and
  WHO assign deaths specifically to malaria.
  [See interpretation and limitations.](#interpretation-limitations)

The sections below describe the data and methods underlying these claims,
followed by the detailed results and sensitivity analyses.

## Analysis structure

### Data selection {#data-and-sample}

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

### Covariate processing {#covariates}

Twenty-seven candidate covariates are evaluated before imputation. Variables
with at most 5% missingness are retained and singly imputed using the
within-country median, with the overall median as a fallback. Variables with
more than 5% missingness are flagged and excluded.

The 19 included covariates are:

| Level | Included covariates |
|---|---|
| Survey-region | Urban residence, regional DTP3 coverage, measles vaccination, facility delivery, maternal education, exclusive breastfeeding, birth interval under 24 months, maternal age at first birth, improved water, improved sanitation and household electricity |
| National, exact survey year | WUENIC Hib3, PCV-completion and final-dose rotavirus coverage |
| National, nearest year | DTP3 coverage, log GDP per capita, health expenditure as a share of GDP, log health expenditure per capita and electricity access |

The eight excluded variables are direct regional pentavalent-dose-3 coverage
(55.1% missing), PCV-dose-3 coverage (57.1%), completion of the
survey-specific rotavirus schedule (59.0%), household wealth (5.56%),
stunting (25.2%), underweight (27.5%), wasting (27.5%) and political stability
(100%).

The direct vaccine measures use standard DHS vaccination recodes:
`h51`–`h53` for pentavalent, `h54`–`h56` for pneumococcal vaccine and
`h57`–`h59` for rotavirus. They estimate weighted coverage among living
children aged 12–23 months from vaccination-card or maternal-report data.
These questions measure receipt, not the national introduction date, and are
concentrated in newer surveys. Missing direct-DHS fields are not coded as zero.

The national series come from UNICEF's 2025-revision WUENIC estimates:
`IM_HIB3`, `IM_PCVC` and `IM_ROTAC`. They are joined by ISO3 and exact survey
year. Reported zeros are preserved; years before a country's first estimate
and countries with no series are coded zero as not-yet-introduced, with a
separate status flag. Gaps after a series starts remain missing. All 936 rows
have values for the three national series, so they enter the ridge block
without imputation.

The local UNICEF dataflow does not contain child HIV prevalence. Its only HIV
epidemiology series is a count of 10–19-year-olds living with HIV for
2010–2024, so it is not relabelled or merged as under-five prevalence. UNAIDS
AIDSinfo remains the appropriate source for country-year paediatric prevalence.
Historical subnational SMC coverage is not currently available as one public
machine-readable panel. The source choices are documented in
[`R_dhs/COVARIATE_SOURCES.md`](R_dhs/COVARIATE_SOURCES.md).

Proportions are logit transformed, continuous variables remain on their
declared scale, and all included variables are standardised. They enter the
model as one ridge-penalised matrix block; the amount of shrinkage is estimated
from the data.

### Main outcome model {#model-comparison}

For each survey region, the mortality rate is converted to an approximate death
count using the DHS birth exposure. Models use:

- a negative-binomial likelihood with log link;
- `offset(log(exposure))`;
- a ridge penalty on the full eligible covariate block;
- country random intercepts;
- country-specific random linear PfPR slopes.

Five models are compared using maximum likelihood and AIC on the same 921
observations:

| Model | Prevalence and time specification | AIC | ΔAIC |
|---|---|---:|---:|
| Spline without interaction | `s(PfPR) + s(year)` | 7112.12 | 0.00 |
| Spline with `ti` interaction | `s(PfPR) + s(year) + ti(PfPR, year)` | 7113.93 | 1.81 |
| Linear without interaction | `PfPR + s(year)` | 7120.28 | 8.16 |
| Full tensor surface | `te(PfPR, year)` | 7123.78 | 11.66 |
| Linear with interaction | `PfPR + PfPR × year + s(year)` | 7128.32 | 16.20 |

The full tensor surface incorporates the PfPR main effect, year main effect and
their interaction in one term and uses 11.45 effective degrees of freedom. It
is more flexible than the selected additive `s(PfPR) + s(year)` decomposition
but fits less well by AIC.

The selected specification is refitted with REML for post-neonatal mortality,
all-under-five mortality and neonatal mortality. The linear no-interaction
model is also refitted for all three outcomes to give an interpretable
percentage change per 10 PfPR percentage points.

The same five-model comparison is also run independently for neonatal mortality
as a diagnostic. The primary negative-control comparison nevertheless retains
the post-neonatal-selected specification to avoid choosing a different
functional form based on the negative-control outcome.

## Main results

### DHS/MIS association {#association-results}

From the ridge-linear summary models, each 10 percentage-point increase in
PfPR2-10 is associated with:

| Outcome | Change in mortality | 95% CI | p-value |
|---|---:|---:|---:|
| Post-neonatal | +8.75% | +5.77% to +11.83% | <0.001 |
| All under-five | +6.41% | +4.12% to +8.76% | <0.001 |
| Neonatal negative control | +0.98% | −1.16% to +3.17% | 0.373 |

The neonatal result reduces, but does not eliminate, concern that the main
association reflects broad differences in child survival between
higher-prevalence and lower-prevalence settings.

### Neonatal model selection {#neonatal-selection}

The additive spline without a time interaction also has the lowest neonatal
AIC, although the preference over the linear no-interaction model is weak:

| Neonatal specification | AIC | ΔAIC |
|---|---:|---:|
| Spline without interaction | 6008.68 | 0.00 |
| Linear without interaction | 6009.64 | 0.97 |
| Linear with interaction | 6010.92 | 2.24 |
| Full tensor surface | 6012.41 | 3.73 |
| Spline with `ti` interaction | 6016.71 | 8.03 |

In the REML refit, the selected neonatal prevalence smooth is effectively
linear (edf 1.00) and remains consistent with no association (`p=0.374`). The
complete post-neonatal and neonatal comparison is in
`results/dhs_rebuild/outcome_specific_aic.csv`.

### Covariate associations

![Ridge-standardized conditional associations for the 19 covariates in the
best-fitting post-neonatal model.](results/dhs_rebuild/figure5_covariate_forest.png){width=92%}

The forest plot reports the percentage change in post-neonatal mortality for a
one-standard-deviation increase in each transformed covariate, conditional on
the complete model. The coefficients and intervals are ridge-shrunk and should
not be interpreted as independent causal effects.

### Attributable fraction {#attributable-fraction}

The selected time-stable spline estimates the following population-average
fractions of post-neonatal mortality attributable to malaria relative to a 1% PfPR
counterfactual:

| PfPR2-10 | Attributable fraction | Model-parameter 95% interval |
|---:|---:|---:|
| 10% | 14.0% | 7.7%–19.8% |
| 30% | 34.0% | 25.5%–41.6% |
| 50% | 40.1% | 30.5%–48.3% |

Country random effects are set to zero for these population-average
predictions and for the national burden estimates.

![Attributable-fraction curves under the alternative time
structures](results/dhs_rebuild/time_surface_af_comparison.png){width=90%}

## Does the time interaction materially affect predictions? {#time-interaction-results}

Yes, particularly for historical trends and the latest aggregate burden, even
though the additive spline is the AIC-selected specification.

At 30% PfPR, the attributable fraction is:

| Model | 2000 | 2014 | 2024 | 2025 |
|---|---:|---:|---:|---:|
| No time interaction | 34.0% | 34.0% | 34.0% | 34.0% |
| Spline `ti` interaction | 24.7% | 34.1% | 40.2% | 40.7% |
| Full `te` surface | 25.8% | 30.2% | 38.8% | 37.7% |

Applied to the same national prevalence and all-cause mortality series, the
estimated post-neonatal malaria burden changes from 2000 to 2024 as follows:

| Model | 2000 deaths | 2024 deaths | Change |
|---|---:|---:|---:|
| No interaction | 1,007,000 | 477,000 | −52.6% |
| Spline `ti` interaction | 761,000 | 568,000 | −25.3% |
| Full `te` surface | 812,000 | 549,000 | −32.4% |

This is the most important sensitivity for the manuscript's temporal bottom
line. The AIC-selected time-stable relationship implies a decline of about
53%. The alternative time-varying relationships imply smaller declines of
about 25%–32%.

![National burden trajectories under the alternative time
structures](results/dhs_rebuild/time_surface_burden_timeseries.png)

For the latest 2025-surface comparison, the estimates are:

| Model | Estimated deaths | Model-parameter 95% interval |
|---|---:|---:|
| No interaction | 477,000 | 342,000–593,000 |
| Spline `ti` interaction | 577,000 | 409,000–723,000 |
| Full `te` surface | 530,000 | 188,000–796,000 |

Thus the spline interaction raises the latest point estimate by about 100,000
deaths, or 21%, compared with the selected no-interaction model. The full `te`
estimate lies between the two.

## Other sensitivity analyses {#sensitivity-results}

### Alternative samples

All sample sensitivities retain a positive nonlinear prevalence relationship:

- Primary sample: AF at 30% PfPR = 34.0%.
- Include regions below 1% PfPR: 33.7%.
- Restrict PfPR to 5–40%: 34.4%.
- Complete cases only: 34.7%.

These choices do not materially change the cross-sectional attributable
fraction.

### Model structure and covariate blocks

- Removing the country random PfPR slope worsens fit by about 13.0 AIC units.
- Regional covariates alone fit about 17.5 AIC units worse than the full
  regional-plus-national block.
- National covariates alone fit very poorly (about 225 AIC units worse).
- Changing the PfPR spline basis to `k=4` improves fit by 0.9 AIC units; `k=8`
  fits about 18.2 units worse. AF at 30% ranges from 32.4% to 34.4% across these
  basis choices.
- Across structural sensitivities, AF at 30% ranges from 32.4% to 44.8%; the
  upper value occurs in the national-covariate-only model, which fits poorly.

### Alternative likelihood

The log-Gaussian rate sensitivity is less precise than the primary
negative-binomial count model:

- ridge-linear post-neonatal estimate: +5.75% per 10 PfPR points
  (95% CI −1.73% to +13.80%; `p=0.136`);
- selected additive-spline prevalence term: `p=0.136`;
- the no-interaction spline remains the selected likelihood sensitivity.

This weakens certainty about the exact effect size and the time interaction,
but the point estimate remains positive.

## Latest national burden comparison

### Aggregate bottom line {#aggregate-burden}

The closest available comparison is not perfectly contemporaneous:

- **Model:** surface evaluated at 2025 using 2024 MAP PfPR, IGME mortality and
  live births for 43 matched sub-Saharan African countries.
- **IHME:** 2025 under-five malaria deaths.
- **WHO:** WMR 2025 estimates for 2024. WHO reports all-age country deaths, so
  the under-five comparator is an approximate 75% share.

| Source | Deaths |
|---|---:|
| Selected time-stable prevalence/all-cause model | **477,000** |
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

### Countries with the largest differences from IHME {#country-differences}

The table below uses the AIC-selected time-stable spline. “Very different” is
defined in the generated results as at least a two-fold difference and an
absolute difference of at least 1,000 deaths.

| Country | Model | IHME 2025 | Model/IHME | Interpretation |
|---|---:|---:|---:|---|
| DR Congo | 108,800 | 68,300 | 1.59 | Largest absolute excess: +40,500 |
| Nigeria | 176,100 | 150,500 | 1.17 | Model higher by about 25,600 |
| Chad | 12,100 | 5,500 | 2.21 | Model substantially higher |
| Sudan | 3,900 | 1,400 | 2.85 | Model substantially higher |
| Uganda | 13,200 | 28,500 | 0.46 | Model lower by about 15,300 |
| Burundi | 3,800 | 15,700 | 0.24 | Model substantially lower |
| Tanzania | 5,000 | 10,700 | 0.47 | Model substantially lower |
| Ethiopia | 3,900 | 9,300 | 0.42 | Model substantially lower |
| Ghana | 3,500 | 7,300 | 0.48 | Model substantially lower |
| Rwanda | 300 | 2,500 | 0.11 | Model substantially lower |

Against the approximate WHO country-level under-five proxy, the clearest
two-fold discrepancies are:

- higher under the model: DR Congo;
- lower under the model: Tanzania, Ethiopia, Ghana, Sudan, Rwanda, Kenya and
  Madagascar.

DR Congo and Nigeria account for most of the model's upward difference from
IHME, but this is largely offset by lower estimates for the countries listed
above. This is why the aggregate model and IHME totals are close despite large
country-level differences.

### Interpretation and limitations {#interpretation-limitations}

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

\enlargethispage{2\baselineskip}

Our current interpretation is that the data provide fairly robust evidence of
a positive cross-sectional association between malaria prevalence and
post-neonatal mortality after adjustment for national vaccine rollout. The
preferred aggregate estimate is now close to IHME, but it is not a validation
against a common estimand. The selected time-stable model implies a decline of
about 53% since 2000; time-varying sensitivities imply smaller declines of
about 25%–32%.

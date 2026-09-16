# Sahel-only sensitivity of the primary MAP models

Retained 378,680 children, 1,230,499 child–age-band records and 15,999 deaths from 196 survey regions in 25 surveys and 7 countries: BFA, GMB, MLI, MRT, NER, SEN, TCD.

## Geographic definition

Include whole survey regions whose cached boundary centroid is at or north of 12°N and west of 36°E, excluding Ethiopia, Eritrea, Somalia and Djibouti. The western/Horn rule follows the project's existing Sahel convention; the requested 12°N replaces the 11°N cutoff only for this sensitivity. All 1,015 primary survey-region keys matched exactly to survey-specific centroid records; no missing locations were dropped. `region_selection.csv` records inclusion for every primary survey region.

This is an approximate geographic proxy for seasonal transmission, not a classification from measured seasonality. Centroids are not child locations or population-weighted locations. Whole regions can cross the cutoff. Nigeria's broad northern-zone centroids are below 12°N, so no Nigerian regions enter; do not interpret this as evidence of absent seasonal transmission in northern Nigeria. The older 11°N monthly mortality analysis is unchanged.

## Model and comparison

Refit seven separate age-band models with primary MAP prevalence, gamma=2, cr PfPR splines (k=5), cr calendar-year splines (k=6), binomial complementary-log-log likelihood, fREML, discrete bam and the fixed full-band-duration offset. Keep the same confounders, globally computed covariate scaling, fixed median child HIV incidence imputation and unweighted likelihood. Each age/subgroup model has its own time curve, confounder effects and survey/country/region random effects. No outcome, exposure or HIV data are rebuilt.

Knots are subgroup quantiles of distinct predictor values, matching the cr construction convention; dimensions are unchanged and knots are recorded in `knots.csv`. This adapts basis locations and re-estimates smoothing parameters for the subgroup. The reference is the saved primary full-sample MAP gamma=2 fit. Prepared-data hashes are checked against its manifest.

Figure caption: Adjusted mortality hazard ratios across MAP PfPR₂–₁₀ in the full primary sample (blue) and Sahel subset defined by regional centroids ≥12°N (orange), by completed-month age band. Both curves are relative to 20% PfPR; ribbons are pointwise 95% model-based intervals conditional on fitted smoothing parameters, fixed exposure and the fixed HIV imputation. Each curve is drawn over its own 2.5th–97.5th percentile exposure range. The hazard-ratio figure uses a logarithmic y-axis; a log-hazard-ratio version is also provided. Covariates and fitting choices match the primary models. The subgroup overlaps the full sample, so these curves are dependent; this comparison does not provide an independent-sample interaction test or establish a seasonal causal effect.

[Hazard-ratio overlay](sahel_vs_full_pfpr_hazard_ratios.png) · [Log-hazard-ratio overlay](sahel_vs_full_pfpr_log_hazard_ratios.png) · [Region selection](region_selection.png)

## Supported 40% to 20% contrasts

Hazard ratios below one imply lower mortality at 20% than at 40%. Intervals are conditional pointwise 95% intervals; compare descriptively, without assuming the fits are independent.

| Age (months) | Full sample HR (95% interval) | Sahel HR (95% interval) | Sahel central PfPR range (%) |
|---|---:|---:|---:|
| <1 | 0.97 (0.94–0.99) | 0.95 (0.89–1.01) | 1.8–67.2 |
| 1-5 | 0.93 (0.89–0.96) | 0.95 (0.87–1.04) | 1.8–67.1 |
| 6-11 | 0.84 (0.80–0.88) | 0.77 (0.69–0.85) | 1.8–67.1 |
| 12-23 | 0.87 (0.81–0.94) | 0.68 (0.61–0.74) | 1.8–65.0 |
| 24-35 | 0.84 (0.78–0.91) | 0.74 (0.66–0.83) | 1.8–64.8 |
| 36-47 | 0.82 (0.76–0.90) | 0.73 (0.63–0.85) | 1.8–64.8 |
| 48-59 | 0.96 (0.86–1.06) | 0.74 (0.65–0.85) | 1.8–64.8 |

The Sahel curves show a stronger positive association above 20% PfPR from six months onward, particularly at 12–23 months; the curves below six months are closer. These are descriptive differences in fitted shapes, not a formal test of effect modification by seasonality.

For 12–23 months, the 40%→20% hazard ratio is 0.68 (0.61–0.74) in the Sahel subset versus 0.87 (0.81–0.94) in the full analysis.

The CSV also records 20% to zero contrasts. Zero is outside observed support in both samples; those contrasts extrapolate and are not shown as supported curve segments.

## Validation and outputs

All seven subgroup models converged, with full rank, finite coefficients/covariances and positive smoothing-Hessian eigenvalues. Every fitted model column was checked against the selected input. Curves were checked against direct link predictions and anchored exactly at HR=1 at 20% PfPR. Input hashes, knots, fit manifests, aggregate sample counts, smoothing EDF and numerical diagnostics are saved alongside this report. Intervals do not account fully for DHS sampling design, HIV/exposure uncertainty or cross-model covariance.

Pure fitting time across seven models: 17.5 seconds (excluding reading/writing).

Reproduce from the project root:

```sh
Rscript R_cbh/sensitivity/sahel/01_fit.R
Rscript R_cbh/sensitivity/sahel/02_report.R
```

Only aggregate tables/figures are published. Individual-level prepared input and fitted model objects remain in the ignored data directory. Primary models, burden estimates and manuscript files are unchanged.

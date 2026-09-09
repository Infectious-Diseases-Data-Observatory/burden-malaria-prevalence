# Exploratory age-band fit

The [national burden stage](../burden/README.md) applies this model's PfPR-to-zero contrasts to country-level 2024 IHME rates and deaths, including the user-selected equal-rate/equal-person-time allocation within ages 2-4 years.

**Active specification, 8 September 2026:** use log child HIV incidence at band entry, with missing/censored rates supplied by the [incidence imputation model](../hiv/README.md). Calendar year now enters through one shared spline across age bands. The earlier incidence model with age-specific time splines remains under `age_band_hiv_incidence_v2/`; the historical prevalence fit also remains available separately.

Run from the project root after completing the dataset build:

```sh
Rscript R_cbh/tests/test_model.R
Rscript R_cbh/tests/test_hiv_incidence.R
Rscript R_cbh/hiv/01_fit_incidence.R --cv
Rscript R_cbh/hiv/02_report_incidence.R
Rscript R_cbh/analysis/01_fit_complete_case.R
Rscript R_cbh/analysis/02_plot_trial.R
Rscript R_cbh/analysis/04_propagate_incidence.R
Rscript R_cbh/analysis/05_report_incidence.R
```

`--prepare-only` writes the exact complete-case input and sample reports without fitting. `--force` replaces matching caches by rebuilding/refitting. No packages are installed and no network requests are made. The implementation uses the installed `mgcv`, `data.table` and `digest` packages.

## Specification

One joint binomial model with complementary log–log link and `offset(log(band_years))`. The seven unordered age bands have separate intercepts, cubic regression splines for PfPR (k=5), and confounder coefficients. A single cubic regression spline for calendar year at entry (k=6) is shared across all age bands. Survey, country-by-age and survey-specific region have random intercepts. Fitting uses discrete `mgcv::bam` with fREML and requests two threads; the currently installed mgcv lacks OpenMP and reports a fallback to one thread.

The trial uses **all complete cases after the child HIV incidence join** for sex, multiple birth, birth order, maternal age at birth, maternal education, wealth quintile, urban residence, log child HIV incidence, log GDP per capita, log health expenditure per capita and political stability. Child HIV incidence is the sole imputed covariate. Vaccine coverage is excluded as requested. It also requires observed outcome, PfPR, time, band width and model grouping variables. It does not require optional unused fields to be observed.

Wealth is an unordered five-level factor, with poorest as the reference; female is the reference sex. Other continuous confounders enter linearly with age-specific coefficients and are standardized using the complete-case sample's means and standard deviations. This is a numerical transformation, not imputation. Scaling parameters are saved. PfPR has age-specific nonlinear effects and calendar year has one shared nonlinear effect; richer confounder functions remain an open modelling decision.

The likelihood is **unweighted**, matching the initial displayed model specification. Raw DHS weights and weights normalized after complete-case selection remain in the saved input for subsequent work. The reported covariance and intervals are approximate conditional model-based quantities, without survey-design, spatial-confounding or exposure-estimation uncertainty. See [open issues](../../docs/OPEN_ANALYSIS_ISSUES.md).

## Outputs

Active model inputs and fitted objects stay under the ignored directory `data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/`. The historical prevalence fit remains under `age_band_complete_case_v1/` and requires `--legacy-prevalence` to select it:

- `complete_case_dataset.rds`: list containing `data`, `scaling`, `selection`, `missing`, `skipped`, the source manifest, specification and cache signature. `data` includes the selected model inputs and pseudonymous sampling/child identifiers. It is research microdata.
- `fit.rds`: list containing `fit` (the `bam` object), elapsed fitting time, captured warnings, specification, signature and R session information. It may contain row-level model data and must remain private.

The cache signature includes the complete dataset manifest, model specification, fitting/preparation code, R version and mgcv version. A failed convergence check preserves the fitted object for diagnosis but does not publish new contrasts. Do not run competing fits into this directory concurrently.

Aggregate reports are written to `results/cbh/age_band_hiv_incidence_shared_time_v3/` (or the legacy directory when explicitly selected):

| File | Content |
|---|---|
| `selection.csv`, `sample_by_country.csv`, `sample_by_age_band.csv` | Denominators and deaths before/after complete-case selection |
| `missing.csv` | Missing counts for each requested covariate among PfPR-available records; counts overlap and must not be summed as unique exclusions |
| `skipped.csv` | Surveys unavailable at the dataset stage |
| `scaling.csv`, `specification.txt` | Model definition and exact covariate transformations |
| `convergence.txt`, `model_summary.txt` | Numerical diagnostics, warnings and model summary |
| `basis_checks.csv` | Seeded basis-dimension checks on up to 10,000 observations; factor smooth checks may be unavailable |
| `pfpr_40_to_20_contrasts.csv` | Age-specific hazard ratios for 40% to 20% PfPR with approximate pointwise 95% model-based intervals and central exposure-support checks |
| `pfpr_40_to_20_by_age.png` | Forest plot of the saved contrasts, generated separately by `02_plot_trial.R` |
| `pfpr_splines_by_age.png`, `pfpr_splines_full_range.png` | Fitted PfPR curves over the central 95% and full observed exposure ranges, on common axes |
| `pfpr_spline_curves.csv` | Aggregate grid predictions, log hazard ratios relative to 20% PfPR, standard errors, pointwise intervals, exposure-support bounds and fitted-model signature |

`02_plot_trial.R` reads the private saved fit locally without refitting and exports only aggregate curve predictions. Spline figures use the installed `ggplot2` and `ragg` packages. For each age band, the curve is `f_a(P) - f_a(20)`; its exponential is the mortality hazard ratio relative to 20% PfPR. Prediction-matrix differences hold all other predictors fixed, and the full fitted coefficient covariance accounts for correlation with the reference prediction. Therefore the contrast and its standard error are exactly zero at 20%. Intervals are pointwise and conditional on smoothing parameters. Exposure percentiles are unweighted, pooled over complete-case records within each band; they do not establish joint confounder overlap. The plotting script verifies direct link predictions, the zero reference contrast and agreement with the inverse of the saved 40-to-20 contrasts.

For band a, the contrast is `exp(f_a(20) - f_a(40))`. The confounder profile and random effects are held fixed, so they cancel in this ratio. It is not a death-probability ratio or a total intervention effect from birth. The code checks the prediction-matrix calculation against differences of predicted links. No sampling of fitting rows is performed; the limited subsample is used only for the basis diagnostic.

Synthetic checks cover complete-input enforcement, factor/scaling setup, seven separate PfPR curves, one calendar-year curve with identical time contrasts across age bands, finite fit outputs, unchanged row counts, covariance-based contrast intervals and a unity ratio for identical exposures.

The [first fit report](../../results/cbh/age_band_complete_case_v1/REPORT.md) records its estimates, convergence and diagnostic flags. The initial fit's PfPR basis-size and smoothing-optimization checks require follow-up before final inference.

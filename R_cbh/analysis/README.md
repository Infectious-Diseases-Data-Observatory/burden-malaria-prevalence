# Exploratory complete-case age-band fit

Run from the project root after completing the dataset build:

```sh
Rscript R_cbh/tests/test_model.R
Rscript R_cbh/analysis/01_fit_complete_case.R
Rscript R_cbh/analysis/02_plot_trial.R
```

`--prepare-only` writes the exact complete-case input and sample reports without fitting. `--force` replaces matching caches by rebuilding/refitting. No packages are installed and no network requests are made. The implementation uses the installed `mgcv`, `data.table` and `digest` packages.

## Specification

One joint binomial model with complementary log–log link and `offset(log(band_years))`. The seven unordered age bands have separate intercepts, cubic regression splines for PfPR (k=5) and entry time (k=6), and confounder coefficients. Survey, country-by-age and survey-specific region have random intercepts. Fitting uses discrete `mgcv::bam` with fREML and requests two threads; the currently installed mgcv lacks OpenMP and reports a fallback to one thread.

The trial uses **all complete cases** for sex, multiple birth, birth order, maternal age at birth, maternal education, wealth quintile, urban residence, log HIV prevalence, log GDP per capita, log health expenditure per capita and political stability. Vaccine coverage is excluded as requested. It also requires observed outcome, PfPR, time, band width and model grouping variables. It does not require optional unused fields to be observed. No outcome/exposure/covariate values are imputed.

Wealth is an unordered five-level factor, with poorest as the reference; female is the reference sex. Other continuous confounders enter linearly with age-specific coefficients and are standardized using the complete-case sample's means and standard deviations. This is a numerical transformation, not imputation. Scaling parameters are saved. Age-specific nonlinear PfPR and time effects are retained as specified; richer confounder functions remain an open modelling decision.

The likelihood is **unweighted**, matching the initial displayed model specification. Raw DHS weights and weights normalized after complete-case selection remain in the saved input for subsequent work. The reported covariance and intervals are approximate conditional model-based quantities, without survey-design, spatial-confounding or exposure-estimation uncertainty. See [open issues](../../docs/OPEN_ANALYSIS_ISSUES.md).

## Outputs

Model inputs and fitted objects stay under the ignored directory `data/derived_cbh/models/age_band_complete_case_v1/`:

- `complete_case_dataset.rds`: list containing `data`, `scaling`, `selection`, `missing`, `skipped`, the source manifest, specification and cache signature. `data` includes the selected model inputs and pseudonymous sampling/child identifiers. It is research microdata.
- `fit.rds`: list containing `fit` (the `bam` object), elapsed fitting time, captured warnings, specification, signature and R session information. It may contain row-level model data and must remain private.

The cache signature includes the complete dataset manifest, model specification, fitting/preparation code, R version and mgcv version. A failed convergence check preserves the fitted object for diagnosis but does not publish new contrasts. Do not run competing fits into this directory concurrently.

Aggregate reports are written to `results/cbh/age_band_complete_case_v1/`:

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

For band a, the contrast is `exp(f_a(20) - f_a(40))`. The confounder profile and random effects are held fixed, so they cancel in this ratio. It is not a death-probability ratio or a total intervention effect from birth. The code checks the prediction-matrix calculation against differences of predicted links. No sampling of fitting rows is performed; the limited subsample is used only for the basis diagnostic.

Synthetic checks cover complete-input enforcement, factor/scaling setup, seven separate PfPR/time curves, finite fit outputs, unchanged row counts, covariance-based contrast intervals and a unity ratio for identical exposures.

The [first fit report](../../results/cbh/age_band_complete_case_v1/REPORT.md) records its estimates, convergence and diagnostic flags. The initial fit's PfPR basis-size and smoothing-optimization checks require follow-up before final inference.

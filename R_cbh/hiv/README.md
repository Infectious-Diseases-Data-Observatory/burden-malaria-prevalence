# Child HIV incidence imputation

This replaces the active mortality model's HIV prevalence adjustment with **child HIV incidence**, as requested on 8 September 2026. Units are new infections per 1,000 uninfected population, both sexes, ages 0–14. The adolescent predictor is the corresponding reported rate for ages 15–19. No population denominator is fetched or calculated.

## Run

From the project root, using the existing installed R packages:

```sh
Rscript R_cbh/tests/test_hiv_incidence.R
Rscript R_cbh/hiv/01_fit_incidence.R --cv
Rscript R_cbh/hiv/02_report_incidence.R
Rscript R_cbh/analysis/01_fit_complete_case.R
Rscript R_cbh/analysis/02_plot_trial.R
Rscript R_cbh/analysis/04_propagate_incidence.R
Rscript R_cbh/analysis/05_report_incidence.R
```

The first fit uses `rstan`, `readxl`, `splines` and existing pipeline utilities. It installs nothing and makes no network requests. `--cv` runs five folds withholding entire countries' child series; their adolescent series remain available. Fits and validation folds are cached with workbook, code and RStan-version signatures. `--force` refits. Failed numerical diagnostics retain the fit for diagnosis and stop publication of that run's downstream outputs.

## Statistical model

Let a[c,t] be log adolescent incidence and y[c,t] log child incidence. The child model is

`y[c,t] = alpha_c + beta_between * mean(a[c,]) + beta_within * (a[c,t] - mean(a[c,])) + s_c(t) + region_effect_c + country_effect_c + e_c[c,t]`.

The adolescent model is

`a[c,t] = alpha_a + s_a(t) + region_effect_a + country_effect_a + e_a[c,t]`.

Separate natural cubic calendar-year spline bases have three interior knots (2006, 2012, 2018), boundaries 2000/2024, centered columns and regularizing Normal(0,2) coefficient priors. Both models have partially pooled UNICEF-region and country intercepts. Within-country residuals are stationary Gaussian AR(1), with separate positive autocorrelations and stationary residual SDs for the two ages. Transitions account for actual gaps in years. Country means of adolescent log incidence are computed inside the model, including latent censored values; between-country and within-country adolescent relationships need not be equal.

Prior scales are explicit in `incidence.stan`: child intercept Normal(0,2), adolescent intercept Normal(-2,3), between/within slopes Normal(1,1), positive region/country SDs half-Normal(0,1) for child and half-Normal(0,2) for adolescent, stationary residual SDs half-Normal(0,1), and AR parameters Beta(3,1). The slope priors permit negative associations. These are regularizing modelling choices, not established biological constants.

Numeric published rates are conditioned on. Values reported as `<0.01` enter as left-censored log rates. A sequential truncated-normal transformation uses latent uniforms and conditional CDF weights; it does not assign 0.005 or 0.01 as observed values. A stable log-probability inverse-normal calculation avoids underflow in very small censoring probabilities. Country-year entries with no usable adolescent series are not imputed from region alone. Regional aggregates and separate-sex rows are excluded.

This is a model of already estimated national series, not a binomial model of sampled infection counts. Published lower/upper bounds are retained in the source audit but **not incorporated as a measurement-error likelihood**: their joint uncertainty across age/year is unavailable. Consequently posterior intervals describe uncertainty under this incidence reconstruction model, conditional on numeric source estimates. They do not represent all UNAIDS estimation uncertainty.

## Prediction and downstream use

Reported numeric child rates are preserved exactly. Censored child rates receive posterior draws below their threshold. Missing child rates receive full predictive draws, including country effects and correlated residual trajectories. An entirely unobserved country's child effect remains a population-distribution draw; it is not fixed to zero. Conditional Gaussian calculations preserve temporal dependence when predicting missing years within an observed child series. Two hundred coherent draws are saved, shared across all DHS records in the same country and entry year.

Nigeria has a numeric adolescent series; Comoros has a censored adolescent series. Both receive child imputations with different amounts of information. Liberia and São Tomé and Príncipe have no usable both-sex adolescent incidence series in this workbook and remain unavailable. Liberia was included in the historical prevalence fit, so changing the HIV definition changes the sample as well as the covariate. The transport assumption is that, conditional on adolescent history, year and region, countries without child estimates follow the same child relationship as countries that report them. Validation probes this assumption but cannot establish it for a country with no child observations.

`analysis/model.R` joins the panel by country and **band-entry calendar year**, then constructs `log_hiv_incidence`. It is an age-specific covariate in the existing seven-band cloglog model. It replaces `log_hiv_prev`; it is not added alongside it, and adolescent incidence does not enter the mortality formula. Other adjustment variables still require complete cases and vaccines remain excluded. This overlay reads the existing validated child-band shards without rebuilding raw DHS histories. The exact incidence-adjusted analysis dataset is saved under `data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/complete_case_dataset.rds`.

The default mortality fit uses posterior-median child incidence values and is a conditional plug-in comparison. `04_propagate_incidence.R` fits ten coherent incidence draws and pools PfPR log-hazard contrasts and curves using within-fit plus between-imputation variance and finite-imputation t intervals. Mortality smoothing parameters are fixed to the median-incidence fit during these refits; HIV scaling is also held fixed. The imputation model receives no DHS mortality outcomes. This is external-covariate uncertainty propagation, not a joint outcome-compatible missing-data model. Ten draws provide exploratory precision; saved Monte Carlo standard errors indicate whether more are needed. Survey-design, source-estimate and mortality smoothing-parameter uncertainty remain outside these intervals.

## Outputs

- `data/derived_cbh/hiv_incidence/`: cached Stan fits, fitted child-rate panel and posterior trajectories.
- `results/cbh/hiv_incidence/`: published country-level source audit, fitted panel, diagnostic/parameter tables and country-held-out validation outputs.
- `results/cbh/age_band_hiv_incidence_shared_time_v3/`: mortality fit using median child incidence, sample accounting, contrasts and curves.
- `results/cbh/age_band_hiv_incidence_shared_time_v3/multiple_imputation/`: conditional-fit diagnostics and pooled contrasts/curves.

Legacy prevalence results remain under `age_band_complete_case_v1`. The fitting and plotting scripts accept `--legacy-prevalence` to explicitly select that specification; it is no longer the default.

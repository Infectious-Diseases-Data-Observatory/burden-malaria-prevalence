# National malaria-attributable mortality in 2024

Run from the repository root after completing the shared-calendar-year mortality model and its ten incidence-imputation refits:

```sh
Rscript R_cbh/burden/01_country_attributable_2024.R
Rscript R_cbh/burden/02_report_country_attributable_2024.R
```

The calculation uses the IHME all-cause export dated **2026-09-09 10-58-22**, restricted to 2024 and both sexes. For country c and model age band g:

```
HR_zero_vs_current = exp(f_g(0) - f_g(PfPR_c))
counterfactual_rate = IHME_rate * HR_zero_vs_current
attributable_rate = IHME_rate * (1 - HR_zero_vs_current)
attributable_deaths = IHME_deaths * (1 - HR_zero_vs_current)
```

Rates are per 100,000 person-years. Annual counts hold the population exposure fixed. The same contrast applies to rates and counts, not directly to conditional death probabilities. The current additive model's calendar-year, confounder and random-intercept terms cancel in the contrast. The common fitted age-specific PfPR curves are transported to each country, including countries outside the mortality model's fitting sample.

The new export directly supplies infancy and ages 12-23 months, while retaining a combined 2-4-year group. **The user selected the same baseline IHME 2-4 rate for ages 24-35, 36-47 and 48-59 months, with deaths/person-time divided equally among them.** Each receives its own age-specific PfPR contrast. This allocation preserves the original 2-4 total exactly.

Early and late neonatal counts are added; their rates are combined using implied person-time (count / rate). Both receive the model's <1-completed-month PfPR effect, with the 0-27-day versus completed-month boundary approximation recorded. Under 1 and Under 5 totals are used for consistency checks only.

Population-weighted national PfPR is extracted locally from `data/pfpr_2to10_africa_2024.tif` and `data/africa_admin0.rds`. Weights use GPW 2020 population density multiplied by raster-cell area and exact polygon overlap. The previous national series omitted cell area from the density weights; both estimates are retained for comparison. Missing MAP is not replaced by zero. Population geography is fixed to 2020 and includes all ages; it is not a 2024 age-specific population surface. This national-mean exposure calculation differs from applying nonlinear effects subnationally and aggregating.

No model is refitted. Compact PfPR coefficients, covariance blocks and smooth objects are cached privately under the active model directory. File hashes invalidate the cache when saved fits change. Every compact-basis prediction is checked against the earlier full-matrix 40%-to-20% contrast for the same fitted model. The ten log-HR estimates are pooled with within-fit covariance and between-imputation variance, using finite-imputation t intervals. Point HRs exponentiate the pooled mean log HR. Bounds condition on source IHME/MAP point estimates and fixed mortality smoothing parameters. They omit source-estimate, survey-design, residual-clustering and allocation uncertainty. Country totals have point estimates only; marginal band interval endpoints are never summed into a total interval.

Signed negative attributable estimates and zero-exposure support flags are retained. A negative value means the fitted association predicts increased mortality under the zero-PfPR contrast; it is not evidence establishing that malaria is protective. The causal interpretation, extrapolation to zero, between-country transport and existing smoothing-stability concern remain unresolved.

Outputs are under `results/cbh/age_band_hiv_incidence_shared_time_v3/country_burden_2024/`:

- `country_age_attributable_2024.csv`: all country-by-band baseline rates/counts, counterfactual estimates, signed attributable estimates, conditional intervals, allocation and exposure-support flags. Missing-country rows remain explicit.
- `country_totals_2024.csv`: national under-five annual death point estimates.
- `national_pfpr_2024.csv`: national prevalence, legacy comparison and raster coverage.
- `ihme_source_2024.csv`, `ihme_disjoint_age_inputs.csv`: public source estimates and nonoverlapping inputs with source uncertainty bounds.
- `individual_log_hazard_contrasts.csv`, `model_pfpr_support.csv`, `provenance.csv`: uncertainty and source auditing.
- `drc_period_life_table_2024.csv`, `drc_malaria_contribution_2024.png`, `REPORT.md`: DRC example and readable findings.

The DRC figure follows Burstein et al. (2018), Figure 1, by showing conditional age-band mortality and survival, with a third panel for annual attributable mortality rates. The period life table uses piecewise constant hazards and a 28-day neonatal boundary, integrating early and late neonatal source rates separately. It describes a synthetic cohort exposed to 2024 mortality conditions; its survival difference is not multiplied by births to calculate annual attributable deaths.

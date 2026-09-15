# National malaria-attributable mortality by year

> The commands below default to historical gamma=1 estimates and also prepare inputs reused by gamma=2. Current primary country results are `map_full` in `results/cbh/map_snow_gamma2_v1/burden/`; do not treat this folder’s defaults as current primary estimates. See [the audit](../../docs/CODE_AUDIT.md).

Scripts accept `--year=YYYY` (2000-2024; default 2024) and `--model=separate|joint` (default **separate**, the primary specification selected on 9 September 2026). Each model/year has a separate output directory. Historical joint-model outputs are preserved. The same fitted PfPR curves and fixed GPW 2020 population geography are used across years.

The primary update uses the **seven existing separate age-band fits and one fixed posterior-median HIV-incidence imputation**, from `data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/sensitivity_single_imputation/age_1.rds` through `age_7.rds`. Each fit has its own calendar-year curve, confounder coefficients and random-effect variances. No mortality model is refitted by the burden scripts. The historical `--model=joint` option uses the joint shared-calendar-year model and ten HIV-incidence uncertainty refits.

Run from the repository root after completing the mortality fits:

```sh
Rscript R_cbh/burden/01_country_attributable.R
Rscript R_cbh/burden/02_report_country_attributable.R
Rscript R_cbh/burden/03_compare_ihme.R

# All requested years, using the separate-age primary by default
for year in 2005 2015 2024; do
  Rscript R_cbh/burden/01_country_attributable.R --year=$year
  Rscript R_cbh/burden/02_report_country_attributable.R --year=$year
  Rscript R_cbh/burden/03_compare_ihme.R --year=$year
done

# Explicit historical model (writes only its separate historical model directory)
Rscript R_cbh/burden/01_country_attributable.R --year=2024 --model=joint
Rscript R_cbh/burden/02_report_country_attributable.R --year=2024 --model=joint
Rscript R_cbh/burden/03_compare_ihme.R --year=2024 --model=joint

# Combined country/year table, changes from previous estimates and summary figure
Rscript R_cbh/burden/04_summarize_primary_update.R
```

The calculation uses the IHME all-cause export dated **2026-09-09 10-58-22**, restricted to the requested year and both sexes. For country c and model age band g:

```
HR_zero_vs_current = exp(f_g(0) - f_g(PfPR_c))
counterfactual_rate = IHME_rate * HR_zero_vs_current
attributable_rate = IHME_rate * (1 - HR_zero_vs_current)
attributable_deaths = IHME_deaths * (1 - HR_zero_vs_current)
```

Rates are per 100,000 person-years. Annual counts hold the population exposure fixed. The same contrast applies to rates and counts, not directly to conditional death probabilities. The current additive model's calendar-year, confounder and random-intercept terms cancel in the contrast. The common fitted age-specific PfPR curves are transported to each country, including countries outside the mortality model's fitting sample.

The new export directly supplies infancy and ages 12-23 months, while retaining a combined 2-4-year group. **The user selected the same baseline IHME 2-4 rate for ages 24-35, 36-47 and 48-59 months, with deaths/person-time divided equally among them.** Each receives its own age-specific PfPR contrast. This allocation preserves the original 2-4 total exactly.

Early and late neonatal counts are added; their rates are combined using implied person-time (count / rate). Both receive the model's <1-completed-month PfPR effect, with the 0-27-day versus completed-month boundary approximation recorded. Under 1 and Under 5 totals are used for consistency checks only.

Population-weighted national PfPR is extracted locally from `data/map_annual/pfpr2_10_YYYY.tif` and `data/africa_admin0.rds`. For 2024, the existing `data/pfpr_2to10_africa_2024.tif` source is retained for exact compatibility. Weights use GPW 2020 population density multiplied by raster-cell area and exact polygon overlap. The previous national series omitted cell area from the density weights; both estimates are retained for comparison. Missing MAP is not replaced by zero. Population geography is fixed to 2020 and includes all ages; it is not a 2024 age-specific population surface. This national-mean exposure calculation differs from applying nonlinear effects subnationally and aggregating.

Compact PfPR coefficients, covariance blocks and smooth objects are cached privately under the selected model directory. File hashes invalidate the cache when saved fits change. For separate-age models, compact log contrasts and variances are checked against full fitted prediction matrices at zero and nonzero exposures, and the 40%-to-20% HR is checked against the saved sensitivity result. Age-band intervals use within-fit coefficient covariance and normal 95% limits on log HR, **conditional on the fixed HIV imputation**. Separate fitting does not imply independent sampling errors across age bands; cross-age covariance is not assumed zero or estimated by this stage. Country totals remain point estimates only, and marginal band interval endpoints are never summed into a total interval.

For the historical joint model only, ten log-HR estimates are pooled with within-fit covariance and between-imputation variance, using finite-imputation t intervals. Point HRs exponentiate the pooled mean log HR. Both methods condition on source IHME/MAP point estimates and fitted mortality smoothing parameters, and omit source-estimate, survey-design, residual-clustering and allocation uncertainty. Comparisons with previous published joint-model totals therefore change both model structure and HIV-imputation treatment; they do not isolate the effect of model structure alone.

Signed negative attributable estimates and zero-exposure support flags are retained. A negative value means the fitted association predicts increased mortality under the zero-PfPR contrast; it is not evidence establishing that malaria is protective. The causal interpretation, extrapolation to zero, between-country transport and existing smoothing-stability concern remain unresolved.

Primary outputs are under `results/cbh/age_band_separate_v1/country_burden_YYYY/`. The historical joint outputs remain under `results/cbh/age_band_hiv_incidence_shared_time_v3/country_burden_YYYY/`. The example filenames below use 2024; each run substitutes its requested year:

- `country_age_attributable_2024.csv`: all country-by-band baseline rates/counts, counterfactual estimates, signed attributable estimates, conditional intervals, allocation and exposure-support flags. Missing-country rows remain explicit.
- `country_totals_2024.csv`: national under-five annual death point estimates.
- `national_pfpr_2024.csv`: national prevalence, legacy comparison and raster coverage.
- `ihme_source_2024.csv`, `ihme_disjoint_age_inputs.csv`: public source estimates and nonoverlapping inputs with source uncertainty bounds.
- `individual_log_hazard_contrasts.csv`, `model_pfpr_support.csv`, `provenance.csv`: uncertainty and source auditing.
- `model_specification.txt`: model identity and HIV-imputation/uncertainty treatment. Model identifiers are also retained in the country and age-band CSVs.
- `drc_period_life_table_2024.csv`, `drc_malaria_contribution_2024.png`, `REPORT.md`: DRC example and readable findings.
- `model_vs_ihme_malaria_2024.csv`, `model_vs_ihme_malaria_2024.png`, `model_vs_ihme_malaria_by_country_2024.png`, `IHME_COMPARISON.md`, `comparison_provenance.csv`: 2024 under-five comparison against the existing IHME malaria death export, including a scatter plot and a chart labelling every matched country. These compare point estimates of all-cause mortality reduction with cause-specific malaria counts. Joint model country-total intervals are not available. Source release/version alignment remains to be verified.

The DRC figure follows Burstein et al. (2018), Figure 1, by showing conditional age-band mortality and survival, with a third panel for annual attributable mortality rates. The period life table uses piecewise constant hazards and a 28-day neonatal boundary, integrating early and late neonatal source rates separately. It describes a synthetic cohort exposed to 2024 mortality conditions; its survival difference is not multiplied by births to calculate annual attributable deaths.

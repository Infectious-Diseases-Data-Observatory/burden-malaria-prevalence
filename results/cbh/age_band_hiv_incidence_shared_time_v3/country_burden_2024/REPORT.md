# National malaria-attributable mortality, 2024

Estimates are available for **42 of 45 countries** in the new IHME export. Countries without usable MAP coverage remain missing: Cape Verde, Lesotho, Sao Tome and Principe.

**Exposure coverage:** 11 estimated countries have MAP covering less than 95% of population weight within the raster footprint. Their estimates apply the covered-area mean to the whole country and should be treated as provisional. Swaziland/Eswatini has only 2.6% coverage; its value is especially poorly representative. Missing pixels are not assumed malaria-free. DR Congo has 96.8% coverage. The 95% threshold is a reporting flag, not an exclusion rule.

Each country's 2024 national PfPR2-10 is evaluated on each fitted age-specific spline. The zero-PfPR hazard ratio is exp[f_g(0)-f_g(P_country)]. Counterfactual mortality rate = IHME rate x HR; attributable rate = IHME rate x (1-HR). Death counts use the same fraction and fixed annual person-time. No adjustment is made to other covariates, calendar year, or country random effects in this contrast: those terms cancel in the current additive model.

The seven PfPR effects come from the latest shared-calendar-year mortality model and its ten HIV-incidence uncertainty refits. Contrasts are pooled on the log-HR scale using within-fit covariance plus between-imputation variance and finite-imputation t intervals. Reported point HRs exponentiate pooled mean log HRs. Bounds condition on IHME and MAP point estimates, fitted smoothing parameters, and the age allocation. They do not include source-estimate, survey-design or residual-clustering uncertainty. National total death counts are point estimates; marginal age-band interval endpoints are not summed into total intervals.

## Source handling and assumptions

- The selected IHME export is the file dated 2026-09-09 10-58-22. Only All causes / Deaths / Both sexes / 2024 is used. Six disjoint source age groups are checked against Under 1 and Under 5 totals; those aggregates are never counted twice. Original rate/count lower and upper bounds are retained in the source output.
- Early and late neonatal deaths are added. Their combined annual rate is total deaths divided by the sum of their implied person-years, not the sum or simple mean of rates. The IHME 0-27-day group uses the model's <1-completed-month PfPR effect, an explicit boundary approximation.
- As authorized, ages 24-35, 36-47 and 48-59 months each use the IHME 2-4-year mortality rate and one-third of its deaths/person-time. Each then receives its own fitted PfPR contrast. These three baseline age estimates are assumed, not separately observed in IHME.
- National exposure is calculated from local 2024 MAP rasters and country polygons, using GPW 2020 density x grid-cell area x polygon overlap as weights. The earlier extraction used density alone. Both values are saved for comparison. This uses fixed all-age 2020 population geography, not a 2024 age-specific population surface. Coverage is conditional on the available raster footprint; missing MAP is never set to zero.
- The national-mean PfPR scenario is the requested approximation; the nonlinear response evaluated at a country mean does not equal a subnational burden aggregation. The fitted PfPR curves are transported to countries outside the mortality fitting sample.
- Negative attributable effects are retained. They mean the fitted association predicts higher mortality at zero PfPR, not protective malaria established by evidence. Zero-PfPR support flags and central-range flags are exported; the model's extrapolation, causal transport and smoothing-stability limitations remain relevant.

## DR Congo

National PfPR: **36.06%**. IHME under-five deaths: **207,594**. Signed model-attributable deaths: **45,859** (**22.1%**), holding the annual population exposure fixed.

All rates below are per 100,000 person-years. Age-band intervals reflect model/HIV-imputation uncertainty only.

| Age | IHME rate | Zero-PfPR rate | Attributable rate (95% interval) | Attributable annual deaths |
|---|---:|---:|---:|---:|
| 0-27 days | 24403.8 | 25948.8 | -1545.0 (-2985.9 to -179.9) | -3614 |
| 1-5 months | 2995.2 | 2875.0 | 120.3 (-75.9 to 303.9) | 1515 |
| 6-11 months | 2463.8 | 1750.0 | 713.7 (502.9 to 901.9) | 10642 |
| 12-23 months | 844.8 | 454.3 | 390.5 (328.7 to 445.0) | 11423 |
| 24-35 months | 598.8 | 246.9 | 351.8 (314.1 to 384.5) | 10053 |
| 36-47 months | 598.8 | 293.9 | 304.9 (252.3 to 349.5) | 8712 |
| 48-59 months | 598.8 | 349.3 | 249.5 (173.4 to 311.9) | 7128 |

![DRC malaria contribution](drc_malaria_contribution_2024.png)

The survival and probability panels describe a synthetic cohort exposed to the 2024 mortality schedule. They use q=1-exp(-H), with IHME neonatal hazards integrated separately over 7 and 21 days, followed by ages 28 days to 6 months, 6-12 months, and annual age intervals. This explicitly handles the source neonatal boundary. The survival product is checked numerically. Neither these probabilities nor their differences are used as annual death-count fractions.

## Outputs

- [All country/age estimates](country_age_attributable_2024.csv)
- [Country totals](country_totals_2024.csv)
- [National prevalence and coverage](national_pfpr_2024.csv)
- [IHME original 2024 rows](ihme_source_2024.csv) and [disjoint age inputs](ihme_disjoint_age_inputs.csv)
- [Individual-fit contrasts](individual_log_hazard_contrasts.csv), [PfPR support](model_pfpr_support.csv)
- [DRC period life table](drc_period_life_table_2024.csv), [provenance](provenance.csv)

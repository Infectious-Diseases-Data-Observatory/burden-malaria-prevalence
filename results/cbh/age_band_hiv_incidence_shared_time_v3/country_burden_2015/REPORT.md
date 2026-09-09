# National malaria-attributable mortality, 2015

Estimates are available for **42 of 45 countries** in the new IHME export. Countries without usable MAP coverage remain missing: Cape Verde, Lesotho, Sao Tome and Principe.

**Exposure coverage:** 11 estimated countries have MAP covering less than 95% of population weight within the raster footprint. Their estimates apply the covered-area mean to the whole country and should be treated as provisional. Swaziland/Eswatini has only 2.6% coverage; its value is especially poorly representative. Missing pixels are not assumed malaria-free. DR Congo has 96.8% coverage. The 95% threshold is a reporting flag, not an exclusion rule.

Each country's 2015 national PfPR2-10 is evaluated on each fitted age-specific spline. The zero-PfPR hazard ratio is exp[f_g(0)-f_g(P_country)]. Counterfactual mortality rate = IHME rate x HR; attributable rate = IHME rate x (1-HR). Death counts use the same fraction and fixed annual person-time. No adjustment is made to other covariates, calendar year, or country random effects in this contrast: those terms cancel in the current additive model.

The seven PfPR effects come from the latest shared-calendar-year mortality model and its ten HIV-incidence uncertainty refits. Contrasts are pooled on the log-HR scale using within-fit covariance plus between-imputation variance and finite-imputation t intervals. Reported point HRs exponentiate pooled mean log HRs. Bounds condition on IHME and MAP point estimates, fitted smoothing parameters, and the age allocation. They do not include source-estimate, survey-design or residual-clustering uncertainty. National total death counts are point estimates; marginal age-band interval endpoints are not summed into total intervals.

## Source handling and assumptions

- The selected IHME export is the file dated 2026-09-09 10-58-22. Only All causes / Deaths / Both sexes / 2015 is used. Six disjoint source age groups are checked against Under 1 and Under 5 totals; those aggregates are never counted twice. Original rate/count lower and upper bounds are retained in the source output.
- Early and late neonatal deaths are added. Their combined annual rate is total deaths divided by the sum of their implied person-years, not the sum or simple mean of rates. The IHME 0-27-day group uses the model's <1-completed-month PfPR effect, an explicit boundary approximation.
- As authorized, ages 24-35, 36-47 and 48-59 months each use the IHME 2-4-year mortality rate and one-third of its deaths/person-time. Each then receives its own fitted PfPR contrast. These three baseline age estimates are assumed, not separately observed in IHME.
- National exposure is calculated from local 2015 MAP rasters and country polygons, using GPW 2020 density x grid-cell area x polygon overlap as weights. The earlier extraction used density alone. Both values are saved for comparison. This uses fixed all-age 2020 population geography, not a 2015 age-specific population surface. Coverage is conditional on the available raster footprint; missing MAP is never set to zero.
- The national-mean PfPR scenario is the requested approximation; the nonlinear response evaluated at a country mean does not equal a subnational burden aggregation. The fitted PfPR curves are transported to countries outside the mortality fitting sample.
- Negative attributable effects are retained. They mean the fitted association predicts higher mortality at zero PfPR, not protective malaria established by evidence. Zero-PfPR support flags and central-range flags are exported; the model's extrapolation, causal transport and smoothing-stability limitations remain relevant.

## DR Congo

National PfPR: **31.29%**. IHME under-five deaths: **248,211**. Signed model-attributable deaths: **57,043** (**23.0%**), holding the annual population exposure fixed.

All rates below are per 100,000 person-years. Age-band intervals reflect model/HIV-imputation uncertainty only.

| Age | IHME rate | Zero-PfPR rate | Attributable rate (95% interval) | Attributable annual deaths |
|---|---:|---:|---:|---:|
| 0-27 days | 29092.8 | 30685.1 | -1592.3 (-3068.6 to -183.7) | -3582 |
| 1-5 months | 3507.4 | 3384.9 | 122.6 (-77.1 to 311.1) | 1499 |
| 6-11 months | 3072.2 | 2236.6 | 835.6 (572.8 to 1070.7) | 11743 |
| 12-23 months | 1128.0 | 617.2 | 510.8 (426.3 to 585.1) | 14058 |
| 24-35 months | 843.6 | 350.6 | 493.0 (438.8 to 540.0) | 12787 |
| 36-47 months | 843.6 | 421.5 | 422.1 (346.3 to 486.3) | 10947 |
| 48-59 months | 843.6 | 473.8 | 369.8 (265.8 to 455.1) | 9591 |

![DRC malaria contribution](drc_malaria_contribution_2015.png)

The survival and probability panels describe a synthetic cohort exposed to the 2015 mortality schedule. They use q=1-exp(-H), with IHME neonatal hazards integrated separately over 7 and 21 days, followed by ages 28 days to 6 months, 6-12 months, and annual age intervals. This explicitly handles the source neonatal boundary. The survival product is checked numerically. Neither these probabilities nor their differences are used as annual death-count fractions.

## Outputs

- [All country/age estimates](country_age_attributable_2015.csv)
- [Country totals](country_totals_2015.csv)
- [National prevalence and coverage](national_pfpr_2015.csv)
- [IHME original 2015 rows](ihme_source_2015.csv) and [disjoint age inputs](ihme_disjoint_age_inputs.csv)
- [Individual-fit contrasts](individual_log_hazard_contrasts.csv), [PfPR support](model_pfpr_support.csv)
- [DRC period life table](drc_period_life_table_2015.csv), [provenance](provenance.csv)

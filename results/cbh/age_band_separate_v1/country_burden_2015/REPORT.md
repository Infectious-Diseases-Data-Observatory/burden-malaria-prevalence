# National malaria-attributable mortality, 2015

**Mortality model:** Separate age-band models; one median HIV imputation (result ID: `age_band_separate_v1`).

Estimates are available for **42 of 45 countries** in the new IHME export. Countries without usable MAP coverage remain missing: Cape Verde, Lesotho, Sao Tome and Principe.

**Exposure coverage:** 11 estimated countries have MAP covering less than 95% of population weight within the raster footprint. Their estimates apply the covered-area mean to the whole country and should be treated as provisional. Swaziland/Eswatini has only 2.6% coverage; its value is especially poorly representative. Missing pixels are not assumed malaria-free. DR Congo has 96.8% coverage. The 95% threshold is a reporting flag, not an exclusion rule.

Each country's 2015 national PfPR2-10 is evaluated on each fitted age-specific spline. The zero-PfPR hazard ratio is exp[f_g(0)-f_g(P_country)]. Counterfactual mortality rate = IHME rate x HR; attributable rate = IHME rate x (1-HR). Death counts use the same fraction and fixed annual person-time. No adjustment is made to other covariates, calendar year, or country random effects in this contrast: those terms cancel in the current additive model.

The seven PfPR effects come from seven separately fitted age-band models, each with its own calendar-year spline, confounder coefficients and random-effect variances. The saved posterior-median child HIV-incidence imputation is held fixed. Age-band intervals use the within-fit coefficient covariance and normal 95% limits on the log-HR scale; HIV-imputation uncertainty is not included. Bounds condition on IHME and MAP point estimates, fitted smoothing parameters, and the age allocation. They do not include source-estimate, survey-design or residual-clustering uncertainty. National total death counts are point estimates; marginal age-band interval endpoints are not summed into total intervals. Separate fitting does not imply independent sampling errors across age bands.

## Source handling and assumptions

- The selected IHME export is the file dated 2026-09-09 10-58-22. Only All causes / Deaths / Both sexes / 2015 is used. Six disjoint source age groups are checked against Under 1 and Under 5 totals; those aggregates are never counted twice. Original rate/count lower and upper bounds are retained in the source output.
- Early and late neonatal deaths are added. Their combined annual rate is total deaths divided by the sum of their implied person-years, not the sum or simple mean of rates. The IHME 0-27-day group uses the model's <1-completed-month PfPR effect, an explicit boundary approximation.
- As authorized, ages 24-35, 36-47 and 48-59 months each use the IHME 2-4-year mortality rate and one-third of its deaths/person-time. Each then receives its own fitted PfPR contrast. These three baseline age estimates are assumed, not separately observed in IHME.
- National exposure is calculated from local 2015 MAP rasters and country polygons, using GPW 2020 density x grid-cell area x polygon overlap as weights. The earlier extraction used density alone. Both values are saved for comparison. This uses fixed all-age 2020 population geography, not a 2015 age-specific population surface. Coverage is conditional on the available raster footprint; missing MAP is never set to zero.
- The national-mean PfPR scenario is the requested approximation; the nonlinear response evaluated at a country mean does not equal a subnational burden aggregation. The fitted PfPR curves are transported to countries outside the mortality fitting sample.
- Negative attributable effects are retained. They mean the fitted association predicts higher mortality at zero PfPR, not protective malaria established by evidence. Zero-PfPR support flags and central-range flags are exported; the model's extrapolation, causal transport and smoothing-stability limitations remain relevant.

## DR Congo

National PfPR: **31.29%**. IHME under-five deaths: **248,211**. Signed model-attributable deaths: **71,147** (**28.7%**), holding the annual population exposure fixed.

All rates below are per 100,000 person-years. Age-band intervals reflect conditional model uncertainty only; HIV imputation held fixed.

| Age | IHME rate | Zero-PfPR rate | Attributable rate (95% interval) | Attributable annual deaths |
|---|---:|---:|---:|---:|
| 0-27 days | 29092.8 | 27970.8 | 1121.9 (-936.0 to 3038.8) | 2524 |
| 1-5 months | 3507.4 | 3105.1 | 402.4 (76.3 to 697.4) | 4922 |
| 6-11 months | 3072.2 | 2130.5 | 941.7 (647.3 to 1200.3) | 13234 |
| 12-23 months | 1128.0 | 583.0 | 545.0 (440.9 to 633.3) | 14999 |
| 24-35 months | 843.6 | 354.2 | 489.4 (420.3 to 547.2) | 12694 |
| 36-47 months | 843.6 | 379.2 | 464.5 (385.0 to 530.1) | 12046 |
| 48-59 months | 843.6 | 430.0 | 413.6 (308.8 to 498.0) | 10728 |

![DRC malaria contribution](drc_malaria_contribution_2015.png)

The survival and probability panels describe a synthetic cohort exposed to the 2015 mortality schedule. They use q=1-exp(-H), with IHME neonatal hazards integrated separately over 7 and 21 days, followed by ages 28 days to 6 months, 6-12 months, and annual age intervals. This explicitly handles the source neonatal boundary. The survival product is checked numerically. Neither these probabilities nor their differences are used as annual death-count fractions.

## Outputs

- [All country/age estimates](country_age_attributable_2015.csv)
- [Country totals](country_totals_2015.csv)
- [National prevalence and coverage](national_pfpr_2015.csv)
- [IHME original 2015 rows](ihme_source_2015.csv) and [disjoint age inputs](ihme_disjoint_age_inputs.csv)
- [Individual-fit contrasts](individual_log_hazard_contrasts.csv), [PfPR support](model_pfpr_support.csv)
- [DRC period life table](drc_period_life_table_2015.csv), [provenance](provenance.csv)

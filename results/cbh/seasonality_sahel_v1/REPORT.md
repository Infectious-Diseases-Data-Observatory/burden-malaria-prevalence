# Sahelian child mortality by calendar month

Regions with a boundary centroid at or above 11 degrees north and west of 36 degrees east, excluding ETH, ERI, SOM, DJI: 247 survey regions in 41 surveys (BEN, BFA, CMR, GIN, GMB, MLI, MRT, NER, NGA, SEN, TCD).

Month of death is derived as month of birth plus age at death in completed months (exact for deaths reported in days, within about half a month for deaths reported in months; deaths reported as one year cannot be placed and those children are excluded). Children with a DHS-imputed birth month are excluded. Exposure is every calendar month in the 60 months before the interview month in which the child was alive and aged 0 to 23 completed months.

## Seasonality by age band

| Band | Surveys | Deaths (weighted) | Child-years (weighted) | Annual rate per 1,000 | Wald chi-square (11 df) | p | AIC gain from month | Peak | Trough | Peak / trough (95% CI) |
|---|---|---|---|---|---|---|---|---|---|---|
| <1 month | 41 | 9,782 | 26,214 | 373.2 | 20.5 | 0.039 | -1.0 | Jun | Oct | 1.28 (1.12 to 1.47) |
| 1-5 months | 41 | 4,025 | 124,917 | 32.2 | 22.6 | 0.020 | 1.0 | Jul | Feb | 1.39 (1.14 to 1.70) |
| 6-11 months | 41 | 3,913 | 145,338 | 26.9 | 20.7 | 0.037 | -1.0 | Sep | Mar | 1.33 (1.10 to 1.62) |
| 12-23 months | 41 | 4,935 | 282,254 | 17.5 | 21.1 | 0.032 | -0.7 | Aug | Mar | 1.41 (1.16 to 1.72) |

## Robustness to age-at-death heaping

Reported ages at death heap at 6, 12 and 18 months. A heaped death's derived month is its birth month plus a round number of months, so it carries the birth-month pattern rather than a death-month pattern. The month-factor model refitted on the unheaped ages only:

| Band | Ages used | Deaths (weighted) | Wald chi-square (11 df) | p | Peak | Trough | Peak / trough (95% CI) |
|---|---|---|---|---|---|---|---|
| 6-11 months | excluding age 6 months | 3,156 | 25.6 | 0.007 | Sep | Mar | 1.44 (1.15 to 1.79) |
| 12-23 months | excluding ages 12 and 18 months | 2,581 | 29.4 | 0.002 | Aug | Jun | 1.53 (1.20 to 1.95) |

## Rate ratio by calendar month against the annual mean

| Band | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec |
|---| ---|---|---|---|---|---|---|---|---|---|---|--- |
| <1 month | 1.04 | 1.00 | 0.95 | 1.02 | 1.07 | 1.11 | 1.00 | 1.05 | 0.94 | 0.87 | 1.02 | 0.96 |
| 1-5 months | 0.96 | 0.83 | 0.95 | 1.07 | 1.05 | 0.96 | 1.15 | 1.07 | 0.94 | 1.10 | 0.87 | 1.11 |
| 6-11 months | 0.91 | 1.03 | 0.86 | 0.98 | 0.98 | 0.98 | 1.02 | 1.12 | 1.15 | 1.13 | 1.03 | 0.87 |
| 12-23 months | 1.01 | 1.03 | 0.83 | 0.89 | 1.07 | 0.91 | 1.04 | 1.17 | 1.07 | 0.98 | 1.10 | 0.95 |
| Births | 1.04 | 0.99 | 1.05 | 1.01 | 0.98 | 0.96 | 0.89 | 1.05 | 1.09 | 1.07 | 0.93 | 0.94 |

## Pooled monthly pattern by regional mean MAP prevalence

Descriptive pooled rates without age adjustment, by tercile-like strata of each region's mean PfPR over the exposure years. Aug-Oct is the mean rate ratio over August to October.

| Band | Stratum | Regions | Deaths (weighted) | Peak | Trough | Peak / trough | Aug-Oct |
|---|---|---|---|---|---|---|---|
| <1 month | PfPR < 10% | 117 | 2,886 | Jun | Oct | 1.55 | 0.93 |
| <1 month | PfPR 10-25% | 53 | 1,652 | Aug | Oct | 1.65 | 0.99 |
| <1 month | PfPR >= 25% | 70 | 5,245 | Jan | Mar | 1.28 | 0.95 |
| 1-5 months | PfPR < 10% | 117 | 961 | Oct | Feb | 1.73 | 1.12 |
| 1-5 months | PfPR 10-25% | 53 | 778 | Dec | Sep | 1.61 | 0.94 |
| 1-5 months | PfPR >= 25% | 70 | 2,286 | Jul | Nov | 1.52 | 1.02 |
| 6-11 months | PfPR < 10% | 117 | 700 | Feb | Mar | 1.56 | 1.06 |
| 6-11 months | PfPR 10-25% | 53 | 740 | Sep | Mar | 1.49 | 1.12 |
| 6-11 months | PfPR >= 25% | 70 | 2,473 | Sep | Jan | 1.76 | 1.19 |
| 12-23 months | PfPR < 10% | 117 | 592 | Oct | Mar | 2.06 | 1.31 |
| 12-23 months | PfPR 10-25% | 53 | 897 | Aug | Oct | 1.74 | 0.97 |
| 12-23 months | PfPR >= 25% | 70 | 3,447 | Aug | Oct | 1.27 | 1.02 |

Deaths reported as one year of age cannot be placed in a calendar month: 4,949 children with such a report are excluded, against 4,935 weighted deaths used in the 12-23 month band, so that band rests on the roughly half of its deaths that were reported in months.

## Children by status

| Status | Children |
|---|---|
| used | 1,027,136 |
| outside_sahel | 632,671 |
| birth_month_imputed | 148,997 |
| region_unmatched | 14,131 |
| death_month_unresolved_year_unit | 4,949 |
| death_after_interview | 253 |

Figure: `sahel_mortality_by_calendar_month.png`. Tables: `month_rates_and_effects.csv`, `cyclic_curves.csv`, `seasonality_tests.csv`, `month_rates_by_prevalence_stratum.csv`, `sahel_regions.csv`, `exclusions.csv`.

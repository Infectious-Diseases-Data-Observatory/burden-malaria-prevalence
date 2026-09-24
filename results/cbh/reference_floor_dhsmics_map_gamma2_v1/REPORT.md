# Burden with a 1% prevalence floor

Counterfactual PfPR[2-10] = min(current, 1%) instead of 0: deaths attributable to *P. falciparum* transmission above 1%, from the saved v7 spline components (no refitting). The 0% floor reproduces the primary exactly. Point estimates only, as in the primary. Deaths occurring at or below 1% transmission are not included; see the note at the end.

## Totals across the 42 countries (primary 0% floor in brackets)

| Year | PfPR-ACM, 1% floor (0%) | IHME | UN IGME | Rate per 1,000 child-years, 1% floor (0%) | IHME rate | UN IGME rate |
|---|---:|---:|---:|---:|---:|---:|
| 2000 | 1,124,878 (1,170,371) | 563,657 | 724,276 | 10.24 (10.66) | 5.13 | 6.59 |
| 2015 | 690,972 (725,613) | 416,122 | 416,815 | 4.36 (4.58) | 2.62 | 2.63 |
| 2024 | 589,119 (620,072) | 428,147 | 440,123 | 3.32 (3.49) | 2.41 | 2.48 |

With the 1% floor the 2024 total is 589,119, 38% above IHME and 34% above UN IGME (primary: 45% and 41%). The rate falls 57.5% from 2000 to 2015 and 23.8% from 2015 to 2024 (primary 57.1% and 23.6%; IHME 48.9% and 8.1%; UN IGME 60.1% and 5.7%). From 2000 to 2024 the rate falls 67.6% (primary 67.2%; IHME 53.0%; UN IGME 62.4%).

Countries in 2024: the model exceeds IHME in 31 of 42 (primary 36); median country ratio 1.38 (primary 1.58); Spearman correlation of rates with IHME 0.91 (primary 0.91).

## 2024 total by floor

| Floor | Deaths | × IHME | × UN IGME | Rate per 1,000 child-years |
|---|---:|---:|---:|---:|
| 0% | 620,072 | 1.45 | 1.41 | 3.49 |
| 1% | 589,119 | 1.38 | 1.34 | 3.32 |
| 2% | 557,570 | 1.30 | 1.27 | 3.14 |
| 5% | 465,565 | 1.09 | 1.06 | 2.62 |

## 2024 deaths by age band, 1% floor (0%)

| Age band | Deaths |
|---|---:|
| <1 | 47,410 (49,853) |
| 1-5 | 51,247 (54,235) |
| 6-11 | 99,178 (104,436) |
| 12-23 | 136,062 (143,183) |
| 24-35 | 96,221 (100,992) |
| 36-47 | 82,213 (86,365) |
| 48-59 | 76,788 (81,009) |

## Attributable fraction within each band, 1% floor (0%)

| Age band | PfPR 10% | 20% | 30% | 40% |
|---|---:|---:|---:|---:|
| <1 | 2.5% (2.8%) | 5.5% (5.8%) | 8.5% (8.8%) | 10.7% (10.9%) |
| 1-5 | 6.9% (7.6%) | 13.1% (13.8%) | 17.5% (18.2%) | 20.4% (21.1%) |
| 6-11 | 16.7% (18.4%) | 29.7% (31.1%) | 37.3% (38.6%) | 41.5% (42.7%) |
| 12-23 | 27.3% (30.0%) | 44.6% (46.7%) | 51.8% (53.5%) | 54.2% (55.8%) |
| 24-35 | 32.6% (35.6%) | 50.8% (53.0%) | 56.4% (58.4%) | 57.6% (59.5%) |
| 36-47 | 26.2% (28.7%) | 42.8% (44.8%) | 49.5% (51.2%) | 51.4% (53.1%) |
| 48-59 | 25.4% (27.9%) | 40.6% (42.6%) | 45.1% (47.0%) | 45.4% (47.2%) |

## Malaria-caused probability of dying before 5 in 2024 (Figure 4 analogue)

Median across the 42 countries 12.5 per 1,000 live births with the 1% floor (primary 13.1); largest SSD 33.8 (primary 35.3).

## Not included

Malaria deaths that would still occur at 1% transmission are excluded by construction. The World Malaria Report's low-transmission method would estimate them as 0.256% of *P. falciparum* cases, with an under-5 share from its quadratic in all-age mortality per 1,000 at risk. That needs a case source (MAP's prevalence-to-incidence curve at 1%, or reported cases), which is not in this repository. Under that method, 10-50 *P. falciparum* cases per 1,000 at risk correspond to roughly 0.04-0.26 under-5 deaths per 1,000 child-years (assuming under-5s are 16% of the population at risk).

Files: `year_totals.csv`, `trend_changes.csv`, `country_year.csv`, `country_year_age.csv`, `agreement_2024.csv`, `country_comparison_2024_floor1.csv`, `deaths_by_age_2024.csv`, `attributable_fraction_by_age.csv`, `under5_probability_2024.csv`. Reproduce: `Rscript R_cbh/sensitivity/reference_floor/01_burden.R`.

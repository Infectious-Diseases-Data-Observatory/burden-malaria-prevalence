# Burden with a 1% prevalence floor

Counterfactual PfPR[2-10] = min(current, 1%) instead of 0: deaths attributable to *P. falciparum* transmission above 1%, from the saved primary spline components (no refitting). The 0% floor reproduces the primary exactly. Point estimates only, as in the primary. Deaths occurring at or below 1% transmission are not included; see the note at the end.

## Totals across the 42 countries (primary 0% floor in brackets)

| Year | PfPR-ACM, 1% floor (0%) | IHME | UN IGME | Rate per 1,000 child-years, 1% floor (0%) | IHME rate | UN IGME rate |
|---|---:|---:|---:|---:|---:|---:|
| 2000 | 1,139,590 (1,185,329) | 563,657 | 724,276 | 10.38 (10.79) | 5.13 | 6.59 |
| 2015 | 701,174 (736,151) | 416,122 | 416,815 | 4.42 (4.64) | 2.62 | 2.63 |
| 2024 | 597,853 (629,136) | 428,147 | 440,123 | 3.37 (3.55) | 2.41 | 2.48 |

With the 1% floor the 2024 total is 597,853, 40% above IHME and 36% above UN IGME (primary: 47% and 43%). The rate falls 57.4% from 2000 to 2015 and 23.8% from 2015 to 2024 (primary 57.0% and 23.6%; IHME 48.9% and 8.1%; UN IGME 60.1% and 5.7%). From 2000 to 2024 the rate falls 67.5% (primary 67.1%; IHME 53.0%; UN IGME 62.4%).

Countries in 2024: the model exceeds IHME in 31 of 42 (primary 36); median country ratio 1.40 (primary 1.61); Spearman correlation of rates with IHME 0.91 (primary 0.91).

## 2024 total by floor

| Floor | Deaths | × IHME | × UN IGME | Rate per 1,000 child-years |
|---|---:|---:|---:|---:|
| 0% | 629,136 | 1.47 | 1.43 | 3.55 |
| 1% | 597,853 | 1.40 | 1.36 | 3.37 |
| 2% | 565,967 | 1.32 | 1.29 | 3.19 |
| 5% | 472,990 | 1.10 | 1.07 | 2.67 |

## 2024 deaths by age band, 1% floor (0%)

| Age band | Deaths |
|---|---:|
| <1 | 51,699 (54,402) |
| 1-5 | 51,321 (54,330) |
| 6-11 | 100,505 (105,813) |
| 12-23 | 136,109 (143,205) |
| 24-35 | 97,759 (102,563) |
| 36-47 | 83,010 (87,170) |
| 48-59 | 77,452 (81,654) |

## Attributable fraction within each band, 1% floor (0%)

| Age band | PfPR 10% | 20% | 30% | 40% |
|---|---:|---:|---:|---:|
| <1 | 2.8% (3.1%) | 6.0% (6.3%) | 9.2% (9.5%) | 11.4% (11.7%) |
| 1-5 | 6.9% (7.7%) | 13.2% (13.9%) | 17.5% (18.2%) | 20.3% (21.0%) |
| 6-11 | 16.9% (18.7%) | 30.1% (31.6%) | 37.8% (39.1%) | 42.0% (43.2%) |
| 12-23 | 27.3% (29.9%) | 44.6% (46.6%) | 51.9% (53.6%) | 54.4% (56.1%) |
| 24-35 | 33.2% (36.3%) | 51.6% (53.8%) | 57.3% (59.3%) | 58.4% (60.4%) |
| 36-47 | 26.4% (29.0%) | 43.2% (45.2%) | 50.0% (51.8%) | 52.1% (53.7%) |
| 48-59 | 25.4% (27.9%) | 40.8% (42.8%) | 45.7% (47.6%) | 46.3% (48.1%) |

## Malaria-caused probability of dying before 5 in 2024 (Figure 4 analogue)

Median across the 42 countries 12.7 per 1,000 live births with the 1% floor (primary 13.3); largest SSD 34.2 (primary 35.7).

## Not included

Malaria deaths that would still occur at 1% transmission are excluded by construction. The World Malaria Report's low-transmission method would estimate them as 0.256% of *P. falciparum* cases, with an under-5 share from its quadratic in all-age mortality per 1,000 at risk. That needs a case source (MAP's prevalence-to-incidence curve at 1%, or reported cases), which is not in this repository. Under that method, 10-50 *P. falciparum* cases per 1,000 at risk correspond to roughly 0.04-0.26 under-5 deaths per 1,000 child-years (assuming under-5s are 16% of the population at risk).

Files: `year_totals.csv`, `trend_changes.csv`, `country_year.csv`, `country_year_age.csv`, `agreement_2024.csv`, `country_comparison_2024_floor1.csv`, `deaths_by_age_2024.csv`, `attributable_fraction_by_age.csv`, `under5_probability_2024.csv`. Reproduce: `Rscript R_cbh/sensitivity/reference_floor/01_burden.R`.

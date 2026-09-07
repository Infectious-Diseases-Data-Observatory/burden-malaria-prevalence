# First full dataset build

Completed locally on 7 September 2026 using the declared input snapshots and R 4.6.0. Only aggregate checks are reported here; generated child records remain under ignored `data/derived_cbh/`. No mortality model was fitted.

## Coverage

| Quantity | Count |
|---|---:|
| Surveys in registry | 124 |
| Surveys built successfully | 120 |
| Surveys unavailable because MAP geography is missing | 4 |
| Processing failures | 0 |
| Eligible child-band rows, before exposure/covariate filtering | 6,431,363 |
| Rows with region and PfPR (`model_ready`) | 6,357,802 |
| Death outcomes in eligible rows | 91,849 |
| Death outcomes in rows with region and PfPR | 90,938 |

The four unavailable surveys are Lesotho LS41FL, LS61FL, LS71FL and LS81FL. Their child histories were not expanded in this build, so the table's eligible counts describe the other 120 surveys only.

| Completed months | Eligible rows | Death outcomes |
|---|---:|---:|
| <1 | 1,144,935 | 33,372 |
| 1–5 | 1,023,016 | 13,695 |
| 6–11 | 978,764 | 12,689 |
| 12–23 | 844,844 | 12,484 |
| 24–35 | 837,491 | 10,069 |
| 36–47 | 811,275 | 6,029 |
| 48–59 | 791,038 | 3,511 |

The five-year entry window contained 994,618 additional band opportunities that were not yet fully observable at interview. Another 418,661 complete-band opportunities fell outside the configured 2000–2024 entry-year range. These are counts of bands, not distinct children. Before expansion, 38 histories had an ambiguous reported death band, 696 had a death age incompatible with interview timing, and 1,810 had an invalid survey weight; those histories were excluded and counted.

## Missingness requiring analysis decisions

There are 73,561 eligible rows without PfPR (1.14%). Of these, 71,904 have unmatched regional geography; the remaining 1,657 lack an exact annual match. The crosswalk contains 28 unmatched source labels across 17 surveys. See `region_crosswalk.csv` and [config/README.md](config/README.md) for unresolved historical groupings, abbreviated/exclusion labels and missing MAP regions. They remain in the audit shards and are excluded by the default reader. Missing exposure is not assigned zero.

Selected candidate-covariate missingness across all eligible rows, before choosing X:

| Variable | Missing |
|---|---:|
| HIV prevalence, ages 0–14 | 12.04% |
| Hib3 coverage | 25.21% |
| PCV coverage | 58.17% |
| Rotavirus coverage | 70.42% |
| Health expenditure per capita | 0.98% |
| Political stability | 2.83% |
| Household wealth quintile | 1.79% |
| Maternal education | 0.05% |
| GDP per capita, sex, PSU | 0% |

Maternal age at birth is missing for four rows. Stratum identifiers are unavailable for all 44,636 eligible rows from UG61FL; its survey design must be reviewed before a stratified bootstrap. Other surveys' stratum identifiers are retained, but their design interpretation is not validated merely by having nonmissing codes.

Vaccine missingness reflects the explicit decision not to turn absent/pre-series estimates into assumed zeros. The country-year audit confirms that HIV is unavailable for Nigeria, Comoros and São Tomé and Príncipe in the cached child-HIV panel; the initial summary omitted the third country. These values are not silently imputed. See the [country-year missingness table](HIV_MISSINGNESS.md). A final adjustment set and missing-data strategy remain to be specified before fitting.

## Checks completed

- All 10 new R files parsed. The synthetic suite passed, covering age/death boundaries, symmetric complete-band eligibility, the exact 60-month entry rule, Gregorian/Ethiopian calendar handling, invalid histories, geographic ambiguity and donor groupings, annual joins, duplicate rejection, missingness, cache invalidation and loading/weight normalization.
- Every saved survey passed unique child-band keys, binary death outcomes, positive predetermined widths, at most one included death per child, no bands after death, entry-year alignment, five-year eligibility and interview-censoring checks.
- The final loader was exercised on ET81FL, ET8AFL and NG7BFL with sex, maternal age and HIV explicitly required. It retained 33,797 and 81,174 Ethiopian rows and reported exclusion of all 195,272 NG7BFL rows because HIV was missing. All seven unordered age-band levels were retained and selected survey weights averaged 1 within each retained survey.
- Completed-manifest hashes matched the code and external snapshots used for the build. Generated outputs were confirmed to be Git-ignored.

This verifies construction and accounting, not causal validity, model fit, cross-survey boundary harmonization, the final confounder set or a survey-design variance estimator. External extraction snapshots retain their upstream assumptions; migrating their ingestion/extraction into the new folder is a later pipeline stage.

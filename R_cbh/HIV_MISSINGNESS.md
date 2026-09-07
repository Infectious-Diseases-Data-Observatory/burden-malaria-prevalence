# HIV missingness by country and band-entry year

Checked against the completed child-band build on 7 September 2026. Counts are unweighted eligible child–age-band records, before filtering on PfPR or other confounders; they are not counts of distinct children.

HIV is missing for 100% of eligible records in Nigeria, Comoros and São Tomé and Príncipe. The cached national HIV panel contains no rows for any of these three countries. The other 33 represented countries have 0% HIV missingness in every represented entry year. Lesotho is not represented because its four surveys lacked MAP inputs and were not built.

This adds São Tomé and Príncipe to the two countries identified in the initial build summary. No HIV values were imputed or changed.

| Band-entry year | Nigeria: missing records | Comoros: missing records | São Tomé and Príncipe: missing records |
|---|---:|---:|---:|
| 2000 | 7,390 | — | — |
| 2001 | 7,191 | — | — |
| 2002 | 5,740 | — | — |
| 2003 | 11,751 | — | 342 |
| 2004 | 35,332 | — | 2,357 |
| 2005 | 37,035 | — | 2,466 |
| 2006 | 43,176 | — | 2,559 |
| 2007 | 39,275 | 697 | 2,406 |
| 2008 | 42,502 | 4,035 | 709 |
| 2009 | 46,816 | 4,137 | — |
| 2010 | 42,606 | 4,114 | — |
| 2011 | 39,619 | 3,914 | — |
| 2012 | 27,861 | 1,101 | — |
| 2013 | 10,263 | — | — |
| 2014 | 45,164 | — | — |
| 2015 | 45,279 | — | — |
| 2016 | 44,559 | — | — |
| 2017 | 41,227 | — | — |
| 2018 | 10,922 | — | — |
| 2019 | 35,272 | — | — |
| 2020 | 38,015 | — | — |
| 2021 | 36,326 | — | — |
| 2022 | 35,073 | — | — |
| 2023 | 16,759 | — | — |
| 2024 | 575 | — | — |
| **Total** | **745,728** | **17,998** | **10,839** |

Every numeric cell above has 100% missingness. A dash means no eligible records for that country–year, not 0% missingness.

Overall: **774,565 / 6,431,363 records (12.04%)** lack HIV. Among records with PfPR and region available: **774,187 / 6,357,802 (12.18%)** lack HIV.

The [complete country–year table](../data/derived_cbh/hiv_missingness_by_country_year.csv) contains all 587 represented country–years across 36 countries, including zero-missingness cells. It reports eligible denominators, missing counts and percentages, and corresponding counts/percentages for the PfPR-available subset. Years absent from this table have no eligible records, rather than an inferred missingness rate.

# Paired adolescent and child HIV values

## Reported incidence rates

Run `Rscript R_cbh/analysis/03_plot_hiv_age_pairs.R --incidence` to generate `adolescent_vs_child_hiv_incidence.png` and `paired_country_year_incidence.csv`. This selects the exact indicator `Estimated incidence rate (new HIV infection per 1,000 uninfected population)` and uses the published rates directly, without fetching or calculating a population denominator. Rates are for both sexes, ages 0–14 and 15–19, matched by country and year; regional aggregates are excluded.

There are 2,121 plotted country-year pairs in 85 countries, 2000–2024; the 2024 snapshot contains 82 countries. Of these pairs, 1,977 contain numeric point estimates in both age groups and 144 contain at least one `<0.01` entry. These limits are shown as open triangles at 0.01 and flagged in the exported table. The dashed line represents equal incidence rates. Nigeria has adolescent incidence data but no paired child incidence series in this workbook.

Missing source strings such as `.` are treated as unavailable, never as zero. Only country-years with entries in both age groups are paired; `incidence_missing_pairs.csv` records paired entries excluded for unavailable or nonpositive plotting values (zero such pairs in this build). Source-key uniqueness and positive finite plotting coordinates are checked. Log scales require no pseudocounts. Incidence describes new infections and is distinct from the existing mortality model's HIV prevalence covariate; the model and its inputs have not been changed.

## Numbers living with HIV

Run `Rscript R_cbh/analysis/03_plot_hiv_age_pairs.R` from the project root to regenerate the figure and country-year table from the existing `data/HIV_Epidemiology_Children_Adolescents_2025.xlsx` workbook. No model is fitted and no values are imputed.

Selection: country entries (excluding regional aggregates), both sexes, indicator `Estimated number of people living with HIV`, paired by ISO3 and year for ages 0–14 and 15–19. There are 2,174 country-year pairs across 87 countries, 2000–2024; the 2024 snapshot has 86 countries. Nigeria has no paired child series and is absent.

The figure shows **counts, not prevalence**. The existing local data include the child population denominator but no identified ages 15–19 denominator. Population size contributes to the association, and the age intervals have different widths. The dashed line represents equal counts, not equal prevalence.

Of the pairs, 1,326 have numeric point estimates in both age groups; 848 have at least one upper-limit entry (`<100`, `<200` or `<500`). Open triangles show the reported upper-limit coordinates, retaining the exact source strings and separate censoring flags in `paired_country_year_counts.csv`. A limit coordinate must not be used as an exact observation in an imputation model. The plot omits source uncertainty intervals for readability. Annual points repeat countries and are not independent observations.

The script checks uniqueness of source country-year-age keys and finite positive plotting coordinates. The two panels use the same logarithmic scales. The table contains published country-level estimates only.

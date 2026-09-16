# Sahelian child mortality by calendar month

Supplementary descriptive analysis, added 16 September 2026. It asks how the all-cause death rate of children under two years varies with the month of the year in Sahelian survey regions, using the same recodes, survey rules, calendar offsets and deterministic geography as the CBH dataset builder. It fits no prevalence model and does not touch the primary outputs.

Run from the project root:

```sh
Rscript R_cbh/seasonality/01_build_cells.R
Rscript R_cbh/seasonality/02_fit_and_plot.R
```

`settings.R` holds the Sahel definition (boundary centroid at or above 11°N, west of 36°E, excluding the Horn of Africa), the 60-month exposure window, the four age bands (<1, 1–5, 6–11, 12–23 completed months) and the output directory `results/cbh/seasonality_sahel_v1/`.

DHS records no date of death. The month of death is birth month plus age at death in completed months: exact for deaths reported in days, within about half a month for deaths reported in completed months, and not recoverable for deaths reported in whole years, so children whose death is reported as "1 year" are excluded and deaths at two years or more lie outside the analysed ages. Children whose birth month was imputed by DHS (flag B10) are excluded. Every valid child contributes a person-month for each calendar month in the window during which it was alive and aged 0–23 completed months; the death month counts as exposure and carries the death.

The second script plots the pooled design-weighted rate by calendar month with leave-one-survey-out jackknife intervals, then fits, for each band separately, a negative-binomial regression of the cell death count on calendar month as a factor, adjusted for age in completed months, with a survey random intercept and log person-months as the offset, plus a cyclic spline version. Month effects are rate ratios against the annual geometric mean. Births by calendar month are drawn for reference because heaped ages at death copy birth seasonality onto derived death months.

Outputs are aggregate only: survey × month × age cells, survey × region × month × band cells, births by month, the region list with centroids and mean MAP prevalence, exclusion counts, model tables, the figure and `REPORT.md`.

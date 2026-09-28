# Tensor product of PfPR and calendar time (no region random intercept)

Seven age-band models without the survey-region random intercept, fitted with bam (fREML, no discretisation) to the collapsed survey-region × entry-month cells (binomial deaths out of children entering the band; the same likelihood as the child-level model). Three specifications per band, otherwise identical (17 covariates, survey and country random intercepts, band-width offset, primary reference knots, gamma = 2):

- **Separate splines:** `s(pfpr_pct, cr, k=5) + s(calendar_year, cr, k=6)` (as the primary).
- **Tensor product:** `te(pfpr_pct, calendar_year, bs = c("cr","cr"), k = c(5, 6))`, which lets the PfPR curve change with calendar time.
- **Interaction test:** separate splines plus `ti(pfpr_pct, calendar_year)`, the pure interaction. The function space contains the separate-spline model, but the penalised fits are not nested: smoothing parameters, including the survey random-effect variance, are re-estimated.

All 21 fits converged with positive-definite smoothing Hessians; none needed a restart.

## Conclusion

Keep the separate-spline (time-constant) PfPR curves as the primary specification and report the tensor product as an exploratory sensitivity. Allowing the PfPR effect to change with calendar time improves AIC in the full-period fits only at <1 and 1–5 months, where the curves barely change. At 24–59 months the fitted low-PfPR gradient steepens over time, and the binned data show the same direction, but the evidence is mixed: in the full-period fits AIC does not favour the interaction and the log-likelihood does not improve, because the interaction takes over between-survey heterogeneity (the survey random-effect EDF falls); it is favoured when the random-effect variances are held at the separate-spline values and at 36–47 months within 2005–2019; it disappears when identified within surveys (survey fixed effects); and it rests on a few dozen deaths in PfPR <5% areas at each end of the period. It therefore cannot be separated from between-survey (country-period) differences. If the steepening were real it would raise the 2024 attributable deaths and shrink the post-2015 decline (below).

## Model fit

| Age (months) | ΔAIC tensor − separate | ΔAIC ti − separate | ti EDF | ti p-value | survey RE EDF: separate / tensor / ti | EDF of PfPR/time terms: separate, tensor |
|---|---:|---:|---:|---:|---:|---:|
| <1 | -16.5 | -18.9 | 4.72 | 0.01168 | 63.9 / 64.5 / 64.1 | 3.9, 6.2 |
| 1-5 | -6.7 | -5.6 | 2.35 | 0.01766 | 28.4 / 28.1 / 28.0 | 5.2, 9.2 |
| 6-11 | -1.0 | 2.9 | 2.63 | 0.28327 | 39.7 / 39.4 / 39.2 | 5.3, 9.2 |
| 12-23 | 4.7 | 2.0 | 1.53 | 0.21662 | 50.4 / 51.0 / 50.2 | 6.5, 9.9 |
| 24-35 | 7.3 | 6.9 | 1.00 | 0.01567 | 33.9 / 32.1 / 32.2 | 4.7, 7.3 |
| 36-47 | 1.3 | 2.3 | 2.31 | 0.00059 | 6.2 / 0.3 / 0.9 | 4.3, 6.7 |
| 48-59 | 1.2 | 2.2 | 2.10 | 0.00374 | 12.3 / 8.0 / 8.1 | 4.2, 6.4 |

AIC is mgcv's (binomial likelihood with the smoothing-corrected degrees of freedom, `aic_df` in [fit_diagnostics.csv](fit_diagnostics.csv); the `loglik` column there is the kernel log-likelihood without the binomial coefficient). Negative ΔAIC favours the more flexible model; |ΔAIC| below about 3 (6–11, 36–47 and 48–59 months) is no discernible difference, and its sign depends on the degrees-of-freedom definition. The ti p-value is mgcv's approximate test conditional on the estimated smoothing parameters; in the older bands it coexists with no gain in log-likelihood because the survey random effect is shrunk when the interaction enters, so it should not be read as stand-alone evidence of effect modification.

## Hazard ratio for PfPR 20% → 0% by calendar year of band entry (evaluated at mid-year)

| Age (months) | Separate (all years) | Tensor mid-2005 | mid-2010 | mid-2015 | mid-2020 |
|---|---:|---:|---:|---:|---:|
| <1 | 0.93 | 0.94 | 0.94 | 0.93 | 0.93 |
| 1-5 | 0.86 | 0.89 | 0.83 | 0.80 | 0.81 |
| 6-11 | 0.68 | 0.70 | 0.69 | 0.66 | 0.67 |
| 12-23 | 0.52 | 0.55 | 0.55 | 0.52 | 0.47 |
| 24-35 | 0.46 | 0.54 | 0.50 | 0.47 | 0.43 |
| 36-47 | 0.55 | 0.69 | 0.60 | 0.53 | 0.46 |
| 48-59 | 0.57 | 0.78 | 0.64 | 0.53 | 0.44 |

## Hazard ratio for PfPR 40% → 20% by calendar year (mid-year)

| Age (months) | Separate (all years) | Tensor mid-2005 | mid-2010 | mid-2015 | mid-2020 |
|---|---:|---:|---:|---:|---:|
| <1 | 0.93 | 0.92 | 0.93 | 0.94 | 0.95 |
| 1-5 | 0.92 | 0.88 | 0.93 | 0.97 | 0.96 |
| 6-11 | 0.83 | 0.83 | 0.86 | 0.85 | 0.81 |
| 12-23 | 0.80 | 0.78 | 0.78 | 0.80 | 0.85 |
| 24-35 | 0.84 | 0.87 | 0.83 | 0.80 | 0.77 |
| 36-47 | 0.84 | 0.87 | 0.83 | 0.79 | 0.75 |
| 48-59 | 0.91 | 0.92 | 0.89 | 0.86 | 0.83 |

## Support and identification checks (24–59 months; [03_checks.R](../../../R_cbh/sensitivity/tensor_pfpr_year/03_checks.R))

Deaths in cells with PfPR below 5%, by period of band entry:

| Age (months) | 2000-04 | 2005-09 | 2010-14 | 2015-19 | 2020-23 |
|---|---:|---:|---:|---:|---:|
| 24-35 | 68 | 139 | 294 | 161 | 50 |
| 36-47 | 44 | 104 | 184 | 108 | 40 |
| 48-59 | 36 |  67 | 147 |  73 | 25 |

Binned adjusted log hazard ratio for PfPR <5% versus 20–40% by period (survey and country random intercepts, calendar spline, 17 covariates), with the deaths in the <5% class:

| Age (months) | 2000-05 | 2006-10 | 2011-15 | 2016-23 |
|---|---:|---:|---:|---:|
| 24-35 | -0.49 (77) | -0.80 (192) | -0.82 (262) | -0.86 (181) |
| 36-47 | -0.35 (48) | -0.66 (133) | -0.78 (167) | -0.78 (132) |
| 48-59 | 0.01 (42) | -0.59 (87) | -0.54 (135) | -0.88 (84) |

Interaction tests (ΔAIC against the matching model without the interaction; HR for PfPR 0% vs 20% at mid-2005 and mid-2020):

| Age (months) | Check | ΔAIC | ti EDF | ti p | survey RE EDF | HR 2005 | HR 2020 |
|---|---|---:|---:|---:|---:|---:|---:|
| 24-35 | ti, free (as in 01_fit.R) | 6.9 | 1.00 | 0.01567 | 32.2 | 0.49 | 0.44 |
| 24-35 | ti, survey and country RE variances fixed at separate fit | -2.2 | 1.00 | 0.01788 | 33.7 | 0.49 | 0.44 |
| 24-35 | ti, survey fixed effects (versus separate with survey fixed effects) | 0.5 | 1.00 | 0.28038 | (fixed effects) | 0.51 | 0.47 |
| 24-35 | tensor, entries 2005-2019 (versus separate, same window) | 2.9 | – | – | 21.0 | 0.54 | 0.40 |
| 24-35 | ti, entries 2005-2019 (versus separate, same window) | -0.2 | 1.00 | 0.04313 | 21.3 | 0.48 | 0.42 |
| 36-47 | ti, free (as in 01_fit.R) | 2.3 | 2.31 | 0.00059 | 0.9 | 0.67 | 0.46 |
| 36-47 | ti, survey and country RE variances fixed at separate fit | -17.9 | 2.27 | 0.00123 | 6.1 | 0.66 | 0.46 |
| 36-47 | ti, survey fixed effects (versus separate with survey fixed effects) | 1.6 | 1.33 | 0.57976 | (fixed effects) | 0.62 | 0.56 |
| 36-47 | tensor, entries 2005-2019 (versus separate, same window) | -0.5 | – | – | 1.3 | 0.60 | 0.47 |
| 36-47 | ti, entries 2005-2019 (versus separate, same window) | -2.2 | 1.00 | 0.03255 | 1.7 | 0.56 | 0.47 |
| 48-59 | ti, free (as in 01_fit.R) | 2.2 | 2.10 | 0.00374 | 8.1 | 0.73 | 0.46 |
| 48-59 | ti, survey and country RE variances fixed at separate fit | -13.9 | 2.07 | 0.00551 | 12.1 | 0.72 | 0.46 |
| 48-59 | ti, survey fixed effects (versus separate with survey fixed effects) | 2.0 | 1.00 | 0.93324 | (fixed effects) | 0.68 | 0.67 |
| 48-59 | tensor, entries 2005-2019 (versus separate, same window) | 6.2 | – | – | 15.3 | 0.70 | 0.45 |
| 48-59 | ti, entries 2005-2019 (versus separate, same window) | 4.3 | 1.54 | 0.31794 | 13.8 | 0.56 | 0.43 |

Full tables: [support_pfpr_class_by_period.csv](support_pfpr_class_by_period.csv), [binned_pfpr_by_period.csv](binned_pfpr_by_period.csv), [interaction_checks.csv](interaction_checks.csv).

## Burden across the 42 countries (IHME all-cause inputs, year-specific hazard ratios)

Calendar time is evaluated at mid-year, held at December 2023 (the last observed band entry and last calendar-year knot) for 2024.

| Model | 2000 | 2015 | 2024 | 2024 × IHME | Rate change 2000–2015 | Rate change 2015–2024 |
|---|---:|---:|---:|---:|---:|---:|
| Primary v9 (with region RE, child-level) | 1,185,329 | 736,151 | 629,136 | 1.47 | -57.0% | -23.6% |
| Separate splines | 1,207,674 | 747,670 | 638,263 | 1.49 | -57.1% | -23.7% |
| Tensor product | 1,059,296 | 794,182 | 752,252 | 1.76 | -48.1% | -15.4% |
| Separate + ti (see note) | 1,209,612 | 786,872 | 769,981 | 1.80 | -54.9% | -12.6% |

2024 under other evaluation times: tensor 755,245 at mid-2024 (extrapolated) and 723,491 with the mid-2019 hazard ratios; ti 783,816 and 681,518.

2024 attributable deaths by age band:

| Age (months) | Separate | Tensor | ti |
|---|---:|---:|---:|
| <1 | 56,927 | 66,641 | 99,743 |
| 1-5 | 54,349 | 63,844 | 72,680 |
| 6-11 | 107,716 | 108,868 | 115,550 |
| 12-23 | 148,198 | 169,517 | 148,936 |
| 24-35 | 102,362 | 113,360 | 109,126 |
| 36-47 | 87,043 | 112,969 | 111,709 |
| 48-59 | 81,667 | 117,054 | 112,236 |

Note: the ti model's neonatal interaction (about 4.6 EDF) oscillates over time around a hazard ratio close to 1; applied to the large neonatal all-cause totals this moves tens of thousands of deaths, so the ti burden row is unstable and is shown only for completeness. The tensor product, which spans essentially the same function space, does not show this. The time-varying models' lower 2000 and higher 2024 totals, and hence their smaller declines, follow directly from the assumed steepening of the PfPR effect and are not separate findings.

PfPR curves of the tensor product by year of band entry (relative to PfPR 20% in the same year, 95% intervals conditional on smoothing parameters; solid within the period's children-weighted 2.5th–97.5th PfPR percentiles, dashed outside), with the time-constant separate-spline curve for reference:

![PfPR curves by year](sfig_pfpr_curves_by_year_tensor.png)

![HR 20% to 0% by year](sfig_hr_20_to_0_by_year.png)

Figure 3 analogue with the tensor-product values (countries, Nigerian states and the annual series, with the primary shown for reference): [figure3/fig3_burden_comparison_tensor.png](figure3/fig3_burden_comparison_tensor.png), caption in [figure3/CAPTION.md](figure3/CAPTION.md) (`04_figure3.R`).

Files: `model_fit_comparison.csv`, `hazard_ratios_by_year.csv`, `hr_20_to_0_by_year.csv`, `pfpr_curves_by_year.csv`, `pfpr_support_by_period.csv`, `year_totals.csv`, `deaths_by_age_2024.csv`, `burden_2024_alternative_times.csv`, `fit_diagnostics.csv`, `smooth_summaries.csv`, `model_formulas.txt`. Reproduce: `Rscript R_cbh/sensitivity/tensor_pfpr_year/01_fit.R`, then `03_checks.R`, `02_report.R` and `04_figure3.R`.

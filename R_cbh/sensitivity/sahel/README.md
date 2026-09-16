# Sahel restriction of the primary PfPR analysis

Run from the project root:

```sh
Rscript R_cbh/sensitivity/sahel/01_fit.R
Rscript R_cbh/sensitivity/sahel/02_report.R
```

`01_fit.R` subsets the saved primary prepared data using exact survey/region keys and the existing boundary-centroid cache. The geographic rule is latitude ≥12°N, longitude <36°E and exclusion of ETH/ERI/SOM/DJI. All primary region keys must match a finite location. The full selection audit is saved; no raw DHS, MAP or HIV rebuilding takes place. The older monthly-mortality analysis in `R_cbh/seasonality` keeps its own 11°N rule.

Seven separate binomial cloglog `bam` models retain the primary adjustment set, global scaling, median HIV incidence imputation, unweighted likelihood, full-band offset, `cr` splines (PfPR k=5; year k=6), fREML and gamma=2. Subgroup knots follow the `cr` convention of quantiles of distinct predictor values and are explicitly supplied/recorded. Fits are cached against input/code/settings signatures; `--force` refits. Numerical diagnostics and model-input equality checks gate completion.

`02_report.R` reads only saved aggregate results. It overlays the subgroup and full primary curves, relative to 20% PfPR, over their respective central 95% exposure ranges. Pointwise conditional intervals do not account for covariance between the nested samples. It also saves support flags, 40%→20% and 20%→0% contrasts, a centroid-selection figure and a report. Zero-PfPR contrasts require extrapolation if flagged outside support.

Results: `results/cbh/sahel_map_gamma2_v1/`. Private fitted objects: `data/derived_cbh/models/sahel_map_gamma2_v1/`. No primary outputs, country burden estimates or Overleaf files are overwritten.

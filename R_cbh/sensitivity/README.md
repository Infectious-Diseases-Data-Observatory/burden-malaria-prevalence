# Single-imputation PfPR sensitivity fits

Run from the project root after the current shared-time mortality model has been fitted:

```sh
Rscript R_cbh/sensitivity/01_fit.R
Rscript R_cbh/sensitivity/02_report.R
```

The saved posterior-median child HIV incidence dataset is held fixed. The base joint fit is reused, and eleven additional `mgcv::bam` models are fitted:

- Seven independent age-band models, with their own PfPR and calendar-time smooths, confounder coefficients and survey/country/region random intercepts.
- Two joint seven-band models for UN M49 Western Africa versus Eastern plus Middle Africa. Southern Africa is excluded from this comparison by default. `--include-south` on both scripts includes it with East/Central and writes a separate output directory.
- Two joint seven-band models split at the median survey year among included surveys. Surveys at or below the median are early; surveys above it are late. Each survey is counted once when computing the median and stays intact when assigning its child-band records. The current cutoff is 2012.

These are three separate sensitivities, not a crossed age-by-geography-by-period analysis. Other covariates and their saved scaling, data eligibility, unweighted likelihood and the band-width offset match the reference. Every new fit re-estimates smoothing parameters. Single-age models replace the country-by-age random intercept with a country random intercept; these identify the same groups within a single age band. Joint subgroup models retain the original formula, including one time spline shared across age bands within that fit.

Curves are plotted as `f_g(P) - f_g(20)`, with pointwise 95% conditional intervals and common axes. Each curve is restricted to its own central 95% observed PfPR range. Full curves and support flags, plus the 40%-to-20% hazard ratio, are saved as aggregate tables. Comparisons with the reference are descriptive because samples overlap; the interval calculation does not estimate covariance between fits.

The single-age and subgroup fits have no HIV imputation-uncertainty propagation. Existing survey-design and smoothing-stability limitations remain. Geography follows [UN M49](https://unstats.un.org/unsd/methodology/m49/overview/), accessed 9 September 2026; exact assignments appear in `survey_groups.csv`.

Outputs:

- Aggregate reports: `results/cbh/age_band_hiv_incidence_shared_time_v3/sensitivity_single_imputation/`.
- Private fitted objects: `data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/sensitivity_single_imputation/` (ignored by Git).
- With `--include-south`, both paths end in `sensitivity_single_imputation_with_south/`.

The fitting script validates convergence, finite coefficients and covariance, unchanged sample size, and agreement between prediction-matrix contrasts and direct link predictions. Identical exposures must give exactly zero log hazard ratio and uncertainty. Cache signatures include the input signature, group assignments, formula, fitting code and R/mgcv versions. `--force` explicitly refits; the report script reads only aggregate files and never refits a model.

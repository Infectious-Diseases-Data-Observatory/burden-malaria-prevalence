# Current primary analysis: regional adjustment

Run from the project root:

```sh
Rscript run_all.R --primary
# Resume with verified fit caches:
Rscript R_cbh/primary/run_regional.R --resume
# Effects, diagnostics and comparison from saved fits:
Rscript R_cbh/primary/run_regional.R --report-only
```

The current output version is `results/cbh/primary_map_regional18_gamma2_v2/`.
Its [report](../../results/cbh/primary_map_regional18_gamma2_v2/REPORT.md) compares
the seven MAP gamma=2 fits with the preserved eleven-variable benchmark in
`results/cbh/primary_map_gamma2_v1/`.

The 18-variable model uses 14 survey-region summaries and four national annual
covariates, as defined by `cbh_regional_spec()` in `covariates/regional.R`.
Urban percentage is included. Hib3, PCV, rotavirus and exclusive breastfeeding
are excluded from both the model and its complete-case requirements.

Stages:

1. `00_prepare_regional.R` assembles existing child-band shards, the audited
   reduced regional overlay and the fixed HIV incidence panel. It verifies
   unique child-band keys and the audit's selection counts, saves new scaling,
   and caches the resulting modelling dataset with source hashes.
2. `01_fit.R` fits seven separate binomial/cloglog models with `gamma=2`,
   PfPR `cr/k=5`, calendar time `cr/k=6`, reference MAP knots, survey/country/region
   random intercepts and a full-band-width offset. Each age has its own
   covariate coefficients, time spline and random-effect variances.
3. `02_effects.R` calculates supported contrasts, zero-PfPR contrasts and
   2005/2015/2024 country estimates using unchanged national MAP/IHME inputs.
4. `04_diagnostics.R` verifies stored fitted inputs and aggregate outcomes.
5. `05_compare_regional.R` produces aggregate CSVs, PNGs and a Markdown report.

The new runner sets `CBH_PRIMARY_VERSION=regional` for its subprocesses.
Shared-stage scripts retain the legacy default for compatibility with the
older paper-reporting workflow. Use the new runner for the current primary
analysis. Direct execution of shared stages requires the same environment
variable. `00_prepare_regional.R` and `05_compare_regional.R` explicitly select
the regional version.

No raw-source extraction, HIV refitting, supplementary fitting or TeX output
is part of this runner. Existing Overleaf figures, annual comparison, Nigerian
state analysis and sensitivity fits remain from the previous iteration until
regenerated separately. The older runner `run.R` reproduces those historical
outputs and includes a TeX table writer; it is not invoked by the current runner.

Model objects and the prepared dataset remain under ignored
`data/derived_cbh/models/primary_map_regional18_gamma2_v2/`. Only aggregate
results, code and provenance are exported. Model cache signatures include data,
formula, knots, code and software versions; invalid caches are refitted.

The comparison reflects both the revised adjustment and the changed sample.
It does not identify which individual covariate or aggregation change caused
an effect difference. Both iterations use the same fixed HIV imputation and
conditional spline uncertainty; uncertainty in exposure, filled covariates and
smoothing parameters is not propagated into the reported intervals.

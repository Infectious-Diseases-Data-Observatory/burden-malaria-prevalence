# Annual Snow prevalence sensitivity

Run from the project root, using installed packages and local files only:

```sh
SNOW_SOURCE=annual_csv SNOW_SKIP_COMPARISONS=1 Rscript R_dhs/51_snow_polygon_prevalence.R
Rscript R_cbh/snow/00_audit_extraction.R
Rscript R_cbh/snow/01_prepare.R
Rscript R_cbh/snow/02_fit.R
Rscript R_cbh/snow/05_check_smoothing.R
Rscript R_cbh/snow/04_verify_saved_fits.R
Rscript R_cbh/snow/03_report.R
```

`annual_csv` is script 51's source option; the actual files are in
`data/snow_prevalence_model/`. The annual posterior mean is in percent, on the
microscopy-equivalent PfPR2–10 scale. This is the supplied annual re-fit of the
Snow survey database, not the published five-year model.

Script 51 now validates the complete polygon-year grid, hashes sources and
survey boundaries before reusing caches, restricts intersections to the same
country, and never carries 2015 forward. It uses installed `exactextractr` for
fractional-cell population sums. Both this method and the previous
`terra::extract(..., exact=TRUE, fun=sum)` correctly account for overlap fractions;
the change of extractor improves speed. Population weights use GPW 2020 density
times cell area and are fixed across years. CSV regional means are exact for
those weights; medians and quantiles are unavailable without aggregating draws.
The independent-polygon SD is only a benchmark, not a universal lower bound.

The audit independently recomputes every region mean from the supplied posterior
draws, validates fractional-cell weighting against an analytic example, checks
country membership and annual keys, and compares with the earlier extraction.
The original annual-draw outputs remain separate from the `_annual_means` outputs.

Preparation reuses the validated child-band shards, without requiring MAP
availability for the Snow analysis. It joins **`period`** from the Snow panel to
the band's **entry year**; the panel's `year` field means survey year and must
not be used as exposure year. Existing reviewed DHS region crosswalks remain in
use. Missing estimates and population coverage below 50% exclude a record;
coverage over 101% is rejected as a geometry problem. Small deviations around
100% are retained and visible in the audit.

The default restricts survey metadata and individual interview years to 2015
or earlier. Full potential bands must finish by the end of 2015. Entry years
remain 2000–2015 and within five years before interview. To retain retrospective
bands completed by 2015 from later surveys, pass `--retain-later-surveys` to
preparation, fitting and reporting; this writes a separate output directory.

Seven separate age-band Snow fits use the same confounders, saved scaling,
fixed median child HIV-incidence imputation, unweighted cloglog likelihood,
full-band-width offset, spline dimensions and random-effect structure as the
primary separate-age model. Every fit re-estimates smoothing parameters.
MAP fits use exactly the same records with both exposures observed. If any
Snow-eligible records lack MAP, seven matched Snow fits are added to isolate
the exposure comparison. Original full-period MAP curves provide context.
`--force` explicitly refits; otherwise input/code/software signatures control
cache reuse.

The final verification rereads each saved fit and checks that the exposure,
outcome, every nuisance predictor and offset match its intended prepared rows.
Pass/fail summaries contain no individual records. Pass the same optional
retrospective flag to this verification when using that variant.

`05_check_smoothing.R` runs only for fits with a nonpositive smoothing Hessian.
It restarts with tighter convergence tolerances and saves the numerical sensitivity
separately. A restart is selected only if it converges with full rank, a positive
smoothing Hessian, a reduced gradient and no worse fREML criterion. Original fits
are preserved; selection does not depend on effect size or significance.
`selected_fit_manifest.csv` explicitly names and hashes the fit supplying every
reported curve. Reporting reads the selected aggregate curves and diagnostics.

Private model data and fitted objects:
`data/derived_cbh/models/age_band_snow_2000_2015_v1/`.
Aggregate audits, selection tables, diagnostics, curves and report:
`results/cbh/age_band_snow_2000_2015_v1/`.
The optional retrospective variant uses
`age_band_snow_2000_2015_retrospective_v1`.

Source prevalence and HIV-imputation uncertainty are held fixed. Intervals are
conditional model intervals and do not account for survey design or smoothing
parameter uncertainty. Source MCMC and sparse-data limitations from the input
README are carried into the report. This stage changes neither the primary
MAP fits nor the national mortality burden estimates.

# Primary MAP pipeline

**17 September specification change:** the commands below reproduce the saved **previous adjustment set**. The revised analysis plan requires regional summaries and the expanded covariates in [../covariates/README.md](../covariates/README.md). Its extraction/availability audit is implemented; fitting and paper results have not yet been updated. Keep these historical outputs separate until a new prepared sample and versioned fitter are promoted.

Run from the project root:

```sh
Rscript run_all.R --primary
```

This forces seven fresh age-band fits, then recalculates contrasts, national mortality for 2005/2015/2024, aggregate diagnostics, primary-only figures, the inclusion flow and the key-results index. It does not rebuild the dataset, refit HIV imputation, extract rasters, or run Snow/sensitivity models.

- `01_fit.R`: unweighted binomial cloglog `bam`, `gamma=2`, `cr` PfPR (`k=5`) and calendar year (`k=6`), fixed full-band-width offset, saved confounder scaling, separate survey/country/region random effects in each band. Prepared data are the fitting input; the committed knot snapshot preserves the promoted primary basis. Fit caches are validated by data/code/knots/software hashes. Tight-tolerance restarts are allowed only for numerical failures, under the same specification.
- `02_effects.R`: validated compact PfPR contrasts, conditional intervals, country-age and country totals, and a DRC synthetic-cohort life table. Reuses the existing national MAP and IHME input columns, with file hashes; no data setup. It does not reuse old attributable estimates when calculating new estimates. The old primary estimates are read only for comparison.
- `04_diagnostics.R`: verifies saved model and prepared-input hashes; aggregates fitted/observed outcomes without exporting individual records. These are in-sample checks, not held-out validation or formal influence diagnostics.
- `03_report.R`: reads only aggregate results to make figures and a report.
- Existing `reporting/` scripts regenerate the inclusion-only flow and link the current primary outputs in `Key results/README.md`.

The prepared sample is `data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/complete_case_dataset.rds`. New fitted objects stay under ignored `data/derived_cbh/models/primary_map_gamma2_v1/`; shareable aggregate outputs are in `results/cbh/primary_map_gamma2_v1/`. Older `map_snow_gamma2_v1` outputs remain unchanged.

```sh
# Continue a failed run, retaining only valid completed fit caches:
Rscript R_cbh/primary/run.R --resume
# Recalculate summaries, diagnostics and figures from current saved primary fits:
Rscript R_cbh/primary/run.R --report-only
```

`reference_knots.csv` was extracted from the hash-verified promoted MAP gamma=2 fits for the first migration. The previous comparison manifest is still read to identify the matching historical fit IDs; its models are not fitting inputs when the committed knot snapshot exists. Reference national inputs live under `results/cbh/age_band_separate_v1/country_burden_<year>/`. These dependencies are explicit and fingerprinted; this is an analysis rerun, not a reconstruction of upstream inputs.

The fixed HIV imputation, unweighted likelihood, conditional covariance, inherited regional MAP extraction issues, and national mean-exposure and IHME age-allocation assumptions remain as documented in `docs/ANALYSIS_PLAN.md` and `docs/CODE_AUDIT.md`. No aggregate confidence intervals are constructed by summing age-band limits or assuming independent age models. Negative contributions and missing-country rows are retained.

## Paper figures

`03_report.R` exports the primary PfPR and country-versus-IHME figures without titles, subtitles or embedded captions. Both country-comparison axes use a true log10 scale and common limits starting at 1,000 deaths; points below the minimum on either axis are outside the plotting window, with their estimates retained in the tables; missing/nonpositive pairs are listed in `burden/log10_plot_exclusions.csv`. The other mortality figures retain their existing scales.

`04_diagnostics.R` exports the aggregate union of included regions per survey. `reporting/03_survey_map.R` combines this with the registry and existing country outlines to regenerate the survey map for the current primary sample (105 surveys/34 countries), without fitting or rebuilding individual records. `reporting/04_paper_figures.R` records the five figure filenames, source hashes and separate captions. These stages are included in the primary runner. The runner writes only inside this project; copying the five PNGs to Overleaf is a separate explicit action and never edits TeX.

Manuscript typography: 20 pt axis titles, 16 pt ticks and legend text, and 18 pt age/year facet labels. Map labels use 4.5 mm text and Figure 4 country labels use 5 mm text. Figure 1 is 13 × 10 inches to accommodate the larger timeline labels and stacked legends.

`Rscript R_cbh/reporting/06_attributable_fraction_by_age.R` generates Figure 3 (11 × 7 inches) and its 28 estimates from the saved primary PfPR components: age bands on x, attributable share of all-cause deaths within each band on y (0–100%), and four lines for 10%, 20%, 30% and 40% PfPR. The fraction is `1 - exp(f_g(0) - f_g(P))`, checked against the saved primary curves. It is not the share of all under-five deaths occurring in each band. Zero prevalence requires extrapolation. Conditional intervals are saved in the CSV; the figure shows points and connecting lines without titles or subtitles. This stage precedes the paper manifest in the runner. Country versus IHME is now Figure 4; TeX figure order remains managed by the author.

## Age-band results table

The primary report includes each band's observed deaths and percentage of the 82,415 observed under-five deaths, PfPR spline EDF, and hazard ratios (conditional 95% intervals) for 40%→20% and 20%→0%. The fitting-time column is omitted. Percentages use largest-remainder rounding at one decimal so the seven displayed shares sum to 100.0%; the CSV also retains unrounded shares.

`Rscript R_cbh/reporting/05_age_band_table.R` regenerates the standalone Markdown, CSV and LaTeX table in `results/cbh/primary_map_gamma2_v1/tables/`. The full reporting stage uses the same generator. Paste `age_band_results.tex` into the manuscript; it uses `booktabs` and has no document preamble. No Overleaf file is changed by these commands.

`Rscript R_cbh/burden/04_annual_comparison.R --audit-only` verifies Figure 5 inputs. Without the flag it calculates annual 2004–2024 national estimates from the saved primary models and raw annual MAP rasters. `Rscript R_cbh/reporting/07_annual_mortality_comparison.R` plots rates for the same 42 countries, with a common annual under-five denominator from IHME all-cause count/rate pairs. These stages are included before the paper manifest in the runner. The UN IGME comparator is the CA-CODE 2026 release from UNICEF's `CME_CAUSE_OF_DEATH` dataflow, not the older all-age WMR proxy. The downloaded input and source metadata are in `data/external/figure5/`; the exact API URL and checksum are preserved in `results/cbh/primary_map_gamma2_v1/annual_comparison/who_source.json`. All scripts use local files, without automatic network requests. Annual totals, country-age estimates, population/exposure coverage, input checks and a full report are in `annual_comparison/`. Precise IHME release alignment remains unverified; small denominator differences are documented.

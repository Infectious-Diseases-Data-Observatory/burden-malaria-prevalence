# Cleanup verification — 15 September 2026

- 133 active and newly archived R files parsed without evaluation: `parse_checks.csv`.
- All 71 preserved R files retain the SHA-256 hashes in `archive_manifest.csv`; the 68 retired script paths are absent. Three old entry-point paths contain replacements.
- 72 literal `source()` / `sys.source()` references across 62 active R files resolve: `literal_source_checks.csv`. Dynamic helper paths were reviewed separately. This is not a claim that all legacy cached-output dependencies have been eliminated; retained input bridges and model-frame dependencies are listed in the component audit.
- `Rscript R_cbh/tests/test_dataset.R`: passed synthetic age/date/eligibility/join/cache/accounting cases.
- `Rscript R_cbh/tests/test_model.R`: passed synthetic historical joint-model helper checks (not a primary subgroup sensitivity).
- `Rscript R_cbh/tests/test_hiv_incidence.R`: passed censoring, country withholding, incidence join and tail-quantile checks.
- `Rscript tests/test_rkey.R` and `Rscript tests/test_region_match.R`: passed the retained legacy geographic-helper tests. The CBH matcher is stricter and is exercised by the dataset tests.
- `Rscript R_cbh/audit/01_verify_saved_primary.R`: passed all seven full-MAP fit/frame/hash/knots/contrast checks and saved national accounting checks. See `primary_fit_checks.csv`, `primary_sample.csv` and `burden_checks.csv`.
- `Rscript run_all.R --help`: lists current explicit commands without starting an analysis.
- `Rscript run_all.R --check-inputs`: completed local prerequisite checking; the four known Lesotho surveys still lack annual MAP availability.
- `Rscript R_cbh/reporting/02_results_index.R`: refreshed the index using explicit MAP gamma=2 selections; no scientific estimate was recomputed.
- Local links in the current plan, audit, new archive README and updated current result/input indices resolve. Archived READMEs preserve their original project-root context; links inside them are historical.
- `git diff --check`: passed for active changes. Original archived files are preserved byte-for-byte, including pre-existing whitespace.

No raw recodes were read during this audit; local prepared model data were used solely for aggregate verification. No network access, installation, raster extraction or DHS/HIV refitting occurred. The root README's pre-existing Snow-script row and unrelated untracked files were preserved.

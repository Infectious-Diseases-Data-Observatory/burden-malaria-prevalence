# Snow prevalence comparison archived on 18 September 2026

The supplementary comparison of MAP with the annual Snow survey-database prevalence estimates was removed from the [analysis plan](../../docs/ANALYSIS_PLAN.md) on 18 September 2026. It is not part of the manuscript. This directory holds the code that implemented it:

- `R_cbh/snow/` — Snow exposure audit, prepared 2000–2015 sample, seven Snow age-band fits, smoothing/penalisation trials, matched MAP–Snow gamma=2 fits and the Snow-based national mortality scenarios, with its original README.
- `R_dhs/51_snow_polygon_prevalence.R` — population-weighted overlay of the Snow polygon estimates onto DHS survey regions.

Files were moved without editing their contents. Original paths, line counts and SHA-256 hashes are in the [archive manifest](../../results/cbh/code_archive_2026_09_18_snow/archive_manifest.csv). The archived code preserves its original project-root paths and is a reference snapshot, not a runnable pipeline from this directory. For historical reproduction, use an isolated checkout at commit `9b4eefd` (the last commit before archiving) with the corresponding authorized input snapshots. Do not copy archived scripts back over the active workflow.

Nothing else was moved or deleted:

- `results/cbh/map_snow_gamma2_v1/` is retained because its `map_full` rows are the **original seven primary MAP gamma=2 fits**, produced by `08_fit_map_comparison_gamma2.R`. `R_cbh/audit/01_verify_saved_primary.R` (run_all `--audit`) and the `legacy` version in `R_cbh/primary/settings.R` still read that directory. The directory name is historical.
- `results/cbh/age_band_snow_2000_2015_v1/` and the `results/dhs_rebuild/snow_*` files are retained as historical outputs. They are no longer linked from the plan or the key/supplementary results indexes.
- `data/snow_prevalence_model/`, `data/snow_region_cache/` and `data/snow_pfpr_by_survey_region*.csv` remain under the ignored `data/` directory.

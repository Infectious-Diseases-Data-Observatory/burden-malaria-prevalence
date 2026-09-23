# MICS surveys in the child age-band analysis

MICS microdata stay in the git-ignored `data/` tree (`data/MICS_Datasets`, `data/MICS_extracted`, `data/derived_mics`); only
code, reviewed configuration and aggregate results are tracked.

`Rscript run_all.R --mics-build` runs `00a_extract.R` and then `00`–`11` in order (`.py` files with `python3`) and stops at the
first failing script. No step uses the network.

> **Rerunning 08 invalidates the primary and every sensitivity fit.** `08_build_child_bands.R` rewrites
> `data/derived_mics/cbh_build/manifest.rds`, with new timestamps, even when every shard is reused from cache. The preparation
> signature of `primary_map_regional17_dhsmics_gamma2_v5` hashes that file, so the next primary run re-prepares the data and
> refits v5, and the subgroup (`subgroups_dhsmics_map_gamma2_v2`), no-nutrition (`nutrition_adjustment_dhsmics_map_gamma2_v2`)
> and imputed-covariate (`primary_map_regional17_dhsmics_imputed_gamma2_v6`) refits go stale. `run_all.R --mics-build` therefore
> skips 08 when the manifest records a completed build; add `--force` to run it, then run `run_all.R --primary` and
> `run_all.R --sensitivities`. The builder writes an incomplete manifest before it starts, so after an interrupted 08 the next
> `--mics-build` runs 08 again. A change by 09 to `data/derived_mics/regional_covariates_wide_mics.csv`, which the v5
> preparation also hashes, has the same effect as rerunning 08, as does any edit of `R_cbh/*.R` or `R_cbh/R/*.R`, which every
> shard signature hashes. 01, 06, 10 and 11 also always run, and rewrite files hashed by the v5 study flow and survey map
> (`survey_inventory.csv`, `survey_registry_mics.csv`, `boundaries/MC_*.rds`, `exclusion_attribution_mics_by_survey.csv`), by
> the subgroup fit signature (registry, `mics_region_centroids.csv`) or by the imputed-covariate provenance (registry,
> `map_pfpr_window_years_mics.csv`); if their bytes change, run at least `Rscript R_cbh/primary/run_regional.R --report-only`
> and then `run_all.R --sensitivities`. `--mics-build` ends by naming every one of these files, and the manifest and overlay,
> whose bytes changed. When 08 is skipped, the existing shards are kept even if 04–07 changed their inputs; use `--force` in
> that case. Do not run the build while a primary or sensitivity job is running.

**Configuration**

- `config/paths.R` defines the one input outside this repository: the Africa admin-1 shapefile from the Snow prevalence
  project. Set the environment variable `MICS_ADMIN1_SHP` to use a copy elsewhere. Sourced by 04, 05 and 06.
- `config/region_overrides.csv` holds the reviewed region-label translations, spellings and merges used by 05.

**Eligibility assessment**

1. `00a_extract.R` extracts the `.sav` members of each download, `data/MICS_Datasets/<survey>/**/*.zip`, into
   `data/MICS_extracted/<survey>/`, flat and with lower-case names (`bh.sav`, `wm.sav`, ...). A zip that wraps another zip
   (Botswana 2000) is opened one level down. A survey whose folder already holds `.sav` files is skipped, so an existing
   extraction is never overwritten. `--list` prints what would be extracted and writes nothing.
2. `00_dump_metadata.R` catalogues every variable name and label.
3. `01_inventory.R` writes `results/mics_inventory/survey_inventory.csv` (birth-history completeness, region, covariate sources).
4. `02_inventory_tables.R` and `03_eligibility_report.py` write `INVENTORY.md` and `ELIGIBILITY.md`.

The inventory covers 93 MICS surveys. 47 have complete birth histories; 46 are built, and Lesotho 2018 (`MC_LSO2018`) is
skipped because Lesotho is malaria free and MAP publishes no prevalence surface for it.

**Integration into the primary analysis (23 September 2026)**

5. `04_boundaries_gnb_caf.R` builds analysis-region polygons for Guinea-Bissau and the Central African Republic from the
   admin-1 shapefile. Capitals are merged with a neighbour.
6. `05_region_map.R` maps each MICS region label to an analysis region (union of boundary polygons), using exact names and the
   reviewed overrides in `config/region_overrides.csv`.
7. `06_polygons_pfpr_registry.R` builds per-survey polygons, extracts annual regional MAP PfPR with the DHS method, and writes the
   MICS registry, builder rules and boundary tables.
8. `07_make_recodes.R` converts each birth history to the DHS Births Recode layout (survival, twins, age at death, weights).
9. `08_build_child_bands.R` runs the unchanged builder into `data/derived_mics/cbh_build/`. The builder's region-override table,
   `data/derived_mics/region_overrides_mics.csv`, is header-only because regions are mapped in 05; 08 creates it only when it is
   missing, and never rewrites an existing copy, whose hash enters every shard signature. `run_all.R --mics-build --force`
   runs 08 without the builder's own `--force`, so shards whose signatures still match are reused;
   `Rscript R_mics/08_build_child_bands.R --force` rebuilds every shard.
10. `09_regional_covariates.R` computes the 13 regional covariates from the microdata, with the WUENIC and within-survey
    fallbacks. It uses every survey the builder marks `built` or `cached`.
11. `10_exclusion_attribution_mics.R` attributes MICS complete-case exclusions for the study flow.
12. `11_region_centroids.R` writes boundary centroids of the MICS analysis regions (same method as the DHS centroid cache), used
    by the subgroup refit.

**Prerequisites** (all local; no step downloads anything)

| Input | Location | Needed by |
|---|---|---|
| MICS downloads | `data/MICS_Datasets/<survey>/**/*.zip` | 00a |
| Raw MICS `.sav` files | `data/MICS_extracted/<survey>/` (written by 00a) | 00, 01, 05, 07, 09 |
| Africa admin-1 shapefile (Snow prevalence project, outside the repository) | `config/paths.R` (`MICS_ADMIN1_SHP`) | 04, 05, 06 |
| DHS boundary polygons | `data/dhs_boundaries/*.rds` | 05, 06 |
| Annual MAP PfPR₂₋₁₀ rasters | `data/map_annual/pfpr2_10_YYYY.tif` | 06 |
| GPW v4 2020 population density | `data/pop/gpw_v4_population_density_rev11_2020_2.5m.tif` | 06 |
| National annual panels and StatCompiler table from the DHS data stage | paths in `R_cbh/00_config.R` (`data/derived_dhs/`, `data/wb_*.csv`) | 08 |
| UNICEF data warehouse snapshot (WUENIC DTP3 and measles fallback) | `data/fusion_GLOBAL_DATAFLOW_UNICEF_1.0_all.csv` | 09 |
| Child HIV incidence panel | `data/derived_cbh/hiv_incidence/child_incidence_country_year.csv` | 10 |
| Hand-curated eligibility table | `results/mics_inventory/eligibility.csv` | 03 |

**Ordering**

- Each step reads earlier outputs: 00 → 01 (variable metadata); 01 → 02, 06, 09, the v5 study flow and survey map
  (`survey_inventory.csv`); 04 → 05, 06; 05 → 06, 07, 09 (`survey_setup.csv`, `region_map.csv`); 06 → 08, 09, 11, the v5
  study flow and survey map, and the subgroup and imputed-covariate refits (registry, rules, PfPR tables, boundaries);
  07 → 08 (recodes); 08 → 09, 10, the v5 preparation and study flow, and the imputed-covariate refit (manifest and build
  ledgers); 09 → 10, the v5 preparation and the imputed-covariate refit (covariate overlay).
- 10 must run before the v5 study flow (`R_cbh/reporting/01_study_flow.R`, a stage of `run_regional.R`), which reads
  `results/mics_inventory/exclusion_attribution_mics_by_survey.csv`.
- 11 must run before the subgroup refit (`R_cbh/sensitivity/subgroups/01_fit.R`), which reads
  `data/derived_mics/mics_region_centroids.csv`.
- 02 and 03 only write `INVENTORY.md` and `ELIGIBILITY.md`; nothing downstream reads them.

Then `Rscript run_all.R --primary` (`R_cbh/primary/run_regional.R`, default version `regional_mics`) prepares, fits and reports
`primary_map_regional17_dhsmics_gamma2_v5`. See `docs/ANALYSIS_PLAN.md` Sections 2.1–2.4 (MICS inputs, geography, records, covariates and inclusion). `Rscript run_all.R --sensitivities` runs
the sensitivity analyses on the same sample: `R_cbh/sensitivity/subgroups/` and `R_cbh/sensitivity/nutrition/` (both default to
DHS+MICS; `--dhs-only` reproduces the DHS-only refits) and `R_cbh/sensitivity/imputation_dhsmics/`.

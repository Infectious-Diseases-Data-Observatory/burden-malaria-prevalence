# MICS surveys in the child age-band analysis

MICS microdata stay in the git-ignored `data/` tree (`data/MICS_Datasets`, `data/MICS_extracted`, `data/derived_mics`); only
code, reviewed configuration and aggregate results are tracked.

**Eligibility assessment**
1. Extract the SPSS files from `data/MICS_Datasets/*.zip` into `data/MICS_extracted/<survey>/`.
2. `00_dump_metadata.R` catalogues every variable name and label.
3. `01_inventory.R` writes `results/mics_inventory/survey_inventory.csv` (birth-history completeness, region, covariate sources).
4. `02_inventory_tables.R` and `03_eligibility_report.py` write `INVENTORY.md` and `ELIGIBILITY.md`.

**Integration into the primary analysis (23 September 2026)**
5. `04_boundaries_gnb_caf.R` builds analysis-region polygons for Guinea-Bissau and the Central African Republic from the local
   Africa admin-1 shapefile in the Snow prevalence project (outside this repository). Capitals are merged with a neighbour.
6. `05_region_map.R` maps each MICS region label to an analysis region (union of boundary polygons), using exact names and the
   reviewed overrides in `config/region_overrides.csv`.
7. `06_polygons_pfpr_registry.R` builds per-survey polygons, extracts annual regional MAP PfPR with the DHS method, and writes the
   MICS registry, builder rules and boundary tables.
8. `07_make_recodes.R` converts each birth history to the DHS Births Recode layout (survival, twins, age at death, weights).
9. `08_build_child_bands.R` runs the unchanged builder into `data/derived_mics/cbh_build/`.
10. `09_regional_covariates.R` computes the 13 regional covariates from the microdata, with the WUENIC and within-survey fallbacks.
11. `10_exclusion_attribution_mics.R` attributes MICS complete-case exclusions for the study flow.

Then `Rscript R_cbh/primary/run_regional.R` (default version `regional_mics`) prepares, fits and reports
`primary_map_regional17_dhsmics_gamma2_v5`. See `docs/ANALYSIS_PLAN.md` section 2.6.

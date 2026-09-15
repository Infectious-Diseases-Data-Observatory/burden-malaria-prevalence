# Retained DHS input bridges and compatibility commands

The primary analysis is the seven-band MAP gamma=2 workflow in [R_cbh](../R_cbh/README.md), specified in the [analysis plan](../docs/ANALYSIS_PLAN.md). The former survey-region and six-band person-time modelling/plotting scripts were moved to [the archive](../archive/2026-09-15-code-audit/R_dhs/). The [original README](../archive/2026-09-15-code-audit/R_dhs/README.md) describes historical methods only.

| Retained code | Purpose / restriction |
|---|---|
| `00_config.R` | Shared extraction/geography helpers. Historical model helpers remain pending input-only separation; do not use its primary-model defaults for CBH. |
| `01_access_dhs_data.R` | Survey/local-file inventory and explicitly requested acquisition. |
| `02_extract_map_pfpr.R`, `41_extract_map_window_years.R` | Producers of the regional MAP snapshots. Survey-year selection, density weighting, year-specific missing coverage and cache provenance require review before a refreshed primary build. |
| `02b_extract_unicef_immunisation.R`, `02c_fetch_statcompiler_covariates.R` | Optional candidate covariate snapshots; vaccines, WASH and wasting are not primary covariates. |
| `03_build_analysis_dataset.R` | Mixed historical panel/input code, retained because it also prepares national covariate snapshots. It is not the primary dataset builder. |
| `51_snow_polygon_prevalence.R` | Supplementary Snow overlay. The specified analysis uses `SNOW_SOURCE=annual_csv SNOW_SKIP_COMPARISONS=1`; the default is historical five-year draws. |
| `11_study_flow.R`, `33_key_results.R` | Compatibility calls to current CBH inclusion-figure and result-index scripts. |
| `run_all.R` | Current command guidance / read-only audit dispatch; no automatic historical fitting. |

No input bridge is automatically run by the CBH dataset builder. Some bridges fetch data when a cache is missing; inspect their documented acquisition behaviour before running them. Current local processing starts with `Rscript R_cbh/01_make_analysis_data.R --check-inputs`.

See the [component audit](../docs/CODE_AUDIT.md) and [file-by-file inventory](../results/cbh/code_audit_2026_09_15/code_inventory.csv) for dispositions and unresolved dependencies.

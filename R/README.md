# Retained legacy input producers

Only four scripts remain here: `00_utils.R`, `01_fetch_data.R`, `03_component3_rdt_microscopy.R` and `11_malaria_deaths_timeseries.R`. They provide unique public-input, measured-parasitaemia/conversion or legacy national-MAP comparison inputs used by retained workflows. They are not current primary analysis commands, and mixed historical fitting branches still need to be separated from data preparation.

The other scripts are preserved in [the dated archive](../archive/2026-09-15-code-audit/R/). See the [audit](../docs/CODE_AUDIT.md) and [analysis plan](../docs/ANALYSIS_PLAN.md). The root runner no longer sources these scripts automatically.

#!/usr/bin/env Rscript
# Regenerate the table alone; no fitting, figure changes or manuscript writes.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/age_band_table.R")
cbh_primary_age_table(cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional_mics"))$out,cbh_config()$age_bands$age_band)
message("Updated age-band table (Markdown, CSV and LaTeX)")

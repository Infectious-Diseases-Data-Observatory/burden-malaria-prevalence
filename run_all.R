# run_all.R — reproduce the whole analysis. Run from the repository root:
#   Rscript run_all.R
# 01 fetches inputs into data/ (idempotent); 02-04 write figures/tables to results/.
# Component 3 runs before Component 2 (it builds the shared per-region prevalence
# table and the RDT->microscopy conversion that Component 2 consumes).

source("R/00_utils.R")             # config + helpers
source("R/01_fetch_data.R")        # inputs -> data/ (needs internet; DHS/IHME logins — see README)
source("R/02_component1_country.R")# Component 1: country-level share vs PfPR (+GDP,+DTP3, outliers)
source("R/03_component3_rdt_microscopy.R") # Component 3: RDT<->microscopy conversion (+ prevalence table)
source("R/04_component2_dhs.R")    # Component 2: DHS multivariable mixed model
message("\nAll components complete. See results/ for figures and tables.")

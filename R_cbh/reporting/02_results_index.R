#!/usr/bin/env Rscript
# Link authoritative current primary outputs; keep historical/supplementary labels.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
root <- cbh_primary_settings()$out
x <- cbh_read_csv(file.path(root,"fit_diagnostics.csv"))
stopifnot(nrow(x)==7L,all(x$series=="map_full"),all(x$gamma==2),all(x$converged),all(x$input_verified))
link <- function(label,path) sprintf("[%s](../%s/%s)",label,root,path)
dir.create("Key results",showWarnings=FALSE)
writeLines(c("# Key results — primary MAP analysis","",
  "The primary analysis uses seven separate age-band models, MAP prevalence, cr PfPR splines and gamma=2. The full sample contains 1,817,912 children, 5,885,022 child-band records and 82,415 deaths from 105 surveys in 34 countries.","",
  link("Primary analysis report and rerun verification","REPORT.md"),"",
  "| Result | Output |","|---|---|",
  paste("| Figure 1: survey map and timing |",link("Figure","survey_map/survey_map_and_timing.png"),"and",link("caption","survey_map/CAPTION.md"),"|"),
  paste("| Inclusion flow |",link("Figure","study_flow/study_flow_diagram.png"),"and",link("caption","study_flow/CAPTION.md"),"|"),
  paste("| Seven PfPR curves |",link("Figure","pfpr_splines.png"),"and",link("estimates","pfpr_curves.csv"),"|"),
  paste("| PfPR 40% to 20% effects |",link("Figure","pfpr_40_to_20.png"),"and",link("contrasts","pfpr_40_to_20_contrasts.csv"),"|"),
  paste("| Diagnostics |",link("Numerical checks","fit_diagnostics.csv"),"and",link("fitted outcomes","fitted_outcome_checks.csv"),"|"),
  paste("| Mortality, 2005/2015/2024 |",link("Country figure","burden/country_deaths_all_years.png"),"and",link("country totals","burden/country_totals.csv"),"|"),
  paste("| Comparison with IHME |",link("Scatter plots","burden/country_vs_ihme.png"),"|"),
  paste("| Age contributions |",link("Figure","burden/deaths_by_age.png"),"and",link("country-age estimates","burden/country_age_estimates.csv"),"|"),
  paste("| DRC synthetic-cohort survival |",link("Figure","burden/drc_survival.png"),"|"),"",
  link("Paper figure manifest and captions","paper_figures/CAPTIONS.md"),"",
  "Snow exposure comparisons are [supplementary](<../Supplementary results/README.md>). The previous MAP gamma=2 outputs remain in `results/cbh/map_snow_gamma2_v1/`, with a numerical comparison in the new report.","",
  "## Historical figures","",
  "Older PNGs in this directory are retained for traceability. The gamma=1 [geography/period/structure figure](pfpr_spline_sensitivity_by_age_geography_period.png) remains a key exploratory result; planned gamma=2 replacements have not been run. See the [original index](../archive/2026-09-15-code-audit/Key%20results/README.md).","",
  "Run `Rscript run_all.R --primary` to refit and regenerate the primary results from prepared inputs. Use `Rscript R_cbh/primary/run.R --report-only` for saved fits. The standalone reporting scripts and old `R_dhs/11_study_flow.R` and `R_dhs/33_key_results.R` wrappers use the same primary selection."
),"Key results/README.md")
message("Updated primary results index")

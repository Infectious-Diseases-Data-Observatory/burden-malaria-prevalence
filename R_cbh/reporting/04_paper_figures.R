#!/usr/bin/env Rscript
# Record paper figure sources/captions locally. Does not access or edit Overleaf.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
root <- cbh_primary_settings()$out
out <- file.path(root,"paper_figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
figures <- data.frame(figure=1:3,
  filename=c("survey_map_and_timing.png","fig2_pfpr_mortality_by_age.png","fig3_country_estimates_vs_ihme.png"),
  source=file.path(root,c("survey_map/survey_map_and_timing.png","pfpr_splines.png","burden/country_vs_ihme.png")))
stopifnot(all(file.exists(figures$source)))
figures$md5 <- vapply(figures$source,cbh_file_hash,"")
cbh_atomic_csv(figures,file.path(out,"manifest.csv"))
writeLines(c("# Main-paper figures and draft captions","",
  "Figures contain no title, subtitle or embedded caption. Axis labels, age/year facet labels and legends remain. File numbering here does not alter the order of figures in TeX.","",
  "## Figure 1 — survey map and timing","",
  "File: `survey_map_and_timing.png`.","",
  readLines(file.path(root,"survey_map/CAPTION.md"))[-c(1,2)],"",
  "## Figure 2 — PfPR and mortality by age","",
  "File: `fig2_pfpr_mortality_by_age.png`.","",
  "Age-specific associations between MAP PfPR[2–10] and all-cause mortality, estimated in seven separate age-band models at gamma=2. Curves show log mortality hazard ratios relative to PfPR=20% across the central 95% of each band's exposure distribution. Shading denotes conditional pointwise 95% intervals. Models use cr splines for PfPR (k=5) and calendar year (k=6), the prespecified adjustment set, and survey/country/region random effects. Intervals condition on fitted smoothing parameters, exposure and the fixed posterior-median child HIV imputation.","",
  "## Figure 3 — country mortality estimates versus IHME","",
  "File: `fig3_country_estimates_vs_ihme.png`.","",
  "Primary-model malaria-attributable deaths before age 5 versus IHME cause-specific malaria deaths, by country in 2005, 2015 and 2024. Both axes use a base-10 logarithmic scale with identical limits starting at 1,000 deaths; each point represents a country and the dashed line denotes equality. The underlying estimates include all 42 estimable countries per year, with no pseudocount; points below 1,000 deaths on either axis are outside the displayed range. The three countries without model estimates remain in the source tables. Model-attributable deaths equal national IHME all-cause deaths multiplied by 1 − exp[f_g(0) − f_g(P_country)], summed over age bands. National MAP prevalence is population weighted using the saved 2020 population weights. Ages 2–4 share the IHME mortality rate and divide deaths/person-time equally. Country totals are point estimates; cross-age covariance and source/design/imputation uncertainty are not propagated. Model-attributable all-cause reductions and IHME cause-specific mortality are different estimands.","",
  "[Source paths and hashes](manifest.csv). Copy only these three images into the manuscript's `figures/` folder. No TeX edits are needed to prepare or copy the images; TeX references/captions are managed separately."),file.path(out,"CAPTIONS.md"))

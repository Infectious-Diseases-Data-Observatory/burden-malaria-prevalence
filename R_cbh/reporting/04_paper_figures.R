#!/usr/bin/env Rscript
# Record paper figure sources/captions locally. Does not access or edit Overleaf.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
root <- cbh_primary_settings()$out
out <- file.path(root,"paper_figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
figures <- data.frame(figure=1:5,
  filename=c("survey_map_and_timing.png","fig2_pfpr_mortality_by_age.png",
    "fig3_malaria_attributable_fraction_by_age.png","fig4_country_estimates_vs_ihme.png",
    "fig5_annual_malaria_mortality.png"),
  source=file.path(root,c("survey_map/survey_map_and_timing.png","pfpr_splines.png",
    "attributable_fraction_by_age.png","burden/country_vs_ihme.png",
    "annual_comparison/fig5_annual_malaria_mortality.png")))
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
  paste0(cbh_paper_model_label(), ": age-specific associations between MAP PfPR[2–10] and all-cause mortality, estimated in seven separate age-band models at gamma=2. Curves show log mortality hazard ratios relative to PfPR=20% across the central 95% of each band's exposure distribution. Shading denotes conditional pointwise 95% intervals. Models use cr splines for PfPR (k=5) and calendar year (k=6), the prespecified adjustment set, and survey/country/region random effects. Intervals condition on fitted smoothing parameters, exposure and the fixed posterior-median child HIV imputation."),"",
  "## Figure 3 — malaria-attributable fraction by age","",
  "File: `fig3_malaria_attributable_fraction_by_age.png`.","",
  paste0(cbh_paper_model_label(), ": estimated share of all-cause deaths attributable to malaria within each completed-month age band, at MAP PfPR[2–10] levels of 10%, 20%, 30% and 40%. Points show 100 × {1 − exp[f_g(0) − f_g(P)]}, using the seven separate primary MAP gamma=2 models; lines connect discrete age-band estimates and do not imply a continuous-age fit. The denominator is all-cause deaths within each band, not all under-five deaths pooled across bands. The hazard-based attributable fraction follows the primary burden convention of holding annual person-time fixed; no country-specific IHME inputs are used here. Zero prevalence lies below observed exposure support, so the counterfactual requires extrapolation. Point estimates are plotted; conditional pointwise intervals and support flags are retained in attributable_fraction_by_age.csv."),"",
  "## Figure 4 — country mortality estimates versus IHME","",
  "File: `fig4_country_estimates_vs_ihme.png`.","",
  paste0(cbh_paper_model_label(), " malaria-attributable deaths before age 5 versus IHME cause-specific malaria deaths, by country in 2005, 2015 and 2024. Both axes use a base-10 logarithmic scale with identical limits starting at 1,000 deaths; each point represents a country and the dashed line denotes equality. The underlying estimates include all 42 estimable countries per year, with no pseudocount; points below 1,000 deaths on either axis are outside the displayed range. The three countries without model estimates remain in the source tables. Model-attributable deaths equal national IHME all-cause deaths multiplied by 1 − exp[f_g(0) − f_g(P_country)], summed over age bands. National MAP prevalence is population weighted using the saved 2020 population weights. Ages 2–4 share the IHME mortality rate and divide deaths/person-time equally. Country totals are point estimates; cross-age covariance and source/design/imputation uncertainty are not propagated. Model-attributable all-cause reductions and IHME cause-specific mortality are different estimands."),"",
  "## Figure 5 — annual under-five malaria mortality","",
  "File: `fig5_annual_malaria_mortality.png`.","",
  readLines(file.path(root,"annual_comparison/CAPTION.md"))[-c(1,2)],"",
  "[Source paths and hashes](manifest.csv). Copy these five images into the manuscript's `figures/` folder. The former fig3_country_estimates_vs_ihme.png is superseded by the Figure 4 filename; retain a legacy copy if existing TeX references still need it. No TeX edits are needed to prepare or copy the images; TeX references/captions are managed separately."),file.path(out,"CAPTIONS.md"))

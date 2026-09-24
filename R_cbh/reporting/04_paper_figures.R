#!/usr/bin/env Rscript
# Record paper figure sources/captions locally. Does not access or edit Overleaf.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
root <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional_mics"))$out
out <- file.path(root,"paper_figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Figure 3 combines the former Figures 4 (2024 panel), 6 and 5 as panels A, B and C.
# Decided 22 September 2026: the log-hazard-by-age curve figure moves to the
# supplement, so the attributable-fraction figure becomes Figure 2 and the
# burden comparison becomes Figure 3.
# Figure 4 (added 24 September 2026): probability of dying before age 5 by country, all causes
# and with malaria transmission removed (R_cbh/reporting/13_under5_death_probability.R).
figures <- data.frame(figure=1:4,
  filename=c("fig1_survey_map_and_timing.png","fig2_malaria_attributable_fraction_by_age.png",
    "fig3_burden_comparison.png","fig4_under5_death_probability.png"),
  source=file.path(root,c("survey_map/survey_map_and_timing.png",
    "attributable_fraction_by_age.png","burden_comparison/fig4_burden_comparison.png",
    "under5_probability/under5_death_probability_2024.png")))
stopifnot(all(file.exists(figures$source)))
figures$md5 <- vapply(figures$source,cbh_file_hash,"")
cbh_atomic_csv(figures,file.path(out,"manifest.csv"))
writeLines(c("# Main-paper figures and draft captions","",
  "Figures contain no title, subtitle or embedded caption. Axis labels, age/year facet labels and legends remain. File numbering here does not alter the order of figures in TeX.","",
  "## Figure 1 — survey map and timing","",
  "File: `fig1_survey_map_and_timing.png` (source `survey_map/survey_map_and_timing.png`).","",
  readLines(file.path(root,"survey_map/CAPTION.md"))[-c(1,2)],"",
  "## Figure 2 — malaria-attributable fraction by age","",
  "File: `fig2_malaria_attributable_fraction_by_age.png`.","",
  paste0(cbh_paper_model_label(), ": estimated share of all-cause deaths attributable to malaria within each completed-month age band, at MAP PfPR[2–10] levels of 10%, 20%, 30% and 40%. Points show 100 × {1 − exp[f_g(0) − f_g(P)]}, using the seven separate primary MAP gamma=2 models; lines connect discrete age-band estimates and do not imply a continuous-age fit. The denominator is all-cause deaths within each band, not all under-five deaths pooled across bands. The hazard-based attributable fraction follows the primary burden convention of holding annual person-time fixed; no country-specific IHME inputs are used here. Zero prevalence lies below observed exposure support, so the counterfactual requires extrapolation. Vertical bars show conditional pointwise 95% intervals, dodged by prevalence level; the interval limits and support flags are retained in attributable_fraction_by_age.csv."),"",
  "## Figure 3 — comparison with IHME and UN IGME as deaths per 1,000 child-years: countries (A), Nigerian states (B) and annual trend (C)","",
  "File: `fig3_burden_comparison.png`. Panel sources: `annual_comparison/country_estimates_2000_2024.csv` (2024 rates), `nigeria_states/state_totals_2024.csv`, `annual_comparison/figure5_data.csv`. The death-count version is the supplementary figure `burden_comparison/sfig_burden_comparison_counts.png` with its own caption (`CAPTION_counts.md`). The separate three-year country scatter (`burden/country_vs_ihme.png`), state figure (`nigeria_states/`) and annual figure (`annual_comparison/`) remain saved results but are no longer paper figures.","",
  readLines(file.path(root,"burden_comparison/CAPTION.md"))[-c(1,2)],"",
  "## Figure 4 — probability of dying before age 5, all causes and with malaria transmission removed","",
  "File: `fig4_under5_death_probability.png` (source `under5_probability/under5_death_probability_2024.png`; values in `under5_probability/under5_death_probability.csv`).","",
  readLines(file.path(root,"under5_probability/CAPTION.md"))[-c(1,2)],"",
  "## Supplementary figure — PfPR and mortality by age (formerly Figure 2)","",
  "File: `sfig_pfpr_mortality_by_age.png` (exported to the manuscript's `Supplementary Figures/` folder; source `pfpr_splines.png`).","",
  paste0(cbh_paper_model_label(), ": age-specific associations between MAP PfPR[2–10] and all-cause mortality, estimated in seven separate age-band models at gamma=2. Curves show log mortality hazard ratios relative to PfPR=20% across the central 95% of each band's exposure distribution. Shading denotes conditional pointwise 95% intervals. Models use cr splines for PfPR (k=5) and calendar year (k=6), 17 survey-region/national annual covariates, including maternal age at first birth, urban percentage, DTP3/measles, delivery, birth interval, WASH, electricity, wasting and stunting, and survey/country/region random effects. Intervals condition on fitted smoothing parameters, exposure, the fixed posterior-median child HIV imputation and filled regional covariates. Sex, multiple births, birth order, Hib3, PCV, rotavirus and exclusive breastfeeding are excluded."),"",
  "## Supplementary figure — count version of Figure 3","",
  "File: `sfig_burden_comparison_counts.png` (exported to the manuscript's `Supplementary Figures/` folder).","",
  readLines(file.path(root,"burden_comparison/CAPTION_counts.md"))[-c(1,2)],"",
  "[Source paths and hashes](manifest.csv). Copy the four main images into the manuscript's `figures/` folder and the supplementary images into `Supplementary Figures/`. Figure 3 is exported as `figures/fig3_burden_comparison.png` (the source file keeps its historical name). The former fig4_country_estimates_vs_ihme.png, fig5_annual_malaria_mortality.png and fig6_nigeria_states_vs_ihme.png are superseded by the combined Figure 3; existing TeX references to them need updating by the author. No TeX edits are needed to prepare or copy the images; TeX references/captions are managed separately."),file.path(out,"CAPTIONS.md"))

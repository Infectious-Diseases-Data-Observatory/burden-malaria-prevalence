#!/usr/bin/env Rscript
# Link authoritative current primary outputs; keep historical/supplementary labels.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
st <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional_mics"));root <- st$out;mics <- isTRUE(st$mics)
x <- cbh_read_csv(file.path(root,"fit_diagnostics.csv"))
sample <- cbh_read_csv(file.path(root,"primary_sample.csv"))
coverage <- cbh_read_csv(file.path(root,"survey_map/survey_coverage.csv"))
stopifnot(nrow(x)==7L,all(x$series=="map_full"),all(x$gamma==2),all(x$converged),all(x$input_verified),
  nrow(coverage)==sample$surveys,all(coverage$type %in% c("DHS","MICS")))
link <- function(label,path) sprintf("[%s](../%s/%s)",label,root,path)
# DHS and MICS sensitivity refits of this primary (the DHS-only versions are history).
sensitivity <- function(label,path) paste0(sprintf("[%s](../%s)",label,path),if(!file.exists(path)) " (not yet written)" else "")
surveys <- if(mics) sprintf("%s surveys (%d DHS and %d MICS)",sample$surveys,sum(coverage$type=="DHS"),sum(coverage$type=="MICS")) else
  sprintf("%s surveys",sample$surveys)
dir.create("Key results",showWarnings=FALSE)
writeLines(c("# Key results — primary MAP analysis","",
  sprintf("The primary analysis uses seven separate age-band MAP models, cr PfPR splines, gamma=2 and 17 regional/annual covariates. The sample contains %s children, %s child-band records and %s deaths from %s in %s countries.",format(sample$distinct_children,big.mark=","),format(sample$records,big.mark=","),format(sample$deaths,big.mark=","),surveys,sample$countries),"",
  link("Current primary figures and tables","RESULTS.md"),"",link("Comparison with the previous iteration","REPORT.md"),"",
  "| Result | Output |","|---|---|",
  paste("| Figure 1: survey map and timing |",link("Figure","survey_map/survey_map_and_timing.png"),"and",link("caption","survey_map/CAPTION.md"),"|"),
  paste("| Figure 2: malaria-attributable share by age |",link("Figure","attributable_fraction_by_age.png"),"and",link("estimates","attributable_fraction_by_age.csv"),"|"),
  paste("| Figure 3: countries (A), Nigerian states (B) and annual trend (C) versus IHME/UN IGME, deaths per 1,000 child-years |",link("Combined figure","burden_comparison/fig4_burden_comparison.png"),"and",link("caption","burden_comparison/CAPTION.md"),"|"),
  paste("| Figure 3A source: country comparison, 2005/2015/2024 |",link("Three-year scatter","burden/country_vs_ihme.png"),"|"),
  paste("| Figure 3B source: Nigerian states, 2024 |",link("Figure","nigeria_states/fig6_nigeria_states_vs_ihme.png"),"and",link("state table and audit","nigeria_states/README.md"),"|"),
  paste("| Figure 3C source: annual mortality, 2000–2024 |",link("Figure","annual_comparison/fig5_annual_malaria_mortality.png"),"and",link("totals and input audit","annual_comparison/README.md"),"|"),
  paste("| Figure 4: probability of dying before age 5 by country, all causes and caused by malaria |",link("Figure","under5_probability/under5_death_probability_2024.png"),"and",link("values","under5_probability/under5_death_probability.csv"),"and",link("caption","under5_probability/CAPTION.md"),"|"),
  paste("| Supplementary: inclusion flow |",link("Figure","study_flow/study_flow_diagram.png"),"and",link("caption","study_flow/CAPTION.md"),"|"),
  paste("| Supplementary: seven PfPR curves |",link("Figure","pfpr_splines.png"),"and",link("estimates","pfpr_curves.csv"),"|"),
  paste("| Supplementary: count version of Figure 3 |",link("Figure","burden_comparison/sfig_burden_comparison_counts.png"),"and",link("caption","burden_comparison/CAPTION_counts.md"),"|"),
  if(mics) c(
    paste("| Supplementary: subgroup refits |",sensitivity("Report",file.path(st$sensitivity_dirs[["subgroups"]],"REPORT.md")),"|"),
    paste("| Supplementary: refit without nutrition covariates |",sensitivity("Contrasts",file.path(st$sensitivity_dirs[["nutrition"]],"CONTRASTS.md")),"|"),
    paste("| Supplementary: imputed covariates |",sensitivity("Report",if(isTRUE(st$liberia)) "results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v8/REPORT.md" else "results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/REPORT.md"),"|")),
  paste("| Age-band results table |",link("Table","tables/age_band_results.md"),"and",link("LaTeX source as text","tables/age_band_results.latex.txt"),"|"),
  paste("| Country, annual and state tables |",link("Table index","tables/README.md"),"|"),
  paste("| PfPR 40% to 20% effects |",link("Figure","pfpr_40_to_20.png"),"and",link("contrasts","pfpr_40_to_20_contrasts.csv"),"|"),
  paste("| Diagnostics |",link("Numerical checks","fit_diagnostics.csv"),"and",link("fitted outcomes","fitted_outcome_checks.csv"),"|"),
  paste("| Mortality, 2005/2015/2024 |",link("Country figure","burden/country_deaths_all_years.png"),"and",link("country totals","burden/country_totals.csv"),"|"),
  paste("| Age contributions |",link("Figure","burden/deaths_by_age.png"),"and",link("country-age estimates","burden/country_age_estimates.csv"),"|"),
  paste("| DRC synthetic-cohort survival |",link("Figure","burden/drc_survival.png"),"|"),"",
  link("Paper figure manifest and captions","paper_figures/CAPTIONS.md"),"",
  paste0("Sensitivity analyses are [supplementary](<../Supplementary results/README.md>). The comparator primary remains in `",st$reference,"/`, with a comparison in the new report",
    if(isTRUE(st$liberia)) "; this sample adds Liberia's 2007, 2013 and 2019-20 DHS to it (child HIV incidence derived from UNAIDS counts)." else if(mics) "; the DHS part of this sample reproduces it exactly." else "."),"",
  "## Historical figures","",
  "Older PNGs in this directory are retained for traceability. The gamma=1 [geography/period/structure figure](pfpr_spline_sensitivity_by_age_geography_period.png) remains a key exploratory result; planned gamma=2 replacements have not been run. See the [original index](../archive/2026-09-15-code-audit/Key%20results/README.md).","",
  paste0("Run `Rscript run_all.R --primary` to refit and regenerate the primary results from prepared inputs. Primary figures and tables follow the current analysis plan. Use `Rscript R_cbh/primary/run_regional.R --report-only` for saved fits. The standalone reporting scripts and old `R_dhs/11_study_flow.R` and `R_dhs/33_key_results.R` wrappers use the same primary selection.",
    if(mics) " The sensitivity refits are run separately from `R_cbh/sensitivity/subgroups/`, `R_cbh/sensitivity/nutrition/` and `R_cbh/sensitivity/imputation_dhsmics/`." else "")
),"Key results/README.md")
message("Updated primary results index")

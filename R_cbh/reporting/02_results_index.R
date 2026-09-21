#!/usr/bin/env Rscript
# Link authoritative current primary outputs; keep historical/supplementary labels.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
root <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional"))$out
x <- cbh_read_csv(file.path(root,"fit_diagnostics.csv"))
sample <- cbh_read_csv(file.path(root,"primary_sample.csv"))
stopifnot(nrow(x)==7L,all(x$series=="map_full"),all(x$gamma==2),all(x$converged),all(x$input_verified))
link <- function(label,path) sprintf("[%s](../%s/%s)",label,root,path)
dir.create("Key results",showWarnings=FALSE)
writeLines(c("# Key results — primary MAP analysis","",
  sprintf("The primary analysis uses seven separate age-band MAP models, cr PfPR splines, gamma=2 and 17 regional/annual covariates. The sample contains %s children, %s child-band records and %s deaths from %s surveys in %s countries.",format(sample$distinct_children,big.mark=","),format(sample$records,big.mark=","),format(sample$deaths,big.mark=","),sample$surveys,sample$countries),"",
  link("Current primary figures and tables","RESULTS.md"),"",link("Comparison with the previous iteration","REPORT.md"),"",
  "| Result | Output |","|---|---|",
  paste("| Figure 1: survey map and timing |",link("Figure","survey_map/survey_map_and_timing.png"),"and",link("caption","survey_map/CAPTION.md"),"|"),
  paste("| Inclusion flow |",link("Figure","study_flow/study_flow_diagram.png"),"and",link("caption","study_flow/CAPTION.md"),"|"),
  paste("| Figure 2: seven PfPR curves |",link("Figure","pfpr_splines.png"),"and",link("estimates","pfpr_curves.csv"),"|"),
  paste("| Figure 3: malaria-attributable share by age |",link("Figure","attributable_fraction_by_age.png"),"and",link("estimates","attributable_fraction_by_age.csv"),"|"),
  paste("| Age-band results table |",link("Table","tables/age_band_results.md"),"and",link("LaTeX source as text","tables/age_band_results.latex.txt"),"|"),
  paste("| Country, annual and state tables |",link("Table index","tables/README.md"),"|"),
  paste("| PfPR 40% to 20% effects |",link("Figure","pfpr_40_to_20.png"),"and",link("contrasts","pfpr_40_to_20_contrasts.csv"),"|"),
  paste("| Diagnostics |",link("Numerical checks","fit_diagnostics.csv"),"and",link("fitted outcomes","fitted_outcome_checks.csv"),"|"),
  paste("| Mortality, 2005/2015/2024 |",link("Country figure","burden/country_deaths_all_years.png"),"and",link("country totals","burden/country_totals.csv"),"|"),
  paste("| Figure 4: countries (A), Nigerian states (B) and annual trend (C) versus IHME/UN IGME, deaths per 1,000 child-years |",link("Combined figure","burden_comparison/fig4_burden_comparison.png"),"and",link("caption","burden_comparison/CAPTION.md"),"|"),
  paste("| Supplementary: count version of Figure 4 |",link("Figure","burden_comparison/sfig_burden_comparison_counts.png"),"and",link("caption","burden_comparison/CAPTION_counts.md"),"|"),
  paste("| Figure 4A source: country comparison, 2005/2015/2024 |",link("Three-year scatter","burden/country_vs_ihme.png"),"|"),
  paste("| Figure 4B source: Nigerian states, 2024 |",link("Figure","nigeria_states/fig6_nigeria_states_vs_ihme.png"),"and",link("state table and audit","nigeria_states/README.md"),"|"),
  paste("| Figure 4C source: annual mortality, 2000–2024 |",link("Figure","annual_comparison/fig5_annual_malaria_mortality.png"),"and",link("totals and input audit","annual_comparison/README.md"),"|"),
  paste("| Age contributions |",link("Figure","burden/deaths_by_age.png"),"and",link("country-age estimates","burden/country_age_estimates.csv"),"|"),
  paste("| DRC synthetic-cohort survival |",link("Figure","burden/drc_survival.png"),"|"),"",
  link("Paper figure manifest and captions","paper_figures/CAPTIONS.md"),"",
  "Sensitivity analyses are [supplementary](<../Supplementary results/README.md>). The previous MAP gamma=2 outputs remain in `results/cbh/primary_map_regional18_gamma2_v2/`, with a comparison in the new report.","",
  "## Historical figures","",
  "Older PNGs in this directory are retained for traceability. The gamma=1 [geography/period/structure figure](pfpr_spline_sensitivity_by_age_geography_period.png) remains a key exploratory result; planned gamma=2 replacements have not been run. See the [original index](../archive/2026-09-15-code-audit/Key%20results/README.md).","",
  "Run `Rscript run_all.R --primary` to refit and regenerate the primary results from prepared inputs. Primary figures and tables follow the current analysis plan; historical sensitivity plots are not relabelled. Use `Rscript R_cbh/primary/run_regional.R --report-only` for saved fits. The standalone reporting scripts and old `R_dhs/11_study_flow.R` and `R_dhs/33_key_results.R` wrappers use the same primary selection."
),"Key results/README.md")
message("Updated primary results index")

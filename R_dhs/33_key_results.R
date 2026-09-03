# =============================================================================
# 33_key_results.R — gather the key sensitivity-analysis figures in one folder.
#
# Copies a curated set of figures from results/dhs_rebuild into "Key results/"
# under descriptive names and writes a README there saying what each one shows
# and which script draws it. Nothing is computed here, so the producing scripts
# must have run first; run_all.R orders them that way.
#
# Output
#   Key results/*.png
#   Key results/README.md
# =============================================================================
source("R_dhs/00_config.R")

KEY_RESULTS_DIR <- file.path(REPO_ROOT, "Key results")
dir.create(KEY_RESULTS_DIR, showWarnings = FALSE)

figure <- function(source, target, script, shows) {
  data.frame(source = source, target = target, script = script, shows = shows,
             stringsAsFactors = FALSE)
}
key_figures <- rbind(
  figure("study_flow_diagram.png",
         "study_flow_diagram.png",
         "11_study_flow.R",
         paste("Data and analysis flow: DHS/MIS recodes, MAP prevalence and the",
               "national covariate series feed the survey-region panel, which is",
               "filtered to the shared-outcome sample and analysed with the",
               "ridge-penalised model, the Bayesian refits, the sensitivity",
               "analyses and the national burden extrapolation. Every count is",
               "read from the pipeline outputs.")),
  figure("figure15_brms_subgroup_curves.png",
         "subgroup_curves_postneonatal.png",
         "31_brms_subgroup_fits.R",
         paste("Dose-response curves for post-neonatal mortality: the",
               "malaria-attributable fraction against MAP PfPR2-10 by era",
               "(before / from 2013), region (West / Central & East) and the",
               "2 x 2 cells. Bayesian (brms) fits with 95% credible intervals;",
               "each curve covers only the prevalence range its own subgroup",
               "observes.")),
  figure("figure13_brms_subgroup_af.png",
         "subgroup_af_10_30_50_postneonatal.png",
         "31_brms_subgroup_fits.R",
         paste("Post-neonatal attributable fraction at 10, 30 and 50%",
               "prevalence in each subgroup, with 95% credible intervals and",
               "the number of surveys and survey regions behind each",
               "estimate.")),
  figure("figure15_brms_subgroup_curves_neonatal.png",
         "subgroup_curves_neonatal.png",
         "31_brms_subgroup_fits.R with BRMS_OUTCOME=neonatal",
         paste("The same dose-response curves for neonatal mortality, the",
               "negative control. Not clipped at zero, so a null effect can",
               "sit on the zero line.")),
  figure("figure13_brms_subgroup_af_neonatal.png",
         "subgroup_af_10_30_50_neonatal.png",
         "31_brms_subgroup_fits.R with BRMS_OUTCOME=neonatal",
         paste("The same 10, 30 and 50% anchors for neonatal mortality, the",
               "negative control.")),
  figure("figure10_survey_map.png",
         "survey_map_and_timing.png",
         "28_survey_map.R",
         paste("Countries contributing to the panel and when each DHS/MIS",
               "survey was fielded; point size gives the number of survey",
               "regions.")),
  figure(sprintf("figure8b_horizon_%d_vs_60.png", CHMORT_PERIOD),
         sprintf("mortality_%d_vs_60_month_window.png", CHMORT_PERIOD),
         "27_horizon_lag_selection.R",
         sprintf(paste(
           "Region-level mortality estimated over the %d-month window (the",
           "primary analysis) against the 60-month window, for post-neonatal",
           "and neonatal mortality, with 95%% delete-one-cluster jackknife",
           "intervals."), CHMORT_PERIOD)),
  figure("map_vs_measured_scatter.png",
         "map_vs_measured_prevalence.png",
         "21_map_vs_measured_prevalence.R",
         paste("MAP modelled PfPR2-10 against parasitaemia measured directly",
               "in the DHS/MIS surveys that tested children 6-59 months (RDT,",
               "or microscopy where available, converted to a",
               "microscopy-equivalent and age-standardised to 2-10 years), by",
               "era of fieldwork.")),
  figure("figure29_person_time_data_and_curves_4m.png",
         "person_time_dose_response_by_age_band.png",
         "45_person_time_data_and_curves.R",
         paste("Person-time design (deaths and child-years by region, 12-month",
               "window and age segment, each window paired with the MAP",
               "prevalence of its own year). Top: observed death rate per 1,000",
               "child-years against prevalence, one point per survey region.",
               "Bottom: fitted dose-response per age band (<1, 1-3, 4-11, 12-23,",
               "24-35, 36-59 months) from one negative-binomial model per band.")),
  figure("figure30_person_time_covariate_forest.png",
         "person_time_covariate_forest.png",
         "46_person_time_covariate_forest.R",
         paste("Ridge-shrunk covariate effects (% change in hazard per SD) by",
               "age band under the person-time model; national covariates are",
               "largely absorbed by the survey intercept.")),
  figure("figure31_person_time_joint_vs_separate.png",
         "person_time_joint_vs_separate.png",
         "47_person_time_joint_model.R",
         paste("Attributable fraction of all-cause mortality against 0% prevalence",
               "at 10, 30 and 50% by age band, from the joint model with",
               "band-specific nuisance terms against the separate per-band fits;",
               "the AIC gap between band-specific and common dose-response",
               "in the subtitle.")),
  figure("figure18_brms_ladder_kfold.png",
         "model_structure_kfold_by_survey.png",
         "35_brms_ladder_kfold.R",
         paste("Survey-grouped 10-fold cross-validation of how prevalence and",
               "calendar time enter the post-neonatal model: additive",
               "s(pfpr10) + s(year_c), the curve shifting linearly in time, and",
               "the full t2 surface. One point per held-out survey, showing",
               "its log predictive density relative to the additive model;",
               "totals, standard errors and stacking weights in the subtitle."))
)

key_figures$source_path <- file.path(RESULTS_DIR, key_figures$source)
missing <- key_figures$source[!file.exists(key_figures$source_path)]
if (length(missing)) {
  stop("Missing figures, run their scripts first: ",
       paste(missing, collapse = ", "))
}
copied <- file.copy(key_figures$source_path,
                    file.path(KEY_RESULTS_DIR, key_figures$target),
                    overwrite = TRUE)
if (!all(copied)) stop("Copy failed for: ",
                       paste(key_figures$target[!copied], collapse = ", "))
key_figures$generated <- format(file.info(key_figures$source_path)$mtime,
                                "%Y-%m-%d %H:%M")

readme <- c(
  "# Key results",
  "",
  paste("Curated figures from the DHS/MIS malaria prevalence and child",
        "mortality analysis. Every file here is a copy of a figure in",
        "`results/dhs_rebuild/`; the source script and the time the figure was",
        "generated are given for each. Regenerate with",
        "`Rscript R_dhs/33_key_results.R` after re-running the producing",
        "scripts (or `Rscript R_dhs/run_all.R`, which does both)."),
  "",
  sprintf(paste("The primary analysis uses the %d-month mortality window and",
                "MAP prevalence in the survey year. Post-neonatal mortality is",
                "the primary outcome; neonatal mortality is the negative",
                "control."), CHMORT_PERIOD),
  "",
  "| File | What it shows | Source | Generated |",
  "|---|---|---|---|",
  sprintf("| `%s` | %s | `%s` (`%s`) | %s |",
          key_figures$target, key_figures$shows, key_figures$script,
          key_figures$source, key_figures$generated)
)
writeLines(readme, file.path(KEY_RESULTS_DIR, "README.md"))
message("Copied ", nrow(key_figures), " figures into ", KEY_RESULTS_DIR,
        " and wrote README.md")

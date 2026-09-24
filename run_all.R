#!/usr/bin/env Rscript
# Explicit routing only: a bare invocation prints the help below and runs nothing (no
# downloads, no long refit). Each mode runs its scripts in order and stops at the first
# failing script, naming it. See docs/ANALYSIS_PLAN.md and R_mics/README.md.
args <- commandArgs(trailingOnly = TRUE)
help <- c(
  "Current primary (23 September 2026): DHS and UNICEF MICS surveys, seven separate MAP age-band models, gamma = 2,",
  "17 regional/annual covariates -> results/cbh/primary_map_regional17_dhsmics_gamma2_v5/",
  "  7,498,459 child-band records, 2,292,089 children, 102,282 deaths, 1,211 survey-regions, 132 surveys",
  "  (95 DHS + 37 MICS; none of the 5 MIS surveys passes complete-case selection), 36 countries.",
  "  The DHS part reproduces the DHS-only v3 (primary_map_regional17_gamma2_v3) exactly.",
  "",
  "Reproduction order (explicit modes; a bare call runs nothing; no mode downloads data):",
  "  0. Rscript run_all.R --check-inputs    local DHS dataset prerequisites (R_cbh/01_make_analysis_data.R --check-inputs)",
  "  1. Rscript run_all.R --mics-build      R_mics/00a_extract.R, then R_mics 00-11 in order (.py via python3).",
  "       Needs the MICS downloads and the inputs listed in R_mics/README.md.",
  "       Skips 08_build_child_bands.R when data/derived_mics/cbh_build/manifest.rds records a completed build; add",
  "       --force to run it (shards whose signatures still match are reused). An interrupted 08 leaves an incomplete",
  "       manifest, and the next --mics-build runs 08 again.",
  "       Rerunning 08 rewrites that manifest, which invalidates v5 and every sensitivity fit: steps 2 and 3 must",
  "       then follow. Do not run the build while a primary or sensitivity job is running.",
  "       09 always rewrites the MICS covariate overlay, which the v5 preparation hashes; 01, 06, 10 and 11 rewrite",
  "       files that the v5 study flow and survey map or the sensitivity fit caches hash. The run ends by naming any",
  "       of these whose bytes changed: a changed manifest or overlay needs steps 2 and 3; any other change needs at",
  "       least run_regional.R --report-only, then step 3.",
  "  2. Rscript run_all.R --primary         v5 via R_cbh/primary/run_regional.R: prepare, fresh fits (about 16 minutes),",
  "       effects, diagnostics, comparison, study flow, survey map, burden, figures, tables, results index, validation.",
  "       Saved fits only:        Rscript R_cbh/primary/run_regional.R --report-only",
  "       DHS-only history (v3):  Rscript R_cbh/primary/run_regional.R --dhs-only",
  "  3. Rscript run_all.R --sensitivities   refits on the v5 sample; fit caches are reused only when their signatures match:",
  "       subgroups (R_cbh/sensitivity/subgroups/)          -> results/cbh/subgroups_dhsmics_map_gamma2_v2/",
  "       no nutrition covariates (sensitivity/nutrition/)   -> results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v2/",
  "       imputed covariates (sensitivity/imputation_dhsmics/) -> results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/",
  "       SMC before/after (sensitivity/smc/; needs data/SMC_rollout/) -> results/cbh/smc_dhsmics_map_gamma2_v1/",
  "       The imputed-covariate fits and multiple-imputation refits take hours; run detached.",
  "  4. python3 R_cbh/reporting/export_paper.py --destination <Overleaf folder> [--copy]",
  "",
  "Upstream DHS data stage (not run by any mode): R_cbh/01_make_analysis_data.R, R_cbh/hiv/01_fit_incidence.R",
  "  [--extended], R_cbh/covariates/run.R, then 11, 12_fetch_nutrition.py (python3; queries the DHS API, reuses cached",
  "  pages), 13 and 14.",
  "Sahel mortality by calendar month (supplementary): R_cbh/seasonality/01_build_cells.R, then 02_fit_and_plot.R.",
  "",
  "Rscript run_all.R --audit   legacy: re-verifies the saved 15 September 2026 primary (results/cbh/map_snow_gamma2_v1),",
  "                            not the current primary.",
  "Old runners: archive/2026-09-15-code-audit/ (historical use only). Code audit: docs/CODE_AUDIT.md"
)
mics <- file.path("R_mics", c("00a_extract.R", "00_dump_metadata.R", "01_inventory.R", "02_inventory_tables.R",
  "03_eligibility_report.py", "04_boundaries_gnb_caf.R", "05_region_map.R", "06_polygons_pfpr_registry.R",
  "07_make_recodes.R", "08_build_child_bands.R", "09_regional_covariates.R", "10_exclusion_attribution_mics.R",
  "11_region_centroids.R"))
sensitivities <- file.path("R_cbh/sensitivity", c("subgroups/01_fit.R", "subgroups/02_report.R",
  "nutrition/01_fit.R", "nutrition/02_report.R", file.path("imputation_dhsmics", c("01_impute_national.R",
  "02_impute_regional.R", "03_prepare.R", "04_fit.R", "05_propagate.R", "06_report.R")),
  file.path("smc", c("01_assign_smc.R", "02_fit.R", "03_report.R", "04_maps.R"))))
modes <- list("--check-inputs" = "R_cbh/01_make_analysis_data.R", "--mics-build" = mics,
  "--primary" = "R_cbh/primary/run_regional.R", "--sensitivities" = sensitivities,
  "--audit" = "R_cbh/audit/01_verify_saved_primary.R")
if (!length(args) || identical(args, "--help")) {
  cat(help, sep = "\n")
} else {
  force <- "--force" %in% args; mode <- setdiff(args, "--force")
  if (length(mode) != 1L || !mode %in% names(modes) || (force && mode != "--mics-build"))
    stop("Use --help, --check-inputs, --mics-build [--force], --primary, --sensitivities or --audit.", call. = FALSE)
  if (!file.exists("R_cbh/primary/run_regional.R")) stop("Run from the project root.", call. = FALSE)
  mics_manifest <- "data/derived_mics/cbh_build/manifest.rds"
  # Files a --mics-build run can rewrite that v5 or its sensitivity caches hash. The first two (08's manifest, 09's
  # covariate overlay) enter the v5 preparation; the rest enter the v5 study flow and survey map (01, 06, 10) or the
  # subgroup and imputed-covariate signatures (06, 11).
  refit_inputs <- c(mics_manifest, "data/derived_mics/regional_covariates_wide_mics.csv")
  watched <- function() tools::md5sum(c(refit_inputs, "results/mics_inventory/survey_inventory.csv",
    "data/derived_mics/survey_registry_mics.csv", "data/derived_mics/map_pfpr_window_years_mics.csv",
    list.files("data/derived_mics/boundaries", "^MC_.*[.]rds$", full.names = TRUE),
    "results/mics_inventory/exclusion_attribution_mics_by_survey.csv", "data/derived_mics/mics_region_centroids.csv"))
  before <- if (mode == "--mics-build") watched()
  for (script in modes[[mode]]) {
    # Skip only a completed build: cbh_build writes an incomplete manifest before it starts.
    if (script == "R_mics/08_build_child_bands.R" && !force && file.exists(mics_manifest) &&
        isTRUE(tryCatch(readRDS(mics_manifest)$complete, error = function(e) FALSE))) {
      message("run_all ", mode, ": skipping ", script, " because ", mics_manifest, " records a completed build. ",
        "Rerunning 08 rewrites that manifest, which invalidates the v5 primary and every sensitivity fit. Add --force ",
        "to run it, then run --primary and --sensitivities.")
      next
    }
    cmd <- if (grepl("[.]py$", script)) Sys.which("python3") else file.path(R.home("bin"), "Rscript")
    if (!nzchar(cmd)) stop("python3 not found; it runs ", script, call. = FALSE)
    trailing <- if (mode == "--check-inputs") "--check-inputs" else character()
    message("run_all ", mode, ": ", script)
    status <- system2(cmd, c(script, trailing))
    if (status != 0L) stop("run_all ", mode, " stopped: ", script, " failed with status ", status, call. = FALSE)
  }
  if (mode == "--mics-build") {
    after <- watched(); keys <- union(names(before), names(after))
    changed <- keys[vapply(keys, function(k) !identical(unname(before[k]), unname(after[k])), NA)]
    if (length(changed)) message("run_all --mics-build: changed ", paste(changed, collapse = ", "), ". ",
      if (any(changed %in% refit_inputs)) paste("The v5 preparation hashes the MICS manifest and covariate overlay:",
        "run --primary and then --sensitivities.")
      else paste("The v5 study flow and survey map or the sensitivity fit caches hash these: run at least",
        "Rscript R_cbh/primary/run_regional.R --report-only, then --sensitivities."))
  }
  message("run_all ", mode, ": complete")
}

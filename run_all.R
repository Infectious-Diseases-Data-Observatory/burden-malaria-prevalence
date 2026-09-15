#!/usr/bin/env Rscript
# Explicit routing only: a bare invocation must not launch superseded models,
# downloads or a long refit. See docs/ANALYSIS_PLAN.md and docs/CODE_AUDIT.md.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || identical(args, "--help")) {
  cat(paste(c(
    "Current primary: seven separate MAP age-band models, gamma=2.",
    "Primary analysis without data setup: Rscript run_all.R --primary",
    "  Forces fresh primary fits, then regenerates effects, diagnostics and figures.",
    "Figures/aggregates from saved primary fits: Rscript R_cbh/primary/run.R --report-only",
    "",
    "Local dataset build: Rscript R_cbh/01_make_analysis_data.R",
    "HIV incidence: Rscript R_cbh/hiv/01_fit_incidence.R (fits the imputation model)",
    "Modelling data only: Rscript R_cbh/analysis/01_fit_complete_case.R --prepare-only",
    "Existing primary fits: results/cbh/primary_map_gamma2_v1/fit_manifest.csv",
    "Current refit implementation: R_cbh/primary/01_fit.R",
    "  Uses the saved prepared sample and fixed MAP knots; no supplementary fitting.",
    "Inclusion figure: Rscript R_cbh/reporting/01_study_flow.R",
    "Results index: Rscript R_cbh/reporting/02_results_index.R",
    "",
    "Rscript run_all.R --check-inputs   checks local dataset prerequisites",
    "Rscript run_all.R --audit          verifies the preserved pre-rerun primary snapshot",
    "",
    "Old runners: archive/2026-09-15-code-audit/ (historical use only).",
    "Full audit and outstanding planned analyses: docs/CODE_AUDIT.md"
  ), collapse = "\n"), "\n")
} else {
  scripts <- c("--primary" = "R_cbh/primary/run.R", "--check-inputs" = "R_cbh/01_make_analysis_data.R", "--audit" = "R_cbh/audit/01_verify_saved_primary.R")
  if (length(args) != 1L || !args %in% names(scripts))
    stop("Use --help, --primary, --check-inputs or --audit; legacy run-all modes have been archived.")
  trailing <- if (args == "--check-inputs") "--check-inputs" else character()
  status <- system2(file.path(R.home("bin"), "Rscript"), c(scripts[[args]], trailing))
  if (status != 0L) stop("Stage failed with status ", status)
}

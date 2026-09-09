#!/usr/bin/env Rscript
# All fits use the same saved median-imputed HIV dataset. No microdata exports.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/sensitivity/model.R")
library(mgcv)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--force", "--include-south"))) stop("Arguments: --force --include-south")
include_south <- "--include-south" %in% args
spec <- cbh_trial_spec()
base_dir <- file.path("data/derived_cbh/models", spec$id)
suffix <- if (include_south) "sensitivity_single_imputation_with_south" else "sensitivity_single_imputation"
private <- file.path(base_dir, suffix)
out <- file.path("results/cbh", spec$id, suffix)
dir.create(private, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
prepared <- readRDS(file.path(base_dir, "complete_case_dataset.rds"))
input_signature <- prepared$signature
scaling <- prepared$scaling
surveys <- prepared$input_manifest$manifest
surveys <- surveys[surveys$survey %in% levels(prepared$data$survey), c("survey", "country", "survey_year")]
cbh_unique(surveys, "survey", "Included survey metadata")
stopifnot(!anyNA(surveys$survey_year), !isTRUE(spec$time_by_age))
cutoff <- median(surveys$survey_year)
surveys$period <- ifelse(surveys$survey_year <= cutoff, "early", "late")
west <- c("BEN", "BFA", "CIV", "GHA", "GIN", "GMB", "MLI", "MRT", "NER", "NGA", "SEN", "SLE", "TGO")
east <- c("BDI", "COM", "ETH", "KEN", "MDG", "MOZ", "MWI", "RWA", "TZA", "UGA", "ZMB", "ZWE")
central <- c("AGO", "CMR", "COD", "COG", "GAB", "TCD")
south <- c("NAM", "SWZ", "ZAF")
stopifnot(setequal(unique(surveys$country), c(west, east, central, south)))
surveys$un_region <- ifelse(surveys$country %in% west, "Western Africa",
  ifelse(surveys$country %in% east, "Eastern Africa",
    ifelse(surveys$country %in% central, "Middle Africa", "Southern Africa")))
surveys$geography <- ifelse(surveys$country %in% west, "west",
  ifelse(surveys$country %in% c(east, central, if (include_south) south), "east_central", "excluded_south"))
vars <- unique(c(all.vars(cbh_trial_formula(spec)), "country"))
d <- prepared$data[vars]
rm(prepared); gc(FALSE)
age <- levels(d$age_band)
survey_index <- match(as.character(d$survey), surveys$survey)
stopifnot(!anyNA(survey_index))
cbh_atomic_csv(surveys, file.path(out, "survey_groups.csv"))
cbh_atomic_csv(scaling, file.path(out, "scaling.csv"))
message("Median included survey year: ", cutoff, "; early <= median, late > median.")
jobs <- data.frame(fit_id = c("reference", paste0("age_", seq_along(age)), "west", "east_central", "early", "late"),
  type = c("reference", rep("single_age", length(age)), rep("joint_subset", 4)),
  age_band = c(NA, age, rep(NA, 4)))
code_files <- c("R_cbh/sensitivity/model.R", "R_cbh/sensitivity/01_fit.R", "R_cbh/analysis/model.R")
code_hash <- vapply(code_files, cbh_file_hash, character(1))
signature <- cbh_hash(list(input_signature, spec, surveys, code_hash, R.version.string, as.character(packageVersion("mgcv"))))
writeLines(c("Single saved posterior-median HIV imputation; no imputation-draw refits.",
  paste("Median survey year:", cutoff, "; early <= median; late > median; one vote per included survey."),
  "Survey years come from the dataset build manifest; all records from each survey stay together.",
  "UN M49 geography: https://unstats.un.org/unsd/methodology/m49/overview/ (accessed 2026-09-09).",
  paste("Southern Africa included in East/Central:", include_south),
  "Geography and period splits are separate sensitivities, not a crossed geography-by-period analysis.",
  "Same complete-case selection, covariate scaling, unweighted cloglog likelihood and full-band offset as the base model.",
  "Every new fit re-estimates smoothing parameters; same cubic basis dimensions, knots adapted to each fitting sample.",
  "Separate age fits have their own time spline and survey/country/region random-effect variances.",
  "Joint subgroup fits retain one time spline shared across ages within that subgroup.",
  "Curves are log hazard ratios relative to 20% PfPR; intervals condition on smoothing parameters and HIV imputation.",
  paste("Input signature:", input_signature), paste("Analysis signature:", signature)), file.path(out, "specification.txt"))
curves <- diagnostics <- counts <- list()
for (k in seq_len(nrow(jobs))) {
  job <- jobs[k, ]; id <- job$fit_id
  use <- if (id == "reference") rep(TRUE, nrow(d)) else if (job$type == "single_age") d$age_band == job$age_band else
    if (id %in% c("west", "east_central")) surveys$geography[survey_index] == id else surveys$period[survey_index] == id
  dd <- droplevels(d[use, , drop = FALSE])
  stopifnot(nrow(dd) > 0, !anyNA(dd), all(dd$band_years > 0))
  if (job$type != "single_age") stopifnot(identical(levels(dd$age_band), age))
  form <- cbh_sensitivity_formula(spec, job$type == "single_age")
  path <- file.path(private, paste0(id, ".rds"))
  fit_signature <- cbh_hash(list(signature, id, deparse(form)))
  message("Starting ", id, ": ", nrow(dd), " records; ", sum(dd$death), " deaths.")
  if (id == "reference") {
    saved <- readRDS(file.path(base_dir, "fit.rds"))
    stopifnot(identical(saved$signature, input_signature), identical(deparse(formula(saved$fit)), deparse(form)))
  } else {
    saved <- if (file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
    if (is.null(saved) || !identical(saved$signature, fit_signature)) {
      warnings <- character()
      set.seed(spec$seed)
      timing <- system.time(fit <- withCallingHandlers(bam(form, data = dd,
        family = binomial(link = "cloglog"), method = "fREML", discrete = TRUE,
        nthreads = spec$nthreads, gc.level = 1, na.action = na.fail,
        control = gam.control(trace = FALSE, maxit = 100)), warning = function(w) {
          warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleWarning")
        }))
      saved <- list(fit = fit, elapsed_seconds = unname(timing["elapsed"]), warnings = warnings,
        signature = fit_signature, input_signature = input_signature, session_info = sessionInfo())
      cbh_atomic_rds(saved, path)
      rm(fit)
    }
  }
  diag <- cbh_sensitivity_diagnostics(saved, dd, id)
  diagnostics[[id]] <- diag
  cbh_atomic_csv(do.call(rbind, diagnostics), file.path(out, "fit_diagnostics.csv"))
  stopifnot(diag$converged, diag$finite_coefficients, diag$finite_covariance,
    length(saved$fit$fitted.values) == nrow(dd))
  curves[[id]] <- cbh_sensitivity_curves(saved$fit, dd, id)
  cbh_atomic_csv(do.call(rbind, curves), file.path(out, "pfpr_curves.csv"))
  tab <- as.data.table(dd)[, .(rows = .N, deaths = sum(death), countries = uniqueN(country), surveys = uniqueN(survey)), by = age_band]
  tab$fit_id <- id; counts[[id]] <- as.data.frame(tab)
  cbh_atomic_csv(do.call(rbind, counts), file.path(out, "sample_by_fit_age.csv"))
  st <- as.data.frame(summary(saved$fit)$s.table)
  st$term <- rownames(st)
  cbh_atomic_csv(st, file.path(out, paste0(id, "_smooth_summary.csv")))
  writeLines(deparse(form), file.path(out, paste0(id, "_formula.txt")))
  message("Completed ", id, " in ", round(diag$elapsed_seconds, 1), " fit seconds.")
  rm(saved, dd); gc(FALSE)
}
stopifnot(length(diagnostics) == 12L)
cbh_atomic_csv(data.frame(file = code_files, md5 = unname(code_hash)), file.path(out, "code_provenance.csv"))
message("All sensitivity fits complete: ", out)

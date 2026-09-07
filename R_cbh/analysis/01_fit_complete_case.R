#!/usr/bin/env Rscript
# Run from project root. Model objects and model input stay in ignored data/.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--force", "--prepare-only"))) stop("Arguments: --force --prepare-only")
spec <- cbh_trial_spec()
input_dir <- file.path(getwd(), "data/derived_cbh")
private_dir <- file.path(input_dir, "models", spec$id)
report_dir <- file.path(getwd(), "results/cbh", spec$id)
dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- readRDS(file.path(input_dir, "manifest.rds"))
stopifnot(isTRUE(manifest$complete))
code_files <- c("R_cbh/analysis/model.R", "R_cbh/analysis/01_fit_complete_case.R", "R_cbh/R/utils.R")
signature <- cbh_hash(list(spec, manifest$manifest, manifest$schema_version,
  vapply(code_files, cbh_file_hash, character(1)), as.character(packageVersion("mgcv")), R.version.string))
data_path <- file.path(private_dir, "complete_case_dataset.rds")
cached <- if (file.exists(data_path) && !"--force" %in% args) readRDS(data_path) else NULL
if (!is.null(cached) && identical(cached$signature, signature)) {
  prepared <- cached; message("Using complete-case dataset with matching provenance.")
} else {
  prepared <- cbh_trial_data(input_dir, spec)
  prepared$signature <- signature; prepared$specification <- spec
  cbh_atomic_rds(prepared, data_path)
}
rm(cached)
d <- prepared$data
message(sprintf("Fitting set: %s rows, %s deaths, %s countries, %s surveys, %s regions.",
  nrow(d), sum(d$death), nlevels(d$country), nlevels(d$survey), nlevels(d$region)))
for (name in c("selection", "missing", "skipped", "scaling")) {
  cbh_atomic_csv(prepared[[name]], file.path(report_dir, paste0(name, ".csv")))
}
band_counts <- aggregate(cbind(rows = rep(1L, nrow(d)), deaths = d$death),
  list(age_band = d$age_band), sum)
cbh_atomic_csv(band_counts, file.path(report_dir, "sample_by_age_band.csv"))
country_counts <- aggregate(prepared$selection[c("eligible_rows", "pfpr_available_rows", "complete_case_rows", "complete_case_deaths")],
  prepared$selection["country"], sum)
cbh_atomic_csv(country_counts, file.path(report_dir, "sample_by_country.csv"))
writeLines(c(paste(deparse(cbh_trial_formula(spec)), collapse = " "),
  "Family: binomial(cloglog); fixed full-band width offset; unweighted conditional likelihood.",
  "All selected complete cases; no vaccine covariates or imputation.",
  "Age-specific linear effects for continuous confounders; wealth is categorical.",
  paste("Input signature:", signature)), file.path(report_dir, "specification.txt"))
if ("--prepare-only" %in% args) quit(save = "no", status = 0)
fit_path <- file.path(private_dir, "fit.rds")
result <- if (file.exists(fit_path) && !"--force" %in% args) readRDS(fit_path) else NULL
if (is.null(result) || !identical(result$signature, signature)) {
  message("Starting joint seven-band mgcv::bam fit on all complete cases.")
  result <- cbh_trial_fit(d, spec)
  result$signature <- signature; result$specification <- spec
  result$session_info <- utils::sessionInfo()
  cbh_atomic_rds(result, fit_path)
}
fit <- result$fit
diagnostics <- list(converged = fit$converged, rank = fit$rank, coefficients = length(coef(fit)),
  finite_coefficients = all(is.finite(coef(fit))), finite_covariance = all(is.finite(fit$Vp)),
  mgcv_convergence = fit$mgcv.conv, outer_info = fit$outer.info,
  elapsed_seconds = result$elapsed_seconds, warnings = result$warnings)
capture.output(diagnostics, file = file.path(report_dir, "convergence.txt"))
capture.output(summary(fit), file = file.path(report_dir, "model_summary.txt"))
if (identical(fit$converged, FALSE) || !diagnostics$finite_coefficients || !diagnostics$finite_covariance) {
  stop("Fit failed convergence/finite-value checks; inspect convergence.txt. No contrasts published.")
}
contrasts <- cbh_trial_contrasts(fit, d)
cbh_atomic_csv(contrasts, file.path(report_dir, "pfpr_40_to_20_contrasts.csv"))
set.seed(spec$seed)
kcheck <- mgcv::k.check(fit, subsample = 10000, n.rep = 200)
cbh_atomic_csv(data.frame(term = rownames(kcheck), kcheck, row.names = NULL, check.names = FALSE),
               file.path(report_dir, "basis_checks.csv"))
print(contrasts, row.names = FALSE)
message("Completed exploratory complete-case fit. Reports: ", report_dir)

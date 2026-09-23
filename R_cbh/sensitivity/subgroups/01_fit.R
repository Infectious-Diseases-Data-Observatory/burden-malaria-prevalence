#!/usr/bin/env Rscript
# Refit the seven primary age-band models in four subgroups of the fitted
# 17-variable sample: Sahel (>=12N), Eastern Africa, and surveys before/after
# the median survey year. Subsets the saved prepared data; no data/HIV rebuild.
# Default: the DHS and MICS primary (v5). --dhs-only: the DHS-only primary (v3).
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/sensitivity/model.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
source("R_cbh/sensitivity/subgroups/settings.R")
library(mgcv)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(all(args %in% c("--force", "--dhs-only")))
st <- cbh_subgroup_settings(cbh_subgroup_sample(args)); primary <- cbh_primary_settings(st$primary_version)
stopifnot(isTRUE(primary$nutrition))
for (p in c(st$out, st$private)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(d, name) cbh_atomic_csv(d, file.path(st$out, name))
code <- c("R_cbh/sensitivity/subgroups/01_fit.R", "R_cbh/sensitivity/subgroups/settings.R",
          "R_cbh/sensitivity/model.R", "R_cbh/analysis/model.R", "R_cbh/primary/settings.R",
          "R_cbh/primary/specification.R")
files <- c(primary$data, st$centroids, st$registry, code)
provenance <- data.frame(file = files, md5 = vapply(files, cbh_file_hash, ""))
reference <- cbh_read_csv(file.path(primary$out, "fit_manifest.csv"))
stopifnot(nrow(reference) == 7L, all(reference$prepared_data_md5 == provenance$md5[1]))
write_csv(provenance, "fit_input_provenance.csv")
message("Reading the fixed primary prepared sample: ", primary$id)
prepared <- readRDS(primary$data)
input <- as.data.table(prepared$data)
write_csv(prepared$scaling, "covariate_scaling.csv")
input_signature <- prepared$signature
rm(prepared); gc(FALSE)
stopifnot(nrow(input) == primary$expected_records, sum(input$death) == primary$expected_deaths,
          uniqueN(input$region) == primary$expected_regions,
          uniqueN(input, by = c("child_id", "age_band")) == nrow(input))

## ---- survey and region attributes ------------------------------------------------
surveys <- unique(as.data.frame(input[, .(survey = as.character(survey), country = as.character(country), survey_year)]))
cbh_unique(surveys, "survey", "Included surveys")
stopifnot(all(is.finite(surveys$survey_year)))
cutoff <- median(surveys$survey_year)
surveys$period <- ifelse(surveys$survey_year <= cutoff, "early", "late")
surveys$east_africa <- surveys$country %in% st$east_africa
write_csv(surveys[order(surveys$country, surveys$survey_year), ], "survey_groups.csv")
message("Median survey year among ", nrow(surveys), " included surveys: ", cutoff,
        " (", sum(surveys$period == "early"), " early, ", sum(surveys$period == "late"), " late)")

read_all <- function(paths, cols) do.call(rbind, lapply(paths, function(x) cbh_read_csv(x)[cols]))
registry <- read_all(st$registry, c("SurveyId", "iso3", "svkey")); cbh_unique(registry, "SurveyId", "Registry")
centroids <- read_all(st$centroids, c("SurveyId", "iso3", "year", "regkey", "lon", "lat"))
cbh_unique(centroids, c("SurveyId", "regkey"), "Centroids")
j <- match(centroids$SurveyId, registry$SurveyId)
stopifnot(!anyNA(j), all(centroids$iso3 == registry$iso3[j]))
centroids$survey <- registry$svkey[j]
centroids$region <- paste(centroids$iso3, centroids$survey, centroids$regkey, sep = ":")
cbh_unique(centroids, "region", "Survey region centroids")
regions <- as.data.frame(input[, .(records = .N, deaths = sum(death), distinct_children = uniqueN(child_id)),
  by = .(survey = as.character(survey), country = as.character(country), region = as.character(region))])
j <- match(regions$region, centroids$region)
stopifnot(!anyNA(j), all(regions$survey == centroids$survey[j]))
regions$lat <- centroids$lat[j]; regions$lon <- centroids$lon[j]
stopifnot(all(is.finite(regions$lat)), all(is.finite(regions$lon)))
k <- match(regions$survey, surveys$survey)
regions$survey_year <- surveys$survey_year[k]
regions$sahel <- regions$lat >= st$sahel$latitude_min & regions$lon < st$sahel$longitude_max &
  !regions$country %in% st$sahel$exclude_countries
regions$east_africa <- surveys$east_africa[k]
regions$early <- surveys$period[k] == "early"
regions$late <- surveys$period[k] == "late"
write_csv(regions, "region_selection.csv")
subgroups <- c("sahel", "east_africa", "early", "late")
summary_sample <- function(d, subgroup) data.frame(subgroup, records = nrow(d), deaths = sum(d$death),
  distinct_children = uniqueN(d$child_id), surveys = uniqueN(d$survey), countries = uniqueN(d$country),
  survey_regions = uniqueN(d$region), countries_list = paste(sort(unique(as.character(d$country))), collapse = " "))
sample <- rbind(summary_sample(input, "map_full"), do.call(rbind, lapply(subgroups, function(s)
  summary_sample(input[region %in% regions$region[regions[[s]]]], s))))
write_csv(sample, "sample_summary.csv")
print(sample[, setdiff(names(sample), "countries_list")])
writeLines(c(sprintf("median_survey_year,%d", as.integer(cutoff)),
  sprintf("early_definition,survey_year <= %d", as.integer(cutoff)),
  sprintf("late_definition,survey_year > %d", as.integer(cutoff)),
  sprintf("sahel_definition,centroid latitude >= %s and longitude < %s and country not in %s",
    st$sahel$latitude_min, st$sahel$longitude_max, paste(st$sahel$exclude_countries, collapse = "/")),
  sprintf("east_africa_definition,UN M49 Eastern Africa: %s", paste(st$east_africa, collapse = " "))),
  file.path(st$out, "subgroup_definitions.csv"))

## ---- fits ---------------------------------------------------------------------------
form <- cbh_primary_regional_formula(primary)
input <- as.data.frame(input[, unique(c(all.vars(form), "age_band", "region")), with = FALSE])
stopifnot(!anyNA(input))
writeLines(deparse(form), file.path(st$out, "model_formula.txt"))
run <- function(d, knots, initial = NULL, strict = FALSE) {
  warnings <- character(); set.seed(primary$seed)
  control <- if (strict) gam.control(epsilon = 1e-9, mgcv.tol = 1e-9, efs.tol = .001, maxit = 200) else gam.control(maxit = 100)
  timing <- system.time(f <- withCallingHandlers(bam(form, data = d, knots = knots,
    family = binomial(link = "cloglog"), method = "fREML", discrete = TRUE,
    gamma = primary$gamma, select = FALSE, nthreads = primary$nthreads, gc.level = 1,
    na.action = na.fail, coef = initial, control = control),
    warning = function(w) { warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleWarning") }))
  list(fit = f, elapsed_seconds = unname(timing["elapsed"]), warnings = warnings, session_info = sessionInfo())
}
ages <- cbh_config()$age_bands$age_band
widths <- c(1, 5, 6, 12, 12, 12, 12) / 12
diags <- curves <- knots_saved <- manifests <- smooths <- restarts <- list()
for (s in subgroups) {
  sub <- droplevels(input[input$region %in% regions$region[regions[[s]]], , drop = FALSE])
  for (i in seq_along(ages)) {
    age <- ages[i]; id <- paste0(s, "_age_", i)
    d <- droplevels(sub[sub$age_band == age, , drop = FALSE])
    stopifnot(nrow(d) > 0, all(d$death %in% 0:1), max(abs(d$band_years - widths[i])) < 1e-12)
    # Subgroup-specific knots: quantiles of distinct predictor values (mgcv's cr convention).
    knots <- list(pfpr_pct = as.numeric(quantile(unique(d$pfpr_pct), seq(0, 1, length.out = 5))),
                  calendar_year = as.numeric(quantile(unique(d$calendar_year), seq(0, 1, length.out = 6))))
    stopifnot(all(diff(knots$pfpr_pct) > 0), all(diff(knots$calendar_year) > 0))
    sig <- cbh_hash(list(provenance, st, primary$gamma, knots, deparse(form), R.version.string,
                         as.character(packageVersion("mgcv"))))
    path <- file.path(st$private, paste0(id, ".rds"))
    saved <- if (file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
    message("Starting ", id, ": ", nrow(d), " records, ", sum(d$death), " deaths")
    if (is.null(saved) || !identical(saved$signature, sig)) {
      saved <- run(d, knots)
      saved$signature <- sig; saved$input_signature <- input_signature
      saved$gamma <- primary$gamma; saved$knots <- knots
      cbh_atomic_rds(saved, path)
    }
    diag <- cbh_sensitivity_diagnostics(saved, d, id)
    version <- "original"
    if (!diag$converged || !is.finite(diag$min_smoothing_hessian_eigenvalue) || diag$min_smoothing_hessian_eigenvalue <= 0) {
      message("Checking strict numerical restart for ", id)
      spath <- file.path(st$private, paste0(id, "_strict.rds"))
      strict <- if (file.exists(spath) && !"--force" %in% args) readRDS(spath) else NULL
      if (is.null(strict) || !identical(strict$signature, sig)) {
        strict <- run(d, knots, coef(saved$fit), TRUE)
        strict$signature <- sig; strict$input_signature <- input_signature
        strict$gamma <- primary$gamma; strict$knots <- knots
        cbh_atomic_rds(strict, spath)
      }
      sd <- cbh_sensitivity_diagnostics(strict, d, id)
      # Accept the restart when it converges with full rank, finite covariance, a
      # positive smoothing Hessian and a fREML criterion no worse than the original.
      # Unlike the primary fitter, no reduction of the (already tiny) smoothing
      # gradient is required; both gradients are recorded.
      take <- sd$converged && sd$finite_covariance && sd$rank == sd$coefficients &&
        is.finite(sd$min_smoothing_hessian_eigenvalue) && sd$min_smoothing_hessian_eigenvalue > 0 &&
        strict$fit$gcv.ubre <= saved$fit$gcv.ubre
      restarts[[id]] <- data.frame(fit_id = id, subgroup = s, original_min_hessian = diag$min_smoothing_hessian_eigenvalue,
        strict_min_hessian = sd$min_smoothing_hessian_eigenvalue, original_max_gradient = diag$max_smoothing_gradient,
        strict_max_gradient = sd$max_smoothing_gradient, original_fREML = unname(saved$fit$gcv.ubre),
        strict_fREML = unname(strict$fit$gcv.ubre), strict_selected = take)
      write_csv(do.call(rbind, restarts), "restarts.csv")
      if (take) { saved <- strict; diag <- sd; version <- "strict_restart"; path <- spath }
      rm(strict)
    }
    f <- saved$fit
    stopifnot(diag$converged, diag$finite_coefficients, diag$finite_covariance, diag$rank == diag$coefficients,
      diag$min_smoothing_hessian_eigenvalue > 0, all(f$prior.weights == 1), f$family$link == "cloglog",
      inherits(f$smooth[[1]], "cr.smooth"), inherits(f$smooth[[2]], "cr.smooth"),
      identical(f$smooth[[1]]$xp, knots$pfpr_pct), identical(f$smooth[[2]]$xp, knots$calendar_year))
    for (v in names(f$model)) {
      want <- if (v == "offset(log(band_years))") log(d$band_years) else d[[v]]
      stopifnot(isTRUE(all.equal(f$model[[v]], want, check.attributes = FALSE)))
    }
    diag$subgroup <- s; diag$age_band <- age; diag$gamma <- primary$gamma; diag$input_verified <- TRUE
    diags[[id]] <- diag
    cc <- cbh_sensitivity_curves(f, d, id); cc$subgroup <- s; curves[[id]] <- cc
    sm <- as.data.frame(summary(f)$s.table); sm$term <- rownames(sm); sm$subgroup <- s; sm$age_band <- age
    smooths[[id]] <- sm
    knots_saved[[id]] <- do.call(rbind, lapply(names(knots), function(v)
      data.frame(subgroup = s, age_band = age, variable = v, index = seq_along(knots[[v]]), value = knots[[v]])))
    manifests[[id]] <- data.frame(fit_id = id, subgroup = s, age_band = age, model_file = path,
      md5 = cbh_file_hash(path), selected_version = version, prepared_data_md5 = provenance$md5[1])
    write_csv(do.call(rbind, diags), "fit_diagnostics.csv")
    write_csv(do.call(rbind, curves), "pfpr_curves.csv")
    write_csv(do.call(rbind, smooths), "smooth_summaries.csv")
    write_csv(do.call(rbind, knots_saved), "knots.csv")
    write_csv(do.call(rbind, manifests), "fit_manifest.csv")
    message("Completed ", id, "; ", round(diag$elapsed_seconds, 1), " fit seconds; Hessian minimum ",
            signif(diag$min_smoothing_hessian_eigenvalue, 3))
    rm(saved, f, d); gc(FALSE)
  }
  rm(sub); gc(FALSE)
}
stopifnot(length(diags) == 28L, identical(provenance$md5, unname(vapply(files, cbh_file_hash, ""))))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(st$out, "session_info.txt"))
message("All 28 subgroup gamma=2 fits completed")

#!/usr/bin/env Rscript
# Run from the project root. Subset the saved primary input; no data/HIV rebuild.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/sensitivity/model.R")
source("R_cbh/primary/settings.R")
source("R_cbh/sensitivity/sahel/settings.R")
library(mgcv)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(all(args %in% "--force"))
st <- cbh_sahel_settings(); primary <- cbh_primary_settings()
for (p in c(st$out, st$private)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(d, name) cbh_atomic_csv(d, file.path(st$out, name))
code <- c("R_cbh/sensitivity/sahel/01_fit.R", "R_cbh/sensitivity/sahel/settings.R",
          "R_cbh/sensitivity/model.R", "R_cbh/analysis/model.R", "R_cbh/primary/settings.R")
files <- c(primary$data, st$centroids, st$registry, code)
provenance <- data.frame(file = files, md5 = vapply(files, cbh_file_hash, ""))
reference <- cbh_read_csv(file.path(primary$out, "fit_manifest.csv"))
stopifnot(nrow(reference) == 7L, all(reference$prepared_data_md5 == provenance$md5[1]))
write_csv(provenance, "fit_input_provenance.csv")
message("Reading the fixed primary prepared sample")
prepared <- readRDS(primary$data)
input <- as.data.table(prepared$data)
write_csv(prepared$scaling, "covariate_scaling.csv")
input_signature <- prepared$signature
rm(prepared); gc(FALSE)
stopifnot(nrow(input) == 5885022L, sum(input$death) == 82415L,
          uniqueN(input, by = c("child_id", "age_band")) == nrow(input))

# Survey-specific boundary keys, not a country-average latitude or an approximate name join.
registry <- cbh_read_csv(st$registry)
cbh_unique(registry, "SurveyId", "Registry")
centroids <- cbh_read_csv(st$centroids)
cbh_unique(centroids, c("SurveyId", "regkey"), "Centroids")
j <- match(centroids$SurveyId, registry$SurveyId)
stopifnot(!anyNA(j), all(centroids$iso3 == registry$iso3[j]))
centroids$survey <- registry$svkey[j]
centroids$region <- paste(centroids$iso3, centroids$survey, centroids$regkey, sep = ":")
cbh_unique(centroids, "region", "Survey region centroids")
regions <- as.data.frame(input[, .(records = .N, deaths = sum(death),
  distinct_children = uniqueN(child_id)), by = .(survey, country, region)])
j <- match(regions$region, centroids$region)
stopifnot(!anyNA(j), all(as.character(regions$survey) == centroids$survey[j]))
regions$regkey <- centroids$regkey[j]
regions$survey_year <- registry$year[match(regions$survey, registry$svkey)]
regions$lat <- centroids$lat[j]; regions$lon <- centroids$lon[j]
stopifnot(all(is.finite(regions$lat)), all(is.finite(regions$lon)))
regions$selected <- regions$lat >= st$latitude_min & regions$lon < st$longitude_max &
  !regions$country %in% st$exclude_countries
regions$exclusion <- ifelse(regions$selected, "included",
  ifelse(regions$country %in% st$exclude_countries | regions$lon >= st$longitude_max,
         "outside_Sahel_geographic_rule", "centroid_south_of_12N"))
write_csv(regions, "region_selection.csv")
write_csv(regions[regions$selected, ], "sahel_regions.csv")
summary_sample <- function(d, series) data.frame(series, records = nrow(d), deaths = sum(d$death),
  distinct_children = uniqueN(d$child_id), surveys = uniqueN(d$survey),
  countries = uniqueN(d$country), survey_regions = uniqueN(d$region))
sample <- summary_sample(input, "map_full")
input <- droplevels(input[region %in% regions$region[regions$selected]])
sample <- rbind(sample, summary_sample(input, "sahel_12N"))
write_csv(sample, "sample_summary.csv")
write_csv(as.data.frame(input[, .(records = .N, deaths = sum(death),
  distinct_children = uniqueN(child_id), surveys = uniqueN(survey),
  survey_regions = uniqueN(region)), by = country]), "country_sample.csv")
write_csv(as.data.frame(input[, .(records = .N, deaths = sum(death),
  distinct_children = uniqueN(child_id), surveys = uniqueN(survey),
  countries = uniqueN(country), survey_regions = uniqueN(region)), by = age_band]), "age_sample.csv")
print(sample)
form <- cbh_sensitivity_formula(cbh_trial_spec(), single_age = TRUE)
input <- as.data.frame(input[, unique(c(all.vars(form), "age_band")), with = FALSE])
stopifnot(!anyNA(input))
writeLines(deparse(form), file.path(st$out, "model_formula.txt"))

run <- function(d, knots, initial = NULL, strict = FALSE) {
  warnings <- character(); set.seed(primary$seed)
  control <- if (strict) gam.control(epsilon = 1e-9, mgcv.tol = 1e-9, efs.tol = .001,
                                    maxit = 200) else gam.control(maxit = 100)
  timing <- system.time(f <- withCallingHandlers(bam(form, data = d, knots = knots,
    family = binomial(link = "cloglog"), method = "fREML", discrete = TRUE,
    gamma = primary$gamma, select = FALSE, nthreads = primary$nthreads, gc.level = 1,
    na.action = na.fail, coef = initial, control = control),
    warning = function(w) { warnings <<- unique(c(warnings, conditionMessage(w)))
                           invokeRestart("muffleWarning") }))
  list(fit = f, elapsed_seconds = unname(timing["elapsed"]), warnings = warnings,
       session_info = sessionInfo())
}
ages <- cbh_config()$age_bands$age_band
diags <- curves <- knots_saved <- manifests <- smooths <- restarts <- list()
for (i in seq_along(ages)) {
  age <- ages[i]; id <- paste0("sahel_12N_age_", i)
  d <- droplevels(input[input$age_band == age, , drop = FALSE])
  stopifnot(nrow(d) > 0, all(d$death %in% 0:1), all(d$band_years == c(1,5,6,12,12,12,12)[i]/12))
  # Subgroup-specific knots as in the plan: quantiles of distinct predictor values
  # (mgcv's cr convention), recorded explicitly before discrete fitting.
  knots <- list(pfpr_pct = as.numeric(quantile(unique(d$pfpr_pct), seq(0, 1, length.out = 5))),
                calendar_year = as.numeric(quantile(unique(d$calendar_year), seq(0, 1, length.out = 6))))
  stopifnot(all(diff(knots$pfpr_pct) > 0), all(diff(knots$calendar_year) > 0))
  sig <- cbh_hash(list(provenance, st, primary$gamma, knots, deparse(form),
    R.version.string, as.character(packageVersion("mgcv"))))
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
  if (!diag$converged || !is.finite(diag$min_smoothing_hessian_eigenvalue) ||
      diag$min_smoothing_hessian_eigenvalue <= 0) {
    message("Checking strict numerical restart for ", id)
    strict <- run(d, knots, coef(saved$fit), TRUE)
    sd <- cbh_sensitivity_diagnostics(strict, d, id)
    take <- sd$converged && sd$finite_covariance && sd$rank == sd$coefficients &&
      is.finite(sd$min_smoothing_hessian_eigenvalue) && sd$min_smoothing_hessian_eigenvalue > 0 &&
      sd$max_smoothing_gradient < diag$max_smoothing_gradient && strict$fit$gcv.ubre <= saved$fit$gcv.ubre
    restarts[[id]] <- data.frame(fit_id = id, original_min_hessian = diag$min_smoothing_hessian_eigenvalue,
      strict_min_hessian = sd$min_smoothing_hessian_eigenvalue,
      original_fREML = unname(saved$fit$gcv.ubre), strict_fREML = unname(strict$fit$gcv.ubre), selected = take)
    write_csv(do.call(rbind, restarts), "restarts.csv")
    if (take) {
      strict$signature <- sig; strict$input_signature <- input_signature
      strict$gamma <- primary$gamma; strict$knots <- knots
      saved <- strict; diag <- sd; version <- "strict_restart"
      path <- file.path(st$private, paste0(id, "_strict.rds")); cbh_atomic_rds(saved, path)
    }
    rm(strict)
  }
  f <- saved$fit
  stopifnot(diag$converged, diag$finite_coefficients, diag$finite_covariance,
    diag$rank == diag$coefficients, diag$min_smoothing_hessian_eigenvalue > 0,
    all(f$prior.weights == 1), f$family$link == "cloglog",
    inherits(f$smooth[[1]], "cr.smooth"), inherits(f$smooth[[2]], "cr.smooth"),
    identical(f$smooth[[1]]$xp, knots$pfpr_pct), identical(f$smooth[[2]]$xp, knots$calendar_year))
  for (v in names(f$model)) {
    want <- if (v == "offset(log(band_years))") log(d$band_years) else d[[v]]
    stopifnot(isTRUE(all.equal(f$model[[v]], want, check.attributes = FALSE)))
  }
  diag$age_band <- age; diag$gamma <- primary$gamma; diag$input_verified <- TRUE
  diags[[id]] <- diag
  curves[[id]] <- cbh_sensitivity_curves(f, d, id)
  sm <- as.data.frame(summary(f)$s.table); sm$term <- rownames(sm); sm$age_band <- age
  smooths[[id]] <- sm
  knots_saved[[id]] <- do.call(rbind, lapply(names(knots), function(v)
    data.frame(age_band = age, variable = v, index = seq_along(knots[[v]]), value = knots[[v]])))
  manifests[[id]] <- data.frame(fit_id = id, age_band = age, model_file = path,
    md5 = cbh_file_hash(path), selected_version = version, prepared_data_md5 = provenance$md5[1])
  write_csv(do.call(rbind, diags), "fit_diagnostics.csv")
  write_csv(do.call(rbind, curves), "pfpr_curves.csv")
  write_csv(do.call(rbind, smooths), "smooth_summaries.csv")
  write_csv(do.call(rbind, knots_saved), "knots.csv")
  write_csv(do.call(rbind, manifests), "fit_manifest.csv")
  message("Completed ", id, "; ", round(diag$elapsed_seconds, 1), " seconds; Hessian minimum ",
          signif(diag$min_smoothing_hessian_eigenvalue, 3))
  rm(saved, f, d); gc(FALSE)
}
stopifnot(length(diags) == 7L, identical(provenance$md5, unname(vapply(files, cbh_file_hash, ""))))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(st$out, "session_info.txt"))
message("All seven Sahel gamma=2 fits completed")

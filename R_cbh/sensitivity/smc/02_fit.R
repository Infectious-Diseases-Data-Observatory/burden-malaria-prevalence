#!/usr/bin/env Rscript
# SMC before/after analysis: refit the seven primary DHS+MICS age-band models in the
# countries that have introduced SMC, adding an admin-1 SMC term. Formula, covariates
# (full-sample scaling), reference knots, gamma, likelihood, offset and random effects
# are the primary's; only the sample and the SMC term differ. Variants and the switch-on
# rule are described in settings.R.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/sensitivity/model.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
source("R_cbh/sensitivity/smc/settings.R")
library(mgcv); library(data.table)
args <- commandArgs(trailingOnly = TRUE); stopifnot(all(args %in% "--force"))
st <- cbh_smc_settings(); primary <- cbh_primary_settings(st$primary_version)
write_csv <- function(d, name) cbh_atomic_csv(as.data.frame(d), file.path(st$out, name))
paths <- file.path(st$private, c("unit_smc.csv", "unit_smc_coverage.csv", "nigeria_child_state.csv"))
code <- c("R_cbh/sensitivity/smc/02_fit.R", "R_cbh/sensitivity/smc/settings.R", "R_cbh/sensitivity/smc/01_assign_smc.R",
          "R_cbh/sensitivity/model.R", "R_cbh/analysis/model.R", "R_cbh/primary/settings.R", "R_cbh/primary/specification.R")
files <- c(primary$data, primary$knots, paths, code)
provenance <- data.frame(file = files, md5 = vapply(files, cbh_file_hash, ""))
reference <- cbh_read_csv(file.path(primary$out, "fit_manifest.csv"))
stopifnot(nrow(reference) == 7L, all(reference$prepared_data_md5 == provenance$md5[1]))
write_csv(provenance, "fit_input_provenance.csv")
units <- fread(paths[1]); coverage <- fread(paths[2]); cross <- fread(paths[3])
prepared <- readRDS(primary$data); input <- as.data.table(prepared$data); input_signature <- prepared$signature
rm(prepared); gc(FALSE)
stopifnot(nrow(input) == primary$expected_records)

## ---- admin-1 SMC assignment -------------------------------------------------------------------
input <- input[as.character(country) %in% unique(units$iso3)]
input[, iso3 := as.character(country)]
start <- unique(units[, .(iso3, national_start)])
input[, national_start := start$national_start[match(iso3, start$iso3)]]
# Admin-1 unit: the analysis region, or the state for Nigerian records.
input[, unit := as.character(regkey)]
# Nigerian records take their state from the child-to-state crosswalk (01_assign_smc.R).
nga <- input$iso3 == "NGA"
input[nga, unit := cross$regkey[match(as.character(child_id), cross$child_id)]]
if (input[iso3 == "NGA" & is.na(unit) & entry_year >= national_start, .N])
  stop("Nigerian records without a state after the national SMC start")
input[iso3 == "NGA" & is.na(unit), unit := paste0("zone:", regkey)]
j <- match(paste(input$iso3, input$unit), paste(units$iso3, units$regkey))
bad <- is.na(j) & input$entry_year >= input$national_start
if (any(bad)) stop("Admin-1 units without an SMC record after the national start: ",
  paste(unique(paste(input$survey, input$unit)[bad]), collapse = ", "))
input[, `:=`(first_all = units$first_year_main[j], first_confirmed = units$first_year_confirmed[j],
  basis = fifelse(is.na(j), "pre_programme", units$status_basis[j]), earliest = units$earliest_cluster_year[j])]
# Main (confirmed status): national-basis units only while no campaign was recorded in them.
input[, main_keep := basis != "national" | is.na(earliest) | entry_year < earliest]
input[, smc_post := as.integer(basis == "confirmed" & !is.na(first_confirmed) & entry_year >= first_confirmed)]
input[, smc_placebo := as.integer(basis == "confirmed" & !is.na(first_confirmed) & entry_year >= first_confirmed - st$placebo_lead)]
input[, smc_national := as.integer(!is.na(first_all) & entry_year >= first_all)]
cv <- match(paste(input$iso3, input$unit, pmin(input$entry_year, max(coverage$year))), paste(coverage$iso3, coverage$regkey, coverage$year))
input[, smc_coverage := fifelse(basis == "confirmed", coverage$coverage_confirmed[cv], 0)]
input[is.na(smc_coverage), smc_coverage := 0]
stopifnot(all(input$smc_post[input$basis != "confirmed"] == 0L), all(input$smc_national[is.na(j)] == 0L),
          all(input$smc_coverage >= 0 & input$smc_coverage <= 1))
# Post-SMC bands (main) that end before 1 July of the switch-on year cannot overlap a campaign season.
input[, pre_season := smc_post == 1L & calendar_year + band_years <= first_confirmed + 0.5]

## ---- samples ---------------------------------------------------------------------------------------
samples <- list(main = input[main_keep == TRUE], no_smc = input[main_keep == TRUE],
  placebo = input[main_keep == TRUE & smc_post == 0L], national = input, coverage = input[main_keep == TRUE])
term <- c(main = "smc_post", no_smc = NA, placebo = "smc_placebo", national = "smc_national", coverage = "smc_coverage")
switchers <- function(d, v) d[, .(pre = any(get(v) == 0), post = any(get(v) > 0)), by = .(survey, iso3, unit)][pre & post]
describe <- function(v) { d <- samples[[v]]; t <- term[[v]]; sw <- switchers(d, t)
  data.frame(variant = v, records = nrow(d), deaths = sum(d$death), post_records = sum(d[[t]] > 0), post_deaths = sum(d$death[d[[t]] > 0]),
    mean_indicator = mean(d[[t]]), surveys = uniqueN(d$survey), countries = uniqueN(d$iso3), regions = uniqueN(d$region),
    admin1_units = uniqueN(d[, .(survey, iso3, unit)]), switching_units = nrow(sw),
    switching_by_country = paste(sw[, .N, by = iso3][order(iso3), sprintf("%s %d", iso3, N)], collapse = ", "),
    countries_with_post = paste(sort(unique(d$iso3[d[[t]] > 0])), collapse = " "), countries_list = paste(sort(unique(d$iso3)), collapse = " ")) }
summary_rows <- do.call(rbind, lapply(c("main", "placebo", "national", "coverage"), describe))
write_csv(summary_rows, "sample_summary.csv"); print(summary_rows[, c("variant", "records", "deaths", "post_records", "switching_units", "switching_by_country")])
by_country <- input[, .(records = .N, deaths = sum(death), main_records = sum(main_keep), main_post_records = sum(smc_post[main_keep]),
  main_post_deaths = sum(death[main_keep & smc_post == 1L]), national_post_records = sum(smc_national), national_post_deaths = sum(death[smc_national == 1L]),
  surveys = uniqueN(survey), admin1_units = uniqueN(paste(survey, unit)), entry_years = paste(range(entry_year), collapse = "-")), by = iso3][order(iso3)]
write_csv(by_country, "country_summary.csv")
write_csv(input[, .(records = .N, deaths = sum(death), main_kept = sum(main_keep), post_main = sum(smc_post), post_national = sum(smc_national),
  entry_first = min(entry_year), entry_last = max(entry_year)), by = .(iso3, survey, unit, basis, first_confirmed, first_all)][order(iso3, survey, unit)],
  "survey_admin1_assignment.csv")
timing <- input[main_keep == TRUE & smc_post == 1L, .(post_records = .N, pre_season_records = sum(pre_season), pre_season_deaths = sum(death[pre_season])), by = age_band][order(age_band)]
write_csv(timing, "timing_pre_season.csv")
omitted <- fread(file.path(primary$out, "prepared_selection_by_survey.csv"))[country %in% st$omitted_smc_countries & records > 0,
  .(records = sum(records), deaths = sum(deaths), surveys = .N), by = country]
write_csv(omitted[order(country)], "omitted_smc_countries.csv")

## ---- fits ----------------------------------------------------------------------------------------------
base <- cbh_primary_regional_formula(primary)
forms <- lapply(term, function(t) if (is.na(t)) base else update(base, as.formula(paste(". ~ . +", t))))
for (v in names(forms)) writeLines(deparse(forms[[v]]), file.path(st$out, paste0("model_formula_", v, ".txt")))
knot_table <- cbh_read_csv(primary$knots)
run <- function(form, d, knots, initial = NULL, strict = FALSE) {
  warnings <- character(); set.seed(primary$seed)
  control <- if (strict) gam.control(epsilon = 1e-9, mgcv.tol = 1e-9, efs.tol = .001, maxit = 200) else gam.control(maxit = 100)
  timing <- system.time(f <- withCallingHandlers(bam(form, data = d, knots = knots,
    family = binomial(link = "cloglog"), method = "fREML", discrete = TRUE,
    gamma = primary$gamma, select = FALSE, nthreads = primary$nthreads, gc.level = 1,
    na.action = na.fail, coef = initial, control = control),
    warning = function(w) { warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleWarning") }))
  list(fit = f, elapsed_seconds = unname(timing["elapsed"]), warnings = warnings, session_info = sessionInfo())
}
ages <- cbh_config()$age_bands$age_band; widths <- c(1, 5, 6, 12, 12, 12, 12) / 12
keep_cols <- unique(c(unlist(lapply(forms, all.vars)), "age_band", "region", "iso3", "survey", "unit"))
diags <- curves <- effects <- manifests <- restarts <- list()
for (v in st$variants) {
  sub <- as.data.frame(samples[[v]][, keep_cols, with = FALSE]); stopifnot(!anyNA(sub[setdiff(keep_cols, "unit")]))
  for (i in seq_along(ages)) {
    age <- ages[i]; id <- paste0(v, "_age_", i)
    d <- droplevels(sub[sub$age_band == age, , drop = FALSE])
    stopifnot(nrow(d) > 0, all(d$death %in% 0:1), max(abs(d$band_years - widths[i])) < 1e-12)
    knots <- setNames(lapply(c("pfpr_pct", "calendar_year"), function(x) {
      kk <- knot_table[knot_table$age_band == age & knot_table$variable == x, ]; kk$value[order(kk$index)] }), c("pfpr_pct", "calendar_year"))
    stopifnot(length(knots$pfpr_pct) == 5L, length(knots$calendar_year) == 6L, min(d$pfpr_pct) >= min(knots$pfpr_pct))
    sig <- cbh_hash(list(provenance, st[c("id", "majority", "placebo_lead", "variants")], v, primary$gamma, knots, deparse(forms[[v]]),
                         R.version.string, as.character(packageVersion("mgcv"))))
    path <- file.path(st$private, paste0(id, ".rds"))
    saved <- if (file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
    message("Starting ", id, ": ", nrow(d), " records, ", sum(d$death), " deaths")
    if (is.null(saved) || !identical(saved$signature, sig)) {
      saved <- run(forms[[v]], d, knots); saved$signature <- sig; saved$input_signature <- input_signature; saved$knots <- knots
      cbh_atomic_rds(saved, path)
    }
    diag <- cbh_sensitivity_diagnostics(saved, d, id); version <- "original"
    if (!diag$converged || !is.finite(diag$min_smoothing_hessian_eigenvalue) || diag$min_smoothing_hessian_eigenvalue <= 0) {
      message("Checking strict numerical restart for ", id)
      spath <- file.path(st$private, paste0(id, "_strict.rds"))
      strict <- if (file.exists(spath) && !"--force" %in% args) readRDS(spath) else NULL
      if (is.null(strict) || !identical(strict$signature, sig)) {
        strict <- run(forms[[v]], d, knots, coef(saved$fit), TRUE); strict$signature <- sig
        strict$input_signature <- input_signature; strict$knots <- knots; cbh_atomic_rds(strict, spath)
      }
      sd <- cbh_sensitivity_diagnostics(strict, d, id)
      # Subgroup-refit acceptance rule (converged, full rank, finite covariance, positive
      # smoothing Hessian, no smoothing-gradient reduction required), with fREML allowed to be
      # worse by at most 0.001: an original fit whose smallest Hessian eigenvalue is ~0 (a
      # variance at its boundary) and a strict refit within 0.001 are the same optimum.
      take <- sd$converged && sd$finite_covariance && sd$rank == sd$coefficients &&
        is.finite(sd$min_smoothing_hessian_eigenvalue) && sd$min_smoothing_hessian_eigenvalue > 0 &&
        strict$fit$gcv.ubre <= saved$fit$gcv.ubre + 1e-3
      restarts[[id]] <- data.frame(fit_id = id, variant = v, original_min_hessian = diag$min_smoothing_hessian_eigenvalue,
        strict_min_hessian = sd$min_smoothing_hessian_eigenvalue, original_max_gradient = diag$max_smoothing_gradient,
        strict_max_gradient = sd$max_smoothing_gradient, original_fREML = unname(saved$fit$gcv.ubre),
        strict_fREML = unname(strict$fit$gcv.ubre), strict_selected = take)
      write_csv(do.call(rbind, restarts), "restarts.csv")
      if (take) { saved <- strict; diag <- sd; version <- "strict_restart"; path <- spath }
      rm(strict)
    }
    f <- saved$fit
    stopifnot(diag$converged, diag$finite_coefficients, diag$finite_covariance, diag$rank == diag$coefficients,
      diag$min_smoothing_hessian_eigenvalue > 0, f$family$link == "cloglog",
      identical(f$smooth[[1]]$xp, knots$pfpr_pct), identical(f$smooth[[2]]$xp, knots$calendar_year))
    for (nm in names(f$model)) {
      want <- if (nm == "offset(log(band_years))") log(d$band_years) else d[[nm]]
      stopifnot(isTRUE(all.equal(f$model[[nm]], want, check.attributes = FALSE)))
    }
    re_sd <- setNames(vapply(c("survey", "country", "region"), function(g) {
      s <- which(vapply(f$smooth, function(z) identical(z$term, g), TRUE)); if (length(s)) 1 / sqrt(f$sp[s]) else NA_real_ }, 0), paste0("re_sd_", c("survey", "country", "region")))
    diag$variant <- v; diag$age_band <- age; diag$input_verified <- TRUE; diags[[id]] <- cbind(diag, as.list(re_sd))
    cc <- cbh_sensitivity_curves(f, d, id); cc$variant <- v; curves[[id]] <- cc
    if (!is.na(term[[v]])) {
      t <- term[[v]]; ix <- match(t, names(coef(f))); stopifnot(!is.na(ix))
      b <- unname(coef(f)[ix]); se <- sqrt(f$Vp[ix, ix])
      dd <- as.data.table(d); sw <- switchers(dd, t)
      effects[[id]] <- data.frame(variant = v, age_band = age, term = t, estimate = b, standard_error = se,
        hazard_ratio = exp(b), lower_95 = exp(b - 1.96 * se), upper_95 = exp(b + 1.96 * se), p_value = 2 * pnorm(-abs(b / se)),
        records = nrow(d), deaths = sum(d$death), post_records = sum(d[[t]] > 0), post_deaths = sum(d$death[d[[t]] > 0]),
        switching_units = nrow(sw))
    }
    manifests[[id]] <- data.frame(fit_id = id, variant = v, age_band = age, model_file = path, md5 = cbh_file_hash(path),
      selected_version = version, prepared_data_md5 = provenance$md5[1])
    write_csv(do.call(rbind, diags), "fit_diagnostics.csv"); write_csv(do.call(rbind, curves), "pfpr_curves.csv")
    write_csv(do.call(rbind, effects), "smc_effects.csv"); write_csv(do.call(rbind, manifests), "fit_manifest.csv")
    message("Completed ", id, "; ", round(diag$elapsed_seconds, 1), " fit seconds; Hessian minimum ", signif(diag$min_smoothing_hessian_eigenvalue, 3))
    rm(saved, f, d); gc(FALSE)
  }
  rm(sub); gc(FALSE)
}
stopifnot(length(diags) == 7L * length(st$variants), identical(provenance$md5, unname(vapply(files, cbh_file_hash, ""))))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(st$out, "session_info.txt"))
message("All SMC fits completed")

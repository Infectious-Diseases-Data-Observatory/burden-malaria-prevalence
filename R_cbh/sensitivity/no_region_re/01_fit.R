#!/usr/bin/env Rscript
# Primary model without the survey-region random intercept (added 24 September 2026, at the
# user's request). Identical to the v7 primary (same prepared data, formula, reference knots,
# gamma = 2 and bam settings) except that s(region, bs = "re") is dropped; survey and country
# random intercepts are kept. Self-contained so that the v7 fits and their code provenance are
# untouched. Fits are cached by signature; --force refits.
source("R_cbh/load_pipeline.R")
source("R_cbh/sensitivity/model.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
suppressPackageStartupMessages(library(mgcv))
source("R_cbh/sensitivity/no_region_re/settings.R")
args <- commandArgs(trailingOnly = TRUE); stopifnot(all(args %in% "--force"))
st <- cbh_no_region_settings(); base <- st$base
for (p in c(st$out, st$private)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ages <- cbh_config()$age_bands$age_band

base_form <- cbh_primary_regional_formula(base)
form <- as.formula(gsub(" + s(region, bs = \"re\")", "", gsub("\\s+", " ", deparse1(base_form)), fixed = TRUE))
stopifnot(!"region" %in% all.vars(form), setequal(setdiff(all.vars(base_form), "region"), all.vars(form)))
writeLines(deparse(form), file.path(st$out, "model_formula.txt"))

v7_manifest <- cbh_read_csv(file.path(base$out, "fit_manifest.csv"))
data_hash <- cbh_file_hash(base$data)
stopifnot(identical(unique(v7_manifest$prepared_data_md5), data_hash))
prepared <- readRDS(base$data); input <- prepared$data; rm(prepared); gc(FALSE)
stopifnot(nrow(input) == base$expected_records, sum(input$death) == base$expected_deaths)
input <- input[unique(c(all.vars(form), "age_band"))]
knot_table <- cbh_read_csv(base$knots)
code <- c("R_cbh/sensitivity/no_region_re/01_fit.R", "R_cbh/sensitivity/no_region_re/settings.R", base$knots, "R_cbh/primary/specification.R")
cbh_atomic_csv(data.frame(file = c(base$data, code), md5 = c(data_hash, vapply(code, cbh_file_hash, ""))), file.path(st$out, "fit_input_provenance.csv"))

diags <- summaries <- components <- manifest <- restarts <- list()
for (a in ages) {
  i <- match(a, ages); fit_id <- paste0("no_region_age_", i)
  d <- droplevels(input[input$age_band == a, , drop = FALSE])
  kt <- knot_table[knot_table$age_band == a, ]
  knots <- lapply(c(pfpr_pct = "pfpr_pct", calendar_year = "calendar_year"), function(v) { z <- kt[kt$variable == v, ]; z$value[order(z$index)] })
  sig <- cbh_hash(list(fit_id, knots, deparse(form), data_hash, R.version.string, as.character(packageVersion("mgcv"))))
  path <- file.path(st$private, paste0(fit_id, ".rds"))
  saved <- if (file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
  if (is.null(saved) || !identical(saved$signature, sig)) {
    message("Fitting ", fit_id, " (", a, "): ", nrow(d), " rows, ", sum(d$death), " deaths")
    warnings <- character(); set.seed(20260907L)
    timing <- system.time(f <- withCallingHandlers(bam(form, data = d, knots = knots, family = binomial(link = "cloglog"),
      method = "fREML", discrete = TRUE, gamma = 2, select = FALSE, nthreads = 2, gc.level = 1, na.action = na.fail,
      control = gam.control(maxit = 100)), warning = function(w) { warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleWarning") }))
    saved <- list(fit = f, elapsed_seconds = unname(timing["elapsed"]), warnings = warnings, signature = sig, knots = knots, gamma = 2)
    cbh_atomic_rds(saved, path)
  }
  diag <- cbh_sensitivity_diagnostics(saved, d, fit_id); version <- "original"
  # Same restart rule as R_cbh/primary/01_fit.R: a fit that has not converged or whose smoothing
  # Hessian is not positive definite is refitted from its coefficients with tighter tolerances,
  # and the restart is kept only if it converges, is positive definite and improves fREML.
  if (!diag$converged || !is.finite(diag$min_smoothing_hessian_eigenvalue) || diag$min_smoothing_hessian_eigenvalue <= 0) {
    spath <- file.path(st$private, paste0(fit_id, "_strict.rds")); ssig <- cbh_hash(list(sig, cbh_file_hash(path), "strict"))
    strict <- if (file.exists(spath) && !"--force" %in% args) readRDS(spath) else NULL
    if (is.null(strict) || !identical(strict$signature, ssig)) {
      message("Strict restart ", fit_id)
      warnings <- character(); set.seed(20260907L)
      timing <- system.time(f2 <- withCallingHandlers(bam(form, data = d, knots = knots, family = binomial(link = "cloglog"),
        method = "fREML", discrete = TRUE, gamma = 2, select = FALSE, nthreads = 2, gc.level = 1, na.action = na.fail, coef = coef(saved$fit),
        control = gam.control(epsilon = 1e-9, mgcv.tol = 1e-9, efs.tol = .001, maxit = 200)),
        warning = function(w) { warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleWarning") }))
      strict <- list(fit = f2, elapsed_seconds = unname(timing["elapsed"]), warnings = warnings, signature = ssig, knots = knots, gamma = 2)
      cbh_atomic_rds(strict, spath); rm(f2)
    }
    sd <- cbh_sensitivity_diagnostics(strict, d, fit_id)
    take <- sd$converged && sd$finite_covariance && sd$rank == sd$coefficients && is.finite(sd$min_smoothing_hessian_eigenvalue) &&
      sd$min_smoothing_hessian_eigenvalue > 0 && strict$fit$gcv.ubre <= saved$fit$gcv.ubre + 1e-3
    restarts[[a]] <- data.frame(fit_id, original_min_hessian = diag$min_smoothing_hessian_eigenvalue, strict_min_hessian = sd$min_smoothing_hessian_eigenvalue,
      original_fREML = unname(saved$fit$gcv.ubre), strict_fREML = unname(strict$fit$gcv.ubre), strict_selected = take)
    cbh_atomic_csv(do.call(rbind, restarts), file.path(st$out, "restarts.csv"))
    if (take) { saved <- strict; diag <- sd; path <- spath; version <- "strict_restart" }
  }
  f <- saved$fit; diag$age_band <- a; diag$selected_version <- version
  diag$aic <- AIC(f); diag$edf_total <- sum(f$edf); diag$fREML <- unname(f$gcv.ubre)
  stopifnot(diag$converged, diag$finite_covariance, diag$rank == diag$coefficients, identical(f$smooth[[1]]$xp, knots$pfpr_pct))
  diags[[a]] <- diag
  s_tab <- as.data.frame(summary(f)$s.table); s_tab$term <- rownames(s_tab); s_tab$age_band <- a; summaries[[a]] <- s_tab
  s <- f$smooth[[1]]; ix <- s$first.para:s$last.para
  components[[fit_id]] <- list(fit_id = fit_id, age_band = a, smooth = s, coef = coef(f)[ix], covariance = f$Vp[ix, ix],
    model_md5 = cbh_file_hash(path), support = quantile(d$pfpr_pct, c(0, .025, .975, 1), names = FALSE))
  manifest[[a]] <- data.frame(fit_id, age_band = a, model_file = path, md5 = cbh_file_hash(path), selected_version = version, prepared_data_md5 = data_hash)
  cbh_atomic_csv(do.call(rbind, diags), file.path(st$out, "fit_diagnostics.csv"))
  cbh_atomic_csv(do.call(rbind, summaries), file.path(st$out, "smooth_summaries.csv"))
  cbh_atomic_csv(do.call(rbind, manifest), file.path(st$out, "fit_manifest.csv"))
  rm(f, saved); if (exists("strict")) rm(strict); gc(FALSE)
}
cbh_atomic_rds(components, file.path(st$private, "pfpr_components.rds"))
message("No-region-RE fits complete: ", st$out)

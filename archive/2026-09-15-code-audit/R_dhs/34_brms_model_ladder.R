# =============================================================================
# 34_brms_model_ladder.R — do the data support a prevalence-by-time surface?
#
# Three models for post-neonatal mortality, identical apart from how prevalence
# and calendar time enter, fitted with shared priors and compared by PSIS-LOO:
#
#   A  additive         s(pfpr10) + s(year_c)
#   C  linear-in-time   A + s(pfpr10, by = year_c, pc = AF_REFERENCE / 10)
#                       the prevalence curve changes linearly with calendar
#                       time; nests A exactly (the by-smooth at zero IS A)
#   D  full surface     t2(pfpr10, year_c)             (script 30's model)
#
# Why C carries a point constraint. With a numeric `by` variable mgcv drops the
# centring constraint, so year_c * g(pfpr10) contains year_c times the constant
# part of g, which duplicates the linear year term already inside s(year_c)
# (design rank 18 of 19 without it). Pinning g to zero at the 1% reference
# prevalence removes that direction and makes g the per-year change in the
# prevalence effect RELATIVE TO THE REFERENCE, which is the attributable
# fraction's own scale.
#
# brms has no te() or ti(), so the mgcv ladder's "main effects plus a pure
# interaction" cannot be reproduced. D is not nested above A or C: each of its
# three penalised blocks bundles a main effect with its linear-in-the-other-
# variable modulation under one smoothness parameter. The comparison is
# therefore predictive (which model expects new region-years better), not a
# test of one term.
#
# Usage
#   Rscript R_dhs/34_brms_model_ladder.R                  fit or load all, compare
#   LADDER_MODEL=C Rscript R_dhs/34_brms_model_ladder.R   fit and cache one model
#
# Outputs
#   results/dhs_rebuild/brms_ladder_loo.csv
#   results/dhs_rebuild/brms_ladder_af_anchors.csv
#   results/dhs_rebuild/brms_ladder_time_contrast.csv
#   results/dhs_rebuild/figure17_brms_ladder.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("brms", "ggplot2", "mgcv", "loo"))

CHAINS <- 4L
ITERATIONS <- 2000L
WARMUP <- 1000L
ADAPT_DELTA <- 0.99
RIDGE_PRIOR_SCALE <- 1
ANCHORS <- c(10, 30, 50)
CONTRAST_YEARS <- c(2005L, 2020L)
CACHE <- file.path(DATA_DIR, "brms_ladder_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

bundle <- readRDS(MODEL_BUNDLE_RDS)
year_center <- unique(bundle$year_center)[1]
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
catalog$included_in_main <- as.logical(catalog$included_in_main)

## ---- the same panel script 30 fits --------------------------------------------
analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
analysis <- analysis[is.finite(analysis$postneonatal_mortality) &
                       analysis$postneonatal_mortality > 0 &
                       is.finite(analysis$exposure) & analysis$exposure > 0 &
                       is.finite(analysis$pfpr10) & is.finite(analysis$year_c),
                     , drop = FALSE]
ridge <- make_ridge_matrix(analysis, catalog, bundle$preprocessing)
covariates <- as.data.frame(ridge$matrix)
names(covariates) <- make.names(colnames(ridge$matrix))
model_data <- data.frame(
  deaths = round(analysis$postneonatal_mortality / 1000 * analysis$exposure),
  log_exposure = log(analysis$exposure),
  pfpr10 = analysis$pfpr10,
  year_c = analysis$year_c,
  country = factor(analysis$iso3),
  svkey = analysis$svkey,
  stringsAsFactors = FALSE
)
model_data <- cbind(model_data, covariates)
message("Model ladder: ", nrow(model_data), " region-years, ",
        nlevels(model_data$country), " countries, ", ncol(covariates),
        " covariates.")

## ---- the ladder -----------------------------------------------------------------
tail_terms <- paste("+", paste(names(covariates), collapse = " + "),
                    "+ (1 | country) + offset(log_exposure)")
PC <- AF_REFERENCE / 10
ladder <- list(
  A = list(label = "A: additive",
           detail = "s(pfpr10) + s(year_c)",
           formula = paste("deaths ~ s(pfpr10, k = 6) + s(year_c, k = 8)",
                           tail_terms)),
  C = list(label = "C: curve changes linearly in time",
           detail = "A + s(pfpr10, by = year_c)",
           formula = paste0(
             "deaths ~ s(pfpr10, k = 6) + s(year_c, k = 8) + ",
             "s(pfpr10, by = year_c, k = 6, pc = ", PC, ")", tail_terms)),
  D = list(label = "D: full surface",
           detail = "t2(pfpr10, year_c)",
           formula = paste("deaths ~ t2(pfpr10, year_c, k = c(6, 6))",
                           tail_terms))
)

# Shared priors, identical to script 30.
priors <- c(
  brms::prior_string(sprintf("normal(0, %g)", RIDGE_PRIOR_SCALE), class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "sds"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(gamma(0.01, 0.01), class = "shape")
)

fit_model <- function(name) {
  spec <- ladder[[name]]
  fingerprint <- list(rows = nrow(model_data),
                      countries = nlevels(model_data$country),
                      deaths = sum(model_data$deaths),
                      covariates = names(covariates), formula = spec$formula,
                      chains = CHAINS, iterations = ITERATIONS,
                      adapt_delta = ADAPT_DELTA)
  cache_file <- file.path(CACHE, paste0(name, ".rds"))
  meta_file <- sub("\\.rds$", "_meta.rds", cache_file)
  if (file.exists(cache_file) && file.exists(meta_file) &&
      isTRUE(all.equal(readRDS(meta_file), fingerprint))) {
    message("  ", spec$label, ": reusing cached fit")
    return(readRDS(cache_file))
  }
  message("  ", spec$label, ": sampling")
  fit <- brms::brm(
    formula = brms::bf(spec$formula), data = model_data,
    family = brms::negbinomial(), prior = priors,
    chains = CHAINS, iter = ITERATIONS, warmup = WARMUP,
    cores = min(CHAINS, max(1L, parallel::detectCores() - 1L)),
    control = list(adapt_delta = ADAPT_DELTA, max_treedepth = 12),
    # Every parameter is kept so PSIS-LOO can moment-match any high-k points.
    save_pars = brms::save_pars(all = TRUE),
    seed = 20260828, refresh = 100
  )
  saveRDS(fit, cache_file)
  saveRDS(fingerprint, meta_file)
  fit
}

only <- Sys.getenv("LADDER_MODEL", "")
if (nzchar(only)) {
  if (!only %in% names(ladder)) stop("LADDER_MODEL must be one of A, C, D")
  invisible(fit_model(only))
  message("Cached model ", only, "; run without LADDER_MODEL to compare.")
  quit(save = "no", status = 0)
}

fits <- lapply(names(ladder), fit_model)
names(fits) <- names(ladder)

## ---- sampling diagnostics ---------------------------------------------------------
diagnose <- function(fit) {
  draws <- posterior::summarise_draws(posterior::as_draws_df(fit))
  key <- draws[grepl("^b_|^bs_|^sds_|^sd_|^shape", draws$variable), ]
  nuts <- brms::nuts_params(fit)
  data.frame(
    divergent = sum(nuts$Parameter == "divergent__" & nuts$Value > 0),
    max_rhat = max(key$rhat, na.rm = TRUE),
    min_bulk_ess = min(key$ess_bulk, na.rm = TRUE))
}
diagnostics <- do.call(rbind, lapply(fits, diagnose))

## ---- PSIS-LOO, moment-matching where the importance weights are unstable ---------
message("\nPSIS-LOO")
loos <- lapply(names(fits), function(name) {
  l <- brms::loo(fits[[name]])
  high <- sum(l$diagnostics$pareto_k > 0.7)
  if (high > 0) {
    message("  ", name, ": ", high, " observations with Pareto k > 0.7; ",
            "moment matching")
    l <- brms::loo(fits[[name]], moment_match = TRUE)
  }
  l
})
names(loos) <- names(fits)
comparison <- loo::loo_compare(loos)
weights <- loo::loo_model_weights(loos, method = "stacking")
print(comparison, simplify = FALSE)
message("\nStacking weights:")
print(round(weights, 3))

loo_table <- data.frame(
  model = rownames(comparison),
  label = vapply(rownames(comparison), function(m) ladder[[m]]$label, ""),
  structure = vapply(rownames(comparison), function(m) ladder[[m]]$detail, ""),
  elpd_loo = comparison[, "elpd_loo"], se_elpd_loo = comparison[, "se_elpd_loo"],
  p_loo = comparison[, "p_loo"],
  elpd_diff = comparison[, "elpd_diff"], se_diff = comparison[, "se_diff"],
  pareto_k_over_0.7 = vapply(rownames(comparison), function(m)
    sum(loos[[m]]$diagnostics$pareto_k > 0.7), 1L),
  stacking_weight = as.numeric(weights[rownames(comparison)]),
  divergent = diagnostics[rownames(comparison), "divergent"],
  max_rhat = diagnostics[rownames(comparison), "max_rhat"],
  min_bulk_ess = diagnostics[rownames(comparison), "min_bulk_ess"],
  stringsAsFactors = FALSE
)
rownames(loo_table) <- NULL
write.csv(loo_table, file.path(RESULTS_DIR, "brms_ladder_loo.csv"),
          row.names = FALSE)

## ---- what the choice does to the attributable fraction ---------------------------
af_draws <- function(fit, prevalence, year) {
  frame <- function(p) {
    out <- data.frame(pfpr10 = p / 10, year_c = year - year_center,
                      log_exposure = 0,
                      country = levels(model_data$country)[1])
    for (name in names(covariates)) out[[name]] <- 0
    out
  }
  high <- brms::posterior_linpred(fit, newdata = frame(prevalence),
                                  re_formula = NA)
  low <- brms::posterior_linpred(
    fit, newdata = frame(rep(AF_REFERENCE, length(prevalence))),
    re_formula = NA)
  pmax(1 - exp(-(high - low)), 0)
}
summarise_af <- function(values) {
  data.frame(af = colMeans(values),
             lo = apply(values, 2, stats::quantile, 0.025),
             hi = apply(values, 2, stats::quantile, 0.975))
}

anchor_rows <- list()
contrast_rows <- list()
curve_rows <- list()
for (name in names(fits)) {
  anchors <- cbind(model = name, prevalence = ANCHORS, year = year_center,
                   summarise_af(af_draws(fits[[name]], ANCHORS, year_center)))
  anchor_rows[[name]] <- anchors
  # The substantive question behind the interaction: does the effect differ
  # between an early and a late year? Zero by construction under A.
  early <- af_draws(fits[[name]], ANCHORS, CONTRAST_YEARS[1])
  late <- af_draws(fits[[name]], ANCHORS, CONTRAST_YEARS[2])
  delta <- late - early
  # Under a model with no time interaction the two years give the same linear
  # predictor and delta is zero up to floating-point noise, so a probability of
  # "positive" would be meaningless; report NA there.
  identically_zero <- apply(abs(delta), 2, max) < 1e-8
  prob_positive <- colMeans(delta > 0)
  prob_positive[identically_zero] <- NA_real_
  contrast_rows[[name]] <- data.frame(
    model = name, prevalence = ANCHORS,
    af_early = colMeans(early), af_late = colMeans(late),
    difference = colMeans(delta),
    lo = apply(delta, 2, stats::quantile, 0.025),
    hi = apply(delta, 2, stats::quantile, 0.975),
    prob_positive = prob_positive,
    stringsAsFactors = FALSE)
  for (year in c(CONTRAST_YEARS[1], year_center, CONTRAST_YEARS[2])) {
    grid <- seq(1, 60, by = 1)
    curve_rows[[paste(name, year)]] <- cbind(
      model = name, label = ladder[[name]]$label, year = year,
      prevalence = grid, summarise_af(af_draws(fits[[name]], grid, year)))
  }
}
anchors <- do.call(rbind, anchor_rows)
contrast <- do.call(rbind, contrast_rows)
curves <- do.call(rbind, curve_rows)
write.csv(anchors, file.path(RESULTS_DIR, "brms_ladder_af_anchors.csv"),
          row.names = FALSE)
write.csv(contrast, file.path(RESULTS_DIR, "brms_ladder_time_contrast.csv"),
          row.names = FALSE)

message("\nAttributable fraction at ", year_center, " (%), 95% credible intervals:")
print(transform(anchors, af = round(100 * af, 1), lo = round(100 * lo, 1),
                hi = round(100 * hi, 1))[, c("model", "prevalence", "af", "lo",
                                             "hi")],
      row.names = FALSE)
message("\nAF in ", CONTRAST_YEARS[2], " minus AF in ", CONTRAST_YEARS[1],
        " (percentage points):")
print(transform(contrast, af_early = round(100 * af_early, 1),
                af_late = round(100 * af_late, 1),
                difference = round(100 * difference, 1),
                lo = round(100 * lo, 1), hi = round(100 * hi, 1),
                prob_positive = round(prob_positive, 2)),
      row.names = FALSE)

## ---- figure ---------------------------------------------------------------------------
curves$label <- factor(curves$label,
                       levels = vapply(ladder, `[[`, "", "label"))
curves$year <- factor(curves$year)
plot <- ggplot2::ggplot(curves, ggplot2::aes(prevalence, 100 * af,
                                             colour = year, fill = year)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * lo, ymax = 100 * hi),
                       alpha = 0.15, colour = NA) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::facet_wrap(~label, nrow = 1) +
  ggplot2::scale_colour_manual(values = c("#8FBEDD", "grey45", "#123F63"),
                               name = "Year") +
  ggplot2::scale_fill_manual(values = c("#8FBEDD", "grey45", "#123F63"),
                             guide = "none") +
  ggplot2::labs(
    x = "MAP PfPR2-10 (%)",
    y = "Malaria-attributable fraction of\npost-neonatal mortality (%)",
    title = "How prevalence and calendar time enter: three models, shared priors",
    subtitle = paste0(
      "Bands are 95% credible intervals. LOO: ",
      paste(sprintf("%s elpd %.1f", loo_table$model, loo_table$elpd_loo),
            collapse = "; "),
      "; stacking weights ",
      paste(sprintf("%s %.2f", loo_table$model, loo_table$stacking_weight),
            collapse = ", "))) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(RESULTS_DIR, "figure17_brms_ladder.png"), plot,
                width = 11, height = 4.4, dpi = 200)
message("\nWrote brms_ladder_loo.csv, brms_ladder_af_anchors.csv, ",
        "brms_ladder_time_contrast.csv and figure17_brms_ladder.png")

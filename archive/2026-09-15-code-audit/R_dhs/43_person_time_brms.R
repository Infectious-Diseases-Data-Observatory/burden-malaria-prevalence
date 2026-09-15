# =============================================================================
# 43_person_time_brms.R — the per-age-band person-time model, fitted in Stan.
#
# The mgcv fits in script 42 are the working results; this refits the smooth
# dose-response model for each of the six age bands (4-month split) in brms so
# the attributable fractions carry full posterior uncertainty, including in the
# smoothing of the prevalence and calendar-year terms and the random-effect
# variances, with the same priors as the region-level Bayesian work (scripts
# 30-37): normal(0, 1) on the standardised covariates, half-t on standard
# deviations.
#
#   deaths_eff ~ s(pfpr10, k = 5) + segment + window + s(year_c) + covariates
#                + (1 | country) + (1 | survey) + offset(log person-months)
#
# Effects are attributable fractions against 0% prevalence at 10, 30 and 50%
# (the dose-response is not linear, so no slope is reported). Reads the
# modelling table script 42 writes; fits are cached by a fingerprint of the
# table; run script 42 first. Six fits take an hour or two.
#
# Outputs (results/dhs_rebuild)
#   person_time_brms_effects.csv     posterior attributable fractions, against mgcv
#   figure24_person_time_brms.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("brms", "ggplot2"))

CHAINS <- 4L
ITERATIONS <- 2000L
WARMUP <- 1000L
# adapt_delta by band: the 4-11 month band had 90 divergent transitions at
# 0.95 (the others 2-31), so it takes a smaller step size; the rest keep 0.95
ADAPT_DELTA <- c(default = 0.95, "4-11 months" = 0.99)
adapt_delta_for <- function(band) if (band %in% names(ADAPT_DELTA)) ADAPT_DELTA[[band]] else ADAPT_DELTA[["default"]]
ANCHORS <- c(10, 30, 50)
BANDS <- AGE6B
CACHE <- file.path(DATA_DIR, "person_time_brms_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
MODEL_DATA_CSV <- file.path(RESULTS_DIR, "person_time_model_data.csv")

bundle <- readRDS(MODEL_BUNDLE_RDS)
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
catalog$included_in_main <- as.logical(catalog$included_in_main)
model_data <- read.csv(MODEL_DATA_CSV, stringsAsFactors = FALSE)
ridge <- make_ridge_matrix(model_data, catalog, bundle$preprocessing)
covariates <- as.data.frame(ridge$matrix)
names(covariates) <- make.names(colnames(ridge$matrix))
model_data <- cbind(model_data[, c("svkey", "regkey", "iso3", "age6b", "segment",
                                   "window", "deaths_eff", "log_pm", "pfpr10", "year_c")],
                    covariates)
model_data$band <- factor(model_data$age6b, levels = BANDS)
model_data$country <- factor(model_data$iso3)
model_data$survey <- factor(model_data$svkey)
model_data$window_f <- factor(model_data$window)

priors <- c(
  brms::prior(normal(0, 1), class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "sds"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(gamma(0.01, 0.01), class = "shape"))

fit_band <- function(band) {
  data <- model_data[model_data$band == band, , drop = FALSE]
  data$segment_f <- droplevels(factor(data$segment))
  terms <- c("s(pfpr10, k = 5)", if (nlevels(data$segment_f) > 1) "segment_f", "window_f",
             "s(year_c, k = 8)", names(covariates), "(1 | country)", "(1 | survey)",
             "offset(log_pm)")
  formula <- paste("deaths_eff ~", paste(terms, collapse = " + "))
  fingerprint <- list(rows = nrow(data), deaths = sum(data$deaths_eff), formula = formula,
                      chains = CHAINS, iterations = ITERATIONS, adapt_delta = adapt_delta_for(band))
  stem <- file.path(CACHE, make.names(band))
  if (file.exists(paste0(stem, ".rds")) && file.exists(paste0(stem, "_meta.rds")) &&
      isTRUE(all.equal(readRDS(paste0(stem, "_meta.rds")), fingerprint))) {
    message("  ", band, ": reusing cached fit")
    return(list(fit = readRDS(paste0(stem, ".rds")), data = data))
  }
  message("  ", band, ": sampling (", nrow(data), " cells, adapt_delta ", adapt_delta_for(band), ")")
  started <- Sys.time()
  fit <- brms::brm(brms::bf(formula), data = data, family = brms::negbinomial(),
                   prior = priors, chains = CHAINS, iter = ITERATIONS, warmup = WARMUP,
                   cores = min(CHAINS, max(1L, parallel::detectCores() - 1L)),
                   control = list(adapt_delta = adapt_delta_for(band), max_treedepth = 12),
                   seed = 20260828, refresh = 200)
  message("  ", band, ": done in ",
          round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), " minutes")
  saveRDS(fit, paste0(stem, ".rds"))
  saveRDS(fingerprint, paste0(stem, "_meta.rds"))
  list(fit = fit, data = data)
}

# posterior draws of the log hazard ratio against the reference for an average
# cell: first segment, window 1, year centred, covariates at their means (zero
# after standardisation), random effects excluded
posterior_log_hr <- function(fit, data, p) {
  segment_levels <- levels(data$segment_f)
  frame <- function(values) {
    out <- data.frame(pfpr10 = values / 10,
                      segment_f = factor(segment_levels[1], levels = segment_levels),
                      window_f = factor("1", levels = levels(model_data$window_f)),
                      year_c = 0, log_pm = 0,
                      country = levels(model_data$country)[1], survey = levels(model_data$survey)[1])
    for (v in names(covariates)) out[[v]] <- 0
    out
  }
  eta_p <- brms::posterior_linpred(fit, newdata = frame(p), re_formula = NA)
  eta_0 <- brms::posterior_linpred(fit, newdata = frame(rep(PERSON_TIME_AF_REFERENCE, length(p))), re_formula = NA)
  eta_p - eta_0
}

rows <- list()
for (band in BANDS) {
  fitted <- fit_band(band)
  fit <- fitted$fit
  lhr <- posterior_log_hr(fit, fitted$data, ANCHORS)
  s <- posterior::summarise_draws(posterior::as_draws_df(fit))
  key <- s[grepl("^b_|^bs_|^sds_|^sd_|^shape", s$variable), ]
  nuts <- brms::nuts_params(fit)
  row <- data.frame(age_band = band, cells = nrow(fit$data), adapt_delta = adapt_delta_for(band),
                    divergent = sum(nuts$Parameter == "divergent__" & nuts$Value > 0),
                    max_rhat = max(key$rhat, na.rm = TRUE),
                    min_bulk_ess = min(key$ess_bulk, na.rm = TRUE), stringsAsFactors = FALSE)
  for (i in seq_along(ANCHORS)) {
    a <- 1 - exp(-lhr[, i])
    row[[paste0("af", ANCHORS[i])]] <- mean(a)
    row[[paste0("af", ANCHORS[i], "_lo")]] <- unname(stats::quantile(a, 0.025))
    row[[paste0("af", ANCHORS[i], "_hi")]] <- unname(stats::quantile(a, 0.975))
    row[[paste0("prob_positive", ANCHORS[i])]] <- mean(a > 0)
  }
  rows[[band]] <- row
}
effects <- do.call(rbind, rows)
rownames(effects) <- NULL
write.csv(effects, file.path(RESULTS_DIR, "person_time_brms_effects.csv"), row.names = FALSE)
message(sprintf("\nPosterior attributable fraction (%%) against %d%% prevalence:", PERSON_TIME_AF_REFERENCE))
fmt <- function(m, lo, hi) sprintf("%.1f (%.1f-%.1f)", 100 * m, 100 * lo, 100 * hi)
print(transform(effects, af10 = fmt(af10, af10_lo, af10_hi), af30 = fmt(af30, af30_lo, af30_hi),
                af50 = fmt(af50, af50_lo, af50_hi), max_rhat = round(max_rhat, 3))[
                  , c("age_band", "af10", "af30", "af50", "prob_positive30", "adapt_delta", "divergent", "max_rhat", "min_bulk_ess")],
      row.names = FALSE)

## ---- against mgcv --------------------------------------------------------------------------
long <- function(x, engine) do.call(rbind, lapply(ANCHORS, function(p) data.frame(
  age_band = x$age_band, anchor = paste0(p, "% prevalence"), engine = engine,
  af = x[[paste0("af", p)]], lo = x[[paste0("af", p, "_lo")]], hi = x[[paste0("af", p, "_hi")]],
  stringsAsFactors = FALSE)))
compare <- long(effects, "brms (posterior mean, 95% credible interval)")
mgcv_file <- file.path(RESULTS_DIR, "person_time_effects.csv")
if (file.exists(mgcv_file)) {
  m <- read.csv(mgcv_file, stringsAsFactors = FALSE)
  m <- m[m$grouping == "six bands (4-month split)" & m$age_group %in% BANDS, ]
  names(m)[names(m) == "age_group"] <- "age_band"
  compare <- rbind(compare, long(m, "mgcv (fREML, 95% confidence interval)"))
}
compare$age_band <- factor(compare$age_band, levels = BANDS)
compare$anchor <- factor(compare$anchor, levels = paste0(ANCHORS, "% prevalence"))
plot <- ggplot2::ggplot(compare, ggplot2::aes(age_band, 100 * af, colour = engine)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = 100 * lo, ymax = 100 * hi), width = 0.2,
                         position = ggplot2::position_dodge(width = 0.4)) +
  ggplot2::geom_point(size = 2.4, position = ggplot2::position_dodge(width = 0.4)) +
  ggplot2::facet_wrap(~anchor, nrow = 1) +
  ggplot2::scale_colour_manual(values = c("#B2182B", "#1D6F8B"), name = NULL) +
  ggplot2::labs(x = NULL, y = sprintf("Attributable fraction of all-cause mortality (%%)\nagainst %d%% prevalence",
                                      PERSON_TIME_AF_REFERENCE),
                title = "Attributable fraction by age band: Stan against mgcv",
                subtitle = "One model per age band with a smooth prevalence effect; person-time design") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
ggplot2::ggsave(file.path(RESULTS_DIR, "figure24_person_time_brms.png"), plot,
                width = 11, height = 4.6, dpi = 200)
message("\nWrote person_time_brms_effects.csv and figure24_person_time_brms.png")

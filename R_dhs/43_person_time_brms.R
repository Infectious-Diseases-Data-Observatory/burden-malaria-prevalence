# =============================================================================
# 43_person_time_brms.R — the per-age-group person-time model, fitted in Stan.
#
# The mgcv fits in script 42 are the working results; this refits the primary
# linear-prevalence model for each of the four age groups in brms so the
# attributable fractions carry full posterior uncertainty, including in the
# smoothing of the calendar-year term and the random-effect variances, with the
# same priors as the region-level Bayesian work (scripts 30-37): normal(0, 1) on
# the standardised covariates, half-t on standard deviations.
#
#   deaths_eff ~ pfpr10 + segment + window + s(year_c) + covariates
#                + (1 | country) + (1 | survey) + offset(log person-months)
#
# Reads the modelling table script 42 writes. Fits are cached by a fingerprint
# of the table; run script 42 first.
#
# Outputs (results/dhs_rebuild)
#   person_time_brms_effects.csv     posterior hazard ratios and attributable fractions
#   figure24_person_time_brms.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("brms", "ggplot2"))

CHAINS <- 4L
ITERATIONS <- 2000L
WARMUP <- 1000L
ADAPT_DELTA <- 0.95
ANCHORS <- c(10, 30, 50)
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
model_data <- cbind(model_data[, c("svkey", "regkey", "iso3", "age_group", "segment",
                                   "window", "deaths_eff", "log_pm", "pfpr10", "year_c")],
                    covariates)
model_data$age_group <- factor(model_data$age_group, levels = AGE_GROUPS)
model_data$country <- factor(model_data$iso3)
model_data$survey <- factor(model_data$svkey)
model_data$window_f <- factor(model_data$window)

priors <- c(
  brms::prior(normal(0, 1), class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "sds"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(gamma(0.01, 0.01), class = "shape"))

fit_group <- function(group) {
  data <- model_data[model_data$age_group == group, , drop = FALSE]
  data$segment_f <- droplevels(factor(data$segment))
  terms <- c("pfpr10", if (nlevels(data$segment_f) > 1) "segment_f", "window_f",
             "s(year_c, k = 8)", names(covariates), "(1 | country)", "(1 | survey)",
             "offset(log_pm)")
  formula <- paste("deaths_eff ~", paste(terms, collapse = " + "))
  fingerprint <- list(rows = nrow(data), deaths = sum(data$deaths_eff), formula = formula,
                      chains = CHAINS, iterations = ITERATIONS, adapt_delta = ADAPT_DELTA)
  stem <- file.path(CACHE, make.names(group))
  if (file.exists(paste0(stem, ".rds")) && file.exists(paste0(stem, "_meta.rds")) &&
      isTRUE(all.equal(readRDS(paste0(stem, "_meta.rds")), fingerprint))) {
    message("  ", group, ": reusing cached fit")
    return(readRDS(paste0(stem, ".rds")))
  }
  message("  ", group, ": sampling (", nrow(data), " cells)")
  started <- Sys.time()
  fit <- brms::brm(brms::bf(formula), data = data, family = brms::negbinomial(),
                   prior = priors, chains = CHAINS, iter = ITERATIONS, warmup = WARMUP,
                   cores = min(CHAINS, max(1L, parallel::detectCores() - 1L)),
                   control = list(adapt_delta = ADAPT_DELTA, max_treedepth = 12),
                   seed = 20260828, refresh = 200)
  message("  ", group, ": done in ",
          round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), " minutes")
  saveRDS(fit, paste0(stem, ".rds"))
  saveRDS(fingerprint, paste0(stem, "_meta.rds"))
  fit
}

rows <- list()
for (group in AGE_GROUPS) {
  fit <- fit_group(group)
  draws <- posterior::as_draws_df(fit)
  beta <- draws$b_pfpr10
  s <- posterior::summarise_draws(draws)
  key <- s[grepl("^b_|^bs_|^sds_|^sd_|^shape", s$variable), ]
  nuts <- brms::nuts_params(fit)
  row <- data.frame(age_group = group, cells = nrow(fit$data),
                    pct_per_10 = 100 * (exp(mean(beta)) - 1),
                    pct_lo = 100 * (exp(stats::quantile(beta, 0.025)) - 1),
                    pct_hi = 100 * (exp(stats::quantile(beta, 0.975)) - 1),
                    prob_positive = mean(beta > 0),
                    divergent = sum(nuts$Parameter == "divergent__" & nuts$Value > 0),
                    max_rhat = max(key$rhat, na.rm = TRUE),
                    min_bulk_ess = min(key$ess_bulk, na.rm = TRUE), stringsAsFactors = FALSE)
  for (p in ANCHORS) {
    a <- 1 - exp(-beta * (p - AF_REFERENCE) / 10)
    row[[paste0("af", p)]] <- mean(a)
    row[[paste0("af", p, "_lo")]] <- unname(stats::quantile(a, 0.025))
    row[[paste0("af", p, "_hi")]] <- unname(stats::quantile(a, 0.975))
  }
  rows[[group]] <- row
}
effects <- do.call(rbind, rows)
rownames(effects) <- NULL
write.csv(effects, file.path(RESULTS_DIR, "person_time_brms_effects.csv"), row.names = FALSE)
message("\nPosterior hazard change per +10 PfPR points and attributable fractions:")
print(transform(effects, pct = sprintf("%+.1f (%+.1f to %+.1f)", pct_per_10, pct_lo, pct_hi),
                af10 = sprintf("%.1f (%.1f-%.1f)", 100 * af10, 100 * af10_lo, 100 * af10_hi),
                af30 = sprintf("%.1f (%.1f-%.1f)", 100 * af30, 100 * af30_lo, 100 * af30_hi),
                af50 = sprintf("%.1f (%.1f-%.1f)", 100 * af50, 100 * af50_lo, 100 * af50_hi),
                max_rhat = round(max_rhat, 3))[
                  , c("age_group", "pct", "prob_positive", "af10", "af30", "af50",
                      "divergent", "max_rhat", "min_bulk_ess")], row.names = FALSE)

mgcv_file <- file.path(RESULTS_DIR, "person_time_effects.csv")
compare <- effects[, c("age_group", "pct_per_10", "pct_lo", "pct_hi")]
compare$engine <- "brms (posterior mean, 95% credible interval)"
if (file.exists(mgcv_file)) {
  m <- read.csv(mgcv_file, stringsAsFactors = FALSE)
  m <- m[m$grouping == "four groups" & !is.na(m$beta_per_10), ]
  compare <- rbind(compare, data.frame(
    age_group = m$age_group, pct_per_10 = 100 * (exp(m$beta_per_10) - 1),
    pct_lo = 100 * (exp(m$beta_per_10 - 1.96 * m$se) - 1),
    pct_hi = 100 * (exp(m$beta_per_10 + 1.96 * m$se) - 1),
    engine = "mgcv (fREML, 95% confidence interval)", stringsAsFactors = FALSE))
}
compare$age_group <- factor(compare$age_group, levels = AGE_GROUPS)
plot <- ggplot2::ggplot(compare, ggplot2::aes(age_group, pct_per_10, colour = engine)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pct_lo, ymax = pct_hi), width = 0.2,
                         position = ggplot2::position_dodge(width = 0.4)) +
  ggplot2::geom_point(size = 2.8, position = ggplot2::position_dodge(width = 0.4)) +
  ggplot2::scale_colour_manual(values = c("#B2182B", "#1D6F8B"), name = NULL) +
  ggplot2::labs(x = NULL, y = "Change in mortality hazard per +10 PfPR2-10 points (%)",
                title = "Prevalence effect by age group: Stan against mgcv",
                subtitle = "One model per age group; person-time design") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(RESULTS_DIR, "figure24_person_time_brms.png"), plot,
                width = 8, height = 4.4, dpi = 200)
message("\nWrote person_time_brms_effects.csv and figure24_person_time_brms.png")

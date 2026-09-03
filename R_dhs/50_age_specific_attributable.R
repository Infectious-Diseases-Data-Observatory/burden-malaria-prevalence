# =============================================================================
# 50_age_specific_attributable.R — all-cause mortality attributable to malaria
# by age band as a function of PfPR2-10, from the six age-band person-time
# models (smooth dose-response, 4-month split).
#
# The figure is drawn from the Stan fits of script 43 (posterior mean and 95%
# credible interval from the full posterior, so the smoothing and variance
# parameters carry their uncertainty). The mgcv fits of script 42 give the same
# curves from 1,000 draws of the coefficient sampling distribution and are kept
# in the table under their own engine label; the figure falls back to them if
# the Stan cache is absent.
#
# Two quantities against 0% prevalence (PERSON_TIME_AF_REFERENCE):
#   attributable fraction  AF_a(p) = 1 - rate_a(0) / rate_a(p): the share of
#       all-cause deaths in band a that the model attributes to prevalence p,
#       holding covariates, calendar year and the random effects at their
#       reference values.
#   attributable rate      rate_a(p) - rate_a(0), in deaths per 1,000
#       child-years, where rate_a(p) is the model's predicted all-cause death
#       rate in band a at prevalence p for an average cell (covariates at their
#       means, year centred, first window, random effects at zero). Segments
#       within a band are weighted by their observed person-time. Tabulated only:
#       for the neonatal band it is a null association times a very high
#       baseline rate and its interval is uninformative.
#
# Run scripts 42 and 43 first.
#
# Outputs (results/dhs_rebuild)
#   person_time_age_specific_attributable.csv   engine x band x PfPR grid: AF and attributable rate
#   person_time_age_specific_anchors.csv        the same at 10, 30, 50% prevalence
#   figure36_age_specific_attributable.png      AF curves from the Stan fits, one panel per band
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "MASS", "ggplot2"))

N_DRAWS <- 1000L
GRID <- seq(0, 80, by = 1)
ANCHORS <- c(10, 30, 50)
BRMS_CACHE <- file.path(DATA_DIR, "person_time_brms_cache")
set.seed(20260903)

bundle <- readRDS(file.path(DERIVED_DIR, "person_time_model_bundle.rds"))
fits <- bundle$bands6b_smooth_fits
model_data <- read.csv(file.path(RESULTS_DIR, "person_time_model_data.csv"), stringsAsFactors = FALSE)
country_levels <- levels(factor(model_data$iso3))
survey_levels <- levels(factor(model_data$svkey))

stan_path <- function(band) file.path(BRMS_CACHE, paste0(make.names(band), ".rds"))
have_stan <- all(file.exists(vapply(AGE6B, stan_path, character(1))))
if (have_stan) required_packages(c("brms", "posterior")) else
  message("Stan fits not found in ", BRMS_CACHE, "; run script 43. The figure will use the mgcv fits.")

segment_shares <- function(band) {
  band_data <- model_data[model_data$age6b == band, ]
  levels_s <- levels(droplevels(factor(band_data$segment)))
  share <- tapply(band_data$person_months_n, factor(band_data$segment, levels = levels_s), sum)
  list(levels = levels_s, share = share / sum(share))
}
# summarise draws (rows) x grid (columns) of the band rate and its reference
summarise_rates <- function(band, engine, rate_grid, rate_ref) {
  af <- 1 - rate_ref / rate_grid
  attributable <- (rate_grid - rate_ref) * 12000   # per 1,000 child-years
  rate <- rate_grid * 12000
  q <- function(m, f) apply(m, 2, f)
  lo <- function(v) stats::quantile(v, 0.025); hi <- function(v) stats::quantile(v, 0.975)
  data.frame(engine = engine, age_band = band, pfpr = GRID,
             af = colMeans(af), af_lo = q(af, lo), af_hi = q(af, hi),
             attributable_per_1000cy = colMeans(attributable),
             attributable_lo = q(attributable, lo), attributable_hi = q(attributable, hi),
             allcause_rate_per_1000cy = colMeans(rate), rate_lo = q(rate, lo), rate_hi = q(rate, hi),
             reference_rate_per_1000cy = mean(rate_ref) * 12000,
             stringsAsFactors = FALSE)
}

## ---- mgcv: 1,000 draws of the coefficients -----------------------------------------------------
band_curves_mgcv <- function(band) {
  fit <- fits[[band]]
  segments <- segment_shares(band)
  frame <- function(p, segment) {
    out <- data.frame(pfpr10 = p / 10,
                      segment_f = factor(segment, levels = segments$levels),
                      window_f = factor("1", levels = as.character(1:5)), year_c = 0, log_pm = 0,
                      country = factor(country_levels[1], levels = country_levels),
                      survey = factor(survey_levels[1], levels = survey_levels))
    G <- matrix(0, nrow(out), length(grep("^G", names(coef(fit)))))
    colnames(G) <- sub("^G", "", grep("^G", names(coef(fit)), value = TRUE))
    out$G <- G
    out
  }
  draws <- MASS::mvrnorm(N_DRAWS, coef(fit), vcov(fit))
  zero_re <- function(X) { X[, grep("^s\\(country\\)|^s\\(survey\\)", colnames(X))] <- 0; X }
  rate_grid <- 0; rate_ref <- 0
  for (s in seq_along(segments$levels)) {
    Xg <- zero_re(predict(fit, frame(GRID, segments$levels[s]), type = "lpmatrix", discrete = FALSE))
    Xr <- zero_re(predict(fit, frame(PERSON_TIME_AF_REFERENCE, segments$levels[s]), type = "lpmatrix", discrete = FALSE))
    rate_grid <- rate_grid + segments$share[[s]] * exp(draws %*% t(Xg))
    rate_ref <- rate_ref + segments$share[[s]] * exp(draws %*% t(Xr))
  }
  summarise_rates(band, "mgcv (fREML, coefficient draws)", rate_grid, as.numeric(rate_ref))
}

## ---- Stan: the posterior --------------------------------------------------------------------------------
band_curves_stan <- function(band) {
  fit <- readRDS(stan_path(band))
  segments <- segment_shares(band)
  cells <- sum(model_data$age6b == band)
  if (nrow(fit$data) != cells) warning(band, ": the Stan fit has ", nrow(fit$data), " cells, the modelling table ", cells,
                                       "; rerun script 43")
  fixed <- c("deaths_eff", "pfpr10", "segment_f", "window_f", "year_c", "country", "survey", "log_pm")
  covariate_names <- setdiff(names(fit$data), fixed)
  # a band with a single segment (the neonatal month) has no segment term in
  # the Stan model, so its data carry no segment_f column
  has_segment <- "segment_f" %in% names(fit$data)
  frame <- function(p, segment) {
    out <- data.frame(pfpr10 = p / 10,
                      window_f = factor("1", levels = levels(fit$data$window_f)),
                      year_c = 0, log_pm = 0,
                      country = levels(fit$data$country)[1], survey = levels(fit$data$survey)[1],
                      stringsAsFactors = FALSE)
    if (has_segment) out$segment_f <- factor(segment, levels = levels(fit$data$segment_f))
    for (v in covariate_names) out[[v]] <- 0
    out
  }
  rate_grid <- 0; rate_ref <- 0
  for (s in seq_along(segments$levels)) {
    segment <- segments$levels[s]
    if (has_segment && !segment %in% levels(fit$data$segment_f)) stop(band, ": segment ", segment, " not in the Stan fit")
    eta_grid <- brms::posterior_linpred(fit, newdata = frame(GRID, segment), re_formula = NA)
    eta_ref <- brms::posterior_linpred(fit, newdata = frame(PERSON_TIME_AF_REFERENCE, segment), re_formula = NA)
    rate_grid <- rate_grid + segments$share[[s]] * exp(eta_grid)
    rate_ref <- rate_ref + segments$share[[s]] * exp(as.numeric(eta_ref))
  }
  summarise_rates(band, "Stan (brms, full posterior)", rate_grid, as.numeric(rate_ref))
}

curves <- do.call(rbind, lapply(AGE6B, band_curves_mgcv))
if (have_stan) {
  message("Stan curves from the cached brms fits")
  curves <- rbind(do.call(rbind, lapply(AGE6B, band_curves_stan)), curves)
}
rownames(curves) <- NULL
curves$age_band <- factor(curves$age_band, levels = AGE6B)
write.csv(curves, file.path(RESULTS_DIR, "person_time_age_specific_attributable.csv"), row.names = FALSE)

anchors <- curves[curves$pfpr %in% ANCHORS, ]
write.csv(anchors, file.path(RESULTS_DIR, "person_time_age_specific_anchors.csv"), row.names = FALSE)
main_engine <- if (have_stan) "Stan (brms, full posterior)" else "mgcv (fREML, coefficient draws)"
message(sprintf("Attributable fraction of all-cause mortality (%%) against %d%% prevalence and attributable deaths per 1,000 child-years, by age band (%s):",
                PERSON_TIME_AF_REFERENCE, main_engine))
show <- anchors[anchors$engine == main_engine, ]
print(transform(show,
                af = sprintf("%.1f (%.1f-%.1f)", 100 * af, 100 * af_lo, 100 * af_hi),
                attributable = sprintf("%.1f (%.1f-%.1f)", attributable_per_1000cy, attributable_lo, attributable_hi),
                allcause = sprintf("%.0f", allcause_rate_per_1000cy),
                reference = sprintf("%.0f", reference_rate_per_1000cy))[
                  , c("age_band", "pfpr", "af", "attributable", "allcause", "reference")],
      row.names = FALSE, right = FALSE)
if (have_stan) {
  message("\nStan against mgcv, AF at 30% (%):")
  a30 <- anchors[anchors$pfpr == 30, ]
  print(reshape(transform(a30, est = sprintf("%.1f (%.1f-%.1f)", 100 * af, 100 * af_lo, 100 * af_hi),
                          engine = ifelse(grepl("^Stan", engine), "stan", "mgcv"))[, c("age_band", "engine", "est")],
                idvar = "age_band", timevar = "engine", direction = "wide"), row.names = FALSE, right = FALSE)
}

## ---- figure: one panel per band, two rows of three -----------------------------------------------
palette <- c("#7570B3", "#1B9E77", "#66A61E", "#E6AB02", "#D95F02", "#A6761D")
names(palette) <- AGE6B
observed <- unique(model_data[, c("svkey", "regkey", "window", "pfpr")])
drawn <- curves[curves$engine == main_engine, ]
plot <- ggplot2::ggplot(drawn, ggplot2::aes(pfpr, 100 * af)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * af_lo, ymax = 100 * af_hi, fill = age_band), alpha = 0.25, colour = NA) +
  ggplot2::geom_line(ggplot2::aes(colour = age_band), linewidth = 1) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55") +
  ggplot2::geom_rug(data = observed[observed$pfpr <= max(GRID), ], ggplot2::aes(x = pfpr),
                    inherit.aes = FALSE, alpha = 0.03, sides = "b") +
  ggplot2::facet_wrap(~ age_band, nrow = 2) +
  ggplot2::scale_colour_manual(values = palette, guide = "none") +
  ggplot2::scale_fill_manual(values = palette, guide = "none") +
  ggplot2::labs(x = expression(PfPR[2-10]~"(%)"),
                y = sprintf("Share of all-cause mortality attributable to malaria (%%), against %d%% prevalence",
                            PERSON_TIME_AF_REFERENCE),
                title = "Malaria-attributable fraction of all-cause mortality by age",
                subtitle = if (have_stan)
                  "Stan fit of one negative-binomial model per age band with a smooth prevalence effect.\nPosterior mean and 95% credible interval; rug shows region-window prevalence."
                else
                  "mgcv fit of one negative-binomial model per age band with a smooth prevalence effect.\nRibbons are 95% intervals; rug shows region-window prevalence.") +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(RESULTS_DIR, "figure36_age_specific_attributable.png"), plot,
                width = 11, height = 7, dpi = 200, bg = "white")
message("Wrote person_time_age_specific_attributable.csv, person_time_age_specific_anchors.csv and figure36_age_specific_attributable.png")

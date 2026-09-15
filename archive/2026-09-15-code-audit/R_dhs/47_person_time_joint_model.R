# =============================================================================
# 47_person_time_joint_model.R — one model across the six age bands, with every
# nuisance term allowed to differ by band.
#
# Script 42 fits each band separately. This fits all bands together so the
# prevalence effects can be compared formally, while keeping what made the
# separate fits necessary: nothing that describes a region's level is shared
# across ages. Band-specific terms:
#
#   prevalence         s(pfpr10, by = band, k = 5)
#   baseline hazard    segment                  (segments nest within bands)
#   window (recall)    one dummy per band x window 2-5, unpenalised
#   calendar year      s(year_c, by = band)
#   covariates         six ridge blocks G_b = G x 1[band = b], each with its own
#                      penalty (paraPen), so each band has its own shrinkage
#   country intercept  s(country, by = band, bs = "re")   own variance per band
#   survey intercept   s(survey,  by = band, bs = "re")   own variance per band
#
# The one shared parameter is the negative-binomial dispersion theta, which
# mgcv::nb() cannot make band-specific; the separate fits gave theta from about
# 4 to 16 across bands. Comparisons against the separate fits show what that
# costs. Whether the dose-response differs by age is judged by AIC (fREML fits)
# between the band-specific smooths and a single smooth common to all bands.
# Effects are reported as attributable fractions against 0% prevalence at 10,
# 30 and 50%, not as slopes: the response is not linear.
#
# Outputs (results/dhs_rebuild)
#   person_time_joint_effects.csv     AF at 10/30/50%: joint against separate, by band
#   person_time_joint_comparison.csv  AIC: common vs band-specific dose-response
#   figure31_person_time_joint_vs_separate.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2"))

BAM_THREADS <- max(1L, min(8L, parallel::detectCores() - 2L))
BANDS <- AGE6B
ANCHORS <- c(10, 30, 50)

bundle <- readRDS(MODEL_BUNDLE_RDS)
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
catalog$included_in_main <- as.logical(catalog$included_in_main)
separate <- readRDS(file.path(DERIVED_DIR, "person_time_model_bundle.rds"))

d <- read.csv(file.path(RESULTS_DIR, "person_time_model_data.csv"), stringsAsFactors = FALSE)
d$band <- factor(d$age6b, levels = BANDS)
d$segment_f <- factor(d$segment)
d$window_f <- factor(d$window)
d$country <- factor(d$iso3)
d$survey <- factor(d$svkey)
G <- make_ridge_matrix(d, catalog, bundle$preprocessing)$matrix

# band-specific ridge blocks and window dummies
penalty <- list()
for (b in seq_along(BANDS)) {
  name <- paste0("G", b)
  d[[name]] <- G * (d$band == BANDS[b])
  penalty[[name]] <- list(diag(ncol(G)))
}
W <- do.call(cbind, lapply(seq_along(BANDS), function(b)
  sapply(2:5, function(w) as.numeric(d$band == BANDS[b] & d$window == w))))
colnames(W) <- paste0("b", rep(seq_along(BANDS), each = 4), "_w", rep(2:5, length(BANDS)))
d$W <- W
message("Joint model data: ", nrow(d), " cells, ", nlevels(d$band), " bands, ",
        ncol(G), " covariates x ", nlevels(d$band), " ridge blocks")

base <- paste(
  "segment_f + W + s(year_c, by = band, k = 8) +",
  paste(paste0("G", seq_along(BANDS)), collapse = " + "),
  "+ s(country, by = band, bs = 're') + s(survey, by = band, bs = 're') + offset(log_pm)")
fit_joint <- function(terms) {
  started <- Sys.time()
  f <- as.formula(paste("deaths_eff ~", terms, "+", base))
  fit <- mgcv::bam(f, family = mgcv::nb(), method = "fREML", paraPen = penalty, data = d,
                   discrete = TRUE, nthreads = BAM_THREADS)
  message(sprintf("  fitted %-28s in %.0f s, theta %.2f", terms,
                  as.numeric(difftime(Sys.time(), started, units = "secs")),
                  fit$family$getTheta(TRUE)))
  fit
}

message("\nFitting")
joint_smooth <- fit_joint("s(pfpr10, by = band, k = 5)")
joint_common <- fit_joint("s(pfpr10, k = 5)")

## ---- attributable fractions: joint against separate ---------------------------------------
# log hazard ratio of prevalence p against the reference in band b, for an
# average cell of that band (its first segment, window 1, year centred,
# covariates at their means, random effects at zero)
joint_lhr <- function(fit, band, p) {
  segment <- levels(droplevels(factor(d$segment[d$band == band])))[1]
  frame <- function(values) {
    out <- data.frame(pfpr10 = values / 10, band = factor(band, levels = BANDS),
                      segment_f = factor(segment, levels = levels(d$segment_f)),
                      year_c = 0, log_pm = 0,
                      country = factor(levels(d$country)[1], levels = levels(d$country)),
                      survey = factor(levels(d$survey)[1], levels = levels(d$survey)))
    out$W <- matrix(0, nrow(out), ncol(W), dimnames = list(NULL, colnames(W)))
    for (b in seq_along(BANDS)) out[[paste0("G", b)]] <- matrix(0, nrow(out), ncol(G), dimnames = list(NULL, colnames(G)))
    out
  }
  Xh <- predict(fit, frame(p), type = "lpmatrix", discrete = FALSE)
  Xl <- predict(fit, frame(rep(PERSON_TIME_AF_REFERENCE, length(p))), type = "lpmatrix", discrete = FALSE)
  re <- grep("^s\\(country\\)|^s\\(survey\\)", colnames(Xh))
  Xh[, re] <- 0; Xl[, re] <- 0
  dX <- Xh - Xl
  data.frame(pfpr = p, est = as.numeric(dX %*% coef(fit)), se = sqrt(rowSums((dX %*% vcov(fit)) * dX)))
}
separate_lhr <- function(fit, band, p) {
  band_data <- d[d$band == band, ]
  segment_levels <- levels(droplevels(factor(band_data$segment)))
  frame <- function(values) {
    out <- data.frame(pfpr10 = values / 10,
                      segment_f = factor(segment_levels[1], levels = segment_levels),
                      window_f = factor("1", levels = levels(d$window_f)), year_c = 0, log_pm = 0,
                      country = factor(levels(d$country)[1], levels = levels(d$country)),
                      survey = factor(levels(d$survey)[1], levels = levels(d$survey)))
    out$G <- matrix(0, nrow(out), ncol(G), dimnames = list(NULL, colnames(G)))
    out
  }
  Xh <- predict(fit, frame(p), type = "lpmatrix", discrete = FALSE)
  Xl <- predict(fit, frame(rep(PERSON_TIME_AF_REFERENCE, length(p))), type = "lpmatrix", discrete = FALSE)
  re <- grep("^s\\(country\\)|^s\\(survey\\)", colnames(Xh))
  Xh[, re] <- 0; Xl[, re] <- 0
  dX <- Xh - Xl
  data.frame(pfpr = p, est = as.numeric(dX %*% coef(fit)), se = sqrt(rowSums((dX %*% vcov(fit)) * dX)))
}
af_rows <- function(lhr, band, engine) data.frame(
  age_band = band, engine = engine, pfpr = lhr$pfpr,
  af = 1 - exp(-lhr$est), af_lo = 1 - exp(-(lhr$est - 1.96 * lhr$se)),
  af_hi = 1 - exp(-(lhr$est + 1.96 * lhr$se)), stringsAsFactors = FALSE)
effects <- do.call(rbind, lapply(BANDS, function(b) rbind(
  af_rows(joint_lhr(joint_smooth, b, ANCHORS), b, "joint model, band-specific nuisance terms"),
  af_rows(separate_lhr(separate$bands6b_smooth_fits[[b]], b, ANCHORS), b, "separate model per band"))))
effects$age_band <- factor(effects$age_band, levels = BANDS)
effects <- effects[order(effects$age_band, effects$pfpr, effects$engine), ]
rownames(effects) <- NULL
write.csv(effects, file.path(RESULTS_DIR, "person_time_joint_effects.csv"), row.names = FALSE)

comparison <- data.frame(
  model = c("common dose-response, all bands", "band-specific dose-response"),
  AIC = c(stats::AIC(joint_common), stats::AIC(joint_smooth)),
  edf = c(sum(joint_common$edf), sum(joint_smooth$edf)),
  theta = c(joint_common$family$getTheta(TRUE), joint_smooth$family$getTheta(TRUE)),
  stringsAsFactors = FALSE)
comparison$dAIC <- comparison$AIC - min(comparison$AIC)
write.csv(comparison, file.path(RESULTS_DIR, "person_time_joint_comparison.csv"), row.names = FALSE)

message(sprintf("\nAttributable fraction (%%) against %d%% prevalence, joint against separate:", PERSON_TIME_AF_REFERENCE))
show <- transform(effects, est = sprintf("%.1f (%.1f to %.1f)", 100 * af, 100 * af_lo, 100 * af_hi),
                  column = paste0(ifelse(grepl("^joint", engine), "joint", "separate"), "_", pfpr))
print(reshape(show[, c("age_band", "column", "est")], idvar = "age_band", timevar = "column",
              direction = "wide"), row.names = FALSE, right = FALSE)
message("\nCommon against band-specific dose-response (AIC, fREML fits):")
print(transform(comparison, AIC = round(AIC, 1), edf = round(edf, 1), theta = round(theta, 2),
                dAIC = round(dAIC, 1)), row.names = FALSE)
message("\nBand-specific variance components (SD of the random intercepts) in the joint model:")
# For a "re" smooth the penalty is the identity, so with the NB scale fixed at 1
# the random-intercept standard deviation is 1 / sqrt(smoothing parameter).
sp <- joint_smooth$sp
re <- sp[grep("country|survey", names(sp))]
print(data.frame(term = names(re), sd = round(1 / sqrt(re), 3)), row.names = FALSE)

## ---- figure -------------------------------------------------------------------------------------
effects$anchor <- factor(paste0(effects$pfpr, "% prevalence"), levels = paste0(ANCHORS, "% prevalence"))
plot <- ggplot2::ggplot(effects, ggplot2::aes(age_band, 100 * af, colour = engine)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = 100 * af_lo, ymax = 100 * af_hi), width = 0.2,
                         position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::geom_point(size = 2.4, position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::facet_wrap(~anchor, nrow = 1) +
  ggplot2::scale_colour_manual(values = c("#B2182B", "#1D6F8B"), name = NULL) +
  ggplot2::labs(x = NULL, y = sprintf("Attributable fraction of all-cause mortality (%%)\nagainst %d%% prevalence",
                                      PERSON_TIME_AF_REFERENCE),
                title = "Attributable fraction by age band: joint model against separate fits",
                subtitle = sprintf(paste0("Joint model: band-specific smooth prevalence, segment, window, calendar, covariate-ridge and\n",
                                          "random-intercept terms, one shared NB dispersion. ",
                                          "AIC favours band-specific over a common dose-response by %.0f"),
                                   comparison$AIC[1] - comparison$AIC[2])) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
ggplot2::ggsave(file.path(RESULTS_DIR, "figure31_person_time_joint_vs_separate.png"), plot,
                width = 11, height = 4.8, dpi = 200)
saveRDS(list(joint_smooth = joint_smooth, joint_common = joint_common),
        file.path(DERIVED_DIR, "person_time_joint_bundle.rds"))
message("\nWrote person_time_joint_effects.csv, person_time_joint_comparison.csv and ",
        "figure31_person_time_joint_vs_separate.png")

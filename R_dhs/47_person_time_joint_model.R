# =============================================================================
# 47_person_time_joint_model.R — one model across the six age bands, with every
# nuisance term allowed to differ by band.
#
# Script 42 fits each band separately. This fits all bands together so the
# prevalence effects can be compared formally, while keeping what made the
# separate fits necessary: nothing that describes a region's level is shared
# across ages. Band-specific terms:
#
#   prevalence         pfpr10 : band            (or s(pfpr10, by = band))
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
# costs. The band-specific structure is compared with a common prevalence
# effect by AIC (fREML fits) and by a Wald test of slope equality.
#
# Outputs (results/dhs_rebuild)
#   person_time_joint_effects.csv     joint against separate estimates, by band
#   person_time_joint_comparison.csv  AIC: common vs band-specific prevalence effect
#   figure31_person_time_joint_vs_separate.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2"))

BAM_THREADS <- max(1L, min(8L, parallel::detectCores() - 2L))
BANDS <- AGE6B

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
joint_linear <- fit_joint("pfpr10:band")
joint_common <- fit_joint("pfpr10")
joint_smooth <- fit_joint("s(pfpr10, by = band, k = 5)")
joint_common_smooth <- fit_joint("s(pfpr10, k = 5)")

## ---- prevalence effects: joint against separate ---------------------------------------
tab <- summary(joint_linear)$p.table
rows <- grep("^pfpr10:band", rownames(tab))
joint <- data.frame(
  age_band = sub("^pfpr10:band", "", rownames(tab)[rows]),
  engine = "joint model, band-specific nuisance terms",
  beta = tab[rows, "Estimate"], se = tab[rows, "Std. Error"], stringsAsFactors = FALSE)
sep <- do.call(rbind, lapply(BANDS, function(b) {
  s <- summary(separate$bands6b_fits[[b]])$p.table
  data.frame(age_band = b, engine = "separate model per band",
             beta = s["pfpr10", "Estimate"], se = s["pfpr10", "Std. Error"],
             stringsAsFactors = FALSE)
}))
effects <- rbind(joint, sep)
effects$pct_per_10 <- 100 * (exp(effects$beta) - 1)
effects$pct_lo <- 100 * (exp(effects$beta - 1.96 * effects$se) - 1)
effects$pct_hi <- 100 * (exp(effects$beta + 1.96 * effects$se) - 1)
effects$age_band <- factor(effects$age_band, levels = BANDS)
effects <- effects[order(effects$age_band, effects$engine), ]
rownames(effects) <- NULL
write.csv(effects, file.path(RESULTS_DIR, "person_time_joint_effects.csv"), row.names = FALSE)

# Wald test that the six slopes are equal
V <- vcov(joint_linear)[rows, rows]
beta <- tab[rows, "Estimate"]
C <- cbind(-1, diag(length(beta) - 1))          # contrasts against the first band
wald <- as.numeric(t(C %*% beta) %*% solve(C %*% V %*% t(C)) %*% (C %*% beta))
wald_p <- stats::pchisq(wald, df = length(beta) - 1, lower.tail = FALSE)

comparison <- data.frame(
  model = c("common linear slope", "band-specific linear slopes",
            "common smooth", "band-specific smooths"),
  AIC = c(stats::AIC(joint_common), stats::AIC(joint_linear),
          stats::AIC(joint_common_smooth), stats::AIC(joint_smooth)),
  edf = c(sum(joint_common$edf), sum(joint_linear$edf),
          sum(joint_common_smooth$edf), sum(joint_smooth$edf)),
  theta = c(joint_common$family$getTheta(TRUE), joint_linear$family$getTheta(TRUE),
            joint_common_smooth$family$getTheta(TRUE), joint_smooth$family$getTheta(TRUE)),
  stringsAsFactors = FALSE)
comparison$dAIC <- comparison$AIC - min(comparison$AIC)
write.csv(comparison, file.path(RESULTS_DIR, "person_time_joint_comparison.csv"), row.names = FALSE)

message("\nPrevalence effect per +10 PfPR points (% change in hazard), joint against separate:")
show <- transform(effects, est = sprintf("%+.1f (%+.1f to %+.1f)", pct_per_10, pct_lo, pct_hi))
print(reshape(show[, c("age_band", "engine", "est")], idvar = "age_band", timevar = "engine",
              direction = "wide"), row.names = FALSE, right = FALSE)
message(sprintf("\nWald test of equal slopes across the six bands: chi-square %.1f on %d df, p = %.2g",
                wald, length(beta) - 1, wald_p))
message("\nCommon against band-specific prevalence effects (AIC, fREML fits):")
print(transform(comparison, AIC = round(AIC, 1), edf = round(edf, 1), theta = round(theta, 2),
                dAIC = round(dAIC, 1)), row.names = FALSE)
message("\nBand-specific variance components (SD of the random intercepts) in the joint linear model:")
# For a "re" smooth the penalty is the identity, so with the NB scale fixed at 1
# the random-intercept standard deviation is 1 / sqrt(smoothing parameter).
sp <- joint_linear$sp
re <- sp[grep("country|survey", names(sp))]
print(data.frame(term = names(re), sd = round(1 / sqrt(re), 3)), row.names = FALSE)
cp <- summary(joint_common)$p.table
message(sprintf("Common slope: %+.1f%% (%+.1f to %+.1f) per 10 points",
                100 * (exp(cp["pfpr10", 1]) - 1), 100 * (exp(cp["pfpr10", 1] - 1.96 * cp["pfpr10", 2]) - 1),
                100 * (exp(cp["pfpr10", 1] + 1.96 * cp["pfpr10", 2]) - 1)))

## ---- figure -------------------------------------------------------------------------------------
plot <- ggplot2::ggplot(effects, ggplot2::aes(age_band, pct_per_10, colour = engine)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pct_lo, ymax = pct_hi), width = 0.2,
                         position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::geom_point(size = 2.6, position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::scale_colour_manual(values = c("#B2182B", "#1D6F8B"), name = NULL) +
  ggplot2::labs(x = NULL, y = "Change in mortality hazard per +10 PfPR2-10 points (%)",
                title = "Prevalence effect by age band: joint model against separate fits",
                subtitle = sprintf(paste("Joint model: band-specific segment, window, calendar, covariate-ridge and",
                                         "random-intercept terms, one shared NB dispersion.\nWald test of equal",
                                         "slopes: chi-square %.1f on %d df, p = %.2g"), wald, length(beta) - 1, wald_p)) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(RESULTS_DIR, "figure31_person_time_joint_vs_separate.png"), plot,
                width = 9, height = 4.8, dpi = 200)
saveRDS(list(joint_linear = joint_linear, joint_smooth = joint_smooth),
        file.path(DERIVED_DIR, "person_time_joint_bundle.rds"))
message("\nWrote person_time_joint_effects.csv, person_time_joint_comparison.csv and ",
        "figure31_person_time_joint_vs_separate.png")

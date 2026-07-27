# =============================================================================
# 12_triangulation_rct.R — RCT triangulation of the prevalence -> mortality link
# (Ported from legacy R/10_triangulation_rct.R; now uses the R_dhs model bundle
# rather than the legacy Component-2 fit.)
#
# Premise: the DHS/MIS model predicts the all-cause post-neonatal child MORTALITY
# effect from the PREVALENCE effect, so a larger drop in parasite prevalence
# should map to a larger mortality reduction. We test this against cluster-RCTs
# of insecticide-treated nets/curtains that measured BOTH prevalence and
# mortality. The primary comparison uses the ridge-LINEAR post-neonatal slope
# from the bundle (per +10 PfPR2-10 points); the selected additive spline is
# retained as a secondary reference.
#
# Trial data hand-extracted from the original papers (+ Cochrane where noted);
# see per-row provenance. Mortality is UNDER-5, neonatal-excluded, matching the
# post-neonatal outcome. Exploratory (8 heterogeneous points), not powered.
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("malariaAtlas", "mgcv", "ggplot2", "patchwork"))

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
beta_lin <- unname(coef(bundle$linear_fits$postneonatal$model)["pfpr10"])
spline_model <- bundle$primary_fits$postneonatal$model

## ---- extracted trial/zone data ----------------------------------------------
# prev_* = P. falciparum prevalence (%) by arm at follow-up (control = counterfactual).
rct <- data.frame(
  study   = c(paste0("D'Alessandro z", 1:5), "Habluetzel 1997", "Phillips-Howard 2003", "Nevill 1996", "Binka 1996"),
  country = c(rep("Gambia", 5), "Burkina Faso", "Kenya", "Kenya", "Ghana"),
  comparison = c(rep("untreated net", 5), "no curtain", "no net", "no net", "no net"),
  prev_treated = c(28.2, 22.8, 15.9, 43.2, 71.0, 83, 54, 11.9, NA),
  prev_control = c(36.7, 33.0, 25.8, 53.3, 44.8, 91, 66, 25.1, NA),
  prev_amin    = c(rep(1, 5), 0.5, 0.25, 0.08, NA),
  prev_amax    = c(rep(5, 5), 5, 3, 0.92, NA),
  d_t  = c(4, 13, 18, 35, 66, NA, NA, NA, NA), cy_t = c(1056, 1472, 1908, 1939, 3079, NA, NA, NA, NA),
  d_c  = c(10, 25, 33, 32, 83, NA, NA, NA, NA), cy_c = c(1008, 1374, 1717, 1676, 4632, NA, NA, NA, NA),
  mort_rr_direct    = c(rep(NA, 5), 0.85, 0.846, 0.70, 0.83),
  mort_logse_direct = c(rep(NA, 5), 0.099, 0.044, 0.142, 0.094),
  include = c(rep(TRUE, 5), TRUE, TRUE, TRUE, FALSE),
  note = c(rep("D'Alessandro Table 2/5; comparator=untreated net; parasitaemia & mortality 1-4y", 5),
           "curtains; ~saturated 85-91% baseline (8-pt drop)",
           "spillover dilutes between-arm contrast; per-arm from ter Kuile Fig 1",
           "INFANT (1-11mo) prevalence -> large age-standardisation; Snow AJTMH 1996",
           "no parasite survey -> excluded from prevalence analysis"),
  stringsAsFactors = FALSE)

# mortality RR + SE(logRR): counts where available (D'Alessandro zones), else direct
rct$mort_rr    <- ifelse(!is.na(rct$d_t), (rct$d_t / rct$cy_t) / (rct$d_c / rct$cy_c), rct$mort_rr_direct)
rct$mort_logse <- ifelse(!is.na(rct$d_t), sqrt(1 / rct$d_t + 1 / rct$d_c), rct$mort_logse_direct)

## ---- age-standardise each arm's prevalence to PfPR2-10 (Smith 2007) ---------
std <- function(p, amin, amax) if (is.na(p)) NA_real_ else
  100 * as.numeric(malariaAtlas::convertPrevalence(p / 100, amin, amax, 2, 10))
rct$pfpr_t <- mapply(std, rct$prev_treated, rct$prev_amin, rct$prev_amax)
rct$pfpr_c <- mapply(std, rct$prev_control, rct$prev_amin, rct$prev_amax)
rct$dPfPR  <- rct$pfpr_c - rct$pfpr_t
rct$obs_red <- (1 - rct$mort_rr) * 100

## ---- per-trial predicted mortality reduction (%) from the R_dhs model --------
# Ridge-linear (primary):  RR = exp( beta_lin * (P_treated - P_control)/10 )
# Selected spline (ref):   RR = exp( s(P_treated) - s(P_control) ), population
#                          prediction with country random effects set to zero.
sm_delta <- function(pt, pc) {
  if (is.na(pt) || is.na(pc)) return(NA_real_)
  high <- newdata_at_mean(spline_model, pt / 10, year_c = 0)
  low  <- newdata_at_mean(spline_model, pc / 10, year_c = 0)
  link_difference(spline_model, high, low)$fit
}
rct$pred_linear <- (1 - exp(beta_lin * (rct$pfpr_t - rct$pfpr_c) / 10)) * 100
rct$pred_spline <- (1 - exp(mapply(sm_delta, rct$pfpr_t, rct$pfpr_c))) * 100

d <- rct[rct$include, ]
write.csv(rct, file.path(RESULTS_DIR, "triangulation_rct_data.csv"), row.names = FALSE)

## ---- meta-regression: observed log(mortality RR) ~ prevalence reduction -----
wls <- lm(log(mort_rr) ~ dPfPR, data = d, weights = 1 / mort_logse^2)
cat(sprintf("Trials/zones in fit: %d\n", nrow(d)))
cat(sprintf("WLS: mortality reduction per +1 PfPR2-10 pt reduced = %.3f%% (per +10 pts = %.1f%%); p=%.3f\n",
            (1 - exp(coef(wls)["dPfPR"])) * 100, (1 - exp(coef(wls)["dPfPR"] * 10)) * 100,
            summary(wls)$coefficients["dPfPR", "Pr(>|t|)"]))
cat(sprintf("R_dhs ridge-linear post-neonatal slope: %.1f%% mortality change per +10 PfPR pts\n",
            (exp(beta_lin) - 1) * 100))
print(d[, c("study", "pfpr_c", "pfpr_t", "dPfPR", "obs_red", "pred_spline", "pred_linear")],
      row.names = FALSE, digits = 3)

## ---- figure A: observed mortality reduction vs prevalence reduction ---------
lab <- sub("Phillips-Howard", "P-Howard", sub("D'Alessandro ", "DA ", sub(" 199.| 2003", "", d$study)))
d$nax_A <- c(0, 2.0, -2.0, 0, 4.5, -1.8, 1.8, -2.4); d$nay_A <- c(8, 8, -8, -8, 8, 8, -8, 7)
d$nax_B <- c(-4.5, 4.5, 4.5, 4.5, 4.5, -4.5, 4.5, -5); d$nay_B <- c(0, 5, -5, 0, 0, 0, 2, 0)
d$obs_lo <- (1 - exp(log(d$mort_rr) + 1.96 * d$mort_logse)) * 100
d$obs_hi <- (1 - exp(log(d$mort_rr) - 1.96 * d$mort_logse)) * 100
pA <- ggplot2::ggplot(d, ggplot2::aes(dPfPR, obs_red)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
  ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
  ggplot2::geom_smooth(method = "lm", mapping = ggplot2::aes(weight = 1 / mort_logse^2), se = TRUE,
                       colour = "grey30", fill = "grey85", linewidth = 0.7) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = (1 - exp(log(mort_rr) + 1.96 * mort_logse)) * 100,
                                      ymax = (1 - exp(log(mort_rr) - 1.96 * mort_logse)) * 100),
                         width = 0.7, alpha = 0.35) +
  ggplot2::geom_point(ggplot2::aes(colour = comparison), size = 3) +
  ggplot2::geom_text(ggplot2::aes(x = dPfPR + nax_A, y = obs_red + nay_A, label = lab), size = 3.3) +
  ggplot2::scale_colour_manual(values = c("no net" = "#08519c", "no curtain" = "#41ab5d", "untreated net" = "#d73027"), name = NULL) +
  ggplot2::labs(x = "Prevalence reduction (age-standardised PfPR2-10 points, control - intervention)",
                y = "Observed U5 mortality reduction (%)") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "top",
                 axis.title = ggplot2::element_text(size = 14), axis.text = ggplot2::element_text(size = 12),
                 legend.text = ggplot2::element_text(size = 12))

## ---- figure B: predicted (ridge-linear) vs observed --------------------------
pv <- data.frame(study = d$study, lab = lab, pred = d$pred_linear, obs = d$obs_red, comparison = d$comparison)
xr <- range(pv$pred, 0) + c(-4, 6)
yr <- range(d$obs_lo, d$obs_hi, 0) + c(-6, 10)
pB <- ggplot2::ggplot(pv, ggplot2::aes(pred, obs)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
  ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
  ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
  ggplot2::geom_errorbar(data = d, ggplot2::aes(x = pred_linear, ymin = obs_lo, ymax = obs_hi),
                         width = 0.9, alpha = 0.35, inherit.aes = FALSE) +
  ggplot2::geom_point(ggplot2::aes(colour = comparison), size = 3) +
  ggplot2::geom_text(data = d, ggplot2::aes(x = pred_linear + nax_B, y = obs_red + nay_B, label = lab),
                     size = 3.0, colour = "grey25", inherit.aes = FALSE) +
  ggplot2::scale_colour_manual(values = c("no net" = "#08519c", "no curtain" = "#41ab5d", "untreated net" = "#d73027"), name = NULL) +
  ggplot2::coord_cartesian(xlim = xr, ylim = yr) +
  ggplot2::labs(x = "Predicted U5 mortality reduction (%)", y = "Observed U5 mortality reduction (%)") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "top",
                 axis.title = ggplot2::element_text(size = 14), axis.text = ggplot2::element_text(size = 12),
                 legend.text = ggplot2::element_text(size = 12))

fig <- pA + pB + patchwork::plot_layout(widths = c(1.35, 1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16))
ggplot2::ggsave(file.path(RESULTS_DIR, "figure6_triangulation_rct.png"), fig, width = 14, height = 6.4, dpi = 300)
cat("saved:", file.path(RESULTS_DIR, "figure6_triangulation_rct.png"),
    "+ triangulation_rct_data.csv\n")

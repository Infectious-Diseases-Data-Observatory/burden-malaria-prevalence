# =============================================================================
# 10_triangulation_rct.R — RCT triangulation of the prevalence -> mortality link
#
# Premise (user): the Component 2 (DHS log-rate) model predicts the all-cause
# child MORTALITY effect from the PREVALENCE effect. So a bigger drop in parasite
# prevalence should map to a bigger mortality reduction. We test this against
# cluster-RCTs of insecticide-treated nets/curtains that measured BOTH.
# (Component 1's linear share model is unsuitable here: extrapolated per-arm to
# these prevalences its fitted share exceeds 100%, so it is not used.)
#
# Because Component 2's exposure-response is NON-LINEAR, we age-standardise each
# arm's parasite prevalence to PfPR2-10 and evaluate the model at BOTH the
# counterfactual (control) and intervention levels — not a single slope.
#
# Data hand-extracted from the ORIGINAL papers (+ Cochrane where noted); see the
# per-row provenance. Mortality is UNDER-5 (neonatal-excluded; the surveillance
# ages exclude neonates), matched to the m1mo5y Component-2 outcome.
#   D'Alessandro 1995 : 5 zones (Table 2 U5=1-4y mortality; Table 5 parasitaemia 1-4y)
#   Habluetzel 1997/99: overall (1-59mo mortality; 6-59mo parasitaemia 83% vs 91%)
#   Phillips-Howard'03: overall (1-59mo mortality; <3y parasitaemia 54% vs 66% [ter Kuile])
#   Nevill 1996       : overall (1-59mo mortality; 1-11mo INFANT parasitaemia 11.9% vs 25.1% [Snow AJTMH])
#   Binka 1996        : EXCLUDED (cluster-RCT but measured no parasite prevalence)
# CAVEAT: 8 heterogeneous points; Nevill infant-age & Phillips-Howard spillover
# make their prevalence-axis placement uncertain. Exploratory, not powered.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(malariaAtlas); library(lme4); library(mgcv); library(ggplot2); library(patchwork)})

## ---- extracted trial/zone data ----------------------------------------------
# prev_* = P. falciparum prevalence (%) by arm at follow-up (control = counterfactual).
# For D'Alessandro zones, mortality RR is computed from U5 (1-4y) death counts
# (d_/cy_ = deaths / child-years, treated & control); other trials give RR+SE directly.
rct <- data.frame(
  study   = c(paste0("D'Alessandro z", 1:5), "Habluetzel 1997", "Phillips-Howard 2003", "Nevill 1996", "Binka 1996"),
  country = c(rep("Gambia",5), "Burkina Faso", "Kenya", "Kenya", "Ghana"),
  comparison = c(rep("untreated net",5), "no curtain", "no net", "no net", "no net"),
  prev_treated = c(28.2, 22.8, 15.9, 43.2, 71.0,  83,  54,  11.9, NA),
  prev_control = c(36.7, 33.0, 25.8, 53.3, 44.8,  91,  66,  25.1, NA),
  prev_amin    = c(rep(1,5),   0.5, 0.25, 0.08, NA),          # prevalence survey age window (yrs)
  prev_amax    = c(rep(5,5),   5,   3,    0.92, NA),
  d_t  = c(4,13,18,35,66, NA, NA, NA, NA), cy_t = c(1056,1472,1908,1939,3079, NA,NA,NA,NA),
  d_c  = c(10,25,33,32,83, NA, NA, NA, NA), cy_c = c(1008,1374,1717,1676,4632, NA,NA,NA,NA),
  mort_rr_direct    = c(rep(NA,5), 0.85, 0.846, 0.70, 0.83),  # Cochrane forest RRs (overall trials)
  mort_logse_direct = c(rep(NA,5), 0.099, 0.044, 0.142, 0.094),
  include = c(rep(TRUE,5), TRUE, TRUE, TRUE, FALSE),
  note = c(rep("D'Alessandro Table 2/5; comparator=untreated net; parasitaemia & mortality 1-4y",5),
           "curtains; ~saturated 85-91% baseline (8-pt drop)",
           "spillover dilutes between-arm contrast; per-arm from ter Kuile Fig 1",
           "INFANT (1-11mo) prevalence -> large age-standardisation; Snow AJTMH 1996",
           "no parasite survey -> excluded from prevalence analysis"),
  stringsAsFactors = FALSE)

# mortality RR + SE(logRR): counts where available (D'Alessandro zones), else direct
rct$mort_rr    <- ifelse(!is.na(rct$d_t), (rct$d_t/rct$cy_t)/(rct$d_c/rct$cy_c), rct$mort_rr_direct)
rct$mort_logse <- ifelse(!is.na(rct$d_t), sqrt(1/rct$d_t + 1/rct$d_c),          rct$mort_logse_direct)

## ---- age-standardise each arm's prevalence to PfPR2-10 (Smith 2007) ---------
std <- function(p, amin, amax) if (is.na(p)) NA_real_ else
  100 * as.numeric(malariaAtlas::convertPrevalence(p/100, amin, amax, 2, 10))
rct$pfpr_t <- mapply(std, rct$prev_treated, rct$prev_amin, rct$prev_amax)
rct$pfpr_c <- mapply(std, rct$prev_control, rct$prev_amin, rct$prev_amax)
rct$dPfPR  <- rct$pfpr_c - rct$pfpr_t                        # prevalence reduction (PfPR2-10 pts)
rct$obs_red <- (1 - rct$mort_rr) * 100                       # observed mortality reduction (%)

## ---- Component 2 model fit (neonatal-excluded outcome, matched to trial ages) --
d2 <- read.csv(file.path(RESULTS, "component2_region_data.csv"))
dd <- d2[complete.cases(d2[, c("m1mo5y", C2_COVS, "country")]) & is.finite(d2$exposure) & d2$exposure > 0, ]
dd$deaths <- round(dd$m1mo5y/1000 * dd$exposure); dd$country <- factor(dd$country)
gam2 <- mgcv::gam(deaths ~ s(pfpr10, k = 4) + dtp3 + log_gdp + pct_urban + s(year_c) +
                    s(country, bs = "re") + s(country, pfpr10, bs = "re") + offset(log(exposure)),
                  family = mgcv::nb(), method = "REML", data = dd)
beta_lin <- fixef(fit_c2_lmm(d2, "m1mo5y")$m)["pfpr10"]      # linear LMM slope (per +10 PfPR pts)
# evaluate the fitted s(pfpr10) smooth (partial, log scale) at arbitrary PfPR2-10 %
# via lpmatrix, isolating the population prevalence-smooth basis columns.
.cf <- coef(gam2); .idx <- grep("^s\\(pfpr10\\)", names(.cf))
sm_eval <- function(pfpr_pct) {
  if (is.na(pfpr_pct)) return(NA_real_)
  nd <- data.frame(pfpr10 = pfpr_pct/10, dtp3 = mean(dd$dtp3), log_gdp = mean(dd$log_gdp),
                   pct_urban = mean(dd$pct_urban), year_c = 0, exposure = 1, country = dd$country[1])
  X <- predict(gam2, nd, type = "lpmatrix")
  sum(X[1, .idx] * .cf[.idx])
}

## ---- per-trial predicted mortality reduction (%) ----------------------------
# Component 2 spline (nonlinear, primary): RR = exp( s(P_treated) - s(P_control) )
# Component 2 linear (LMM, reference):      RR = exp( beta * (P_treated - P_control)/10 )
rct$pred_c2_spline <- (1 - exp(mapply(sm_eval, rct$pfpr_t) - mapply(sm_eval, rct$pfpr_c))) * 100
rct$pred_c2_linear <- (1 - exp(beta_lin * (rct$pfpr_t - rct$pfpr_c)/10)) * 100

d <- rct[rct$include, ]
write.csv(rct, file.path(RESULTS, "triangulation_rct_data.csv"), row.names = FALSE)

## ---- meta-regression: observed log(mortality RR) ~ prevalence reduction -----
wls <- lm(log(mort_rr) ~ dPfPR, data = d, weights = 1/mort_logse^2)
cat(sprintf("Trials/zones in fit: %d\n", nrow(d)))
cat(sprintf("WLS: mortality reduction per +1 PfPR2-10 pt reduced = %.3f%% (per +10 pts = %.1f%%); p=%.3f\n",
            (1-exp(coef(wls)["dPfPR"]))*100, (1-exp(coef(wls)["dPfPR"]*10))*100,
            summary(wls)$coefficients["dPfPR","Pr(>|t|)"]))
cat(sprintf("Component 2 LMM slope: %.1f%% mortality change per +10 PfPR pts\n", (exp(beta_lin)-1)*100))
print(d[, c("study","pfpr_c","pfpr_t","dPfPR","obs_red","pred_c2_spline","pred_c2_linear")],
      row.names = FALSE, digits = 3)

## ---- figure A: observed mortality reduction vs prevalence reduction ---------
lab <- sub("Phillips-Howard", "P-Howard", sub("D'Alessandro ", "DA ", sub(" 199.| 2003", "", d$study)))
# hand-tuned label offsets (d row order: DA z1-5, Habluetzel, P-Howard, Nevill)
d$nax_A <- c(0, 2.0, -2.0, 0,  4.5, -1.8, 1.8, -2.4); d$nay_A <- c(8, 8, -8, -8, 8, 8, -8, 7)
d$nax_B <- c(-4.5, 4.5, 4.5, 4.5, 4.5, -4.5, 4.5, -5); d$nay_B <- c(0, 5, -5, 0, 0, 0, 2, 0)  # beside error bars
# 95% CI on observed mortality reduction (from SE of log RR) — used in both panels
d$obs_lo <- (1 - exp(log(d$mort_rr) + 1.96*d$mort_logse)) * 100
d$obs_hi <- (1 - exp(log(d$mort_rr) - 1.96*d$mort_logse)) * 100
pA <- ggplot(d, aes(dPfPR, obs_red)) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_smooth(method = "lm", mapping = aes(weight = 1/mort_logse^2), se = TRUE,
              colour = "grey30", fill = "grey85", linewidth = 0.7) +
  geom_errorbar(aes(ymin = (1-exp(log(mort_rr)+1.96*mort_logse))*100,
                    ymax = (1-exp(log(mort_rr)-1.96*mort_logse))*100), width = 0.7, alpha = 0.35) +
  geom_point(aes(colour = comparison), size = 3) +
  geom_text(aes(x = dPfPR + nax_A, y = obs_red + nay_A, label = lab), size = 2.6) +
  scale_colour_manual(values = c("no net"="#08519c","no curtain"="#41ab5d","untreated net"="#d73027"), name = NULL) +
  labs(x = "Prevalence reduction (age-standardised PfPR2-10 points, control - intervention)",
       y = "Observed U5 mortality reduction (%)",
       title = "RCT dose-response: mortality reduction vs prevalence reduction",
       subtitle = sprintf("%d cluster-RCT points (D'Alessandro 5 zones + 3 trials). Grey = inverse-variance weighted fit.", nrow(d))) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "top")

## ---- figure B: Component 2 predicted vs observed ----------------------------
# single series: Component 2 nonlinear spline (linear LMM retained in the CSV only)
pv <- data.frame(study = d$study, lab = lab, pred = d$pred_c2_spline, obs = d$obs_red, comparison = d$comparison)
xr <- range(pv$pred, 0) + c(-4, 6)                          # x: predictions (tight)
yr <- range(d$obs_lo, d$obs_hi, 0) + c(-6, 10)              # y: observed + 95% CI (wide)
pB <- ggplot(pv, aes(pred, obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_errorbar(data = d, aes(x = pred_c2_spline, ymin = obs_lo, ymax = obs_hi),
                width = 0.9, alpha = 0.35, inherit.aes = FALSE) +
  geom_point(aes(colour = comparison), size = 3) +
  geom_text(data = d, aes(x = pred_c2_spline + nax_B, y = obs_red + nay_B, label = lab),
            size = 2.4, colour = "grey25", inherit.aes = FALSE) +
  scale_colour_manual(values = c("no net"="#08519c","no curtain"="#41ab5d","untreated net"="#d73027"), name = NULL) +
  coord_cartesian(xlim = xr, ylim = yr) +
  labs(x = "Component 2 predicted U5 mortality reduction (%)", y = "Observed U5 mortality reduction (%)",
       title = "Predicted vs observed (Component 2, nonlinear spline)",
       subtitle = "Dotted = identity (model matches trial). Error bars = 95% CI on observed (SE of log RR).") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(RESULTS, "triangulation_rct.png"), pA + pB + patchwork::plot_layout(widths = c(1.35, 1)),
       width = 14, height = 6.4, dpi = 140)
cat("saved: results/triangulation_rct.png + triangulation_rct_data.csv\n")

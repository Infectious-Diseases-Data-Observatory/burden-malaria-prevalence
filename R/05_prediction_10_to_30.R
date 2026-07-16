# =============================================================================
# 05_prediction_10_to_30.R — model prediction: how does TOTAL all-cause child
# mortality change as PfPR2-10 rises 10% -> 30%, holding covariates at their
# means, for each component, anchored to the SAME baseline at 10% (within each
# outcome). Done for both ALL under-5 and 1mo-5y (neonatal-excluded).
#
# Component 2 models total mortality directly (log-rate mixed model).
# Component 1 models malaria's SHARE; total recovered assuming NON-malaria
#   mortality is constant: total = non-malaria/(1 - share).
# NB: the 1mo-5y exposure-response is nonlinear (script 04 spline, edf~5); the
# log-linear prediction here is an approximation over the 10-30% span.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(lme4); library(ggplot2)})

c1 <- read.csv(file.path(RESULTS, "component1_country_data.csv"), stringsAsFactors = FALSE); c1$log_gdp <- log(c1$gdp_pc)
c2 <- read.csv(file.path(RESULTS, "component2_region_data.csv"), stringsAsFactors = FALSE)
covs <- C2_COVS     # match Component 2 primary: no-stunting (00_utils.R)
grid <- seq(10, 30, by = 0.5)

predict_outcome <- function(share_col, mort_col, label) {
  ## Component 1 — share model (shared spec: 00_utils.R)
  m1  <- fit_c1_share(c1, share_col)
  b1  <- coef(m1)["pfpr_pct"]; se1 <- summary(m1)$coefficients["pfpr_pct", "Std. Error"]
  s10 <- predict(m1, data.frame(pfpr_pct = 10, log_gdp = mean(c1$log_gdp), dtp3 = mean(c1$dtp3)))
  s30 <- predict(m1, data.frame(pfpr_pct = 30, log_gdp = mean(c1$log_gdp), dtp3 = mean(c1$dtp3)))
  ## Component 2 — total-mortality model (shared LMM spec: 00_utils.R)
  res2 <- fit_c2_lmm(c2, mort_col, covs); m2 <- res2$m; cc <- res2$dd
  b2 <- fixef(m2)["pfpr10"]; se2 <- sqrt(vcov(m2)["pfpr10", "pfpr10"])
  ndm <- data.frame(pfpr10 = 1, dtp3 = mean(cc$dtp3), log_gdp = mean(cc$log_gdp),
                    pct_urban = mean(cc$pct_urban), year_c = 0)
  B  <- exp(predict(m2, ndm, re.form = NA))                       # common baseline at 10% PfPR
  ## curves over 10 -> 30
  c2f <- B * exp(b2 * (grid - 10) / 10)
  c1t <- function(sl) { sp <- s10 + sl * (grid - 10); B * (1 - s10/100) / (1 - sp/100) }
  out <- rbind(
    data.frame(pfpr = grid, mort = c1t(b1), lo = c1t(b1 - 1.96*se1), hi = c1t(b1 + 1.96*se1),
               component = "Component 1 (share, non-malaria fixed)", outcome = label),
    data.frame(pfpr = grid, mort = c2f, lo = B*exp((b2-1.96*se2)*(grid-10)/10), hi = B*exp((b2+1.96*se2)*(grid-10)/10),
               component = "Component 2 (total-mortality model)", outcome = label))
  chg <- function(v) (v[grid == 30] / v[grid == 10] - 1) * 100
  cat(sprintf("[%s] baseline %.1f/1000 at 10%%; malaria share %.1f%%->%.1f%%. Total 10->30%%: C1 %+.1f%% (->%.1f), C2 %+.1f%% (->%.1f)\n",
              label, B, s10, s30, chg(out$mort[out$component == unique(out$component)[1]]), c1t(b1)[grid==30],
              chg(c2f), c2f[grid == 30]))
  out
}
pred <- rbind(predict_outcome("share_u5", "u5mr", "All under-5 (U5MR)"),
              predict_outcome("share_1mo5y", "m1mo5y", "1 month - 5 years (neonatal excl.)"))
write.csv(pred, file.path(RESULTS, "prediction_10_to_30.csv"), row.names = FALSE)

p <- ggplot(pred, aes(pfpr, mort, colour = component, fill = component)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(10, 30), linetype = "dotted", colour = "grey60") +
  geom_point(data = pred[pred$pfpr %in% c(10, 30), ], size = 2.2) +
  facet_wrap(~outcome, scales = "free_y") +
  scale_colour_manual(values = c("#08519c", "#d73027"), name = NULL) +
  scale_fill_manual(values = c("#08519c", "#d73027"), name = NULL) +
  labs(x = expression(PfPR[2-10]*" (%)"), y = "Total all-cause mortality (per 1,000 live births)",
       title = "Predicted total child mortality vs prevalence, 10% → 30% PfPR2-10",
       subtitle = "Covariates at means; within each panel both curves share the baseline at 10% PfPR. Shaded = 95% CI on the prevalence effect.",
       caption = "C2: direct total-mortality model. C1: malaria share model, total assuming non-malaria mortality fixed. 1mo-5y response is nonlinear (approx.).") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave(file.path(RESULTS, "prediction_10_to_30.png"), p, width = 12, height = 6, dpi = 300)
cat("saved: results/prediction_10_to_30.png + prediction_10_to_30.csv\n")

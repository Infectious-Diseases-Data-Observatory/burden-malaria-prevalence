# =============================================================================
# 05_prediction_10_to_30.R — model prediction: how does TOTAL all-cause under-5
# mortality change as PfPR2-10 rises from 10% to 30%, holding covariates at their
# means, for each component, anchored to the SAME baseline U5MR at 10%.
#
# Component 2 models total U5 mortality directly.
# Component 1 models malaria's SHARE of U5 deaths; total is recovered under the
#   assumption that NON-malaria mortality is constant: total = non-malaria/(1-share).
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(lme4); library(ggplot2)})

## ---- Component 1: share model -> total via constant non-malaria -------------
c1 <- read.csv(file.path(RESULTS, "component1_country_data.csv"), stringsAsFactors = FALSE)
c1$log_gdp <- log(c1$gdp_pc)
m1 <- lm(share_u5 ~ pfpr_pct + log_gdp + dtp3, data = c1)     # malaria % of all-U5 deaths
b1 <- coef(m1)["pfpr_pct"]; se1 <- summary(m1)$coefficients["pfpr_pct", "Std. Error"]
share <- function(p) predict(m1, data.frame(pfpr_pct = p, log_gdp = mean(c1$log_gdp), dtp3 = mean(c1$dtp3)))
s10 <- share(10)

## ---- Component 2: total U5 mortality model ---------------------------------
c2 <- read.csv(file.path(RESULTS, "component2_region_data.csv"), stringsAsFactors = FALSE)
covs <- c("pfpr10", "dtp3", "log_gdp", "pct_urban", "year_c", "stunting")
cc <- c2[complete.cases(c2[, c("u5mr", covs, "country")]), ]
m2 <- lmer(log(u5mr) ~ pfpr10 + dtp3 + log_gdp + pct_urban + year_c + stunting + (1 + pfpr10 || country),
           data = cc, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
b2 <- fixef(m2)["pfpr10"]; se2 <- sqrt(vcov(m2)["pfpr10", "pfpr10"])   # per +10 PfPR pts (pfpr10 units)

# common baseline B = Component 2's predicted U5MR at PfPR 10%, covariates at means
ndm <- data.frame(pfpr10 = 1, dtp3 = mean(cc$dtp3), log_gdp = mean(cc$log_gdp),
                  pct_urban = mean(cc$pct_urban), year_c = 0, stunting = mean(cc$stunting))
B <- exp(predict(m2, ndm, re.form = NA))

## ---- prediction curves over PfPR 10 -> 30 -----------------------------------
grid <- seq(10, 30, by = 0.5)
# Component 2: multiplicative in log-rate; anchored at 10%
c2_fit <- B * exp(b2 * (grid - 10) / 10)
c2_lo  <- B * exp((b2 - 1.96 * se2) * (grid - 10) / 10)
c2_hi  <- B * exp((b2 + 1.96 * se2) * (grid - 10) / 10)
# Component 1: total = non-malaria/(1-share); non-malaria fixed at baseline
tot_c1 <- function(slope) { s_p <- s10 + slope * (grid - 10)
  B * (1 - s10/100) / (1 - s_p/100) }
c1_fit <- tot_c1(b1); c1_lo <- tot_c1(b1 - 1.96 * se1); c1_hi <- tot_c1(b1 + 1.96 * se1)

pred <- rbind(
  data.frame(pfpr = grid, u5mr = c1_fit, lo = c1_lo, hi = c1_hi, component = "Component 1 (share → total, non-malaria fixed)"),
  data.frame(pfpr = grid, u5mr = c2_fit, lo = c2_lo, hi = c2_hi, component = "Component 2 (total mortality model)"))
write.csv(pred, file.path(RESULTS, "prediction_10_to_30.csv"), row.names = FALSE)

chg <- function(v) (v[grid == 30] / v[grid == 10] - 1) * 100
cat(sprintf("Baseline U5MR at 10%% PfPR (common anchor): %.1f per 1,000\n", B))
cat(sprintf("Malaria %%-share of U5 deaths: %.1f%% at 10%% -> %.1f%% at 30%% PfPR\n", s10, share(30)))
cat(sprintf("TOTAL U5MR change, 10%% -> 30%% PfPR2-10:\n"))
cat(sprintf("  Component 1 (non-malaria fixed): %.1f -> %.1f per 1,000  (%+.1f%%)\n", B, c1_fit[grid==30], chg(c1_fit)))
cat(sprintf("  Component 2 (direct model)     : %.1f -> %.1f per 1,000  (%+.1f%%)\n", B, c2_fit[grid==30], chg(c2_fit)))

## ---- figure -----------------------------------------------------------------
p <- ggplot(pred, aes(pfpr, u5mr, colour = component, fill = component)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(10, 30), linetype = "dotted", colour = "grey60") +
  geom_point(data = pred[pred$pfpr %in% c(10, 30), ], size = 2.4) +
  scale_colour_manual(values = c("#08519c", "#d73027"), name = NULL) +
  scale_fill_manual(values = c("#08519c", "#d73027"), name = NULL) +
  labs(x = expression(PfPR[2-10]*" (%)"), y = "Total all-cause under-5 mortality (per 1,000 live births)",
       title = "Predicted total U5 mortality vs prevalence, 10% → 30% PfPR2-10",
       subtitle = sprintf("Both anchored at %.0f per 1,000 at 10%% PfPR, covariates held at means. Shaded = 95%% CI on the prevalence effect.", B),
       caption = "Component 2: direct total-mortality model. Component 1: malaria share model, total recovered assuming non-malaria mortality is fixed.") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave(file.path(RESULTS, "prediction_10_to_30.png"), p, width = 9, height = 6, dpi = 140)
cat("saved: results/prediction_10_to_30.png + prediction_10_to_30.csv\n")

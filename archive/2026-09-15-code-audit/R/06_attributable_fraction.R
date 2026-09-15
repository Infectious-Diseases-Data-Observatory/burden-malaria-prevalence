# =============================================================================
# 06_attributable_fraction.R — malaria-attributable fraction of child mortality
# implied by the PRIMARY Component 2 (Method 2) model.
#
# Primary model = linear-in-prevalence negative-binomial GAM on the analysis
# sample (PfPR2-10 >= 1%); see 00_utils.R::fit_c2_primary and SText 1.
# Attributable fraction is referenced to a 1% PfPR counterfactual:
#   AF(p) = 1 - exp(-beta * (p - 1)/10),   p >= 1%,
# the proportion of child deaths at prevalence p attributable to malaria relative
# to a 1%-prevalence setting (covariates and country effects cancel in the ratio).
# Outcomes: U5MR (5q0) and 1mo-5y (neonatal-excluded, primary outcome).
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})

d      <- read.csv(file.path(RESULTS, "component2_region_data.csv"), stringsAsFactors = FALSE)
grid_p <- seq(C2_MIN_PFPR, 80, by = 0.5)                          # PfPR2-10 (%), from the 1% floor
LAB    <- c(u5mr = "All under-5 (U5MR)", m1mo5y = "1mo-5y (neonatal excluded)")

## ---- attributable fraction from the linear model coefficient ----------------
af_one <- function(outcome) {
  f <- fit_c2_primary(d, outcome); b <- f$beta; se <- f$se
  cat(sprintf("  %-7s: +%.1f%% per +10 PfPR2-10 pts (beta=%.4f, SE=%.4f; n=%d, %d countries)\n",
              outcome, (exp(b) - 1) * 100, b, se, nrow(f$dd), nlevels(f$dd$country)))
  data.frame(pfpr = grid_p, af = af_c2(b, grid_p),
             lo = af_c2(b - 1.96 * se, grid_p), hi = af_c2(b + 1.96 * se, grid_p),
             outcome = LAB[[outcome]])
}
cat("Method 2 primary model (linear nb-GAM, PfPR2-10 >= 1%):\n")
af <- do.call(rbind, lapply(c("u5mr", "m1mo5y"), af_one))
af$outcome <- factor(af$outcome, levels = LAB)
write.csv(af, file.path(RESULTS, "attributable_fraction.csv"), row.names = FALSE)

## ---- figure -----------------------------------------------------------------
p <- ggplot(af, aes(pfpr, 100 * af)) +
  geom_ribbon(aes(ymin = 100 * lo, ymax = 100 * hi), fill = "#d73027", alpha = 0.15) +
  geom_line(colour = "#d73027", linewidth = 1) +
  facet_wrap(~outcome) +
  scale_x_continuous(breaks = seq(0, 80, 20)) +
  labs(x = expression("Age-standardised "*italic(Pf)*"PR"[2-10]*" (%)"),
       y = "Malaria-attributable fraction of child deaths (%)",
       title = "Model-implied malaria-attributable fraction of child mortality (Method 2, primary model)",
       subtitle = "AF(p) = 1 - exp(-β(p-1)/10): proportion of child deaths attributable to malaria vs a 1% PfPR counterfactual.",
       caption = "Linear nb-GAM on PfPR2-10 >= 1% (467 survey-regions, 22 countries). Shaded = 95% CI (prevalence-coefficient sampling error only).") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(file.path(RESULTS, "attributable_fraction.png"), p, width = 11, height = 5.5, dpi = 300)

## ---- console summary at round prevalences -----------------------------------
show <- af[af$pfpr %in% c(10, 20, 30, 40, 50), c("outcome", "pfpr", "af", "lo", "hi")]
show[, c("af", "lo", "hi")] <- round(100 * show[, c("af", "lo", "hi")], 1)
cat("\nMalaria-attributable fraction of child deaths (%), referenced to 1% PfPR:\n"); print(show, row.names = FALSE)
cat("\nsaved: results/attributable_fraction.png + attributable_fraction.csv\n")

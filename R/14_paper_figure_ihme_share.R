# =============================================================================
# 14_paper_figure_ihme_share.R — PAPER FIGURE (Method 1 / IHME deconvolution)
# Estimated proportion of all-cause POST-NEONATAL child deaths that IHME/GBD
# attributes to malaria, as a function of national age-standardised PfPR2-10,
# fitted with a penalised SPLINE and adjusted for log GDP per capita and DTP3
# coverage (same covariates as the linear Method-1 model). Country points +
# spline fit + 95% CI; shows the plateau at high prevalence. No header text.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})

c1 <- read.csv(file.path(RESULTS, "component1_country_data.csv"), stringsAsFactors = FALSE)
c1$log_gdp <- log(c1$gdp_pc)
c1 <- c1[is.finite(c1$share_1mo5y) & is.finite(c1$pfpr_pct) &
           is.finite(c1$log_gdp) & is.finite(c1$dtp3), ]
c1 = c1[c1$pfpr_pct>0.5, ]

## spline model: malaria's share of post-neonatal deaths ~ s(PfPR) + covariates
m <- mgcv::gam(share_1mo5y ~ s(pfpr_pct) + log_gdp + dtp3, data = c1, method = "REML")
edf <- summary(m)$s.table["s(pfpr_pct)", "edf"]

g <- data.frame(pfpr_pct = exp(seq(log(min(c1$pfpr_pct)), log(max(c1$pfpr_pct)), length.out = 200)),
                log_gdp = mean(c1$log_gdp), dtp3 = mean(c1$dtp3))
pr <- predict(m, g, se.fit = TRUE)
g$fit <- pr$fit; g$lo <- pr$fit - 1.96 * pr$se.fit; g$hi <- pr$fit + 1.96 * pr$se.fit

BLU <- "#2c7fb8"
p <- ggplot() +
  geom_ribbon(data = g, aes(pfpr_pct, ymin = lo, ymax = hi), fill = BLU, alpha = 0.18) +
  geom_line(data = g, aes(pfpr_pct, fit), colour = BLU, linewidth = 1.2) +
  geom_point(data = c1, aes(pfpr_pct, share_1mo5y), colour = "grey35", alpha = 0.7, size = 2.2) +
  geom_hline(yintercept = 0, colour = "grey85") +
  scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10, 30)) +
  labs(x = expression("National age-standardised "*italic(Pf)*"PR"[2-10]*" (%), log scale"),
       y = "Malaria's share of post-neonatal\nchild deaths, IHME/GBD (%)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14), axis.text = element_text(size = 13))

OUT <- file.path(RESULTS, "fig_ihme_share_spline.png")
ggsave(OUT, p, width = 8, height = 5.6, dpi = 320)
cat(sprintf("saved: %s\n", OUT))
at <- function(x) as.numeric(predict(m, data.frame(pfpr_pct = x, log_gdp = mean(c1$log_gdp), dtp3 = mean(c1$dtp3))))
cat(sprintf("n=%d countries; s(PfPR) edf=%.2f; fitted share at PfPR 10/20/30/36%% = %.0f/%.0f/%.0f/%.0f%%\n",
            nrow(c1), edf, at(10), at(20), at(30), at(max(c1$pfpr_pct))))
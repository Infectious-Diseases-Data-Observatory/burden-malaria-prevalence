# =============================================================================
# 12_paper_figure1.R — PAPER FIGURE 1 (two-panel; PRIMARY model M1)
# Primary Method-2 model: linear-in-prevalence negative-binomial GAM on the
# analysis sample (PfPR2-10 >= 1%); see 00_utils.R::fit_c2_primary and SText 1.
#
# A: all analysis-sample survey-regions (points sized by birth-exposure) on
#    log-log axes, with the population-level mean fit (fixed effects only;
#    country random effects excluded).
# B: malaria-attributable FRACTION of post-neonatal mortality implied by the
#    same fit, referenced to a 1% PfPR counterfactual:  1 - exp(-beta(p-1)/10).
# No panel titles/subtitles; panels tagged A / B; grid lines on.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2); library(patchwork)})

d <- read.csv(file.path(RESULTS, "component2_region_data.csv"), stringsAsFactors = FALSE)
f <- fit_c2_primary(d, "m1mo5y")            # M1: linear nb-GAM, PfPR2-10 >= 1%
m <- f$m; dd <- f$dd; beta <- f$beta; se <- f$se
re_terms <- grep("country", vapply(m$smooth, function(s) s$label, ""), value = TRUE)

## ---- population-level mean fit (country RE excluded) ------------------------
grid <- data.frame(pfpr2_10 = exp(seq(log(C2_MIN_PFPR), log(max(dd$pfpr2_10)), length.out = 200)),
                   dtp3 = mean(dd$dtp3), log_gdp = mean(dd$log_gdp),
                   pct_urban = mean(dd$pct_urban), year_c = 0, exposure = 1,
                   country = dd$country[1])
grid$pfpr10 <- grid$pfpr2_10 / 10
pr <- predict(m, grid, type = "link", se.fit = TRUE, exclude = re_terms)
grid$rate <- 1000 * exp(pr$fit)
grid$lo   <- 1000 * exp(pr$fit - 1.96 * pr$se.fit)
grid$hi   <- 1000 * exp(pr$fit + 1.96 * pr$se.fit)
# attributable fraction referenced to 1% (closed form; coefficient CI)
grid$af    <- 100 * af_c2(beta, grid$pfpr2_10)
grid$af_lo <- 100 * af_c2(beta - 1.96 * se, grid$pfpr2_10)
grid$af_hi <- 100 * af_c2(beta + 1.96 * se, grid$pfpr2_10)

GRID_MINOR <- element_line(colour = "grey92", linewidth = 0.3)
GRID_MAJOR <- element_line(colour = "grey85", linewidth = 0.4)
XSC  <- scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 80), expand = expansion(mult = c(0.02, 0.03)))
XLAB <- expression("Age-standardised "*italic(Pf)*"PR"[2-10]*" (%), log scale")
BLU  <- "#08519c"

## ---- panel A: scatter + mean fit -------------------------------------------
pA <- ggplot() +
  geom_point(data = dd, aes(pfpr2_10, m1mo5y, size = exposure),
             colour = "grey45", alpha = 0.40) +
  geom_ribbon(data = grid, aes(pfpr2_10, ymin = lo, ymax = hi), fill = BLU, alpha = 0.18) +
  geom_line(data = grid, aes(pfpr2_10, rate), colour = BLU, linewidth = 1.1) +
  scale_size_area(max_size = 6, name = "Births in\nsurvey-region",
                  breaks = c(500, 2000, 5000), labels = scales::comma) +
  XSC + scale_y_log10(breaks = c(10, 20, 50, 100, 200)) + coord_cartesian(ylim = c(10, 200)) +
  labs(x = XLAB,
       y = "All-cause post-neonatal mortality\n(1 month-5 years, per 1,000 live births, log scale)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        panel.grid.minor = GRID_MINOR, panel.grid.major = GRID_MAJOR,
        legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
        legend.key.size = unit(4, "mm"), legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

## ---- panel B: attributable share vs 1% counterfactual ----------------------
pB <- ggplot(grid, aes(pfpr2_10, af)) +
  geom_ribbon(aes(ymin = af_lo, ymax = af_hi), fill = "#d73027", alpha = 0.16) +
  geom_line(colour = "#d73027", linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
  XSC +
  labs(x = XLAB,
       y = "Share of post-neonatal deaths attributable\nto malaria (%, vs 1% PfPR counterfactual)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        panel.grid.minor = GRID_MINOR, panel.grid.major = GRID_MAJOR)

fig <- pA + pB + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

OUT <- file.path(RESULTS, "fig1_survey_gam_postneonatal.png")
ggsave(OUT, fig, width = 13, height = 5.8, dpi = 320)
cat(sprintf("saved: %s\n", OUT))
cat(sprintf("M1 primary: n=%d survey-regions, %d countries; +%.1f%% per +10 PfPR2-10 pts (95%% CI %+.1f to %+.1f)\n",
            nrow(dd), nlevels(dd$country), (exp(beta)-1)*100, (exp(beta-1.96*se)-1)*100, (exp(beta+1.96*se)-1)*100))
for (p in c(10, 40, 80)) cat(sprintf("  attributable share at PfPR2-10=%d%%: %.0f%% (95%% CI %.0f-%.0f)\n",
    p, 100*af_c2(beta,p), 100*af_c2(beta-1.96*se,p), 100*af_c2(beta+1.96*se,p)))

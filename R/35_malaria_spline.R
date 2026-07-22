# =============================================================================
# 35_malaria_spline.R — shape of the fitted GAM malaria spline s(pfpr10), for
# POST-NEONATAL (m_comb, R/25) vs NEONATAL (n_comb, R/26), on the log-mortality
# partial-effect scale (centred, as mgcv fits it) with 95% CI. Both models share
# the same 847-region-year sample and covariate set, so the centred splines are
# directly comparable. Rug = observed region-year prevalence.
# =============================================================================
source("R/00_utils.R"); suppressMessages({library(mgcv); library(ggplot2)})
mpn <- readRDS(file.path(RESULTS, "combined_models.rds"))$m_comb   # post-neonatal
mnn <- readRDS(file.path(RESULTS, "neonatal_models.rds"))$n_comb   # neonatal
grid <- data.frame(pfpr2_10 = seq(1, max(mpn$model$pfpr10)*10, length.out = 250)); grid$pfpr10 <- grid$pfpr2_10/10

sp <- function(mod, lab) { mf <- mod$model
  nd <- data.frame(pfpr10 = grid$pfpr10, year_c = 0, dtp3_reg = mean(mf$dtp3_reg), facility = mean(mf$facility),
    educ_yrs = mean(mf$educ_yrs), wealth_q = mean(mf$wealth_q), log_gdp = 0, pct_urban = mean(mf$pct_urban),
    exposure = 1, country = levels(mf$country)[1])
  pr <- predict(mod, nd, type = "terms", terms = "s(pfpr10)", se.fit = TRUE)
  data.frame(pfpr2_10 = grid$pfpr2_10, outcome = lab, fit = as.numeric(pr$fit),
             lo = as.numeric(pr$fit) - 1.96*as.numeric(pr$se.fit),
             hi = as.numeric(pr$fit) + 1.96*as.numeric(pr$se.fit)) }
L <- rbind(sp(mpn, "Post-neonatal (1mo-5y)"), sp(mnn, "Neonatal (<1mo)"))
cols <- c("Post-neonatal (1mo-5y)" = "#08519c", "Neonatal (<1mo)" = "#d73027")

p <- ggplot(L, aes(pfpr2_10, fit, colour = outcome, fill = outcome)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_rug(data = data.frame(x = mpn$model$pfpr10*10), aes(x = x), inherit.aes = FALSE, alpha = 0.22, sides = "b") +
  scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
  labs(x = expression(italic(Pf)*"PR"[2-10]~"(%)"),
       y = "s(PfPR) partial effect on mortality (log scale)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave(file.path(RESULTS, "malaria_spline_shape.png"), p, width = 8, height = 5.8, dpi = 300)
cat(sprintf("post-neonatal s(pfpr10): edf=%.2f p=%.1e\n", summary(mpn)$s.table["s(pfpr10)","edf"], summary(mpn)$s.table["s(pfpr10)","p-value"]))
cat(sprintf("neonatal      s(pfpr10): edf=%.2f p=%.3f\n", summary(mnn)$s.table["s(pfpr10)","edf"], summary(mnn)$s.table["s(pfpr10)","p-value"]))
cat("saved: results/malaria_spline_shape.png\n")

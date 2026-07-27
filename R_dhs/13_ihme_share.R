# =============================================================================
# 13_ihme_share.R — manuscript figure: the share of all-cause POST-NEONATAL
# child deaths that IHME/GBD attributes to malaria, as a function of national
# age-standardised PfPR2-10, fitted with a penalised spline and adjusted for
# log GDP per capita and DTP3 coverage. (Ported from legacy R/14.)
#
# This is a Method-1 (IHME/GBD deconvolution) comparison figure, distinct from
# the DHS/MIS survey model. It reads the country-level panel produced by the
# Method-1 script R/02_component1_country.R (results/component1_country_data.csv,
# a tracked artifact); it is skipped with a message if that panel is absent.
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2"))

C1_PANEL <- file.path(REPO_ROOT, "results", "component1_country_data.csv")
if (!file.exists(C1_PANEL)) {
  message(
    "Method-1 country panel not found (", C1_PANEL, "); skipping the IHME ",
    "attributable-share figure. Generate it with R/02_component1_country.R."
  )
  quit(save = "no", status = 0)
}

c1 <- read.csv(C1_PANEL, stringsAsFactors = FALSE)
c1$log_gdp <- log(c1$gdp_pc)
c1 <- c1[is.finite(c1$share_1mo5y) & is.finite(c1$pfpr_pct) &
           is.finite(c1$log_gdp) & is.finite(c1$dtp3) & c1$pfpr_pct > 0.5, ]

## spline model: malaria's share of post-neonatal deaths ~ s(PfPR) + covariates
m <- mgcv::gam(share_1mo5y ~ s(pfpr_pct) + log_gdp + dtp3, data = c1, method = "REML")
edf <- summary(m)$s.table["s(pfpr_pct)", "edf"]

g <- data.frame(
  pfpr_pct = exp(seq(log(min(c1$pfpr_pct)), log(max(c1$pfpr_pct)), length.out = 200)),
  log_gdp = mean(c1$log_gdp), dtp3 = mean(c1$dtp3)
)
pr <- predict(m, g, se.fit = TRUE)
g$fit <- pr$fit; g$lo <- pr$fit - 1.96 * pr$se.fit; g$hi <- pr$fit + 1.96 * pr$se.fit

BLU <- "#2c7fb8"
p <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(data = g, ggplot2::aes(pfpr_pct, ymin = lo, ymax = hi), fill = BLU, alpha = 0.18) +
  ggplot2::geom_line(data = g, ggplot2::aes(pfpr_pct, fit), colour = BLU, linewidth = 1.2) +
  ggplot2::geom_point(data = c1, ggplot2::aes(pfpr_pct, share_1mo5y), colour = "grey35", alpha = 0.7, size = 2.2) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
  ggplot2::scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10, 30)) +
  ggplot2::labs(x = expression("National age-standardised " * italic(Pf) * "PR"[2-10] * " (%), log scale"),
                y = "Malaria's share of post-neonatal\nchild deaths, IHME/GBD (%)") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 axis.title = ggplot2::element_text(size = 14), axis.text = ggplot2::element_text(size = 13))

out <- file.path(RESULTS_DIR, "figure7_ihme_share.png")
ggplot2::ggsave(out, p, width = 8, height = 5.6, dpi = 320)
write.csv(
  c1[, c("iso3", "pfpr_pct", "share_1mo5y", "gdp_pc", "dtp3")],
  file.path(RESULTS_DIR, "ihme_share_country_data.csv"), row.names = FALSE
)
at <- function(x) as.numeric(predict(m, data.frame(pfpr_pct = x, log_gdp = mean(c1$log_gdp), dtp3 = mean(c1$dtp3))))
cat(sprintf("saved: %s\n", out))
cat(sprintf("n=%d countries; s(PfPR) edf=%.2f; fitted share at PfPR 10/20/30/max = %.0f/%.0f/%.0f/%.0f%%\n",
            nrow(c1), edf, at(10), at(20), at(30), at(max(c1$pfpr_pct))))

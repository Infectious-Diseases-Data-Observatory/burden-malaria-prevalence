# =============================================================================
# 28_country_af_curves.R — malaria-attributable SHARE of post-neonatal deaths as
# a function of PfPR2-10, one line per country, from the best-fitting combined
# model (m_comb, R/25). Counterfactual = 1% PfPR.
#
# Each country's response = population nonlinear smooth s(pfpr10) + its own linear
# random-slope tilt b_c from s(country,pfpr10,bs="re"). Referenced to 1% the AF is:
#   AF_c(p) = 1 - exp(-[ (f(p) - f(1)) + b_c*(p-1)/10 ]),
# with f() the population smooth (intercept, covariates, country intercept, year
# all cancel). Each curve is drawn only over that country's observed PfPR range
# (no extrapolation of its linear tilt beyond data). Share floored at 0.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})
m <- readRDS(file.path(RESULTS, "combined_models.rds"))$m_comb

## ---- population smooth difference f(p) - f(1) -------------------------------
# isolate the s(pfpr10) contribution via type="terms"; the smooth's centering
# constant cancels in the p-vs-1 difference, so covariates/country are irrelevant.
mf <- m$model
nd <- function(p10) data.frame(pfpr10 = p10, year_c = 0, dtp3_reg = mean(mf$dtp3_reg),
  facility = mean(mf$facility), educ_yrs = mean(mf$educ_yrs), wealth_q = mean(mf$wealth_q),
  log_gdp = mean(mf$log_gdp), pct_urban = mean(mf$pct_urban), exposure = 1,
  country = levels(mf$country)[1])
sterm <- function(p10) as.numeric(predict(m, nd(p10), type = "terms", terms = "s(pfpr10)"))
e1 <- sterm(0.1)
fpop <- function(p) sterm(p/10) - e1

## ---- country random-slope deviations, keyed by iso3 -------------------------
sm  <- m$smooth[[which(vapply(m$smooth, function(s) s$label, "") == "s(country,pfpr10)")]]
b_c <- setNames(as.numeric(coef(m)[sm$first.para:sm$last.para]), levels(mf$country))
pmax_obs <- tapply(mf$pfpr10 * 10, mf$country, max)                 # each country's observed max PfPR (%)

## ---- assemble one curve per country over its observed range -----------------
cur <- do.call(rbind, lapply(names(b_c), function(iso) {
  hi <- max(2, pmax_obs[[iso]]); ps <- seq(1, hi, length.out = 60)
  data.frame(iso3 = iso, pfpr = ps,
             af = pmax(1 - exp(-(fpop(ps) + b_c[[iso]] * (ps - 1) / 10)), 0),
             slope_dev = b_c[[iso]])
}))
ends <- do.call(rbind, by(cur, cur$iso3, function(x) x[which.max(x$pfpr), ]))
popc <- data.frame(pfpr = seq(1, max(pmax_obs), length.out = 80))
popc$af <- pmax(1 - exp(-fpop(popc$pfpr)), 0)                       # population mean (b_c = 0)

## ---- plot -------------------------------------------------------------------
p <- ggplot(cur, aes(pfpr, 100 * af, group = iso3, colour = slope_dev)) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  geom_line(data = popc, aes(pfpr, 100 * af), inherit.aes = FALSE,
            colour = "black", linewidth = 1.4, linetype = "22") +
  geom_text(data = ends, aes(pfpr, 100 * af, label = iso3), size = 2.6, hjust = -0.15,
            show.legend = FALSE, check_overlap = TRUE) +
  scale_colour_viridis_c(option = "C", end = 0.92,
                         name = "Country slope\ndeviation (log, /+10pts)") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.08))) +
  labs(x = expression(PfPR[2-10]~"(%)"),
       y = "Malaria-attributable share of post-neonatal deaths (%)  [vs 1% counterfactual]") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")
ggsave(file.path(RESULTS, "country_af_curves_postneonatal.png"), p, width = 9, height = 6.5, dpi = 300)

## ---- log-log variant --------------------------------------------------------
pll <- ggplot(subset(cur, af >= 0.005), aes(pfpr, 100 * af, group = iso3, colour = slope_dev)) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  geom_line(data = subset(popc, af >= 0.005), aes(pfpr, 100 * af), inherit.aes = FALSE,
            colour = "black", linewidth = 1.4, linetype = "22") +
  geom_text(data = ends, aes(pfpr, 100 * af, label = iso3), size = 2.6, hjust = -0.15,
            show.legend = FALSE, check_overlap = TRUE) +
  scale_colour_viridis_c(option = "C", end = 0.92, name = "Country slope\ndeviation (log, /+10pts)") +
  scale_x_log10(breaks = c(1,2,5,10,20,40,60,80)) + scale_y_log10(breaks = c(1,2,5,10,20,40,60)) +
  annotation_logticks(sides = "bl") +
  labs(x = expression(PfPR[2-10]~"(%)"),
       y = "Malaria-attributable share of post-neonatal deaths (%)  [vs 1% counterfactual]") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), legend.position = "right")
ggsave(file.path(RESULTS, "country_af_curves_postneonatal_loglog.png"), pll, width = 9, height = 6.5, dpi = 300)

## ---- console: spread at a common prevalence ---------------------------------
af_at <- function(pp) sort(setNames(vapply(names(b_c),
  function(iso) 100 * max(1 - exp(-(fpop(pp) + b_c[[iso]] * (pp - 1)/10)), 0), 0), names(b_c)))
cat(sprintf("Population AF at PfPR=30%%: %.1f%%\n", 100 * (1 - exp(-fpop(30)))))
a30 <- af_at(30)
cat(sprintf("Country AF at PfPR=30%% (of %d countries): range %.1f%% (%s) to %.1f%% (%s), median %.1f%%\n",
            length(a30), a30[1], names(a30)[1], a30[length(a30)], names(a30)[length(a30)], median(a30)))
cat("dashed black line = population mean; black-dashed reference in figure.\n")
cat("saved: results/country_af_curves_postneonatal.png\n")

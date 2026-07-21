# =============================================================================
# 20_te_pfpr_year.R — full tensor-product te(PfPR2-10, calendar year) on the
# EXPANDED 2000-2024 panel (R/18). Lets the prevalence->mortality relationship
# bend freely in BOTH prevalence and year, instead of the linear-in-PfPR /
# linear-interaction forms. Compares fit (AIC) against:
#   m_lin : pfpr10 + country random slope           (no year interaction)
#   m_int : pfpr10 + pfpr10:year_c                   (linear interaction)
#   m_add : te-marginals only, s(pfpr10)+s(year_c)   (additive, NO interaction)
#   m_te  : te(pfpr10, year_c)                        (full 2D surface)
# Country random intercept + random slope retained throughout; nb GAM, REML.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})

d <- read.csv(file.path(RESULTS, "component2_region_data_expanded.csv"), stringsAsFactors = FALSE)
mods <- readRDS(file.path(RESULTS, "expanded_models.rds")); center <- mods$center

## ---- reconstruct the fit sample exactly as R/18 (PfPR2-10 >= 1%) -------------
need <- c("m1mo5y","pfpr10","dtp3","log_gdp","pct_urban","year_c","iso3","svkey")
fitd <- d[complete.cases(d[, need]) & is.finite(d$exposure) & d$exposure > 0 & d$pfpr2_10 >= 1, ]
fitd$deaths  <- round(fitd$m1mo5y / 1000 * fitd$exposure)
fitd$country <- factor(fitd$iso3); fitd$svkey <- factor(fitd$svkey)
cat(sprintf("fit sample: %d region-years, %d countries, years %d-%d\n",
            nrow(fitd), nlevels(fitd$country), min(fitd$year), max(fitd$year)))

RE  <- "s(country,bs=\"re\") + s(country,pfpr10,bs=\"re\") + offset(log(exposure))"
COV <- "dtp3 + log_gdp + pct_urban"
f   <- function(core) as.formula(paste("deaths ~", core, "+", COV, "+", RE))

m_lin <- mods$m_lin; m_int <- mods$m_int                       # reuse the R/18 fits
m_add <- gam(f("s(pfpr10) + s(year_c)"),          family = nb(), method = "REML", data = fitd)
m_te  <- gam(f("te(pfpr10, year_c, k = c(6, 5))"), family = nb(), method = "REML", data = fitd)

## ---- fit comparison ---------------------------------------------------------
edf_te <- sum(m_te$edf[grep("te\\(pfpr10,year_c\\)", names(m_te$edf))])
tab <- data.frame(
  model = c("m_lin (no interaction)", "m_int (linear pfpr x year)",
            "m_add (s(pfpr)+s(year), additive)", "m_te (full tensor)"),
  edf   = round(c(sum(m_lin$edf), sum(m_int$edf), sum(m_add$edf), sum(m_te$edf)), 1),
  AIC   = round(c(AIC(m_lin), AIC(m_int), AIC(m_add), AIC(m_te)), 1))
tab$dAIC <- round(tab$AIC - min(tab$AIC), 1)
cat("\n--- fit comparison (lower AIC = better) ---\n"); print(tab, row.names = FALSE)
cat(sprintf("\nte(pfpr10,year_c): edf = %.1f\n", edf_te))
st <- summary(m_te)$s.table; print(round(st[grep("te\\(", rownames(st)), , drop = FALSE], 4))

## ---- AF vs 1% implied by the te surface, by year ----------------------------
# AF_year(p) = 1 - exp(-(eta_te(p,year) - eta_te(1,year))); covariates/country cancel.
af_te <- function(p, yr) {
  nd0 <- data.frame(pfpr10 = 0.1, year_c = yr - center, dtp3 = 0, log_gdp = 0,
                    pct_urban = 0, exposure = 1, country = fitd$country[1])
  nd1 <- nd0; nd1$pfpr10 <- p / 10
  ex  <- c("s(country)", "s(country,pfpr10)")
  e0  <- predict(m_te, nd0, type = "link", exclude = ex)
  e1  <- predict(m_te, nd1, type = "link", exclude = ex)
  as.numeric(1 - exp(-(e1 - e0)))
}
cat("\n--- te-implied attributable fraction vs 1% at PfPR2-10 = 30% ---\n")
for (yr in c(2000, 2012, 2024)) cat(sprintf("  %d: AF(30%%) = %.1f%%\n", yr, 100 * af_te(30, yr)))

## ---- figure: te prevalence-response at representative years -----------------
yrs <- c(2005, 2012, 2019, 2024)
grid <- do.call(rbind, lapply(yrs, function(yr) {
  g <- data.frame(pfpr2_10 = seq(1, 70, 0.5)); g$pfpr10 <- g$pfpr2_10 / 10
  g$year_c <- yr - center; g$dtp3 <- mean(fitd$dtp3); g$log_gdp <- mean(fitd$log_gdp)
  g$pct_urban <- mean(fitd$pct_urban); g$exposure <- 1000; g$country <- fitd$country[1]
  g$rate <- predict(m_te, g, type = "response", exclude = c("s(country)", "s(country,pfpr10)"))
  g$year <- factor(yr); g
}))
p <- ggplot(grid, aes(pfpr2_10, rate, colour = year)) +
  geom_point(data = fitd, aes(pfpr2_10, m1mo5y), inherit.aes = FALSE,
             alpha = 0.15, size = 0.8, colour = "grey40") +
  geom_line(linewidth = 1) +
  scale_colour_viridis_d(end = 0.9, name = "Year") +
  labs(x = expression(PfPR[2-10]~"(%)"),
       y = "Predicted post-neonatal deaths per 1000 (te surface)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave(file.path(RESULTS, "te_pfpr_year_expanded.png"), p, width = 7.5, height = 6, dpi = 300)
cat("\nsaved: results/te_pfpr_year_expanded.png\n")
saveRDS(list(m_add = m_add, m_te = m_te, center = center, tab = tab),
        file.path(RESULTS, "te_models.rds"))

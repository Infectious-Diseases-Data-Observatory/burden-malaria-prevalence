# =============================================================================
# 22_temporal_te.R — temporal malaria burden under the FULL TENSOR model
# (m_te from R/20): te(pfpr10, year_c) lets the prevalence response bend in
# BOTH prevalence and calendar year. Plotted alongside m_add (time-stable
# response, R/21) so the effect of allowing a time-varying surface is visible.
#
# Effects used (same convention as R/19 & R/21): the AF vs a 1% counterfactual
# uses the POPULATION (fixed/mean) surface PLUS the country RANDOM SLOPES
# [s(country,pfpr10,bs="re")]. Country random INTERCEPTS, covariates and the
# year main effect all cancel in the within-year ratio, so they do not enter.
#   AF_c(p,yr) = 1 - exp(-[ (g(p,yr) - g(1,yr)) + b_c*(p-1)/10 ]).
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})
mods <- readRDS(file.path(RESULTS, "te_models.rds")); mt <- mods$m_te; center <- mods$center

## ---- country random slopes from the te fit, keyed by iso3 -------------------
cf <- coef(mt)
sm <- mt$smooth[[which(vapply(mt$smooth, function(s) s$label, "") == "s(country,pfpr10)")]]
bsl <- setNames(as.numeric(cf[sm$first.para:sm$last.para]), levels(mt$model$country))

## ---- te population surface AF: g(p,yr) - g(1,yr), by year -------------------
EX  <- c("s(country)", "s(country,pfpr10)")
nd  <- function(p10, yr) data.frame(pfpr10 = p10, year_c = yr - center,
  dtp3 = mean(mt$model$dtp3), log_gdp = mean(mt$model$log_gdp),
  pct_urban = mean(mt$model$pct_urban), exposure = 1, country = levels(mt$model$country)[1])
# vectorised over (pfpr, year) pairs
gdiff <- function(pfpr, year) {
  hi <- data.frame(pfpr10 = pfpr/10, year_c = year - center, dtp3 = mean(mt$model$dtp3),
    log_gdp = mean(mt$model$log_gdp), pct_urban = mean(mt$model$pct_urban),
    exposure = 1, country = levels(mt$model$country)[1])
  lo <- hi; lo$pfpr10 <- 0.1
  as.numeric(predict(mt, hi, type = "link", exclude = EX) -
             predict(mt, lo, type = "link", exclude = EX))
}

## ---- malaria deaths per country-year ----------------------------------------
ts  <- read.csv(file.path(DATA, "wb_mortality_timeseries.csv"))
pcy <- read.csv(file.path(DATA, "pfpr_by_country_year.csv"))
d <- merge(ts[, c("iso3","year","allcause_1mo5y")], pcy, by = c("iso3","year"))
bc <- unname(bsl[d$iso3]); bc[is.na(bc)] <- 0
d$af  <- pmax(ifelse(d$pfpr_pct >= 1,
                     1 - exp(-(gdiff(d$pfpr_pct, d$year) + bc * (d$pfpr_pct - 1)/10)), 0), 0)
d$mal <- d$af * d$allcause_1mo5y
agg <- aggregate(mal ~ year, d, sum)
o <- function(y) agg$mal[agg$year == y]
cat(sprintf("Post-neonatal malaria deaths (m_te): 2000=%.0f, 2012=%.0f, 2024=%.0f  (2000->2024 %+.0f%%)\n",
            o(2000), o(2012), o(2024), 100 * (o(2024)/o(2000) - 1)))
cat("te population AF vs 1% at PfPR2-10=30%: ")
cat(sprintf("2000=%.1f%%, 2012=%.1f%%, 2024=%.1f%%\n",
            100*(1-exp(-gdiff(30,2000))), 100*(1-exp(-gdiff(30,2012))), 100*(1-exp(-gdiff(30,2024)))))

## ---- assemble comparison figure: te + m_add + WHO + IHME --------------------
add <- read.csv(file.path(RESULTS, "malaria_deaths_timeseries_add.csv"))   # m_add trajectory (R/21)
U5 <- 0.75
who_af <- data.frame(year = 2000:2024,
  point = c(804,815,786,758,751,712,723,704,666,673,650,618,578,556,550,550,546,548,552,545,598,577,573,567,579)*1000)
ihme <- tryCatch(read.csv(file.path(RESULTS, "ihme_u5_deaths_ssa_timeseries.csv")), error = function(e) NULL)
TE <- "te(PfPR, year) — time-varying response"; AD <- "s(PfPR)+s(year) — time-stable response"
WHOL <- "WHO — African-region under-5"; IHL <- "IHME/GBD — SSA under-5"
pl <- rbind(data.frame(year = agg$year, series = TE, deaths = agg$mal),
            data.frame(year = add$year, series = AD, deaths = add$mal),
            data.frame(year = who_af$year, series = WHOL, deaths = who_af$point*U5))
if (!is.null(ihme)) pl <- rbind(pl, data.frame(year = ihme$year, series = IHL, deaths = ihme$point))
pl <- pl[pl$year <= 2024, ]; lev <- c(TE, AD, WHOL, IHL); pl$series <- factor(pl$series, levels = lev)
p <- ggplot(pl, aes(year, deaths/1000, colour = series)) +
  geom_line(aes(linetype = series), linewidth = 1) + geom_point(size = 1.1) +
  scale_colour_manual(values = setNames(c("#d73027","#fc8d59","grey35","#238b45"), lev), name = NULL) +
  scale_linetype_manual(values = setNames(c("solid","solid","22","44"), lev), name = NULL) +
  scale_x_continuous(breaks = c(2000,2005,2010,2015,2020,2024)) + expand_limits(y = 0) +
  labs(x = NULL, y = "Malaria-attributable child deaths (thousands/yr)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.title = element_text(size = 13), axis.text = element_text(size = 11), legend.text = element_text(size = 9))
ggsave(file.path(RESULTS, "malaria_deaths_timeseries_te.png"), p, width = 11.5, height = 6.6, dpi = 300)
write.csv(agg, file.path(RESULTS, "malaria_deaths_timeseries_te.csv"), row.names = FALSE)
cat("saved: results/malaria_deaths_timeseries_te.{png,csv}\n")

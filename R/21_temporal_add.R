# =============================================================================
# 21_temporal_add.R — temporal malaria burden under the BEST-FITTING model
# (m_add from R/20): smooth, NONLINEAR, but TIME-STABLE prevalence response
#   deaths ~ s(pfpr10) + s(year_c) + covars + country RE + country random slope.
# No prevalence x year interaction (AIC-preferred over m_int and m_te).
#
# AF vs a 1% counterfactual depends only on the prevalence terms (year main
# effect, covariates, country intercept all cancel within a year):
#   AF_c(p) = 1 - exp(-[ (f(p) - f(1)) + b_c*(p-1)/10 ]),
# where f() is the population smooth s(pfpr10) and b_c the country random slope.
# Because f() has no year term, the AF at a given prevalence is the SAME every
# year -> the trajectory is driven purely by falling prevalence + falling
# all-cause mortality. Contrast R/19 (steepening slope -> flat burden).
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})
mods <- readRDS(file.path(RESULTS, "te_models.rds")); ma <- mods$m_add; center <- mods$center

## ---- country random slopes (linear deviations), keyed by iso3 ---------------
cf <- coef(ma)
sm <- ma$smooth[[which(vapply(ma$smooth, function(s) s$label, "") == "s(country,pfpr10)")]]
bsl <- setNames(as.numeric(cf[sm$first.para:sm$last.para]), levels(ma$model$country))

## ---- population smooth AF: f(p) - f(1), referenced to 1% --------------------
# predict at prevalence p vs 1% with everything else held equal & country
# excluded; the difference isolates the nonlinear population prevalence term.
base_nd <- function(p10) data.frame(pfpr10 = p10, year_c = 0, dtp3 = mean(ma$model$dtp3),
  log_gdp = mean(ma$model$log_gdp), pct_urban = mean(ma$model$pct_urban),
  exposure = 1, country = levels(ma$model$country)[1])
EX  <- c("s(country)", "s(country,pfpr10)")
e1  <- predict(ma, base_nd(0.1), type = "link", exclude = EX)          # f(1%)
fpop <- function(p) as.numeric(predict(ma, base_nd(p/10), type = "link", exclude = EX) - e1)

af_c <- function(iso3, pfpr) {                                         # country-specific, time-invariant
  bc <- unname(bsl[iso3]); bc[is.na(bc)] <- 0
  pmax(ifelse(pfpr >= 1, 1 - exp(-(fpop(pfpr) + bc * (pfpr - 1) / 10)), 0), 0)
}

## ---- population AF at prevalence anchors (same every year) ------------------
cat("Population AF vs 1% (time-invariant, country slope = 0):\n")
for (p in c(5, 10, 20, 30, 50)) cat(sprintf("  PfPR2-10=%2d%%: AF=%.1f%%\n", p, 100 * (1 - exp(-fpop(p)))))

## ---- malaria deaths per country-year ----------------------------------------
ts  <- read.csv(file.path(DATA, "wb_mortality_timeseries.csv"))
pcy <- read.csv(file.path(DATA, "pfpr_by_country_year.csv"))
d <- merge(ts[, c("iso3","year","allcause_1mo5y")], pcy, by = c("iso3","year"))
d$af  <- mapply(af_c, d$iso3, d$pfpr_pct)
d$mal <- d$af * d$allcause_1mo5y
agg <- aggregate(mal ~ year, d, sum)
o <- function(y) agg$mal[agg$year == y]
cat(sprintf("\nPost-neonatal malaria deaths (m_add): 2000=%.0f, 2012=%.0f, 2024=%.0f  (2000->2024 %+.0f%%)\n",
            o(2000), o(2012), o(2024), 100 * (o(2024) / o(2000) - 1)))

## ---- reference series (WHO WMR 2025 African-region; IHME/GBD SSA) ------------
U5 <- 0.75
who_af <- data.frame(year = 2000:2024,
  point = c(804,815,786,758,751,712,723,704,666,673,650,618,578,556,550,550,546,548,552,545,598,577,573,567,579)*1000)
ihme <- tryCatch(read.csv(file.path(RESULTS, "ihme_u5_deaths_ssa_timeseries.csv")), error = function(e) NULL)
OPN <- "Prevalence/all-cause mortality (m_add: time-stable response)"
WHOL <- "WHO — African-region under-5"; IHL <- "IHME/GBD — SSA under-5"
pl <- rbind(data.frame(year = agg$year, series = OPN, deaths = agg$mal),
            data.frame(year = who_af$year, series = WHOL, deaths = who_af$point * U5))
if (!is.null(ihme)) pl <- rbind(pl, data.frame(year = ihme$year, series = IHL, deaths = ihme$point))
pl <- pl[pl$year <= 2024, ]; lev <- c(OPN, WHOL, IHL); pl$series <- factor(pl$series, levels = lev)
p <- ggplot(pl, aes(year, deaths / 1000, colour = series)) +
  geom_line(aes(linetype = series), linewidth = 1) + geom_point(size = 1.1) +
  scale_colour_manual(values = setNames(c("#d73027","grey35","#238b45"), lev), name = NULL) +
  scale_linetype_manual(values = setNames(c("solid","22","44"), lev), name = NULL) +
  scale_x_continuous(breaks = c(2000,2005,2010,2015,2020,2024)) + expand_limits(y = 0) +
  labs(x = NULL, y = "Malaria-attributable child deaths (thousands/yr)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.title = element_text(size = 13), axis.text = element_text(size = 11), legend.text = element_text(size = 9))
ggsave(file.path(RESULTS, "malaria_deaths_timeseries_add.png"), p, width = 11.5, height = 6.4, dpi = 300)
write.csv(agg, file.path(RESULTS, "malaria_deaths_timeseries_add.csv"), row.names = FALSE)
cat("saved: results/malaria_deaths_timeseries_add.{png,csv}\n")

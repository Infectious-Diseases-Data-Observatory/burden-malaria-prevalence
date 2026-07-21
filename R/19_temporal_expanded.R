# =============================================================================
# 19_temporal_expanded.R — temporal malaria burden from the EXPANDED 2000-2024
# panel model (R/18). Uses country random slopes + prevalence x year interaction,
# now estimated from real data back to 2000 (no pre-2009 slope cap). AF referenced
# to a 1% counterfactual: beta_ct = beta_pop + b_country + gamma*year_c.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})
mods <- readRDS(file.path(RESULTS, "expanded_models.rds")); mi <- mods$m_int; center <- mods$center
cf <- coef(mi); b_pop <- unname(cf["pfpr10"]); g <- unname(cf["pfpr10:year_c"])
sm <- mi$smooth[[which(vapply(mi$smooth, function(s) s$label, "") == "s(country,pfpr10)")]]
bsl <- setNames(as.numeric(cf[sm$first.para:sm$last.para]), levels(mi$model$country))

ts  <- read.csv(file.path(DATA, "wb_mortality_timeseries.csv"))
pcy <- read.csv(file.path(DATA, "pfpr_by_country_year.csv"))
d <- merge(ts[, c("iso3","year","allcause_1mo5y")], pcy, by = c("iso3","year"))
d$year_c <- d$year - center                              # no cap: panel covers 2000-2024
d$bc <- unname(bsl[d$iso3]); d$bc[is.na(d$bc)] <- 0
d$slope <- b_pop + d$bc + g * d$year_c
d$af  <- pmax(ifelse(d$pfpr_pct >= 1, 1 - exp(-d$slope * (d$pfpr_pct - 1)/10), 0), 0)
d$mal <- d$af * d$allcause_1mo5y
agg <- aggregate(mal ~ year, d, sum)
sl <- function(y) (exp(b_pop + g * (y - center)) - 1) * 100
cat(sprintf("expanded-model population slope: 2000=%+.1f%%, 2012=%+.1f%%, 2024=%+.1f%% per +10 pts (interaction %.5f/yr)\n",
            sl(2000), sl(2012), sl(2024), g))
o <- function(y) agg$mal[agg$year == y]
cat(sprintf("Post-neonatal malaria deaths: 2000=%.0f, 2012=%.0f, 2024=%.0f  (2000->2024 %+.0f%%)\n",
            o(2000), o(2012), o(2024), 100*(o(2024)/o(2000)-1)))

U5 <- 0.75
who_af <- data.frame(year = 2000:2024,
  point = c(804,815,786,758,751,712,723,704,666,673,650,618,578,556,550,550,546,548,552,545,598,577,573,567,579)*1000)
ihme <- tryCatch(read.csv(file.path(RESULTS, "ihme_u5_deaths_ssa_timeseries.csv")), error = function(e) NULL)
OPN <- "Prevalence/all-cause mortality (expanded 2000-2024 panel)"
WHOL <- "WHO — African-region under-5"; IHL <- "IHME/GBD — SSA under-5"
pl <- rbind(data.frame(year = agg$year, series = OPN, deaths = agg$mal),
            data.frame(year = who_af$year, series = WHOL, deaths = who_af$point*U5))
if (!is.null(ihme)) pl <- rbind(pl, data.frame(year = ihme$year, series = IHL, deaths = ihme$point))
pl <- pl[pl$year <= 2024, ]; lev <- c(OPN, WHOL, IHL); pl$series <- factor(pl$series, levels = lev)
p <- ggplot(pl, aes(year, deaths/1000, colour = series)) +
  geom_line(aes(linetype = series), linewidth = 1) + geom_point(size = 1.1) +
  scale_colour_manual(values = setNames(c("#d73027","grey35","#238b45"), lev), name = NULL) +
  scale_linetype_manual(values = setNames(c("solid","22","44"), lev), name = NULL) +
  scale_x_continuous(breaks = c(2000,2005,2010,2015,2020,2024)) + expand_limits(y = 0) +
  labs(x = NULL, y = "Malaria-attributable child deaths (thousands/yr)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.title = element_text(size = 13), axis.text = element_text(size = 11), legend.text = element_text(size = 10))
ggsave(file.path(RESULTS, "malaria_deaths_timeseries_expanded.png"), p, width = 11.5, height = 6.4, dpi = 300)
cat("saved: results/malaria_deaths_timeseries_expanded.png\n")

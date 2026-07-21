# =============================================================================
# 11_malaria_deaths_timeseries.R — "our" longitudinal estimate of malaria-
# attributable child deaths in sub-Saharan Africa, 2000-2024, via the Component 2
# GAM. For each country-year: malaria deaths = AF(PfPR2-10) x all-cause deaths,
# where AF is the Component 2 negative-binomial GAM attributable fraction
# (1 - exp(-(s(PfPR) - s(0)))). Done for ALL under-5 and for 1mo-5y (neonatal-
# excluded); summed across countries and plotted over time. This is the version
# to be compared later against IHME/WHO longitudinal estimates.
#
# Inputs (cached in data/; heavy, so NOT part of run_all.R):
#   * wb_mortality_timeseries.csv  — WB U5MR/NNMR/CBR/pop -> all-cause deaths/yr
#   * pfpr_by_country_year.csv     — MAP PfPR2-10 pop-weighted per country-year
#     (25 annual MAP rasters, ~350 MB download, aggregated with GPW + admin-0)
# NOTE: population weights use the single GPW-2020 surface for every year (the
# pop-weighted national PfPR is insensitive to this; documented approximation).
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(malariaAtlas); library(terra); library(sf); library(ggplot2)})

## ---- 1. World Bank all-cause child mortality time series (cache) ------------
TS <- file.path(DATA, "wb_mortality_timeseries.csv")
if (!file.exists(TS)) {
  message("Fetching World Bank child-mortality time series 2000-2024 ...")
  ind <- c(u5mr = "SH.DYN.MORT", nnmr = "SH.DYN.NMRT", cbr = "SP.DYN.CBRT.IN", pop = "SP.POP.TOTL")
  tsd <- NULL
  for (nm in names(ind)) { x <- wb_fetch(ind[[nm]], "2000:2024")[, c("iso3","year","value")]; names(x)[3] <- nm
    tsd <- if (is.null(tsd)) x else merge(tsd, x, by = c("iso3","year"), all = TRUE) }
  tsd$region <- countrycode::countrycode(tsd$iso3, "iso3c", "region", warn = FALSE)
  tsd <- tsd[!is.na(tsd$region) & tsd$region == "Sub-Saharan Africa", ]
  tsd$births         <- tsd$cbr/1000 * tsd$pop
  tsd$allcause_u5    <- tsd$u5mr/1000 * tsd$births
  tsd$allcause_1mo5y <- (tsd$u5mr - tsd$nnmr)/1000 * tsd$births
  write.csv(tsd[complete.cases(tsd[, c("u5mr","nnmr","births")]), ], TS, row.names = FALSE)
}
ts <- read.csv(TS, stringsAsFactors = FALSE)

## ---- 2. MAP national PfPR2-10 per country-year (cache) ----------------------
PCY <- file.path(DATA, "pfpr_by_country_year.csv")
if (!file.exists(PCY)) {
  message("Downloading + aggregating MAP annual PfPR2-10 rasters 2000-2024 (~350 MB) ...")
  dl <- file.path(tempdir(), "map_years"); dir.create(dl, showWarnings = FALSE)
  ext <- matrix(c(-18, -35, 52, 38), nrow = 2)                       # Africa bbox
  invisible(tryCatch(getRaster(dataset_id = "Malaria__202508_Global_Pf_Parasite_Rate",
                               year = 2000:2024, extent = ext, file_path = dl),
                     error = function(e) message("note: ", conditionMessage(e))))
  tifs   <- list.files(dl, pattern = "\\.tiff?$", full.names = TRUE)
  pf     <- terra::rast(tifs)
  mean_l <- pf[[grep("_1$", names(pf))]]                             # band 1 = mean estimate, per year
  yrs    <- as.integer(sub(".*_(\\d{4})-.*", "\\1", names(mean_l)))
  den    <- terra::resample(terra::crop(terra::rast(GPW_TIF), mean_l[[1]]), mean_l[[1]], method = "bilinear")
  v      <- terra::makeValid(terra::vect(sf::st_make_valid(readRDS(file.path(DATA, "africa_admin0.rds")))))
  iso    <- terra::values(v)$iso
  rows <- lapply(seq_len(terra::nlyr(mean_l)), function(i) {
    pfy <- mean_l[[i]]; w <- terra::mask(den, pfy)                   # pop weight where PfPR defined
    an <- terra::extract(pfy * w, v, fun = sum, na.rm = TRUE, ID = FALSE)[[1]]
    aw <- terra::extract(w,       v, fun = sum, na.rm = TRUE, ID = FALSE)[[1]]
    data.frame(iso3 = iso, year = yrs[i], pfpr_pct = 100 * an / aw) })
  pcy <- do.call(rbind, rows); pcy <- pcy[is.finite(pcy$pfpr_pct), ]
  write.csv(pcy, PCY, row.names = FALSE)
}
pcy <- read.csv(PCY, stringsAsFactors = FALSE)

## ---- 3. Component 2 GAM attributable fraction AF(PfPR2-10), by outcome ------
d2 <- read.csv(file.path(RESULTS, "component2_region_data.csv"), stringsAsFactors = FALSE)
af_fun <- function(oc) {                                             # primary linear-model AF + coefficient CI
  f <- fit_c2_primary(d2, oc); b <- f$beta; se <- f$se               # M1: linear nb-GAM, PfPR2-10 >= 1%
  list(af = function(p) af_c2(b, p),                                 # AF referenced to 1% PfPR counterfactual
       lo = function(p) af_c2(b - 1.96 * se, p),
       hi = function(p) af_c2(b + 1.96 * se, p)) }
AFU <- af_fun("u5mr"); AFP <- af_fun("m1mo5y")

## ---- 4. malaria deaths per country-year, both denominators ------------------
d <- merge(ts[, c("iso3","year","allcause_u5","allcause_1mo5y")], pcy, by = c("iso3","year"))
d$mal_u5    <- AFU$af(d$pfpr_pct)  * d$allcause_u5;    d$mal_u5_lo <- AFU$lo(d$pfpr_pct)*d$allcause_u5;    d$mal_u5_hi <- AFU$hi(d$pfpr_pct)*d$allcause_u5
d$mal_1mo5y <- AFP$af(d$pfpr_pct)  * d$allcause_1mo5y; d$mal_pn_lo <- AFP$lo(d$pfpr_pct)*d$allcause_1mo5y; d$mal_pn_hi <- AFP$hi(d$pfpr_pct)*d$allcause_1mo5y
write.csv(d, file.path(RESULTS, "malaria_deaths_timeseries_by_country.csv"), row.names = FALSE)

agg <- aggregate(cbind(mal_u5, mal_u5_lo, mal_u5_hi, mal_1mo5y, mal_pn_lo, mal_pn_hi) ~ year, d, sum)
tot <- rbind(
  data.frame(year = agg$year, outcome = "All under-5",             deaths = agg$mal_u5,    lo = agg$mal_u5_lo, hi = agg$mal_u5_hi),
  data.frame(year = agg$year, outcome = "1mo-5y (neonatal excl.)", deaths = agg$mal_1mo5y, lo = agg$mal_pn_lo, hi = agg$mal_pn_hi))
write.csv(tot, file.path(RESULTS, "malaria_deaths_timeseries_total.csv"), row.names = FALSE)
cat(sprintf("countries summed: %d | years %d-%d\n", length(unique(d$iso3)), min(d$year), max(d$year)))
cat(sprintf("All-U5 malaria deaths: 2000=%.0f, 2024=%.0f (%+.0f%%)\n",
    agg$mal_u5[agg$year==2000], agg$mal_u5[agg$year==2024], 100*(agg$mal_u5[agg$year==2024]/agg$mal_u5[agg$year==2000]-1)))

## ---- 5. WHO reference (WMR 2025) --------------------------------------------
# Table 2.1 = GLOBAL, all-ages deaths (context, saved). Table 2.4 = WHO AFRICAN
# REGION, all-ages deaths; x the ~75% under-5 share (WMR: "just over 75% of all
# deaths in the region are of children aged under 5") = like-for-like with ours.
who_glob <- data.frame(year = 2000:2024,
  point = c(864,873,840,811,806,767,771,747,708,715,693,655,610,583,579,578,576,574,575,567,621,601,598,598,610)*1000,
  lo    = c(833,840,809,781,770,734,739,717,679,683,659,625,583,554,546,543,542,540,536,527,575,558,554,550,561)*1000,
  hi    = c(904,916,881,854,863,815,818,790,745,762,744,696,651,625,632,635,634,638,649,649,736,716,722,725,738)*1000)
write.csv(who_glob, file.path(RESULTS, "who_wmr2025_global_deaths.csv"), row.names = FALSE)
who_af <- data.frame(year = 2000:2024,
  point = c(804,815,786,758,751,712,723,704,666,673,650,618,578,556,550,550,546,548,552,545,598,577,573,567,579)*1000,
  lo    = c(781,789,761,733,722,686,695,677,641,646,620,592,552,526,518,517,513,513,514,505,553,535,530,520,531)*1000,
  hi    = c(834,851,819,792,803,755,766,742,700,716,698,658,616,598,601,605,602,608,624,626,711,691,698,694,706)*1000)
write.csv(who_af, file.path(RESULTS, "who_wmr2025_africa_deaths.csv"), row.names = FALSE)
U5 <- 0.75                                             # constant U5 share (WMR); flagged approximation

## ---- IHME/GBD U5 malaria deaths, SSA (user-supplied export; access-controlled) ----
ihme_f <- list.files(DATA, pattern = "export_trigger.*Sub-Saharan|Data Explorer.*Sub-Saharan", full.names = TRUE)
ihme <- NULL
if (length(ihme_f)) {
  xi <- read.csv(ihme_f[1], stringsAsFactors = FALSE, check.names = FALSE)
  xi <- xi[xi$Measure == "Deaths" & xi$Age == "Under 5" & xi$Unit == "Number" & xi$Location == "Sub-Saharan Africa", ]
  ihme <- data.frame(year = xi$Year, point = xi$Value, lo = xi$Lower, hi = xi$Upper)
  ihme <- ihme[order(ihme$year), ]
  write.csv(ihme, file.path(RESULTS, "ihme_u5_deaths_ssa_timeseries.csv"), row.names = FALSE)
} else message("IHME SSA export not found in data/ — plotting ours + WHO only.")

## ---- 6. plot: our SSA under-5 vs WHO & IHME (like-for-like under-5) ----------
OPN  <- "Prevalence/all-cause mortality method"          # our post-neonatal estimate
WHOL <- "WHO — African-region under-5 (WMR 2025)"; IHL <- "IHME/GBD — SSA under-5"
tot_pn <- tot[tot$outcome == "1mo-5y (neonatal excl.)", ]  # post-neonatal only (drop all-U5 series)
pl <- rbind(
  data.frame(year = tot_pn$year, series = OPN, deaths = tot_pn$deaths, lo = tot_pn$lo, hi = tot_pn$hi),
  data.frame(year = who_af$year, series = WHOL, deaths = who_af$point*U5, lo = who_af$lo*U5, hi = who_af$hi*U5))
if (!is.null(ihme)) pl <- rbind(pl, data.frame(year = ihme$year, series = IHL, deaths = ihme$point, lo = ihme$lo, hi = ihme$hi))
lev <- c(OPN, WHOL, IHL); pl$series <- factor(pl$series, levels = lev)
cols <- setNames(c("#d73027","grey35","#238b45"), lev)
lty  <- setNames(c("solid","22","44"), lev)
pl <- pl[pl$year <= 2024, ]                              # stop the series at 2024 (drop any 2025 point)

p <- ggplot(pl, aes(year, deaths/1000, colour = series)) +
  geom_line(aes(linetype = series), linewidth = 1) + geom_point(size = 1.1) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_linetype_manual(values = lty, name = NULL) +
  scale_x_continuous(breaks = c(2000, 2005, 2010, 2015, 2020, 2024),
                     expand = expansion(mult = c(0.01, 0.02))) + expand_limits(y = 0) +
  labs(x = NULL, y = "Malaria-attributable child deaths (thousands/yr)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.title = element_text(size = 14), axis.text = element_text(size = 12),
        legend.text = element_text(size = 12))
ggsave(file.path(RESULTS, "malaria_deaths_timeseries.png"), p, width = 11.5, height = 6.4, dpi = 300)
if (!is.null(ihme)) cat(sprintf("IHME SSA U5: 2000=%.0f, 2024=%.0f (%+.0f%%)\n",
    ihme$point[ihme$year==2000], ihme$point[ihme$year==2024], 100*(ihme$point[ihme$year==2024]/ihme$point[ihme$year==2000]-1)))
cat("saved: results/malaria_deaths_timeseries.png + timeseries CSVs + who_wmr2025_{global,africa}_deaths.csv + ihme_u5_deaths_ssa_timeseries.csv\n")

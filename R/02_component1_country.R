# =============================================================================
# 02_component1_country.R — COUNTRY LEVEL
# Malaria's share of all-cause under-5 deaths (all U5, and neonatal-excluded)
# vs national PfPR2-10, plus a multivariable model adding GDP per capita and
# DTP3 coverage, and identification of outlier countries where malaria is a
# higher % of child deaths than prevalence predicts.
#
# Numerator : IHME/GBD under-5 malaria deaths (2025 export).
# Denominator: IGME all-cause U5 deaths = U5MR/1000 x live births (and the
#              neonatal-excluded 1mo-5y version = (U5MR-NNMR)/1000 x births).
# Exposure  : population-weighted national PfPR2-10 (MAP, 2024).
# Covariates: log GDP per capita, DTP3 coverage (World Bank, latest year).
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(malariaAtlas); library(ggplot2); library(patchwork)})

## ---- national PfPR2-10 per country (cached; raster op is slow) --------------
PFPR_CSV <- file.path(DATA, "pfpr_by_country.csv")
if (!file.exists(PFPR_CSV)) {
  message("Computing population-weighted national PfPR2-10 ...")
  a <- sf::st_make_valid(readRDS(file.path(DATA, "africa_admin0.rds")))
  v <- terra::makeValid(terra::vect(a))
  a$pfpr_pct <- pfpr_by_admin(read_pfpr(), v)
  write.csv(data.frame(iso3 = a$iso, pfpr_pct = a$pfpr_pct)[is.finite(a$pfpr_pct), ], PFPR_CSV, row.names = FALSE)
}
pfpr <- read.csv(PFPR_CSV, stringsAsFactors = FALSE)

## ---- assemble country table -------------------------------------------------
ihme   <- parse_ihme(file.path(DATA, "ihme_malaria_u5_deaths_by_country.csv"))   # iso3, deaths, rate, year
igme   <- read.csv(file.path(DATA, "igme_mortality_by_country.csv"), stringsAsFactors = FALSE)
births <- read.csv(file.path(DATA, "wb_livebirths_by_country.csv"), stringsAsFactors = FALSE)[, c("iso3", "births")]
latest <- function(f, col) { d <- read.csv(f, stringsAsFactors = FALSE); d <- d[order(d$iso3, -d$year), ]
                             d <- d[!duplicated(d$iso3), c("iso3", col)]; d }
gdp    <- latest(file.path(DATA, "wb_gdp_pc.csv"), "gdp_pc")
dtp3   <- latest(file.path(DATA, "wb_dtp3.csv"),   "dtp3")

m <- merge(ihme[, c("iso3", "deaths")], igme[, c("iso3", "u5mr", "m_1mo_5y")], by = "iso3")
m <- merge(m, births, by = "iso3"); m <- merge(m, pfpr, by = "iso3")
m <- merge(m, gdp, by = "iso3");    m <- merge(m, dtp3, by = "iso3")
m$region  <- countrycode::countrycode(m$iso3, "iso3c", "region", warn = FALSE)     # World Bank region
m <- m[!is.na(m$region) & m$region == "Sub-Saharan Africa", ]
m$country <- countrycode::countrycode(m$iso3, "iso3c", "country.name", warn = FALSE)

# malaria as % of all-cause U5 deaths (two age windows); covariates
m$share_u5     <- 100 * m$deaths / (m$u5mr    / 1000 * m$births)   # all under-5
m$share_1mo5y  <- 100 * m$deaths / (m$m_1mo_5y / 1000 * m$births)  # neonatal excluded
m$log_gdp      <- log(m$gdp_pc)
m <- m[is.finite(m$share_u5) & is.finite(m$share_1mo5y) & is.finite(m$log_gdp) & is.finite(m$dtp3) & is.finite(m$pfpr_pct), ]
cat(sprintf("Component 1: %d SSA countries (IHME year %s, MAP PfPR %d)\n", nrow(m), unique(ihme$year), MAP_YEAR))

## ---- models -----------------------------------------------------------------
# "prevalence alone" (defines outliers) and the multivariable adjustment.
coefs <- list(); m$resid_u5 <- m$resid_1mo5y <- NA_real_
for (out in c("share_u5", "share_1mo5y")) {
  uni  <- lm(reformulate("pfpr_pct", out), data = m)                       # prevalence only
  mult <- lm(reformulate(c("pfpr_pct", "log_gdp", "dtp3"), out), data = m) # + GDP + DTP3
  m[[paste0("resid_", sub("share_", "", out))]] <- rstudent(uni)           # outlier metric = prevalence-only residual
  ct <- as.data.frame(summary(mult)$coefficients); ct$term <- rownames(ct); ct$outcome <- out; ct$model <- "multivariable"
  cu <- as.data.frame(summary(uni)$coefficients);  cu$term <- rownames(cu); cu$outcome <- out; cu$model <- "prevalence_only"
  coefs[[out]] <- rbind(cu, ct)
  cat(sprintf("  %-11s: PfPR effect (prev-only) = %.2f %%/PfPR-pt (R2=%.2f);  adj. for GDP+DTP3 = %.2f (R2=%.2f)\n",
              out, coef(uni)["pfpr_pct"], summary(uni)$r.squared, coef(mult)["pfpr_pct"], summary(mult)$r.squared))
}
coef_tab <- do.call(rbind, coefs); names(coef_tab)[1:4] <- c("estimate","std_error","t_value","p_value")
write.csv(coef_tab[, c("outcome","model","term","estimate","std_error","p_value")],
          file.path(RESULTS, "component1_model_coefficients.csv"), row.names = FALSE)

## ---- outliers (malaria a higher % of U5 deaths than prevalence predicts) ----
m$outlier <- m$resid_1mo5y > 2                                            # studentized resid > 2 (prevalence-only, 1mo-5y)
out_tab <- m[order(-m$resid_1mo5y), c("iso3","country","pfpr_pct","share_u5","share_1mo5y","resid_1mo5y","resid_u5","log_gdp","dtp3")]
write.csv(out_tab, file.path(RESULTS, "component1_outliers.csv"), row.names = FALSE)
write.csv(m[, c("iso3","country","pfpr_pct","deaths","u5mr","m_1mo_5y","births","gdp_pc","dtp3","share_u5","share_1mo5y")],
          file.path(RESULTS, "component1_country_data.csv"), row.names = FALSE)
cat("outliers (resid>2, 1mo-5y): ", paste(m$country[m$outlier], collapse = ", "), "\n")

## ---- figure -----------------------------------------------------------------
mk <- function(yv, ylab) ggplot(m, aes(pfpr_pct, .data[[yv]])) +
  geom_smooth(method = "lm", se = TRUE, colour = "#2c7fb8", fill = "#c7e0ee", linewidth = 0.6) +
  geom_point(aes(colour = outlier), size = 2.2) +
  geom_text(data = m[m$outlier, ], aes(label = country), size = 2.6, vjust = -0.8, colour = "#d73027") +
  scale_colour_manual(values = c("FALSE" = "#08519c", "TRUE" = "#d73027"), guide = "none") +
  labs(x = expression("National "*PfPR[2-10]*" (%, MAP 2024)"), y = ylab) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
p <- mk("share_1mo5y", "Malaria % of 1mo-5y deaths (neonatal excl.)") +
     mk("share_u5",    "Malaria % of all U5 deaths") +
     plot_annotation(
       title = "Component 1 — Country level: malaria's share of child deaths vs prevalence",
       subtitle = "Line = fit on prevalence alone; red = outliers (studentized resid > 2) where malaria exceeds the prevalence prediction. Multivariable model adjusts for GDP p.c. + DTP3.",
       caption = "IHME/GBD U5 malaria deaths (2025); IGME all-cause U5 mortality; MAP PfPR2-10 (2024); World Bank GDP p.c. & DTP3.")
ggsave(file.path(RESULTS, "component1_share_vs_pfpr.png"), p, width = 13, height = 6.5, dpi = 140)
cat("saved: results/component1_share_vs_pfpr.png + component1_{model_coefficients,outliers,country_data}.csv\n")

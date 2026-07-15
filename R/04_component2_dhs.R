# =============================================================================
# 04_component2_dhs.R — DHS SUBNATIONAL + NATIONAL
# Association between child malaria prevalence and all-cause under-5 mortality
# across all usable SSA DHS/MIS survey-regions, adjusted for DTP3, GDP per
# capita, % urban, calendar year and stunting, via a mixed-effects model with
# country random slopes. Outcome on the log scale so the malaria coefficient
# reads as a % change in mortality. Fit for U5MR and for 1mo-5y (neonatal excl.).
#
# Exposure = microscopy-equivalent prevalence (RDT->microscopy via Component 3
# conversion where microscopy missing), age-standardised 0.5-5y -> PfPR2-10.
# Prevalence + covariate table comes from Component 3; mortality from DHS.rates.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(DHS.rates); library(malariaAtlas); library(lme4); library(ggplot2)})
DHS <- file.path(DATA, "dhs")

## ---- inputs: per-region prevalence/covariates + RDT->micro conversion -------
prev <- read.csv(file.path(DATA, "dhs_prevalence_by_region.csv"), stringsAsFactors = FALSE)
prev$svkey <- paste0(substr(prev$survey, 1, 2), substr(prev$survey, 5, 8))   # country+phase, links PR<->BR
conv <- read.csv(file.path(RESULTS, "rdt_microscopy_conversion.csv"), stringsAsFactors = FALSE)
k <- conv$slope[conv$model == "through_origin"]                              # microscopy = k x RDT

## ---- assemble modelling table (cached: the chmort loop is slow) -------------
REGION_DATA <- file.path(RESULTS, "component2_region_data.csv")
if (file.exists(REGION_DATA)) {
  d <- read.csv(REGION_DATA, stringsAsFactors = FALSE)
} else {
  message("Computing region mortality (DHS.rates::chmort) ...")
  bases <- toupper(sub("[.]rds$", "", list.files(DHS, pattern = "rds$")))
  mrows <- list()
  for (b in bases[substr(bases, 3, 4) == "BR"]) {
    cc <- substr(b, 1, 2); rvb <- if (cc == "NG") "sstate" else "v024"
    mo <- tryCatch(mort_by_region(readRDS(file.path(DHS, paste0(b, ".rds"))), rvb), error = function(e) NULL)
    if (is.null(mo)) next
    mo$svkey <- paste0(cc, substr(b, 5, 8)); mrows[[b]] <- mo
  }
  d <- merge(prev, do.call(rbind, mrows), by = c("svkey", "regkey"))

  # GDP p.c. and DTP3 matched to each survey's country & nearest year
  gdp <- read.csv(file.path(DATA, "wb_gdp_pc.csv"), stringsAsFactors = FALSE)
  dtp <- read.csv(file.path(DATA, "wb_dtp3.csv"),   stringsAsFactors = FALSE)
  nearest <- function(panel, iso, yr, col) { s <- panel[panel$iso3 == iso & is.finite(panel[[col]]), ]
    if (!nrow(s)) return(NA_real_); s[[col]][which.min(abs(s$year - yr))] }
  d$log_gdp <- log(mapply(function(i, y) nearest(gdp, i, y, "gdp_pc"), d$iso3, d$year))
  d$dtp3    <-     mapply(function(i, y) nearest(dtp, i, y, "dtp3"),   d$iso3, d$year)

  # exposure: microscopy-equivalent -> age-standardised PfPR2-10 (%)
  d$mic_eq   <- ifelse(is.finite(d$mic), d$mic, k * d$rdt)
  d$pfpr2_10 <- 100 * to_pfpr210(pmin(pmax(d$mic_eq, 0), 100) / 100)
  d$pfpr10   <- d$pfpr2_10 / 10                                 # per +10 PfPR2-10 points
  d$year_c   <- d$year - round(mean(d$year, na.rm = TRUE))

  # drop chmort failures (U5MR ~ 0) and non-finite exposure
  d <- d[is.finite(d$u5mr) & d$u5mr > 5 & is.finite(d$m1mo5y) & d$m1mo5y > 0 & is.finite(d$pfpr2_10), ]
  write.csv(d, REGION_DATA, row.names = FALSE)
}
cat(sprintf("Component 2: %d survey-regions, %d surveys, %d countries\n",
            nrow(d), length(unique(d$svkey)), length(unique(d$iso3))))

## ---- mixed models (log outcome; country random slope on prevalence) ---------
covs <- c("pfpr10", "dtp3", "log_gdp", "pct_urban", "year_c", "stunting")
fit_one <- function(outcome, dat = d) {
  dd <- dat[complete.cases(dat[, c(outcome, covs, "country", "svkey")]), ]
  f  <- as.formula(paste0("log(", outcome, ") ~ ",   # || = uncorrelated country intercept & slope
                          paste(c(covs, "(1 + pfpr10 || country)", "(1 | svkey)"), collapse = " + ")))
  m  <- lmer(f, data = dd, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  list(m = m, dd = dd)
}
summarise <- function(res, outcome) {
  m <- res$m; fe <- fixef(m); se <- sqrt(diag(vcov(m)))
  tab <- data.frame(outcome = outcome, term = names(fe), beta = as.numeric(fe), se = as.numeric(se),
                    pct_change_per_unit = (exp(fe) - 1) * 100,
                    p = 2 * pnorm(-abs(fe / se)), stringsAsFactors = FALSE)
  b <- fe["pfpr10"]; s <- se["pfpr10"]
  cat(sprintf("  %-7s: %+.1f%% U5 mortality per +10 PfPR2-10 pts (95%% CI %+.1f to %+.1f); n=%d, %d countries; singular=%s\n",
              outcome, (exp(b) - 1) * 100, (exp(b - 1.96 * s) - 1) * 100, (exp(b + 1.96 * s) - 1) * 100,
              nrow(res$dd), length(unique(res$dd$country)), isSingular(m)))
  tab
}
cat("MAIN model (all survey-regions):\n")
r_u5 <- fit_one("u5mr"); r_pn <- fit_one("m1mo5y")
coef_full <- rbind(summarise(r_u5, "u5mr"), summarise(r_pn, "m1mo5y")); coef_full$sample <- "full"

# sensitivity: restrict to mid-transmission regions (PfPR2-10 in [5, 50]%)
d_s <- d[d$pfpr2_10 >= 5 & d$pfpr2_10 <= 50, ]
cat(sprintf("\nSENSITIVITY — PfPR2-10 in [5,50]%% (%d of %d survey-regions):\n", nrow(d_s), nrow(d)))
s_u5 <- fit_one("u5mr", d_s); s_pn <- fit_one("m1mo5y", d_s)
coef_sens <- rbind(summarise(s_u5, "u5mr"), summarise(s_pn, "m1mo5y")); coef_sens$sample <- "pfpr_5_50"

write.csv(rbind(coef_full, coef_sens), file.path(RESULTS, "component2_model_coefficients.csv"), row.names = FALSE)

# country-specific prevalence slopes (per +10 pts) from the U5MR model
cc <- coef(r_u5$m)$country
csd <- data.frame(country = rownames(cc), pct_change_per10 = (exp(cc[, "pfpr10"]) - 1) * 100)
write.csv(csd[order(-csd$pct_change_per10), ], file.path(RESULTS, "component2_country_slopes.csv"), row.names = FALSE)

## ---- figure: adjusted country-specific prevalence effects -------------------
dd <- r_u5$dd; cf <- coef(r_u5$m)$country; cf$country <- rownames(cf)
xr <- range(dd$pfpr2_10); fe <- fixef(r_u5$m)
lines <- do.call(rbind, lapply(seq_len(nrow(cf)), function(i) data.frame(
  country = cf$country[i], x = xr,
  y = exp(cf[i, "(Intercept)"] + cf[i, "pfpr10"] * (xr / 10) +
          fe["dtp3"] * mean(dd$dtp3) + fe["log_gdp"] * mean(dd$log_gdp) +
          fe["pct_urban"] * mean(dd$pct_urban) + fe["stunting"] * mean(dd$stunting)))))
p <- ggplot() +
  geom_point(data = dd, aes(pfpr2_10, u5mr, colour = country), size = 1.1, alpha = 0.35) +
  geom_line(data = lines, aes(x, y, colour = country, group = country), linewidth = 0.6, alpha = 0.9) +
  annotate("text", x = xr[1], y = max(dd$u5mr), hjust = 0, vjust = 1, fontface = "bold", size = 3.4,
           label = sprintf("adjusted: %+.1f%% U5 mortality per +10 PfPR2-10 pts",
                           (exp(fixef(r_u5$m)["pfpr10"]) - 1) * 100)) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2, alpha = 1, linewidth = 1))) +
  labs(x = expression("Age-standardised "*PfPR[2-10]*" (%), microscopy-equivalent"),
       y = "All-cause under-5 mortality, 5q0 (per 1,000 live births)",
       title = "Component 2 — DHS: adjusted malaria-prevalence effect on child mortality",
       subtitle = "Lines = country-specific fitted U5MR vs prevalence at mean covariates (mixed model, country random slopes).",
       caption = "DHS.rates 5q0 & microscopy-equivalent PfPR2-10; adjusted for DTP3, GDP p.c., % urban, year, stunting.") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
    legend.key.size = unit(3.3, "mm"), legend.text = element_text(size = 6), legend.title = element_text(size = 8))
ggsave(file.path(RESULTS, "component2_country_slopes.png"), p, width = 12, height = 7.5, dpi = 140)
cat("saved: results/component2_country_slopes.png + component2_{model_coefficients,country_slopes,region_data}.csv\n")

# =============================================================================
# 13_ihme_share.R — Method-1 (IHME/GBD) comparison figure.
#
# Malaria's share of all-cause POST-NEONATAL under-5 deaths implied by the IHME
# estimates, by country-year (2000-2024), as a function of national PfPR2-10.
# Adjusted for the SAME national covariate block as the DHS/MIS survey analysis:
# WUENIC Hib3, PCV-completion and rotavirus coverage, child HIV prevalence
# (log), log GDP per capita and log health expenditure per capita; plus a
# country random effect and a calendar-year term for the panel structure.
#
#   Numerator  : IHME/GBD under-5 malaria deaths (Data Explorer export,
#                data/ihme_malaria_u5_deaths_by_country_year.csv).
#   Denominator: IHME/GBD all-cause POST-NEONATAL deaths = under-5 minus neonatal
#                (0-6 day + 7-27 day), from
#                data/ihme_allcause_by_age_country_year.csv.
# Both numerator and denominator are IHME estimates. IHME assigns ~no malaria
# deaths to neonates, so the under-5 malaria count is used for the numerator.
# This is a Method-1 comparison, independent of the DHS survey model.
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2", "countrycode"))

IHME_CSV <- file.path(DATA_DIR, "ihme_malaria_u5_deaths_by_country_year.csv")
if (!file.exists(IHME_CSV)) {
  message("IHME country-year file not found (", IHME_CSV, "); skipping the IHME share figure.")
  quit(save = "no", status = 0)
}

## ---- IHME under-5 malaria deaths (numerator), by country-year ---------------
ih <- read.csv(IHME_CSV, stringsAsFactors = FALSE, check.names = FALSE)
ih <- ih[ih$Measure == "Deaths" & ih$Unit == "Number" & ih$Age == "Under 5", ]
ih$iso3 <- countrycode::countrycode(ih$Location, "country.name", "iso3c", warn = FALSE)
ih$year <- as.integer(ih$Year)
ih$malaria_deaths <- suppressWarnings(as.numeric(ih$Value))
ih <- ih[!is.na(ih$iso3) & is.finite(ih$year) & is.finite(ih$malaria_deaths),
         c("iso3", "year", "malaria_deaths"), drop = FALSE]

## ---- IHME all-cause POST-NEONATAL deaths (denominator), by country-year -----
AC_CSV <- file.path(DATA_DIR, "ihme_allcause_by_age_country_year.csv")
if (!file.exists(AC_CSV)) {
  message("IHME all-cause file not found (", AC_CSV, "); skipping the IHME share figure.")
  quit(save = "no", status = 0)
}
ac <- read.csv(AC_CSV, stringsAsFactors = FALSE, check.names = FALSE)
ac <- ac[ac$Measure == "Deaths" & ac$Unit == "Number", ]
ac$iso3 <- countrycode::countrycode(ac$Location, "country.name", "iso3c", warn = FALSE)
ac$year <- suppressWarnings(as.integer(ac$Year))
ac$val <- suppressWarnings(as.numeric(ac$Value))
ac <- ac[!is.na(ac$iso3) & is.finite(ac$year) & is.finite(ac$val), , drop = FALSE]
ackey <- paste(ac$iso3, ac$year, sep = "|")
u5  <- tapply(ac$val[ac$Age == "Under 5"], ackey[ac$Age == "Under 5"], sum)
neo <- tapply(ac$val[grepl("neonatal", ac$Age)], ackey[grepl("neonatal", ac$Age)], sum)
den <- data.frame(key = names(u5), allcause_u5 = as.numeric(u5), stringsAsFactors = FALSE)
den$neonatal <- as.numeric(neo[den$key])
den$allcause_pn <- den$allcause_u5 - den$neonatal        # all-cause post-neonatal (IHME)
den$iso3 <- sub("\\|.*$", "", den$key)
den$year <- as.integer(sub("^.*\\|", "", den$key))
den <- den[is.finite(den$allcause_pn) & den$allcause_pn > 0, , drop = FALSE]

m <- merge(ih, den[, c("iso3", "year", "allcause_pn")], by = c("iso3", "year"))
m$share_pn <- 100 * m$malaria_deaths / m$allcause_pn     # malaria % of post-neonatal deaths (IHME/IHME)

## ---- national covariates (same block as the DHS analysis) -------------------
kk <- function(a, b) paste(a, b, sep = "|")
pf  <- read.csv(file.path(DATA_DIR, "pfpr_by_country_year.csv"), stringsAsFactors = FALSE)
imm <- read.csv(UNICEF_IMMUNISATION_CSV, stringsAsFactors = FALSE)[
  , c("iso3", "year", "hib3_wuenic", "pcv3_wuenic", "rotac_wuenic")]
hiv <- read.csv(HIV_PANEL_CSV, stringsAsFactors = FALSE)[, c("iso3", "year", "hiv_prev")]
gdp <- read.csv(file.path(DATA_DIR, "wb_gdp_pc.csv"), stringsAsFactors = FALSE)
hex <- read.csv(file.path(DATA_DIR, "wb_hexp_pc.csv"), stringsAsFactors = FALSE)

m <- merge(m, pf[, c("iso3", "year", "pfpr_pct")], by = c("iso3", "year"))
m <- merge(m, imm, by = c("iso3", "year"), all.x = TRUE)
m$hiv_prev <- hiv$hiv_prev[match(kk(m$iso3, m$year), kk(hiv$iso3, hiv$year))]
m$gdp_pc  <- mapply(function(i, y) nearest_panel_value(gdp, i, y, "gdp_pc"), m$iso3, m$year)
m$hexp_pc <- mapply(function(i, y) nearest_panel_value(hex, i, y, "hexp_pc"), m$iso3, m$year)
m$log_gdp <- log(m$gdp_pc)
m$log_hexp_pc <- log(m$hexp_pc)
m$log_hiv_prev <- log(m$hiv_prev)

# single-impute log child HIV prevalence for countries with no UNAIDS 0-14
# series, exactly as in the DHS analysis, with a flag
m$hiv_imputed <- !is.finite(m$log_hiv_prev)
m$log_hiv_prev[m$hiv_imputed] <- median(m$log_hiv_prev[!m$hiv_imputed], na.rm = TRUE)

m$year_c <- m$year - 2012
keep <- is.finite(m$share_pn) & m$share_pn > 0 & is.finite(m$pfpr_pct) & m$pfpr_pct > 0.5 &
  is.finite(m$hib3_wuenic) & is.finite(m$pcv3_wuenic) & is.finite(m$rotac_wuenic) &
  is.finite(m$log_gdp) & is.finite(m$log_hexp_pc)
m <- m[keep, , drop = FALSE]
m$iso3 <- factor(m$iso3)
cat(sprintf(
  "IHME-share panel: %d country-years, %d countries, %d-%d (HIV imputed for %d).\n",
  nrow(m), nlevels(m$iso3), min(m$year), max(m$year), sum(m$hiv_imputed)))

## ---- model: share ~ s(PfPR) + DHS national covariates + year + country RE ---
fit <- mgcv::gam(
  share_pn ~ s(pfpr_pct) + hib3_wuenic + pcv3_wuenic + rotac_wuenic +
    log_hiv_prev + log_gdp + log_hexp_pc + year_c + s(iso3, bs = "re"),
  data = m, method = "REML")
sm <- summary(fit)
cat(sprintf("s(PfPR2-10): edf=%.2f, p=%.3g; adj R2=%.2f; deviance explained=%.0f%%\n",
            sm$s.table["s(pfpr_pct)", "edf"], sm$s.table["s(pfpr_pct)", "p-value"],
            sm$r.sq, 100 * sm$dev.expl))
print(round(sm$p.table, 4))

## ---- variance explained: nested models (RE absorbs stable country baselines)-
dv <- function(f) 100 * summary(f)$dev.expl
mp  <- mgcv::gam(share_pn ~ s(pfpr_pct), data = m, method = "REML")
mpc <- mgcv::gam(share_pn ~ s(pfpr_pct) + hib3_wuenic + pcv3_wuenic + rotac_wuenic +
                   log_hiv_prev + log_gdp + log_hexp_pc, data = m, method = "REML")
cat(sprintf("Deviance explained: PfPR only=%.0f%%; +national covariates=%.0f%%; full (+year+country RE)=%.0f%%\n",
            dv(mp), dv(mpc), dv(fit)))

## ---- population share-vs-prevalence curve (country RE excluded) --------------
g <- data.frame(
  pfpr_pct = exp(seq(log(min(m$pfpr_pct)), log(max(m$pfpr_pct)), length.out = 200)),
  hib3_wuenic = mean(m$hib3_wuenic), pcv3_wuenic = mean(m$pcv3_wuenic),
  rotac_wuenic = mean(m$rotac_wuenic), log_hiv_prev = mean(m$log_hiv_prev),
  log_gdp = mean(m$log_gdp), log_hexp_pc = mean(m$log_hexp_pc),
  year_c = 0, iso3 = m$iso3[1])
pr <- predict(fit, g, se.fit = TRUE, exclude = "s(iso3)")
g$fit <- pr$fit; g$lo <- pr$fit - 1.96 * pr$se.fit; g$hi <- pr$fit + 1.96 * pr$se.fit
at <- function(x) { nd <- g[1, ]; nd$pfpr_pct <- x
  as.numeric(predict(fit, nd, exclude = "s(iso3)")) }
cat(sprintf("Fitted GBD share at PfPR 10/20/30/40%% = %.0f/%.0f/%.0f/%.0f%%\n",
            at(10), at(20), at(30), at(40)))

BLU <- "#2c7fb8"
p <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(data = g, ggplot2::aes(pfpr_pct, ymin = lo, ymax = hi), fill = BLU, alpha = 0.18) +
  ggplot2::geom_line(data = g, ggplot2::aes(pfpr_pct, fit), colour = BLU, linewidth = 1.2) +
  ggplot2::geom_point(data = m, ggplot2::aes(pfpr_pct, share_pn), colour = "grey35", alpha = 0.45, size = 1.6) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
  ggplot2::scale_x_log10(breaks = c(1, 3, 10, 30)) +
  ggplot2::labs(x = expression("National " * italic(Pf) * "PR"[2-10] * " (%), log scale"),
                y = "Malaria's share of post-neonatal\nchild deaths, IHME/GBD (%)") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 axis.title = ggplot2::element_text(size = 14), axis.text = ggplot2::element_text(size = 13))
ggplot2::ggsave(file.path(RESULTS_DIR, "figure7_ihme_share.png"), p, width = 8, height = 5.6, dpi = 320)

ct <- as.data.frame(sm$p.table); ct$term <- rownames(ct)
write.csv(ct[, c("term", "Estimate", "Std. Error", "Pr(>|t|)")],
          file.path(RESULTS_DIR, "ihme_share_model_coefficients.csv"), row.names = FALSE)
write.csv(
  m[, c("iso3", "year", "pfpr_pct", "share_pn", "malaria_deaths", "allcause_pn",
        "hib3_wuenic", "pcv3_wuenic", "rotac_wuenic", "log_hiv_prev", "log_gdp",
        "log_hexp_pc", "hiv_imputed")],
  file.path(RESULTS_DIR, "ihme_share_country_year_data.csv"), row.names = FALSE)
cat("saved: figure7_ihme_share.png + ihme_share_model_coefficients.csv + data\n")

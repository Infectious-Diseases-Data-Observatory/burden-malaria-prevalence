# =============================================================================
# 15_ihme_share_time.R — Method 1 (IHME deconvolution) as a COUNTRY-YEAR PANEL
# Malaria's share of all-cause POST-NEONATAL under-5 deaths from IHME/GBD 2023
# (both numerator and denominator from GBD), 2013-2023, merged with year-specific
# MAP population-weighted PfPR2-10. Spline model adjusted for log GDP p.c. and
# DTP3, with a country random effect; a prevalence x calendar-year interaction
# tests whether the prevalence-share relationship has changed over time.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})

## ---- IHME GBD deaths -> post-neonatal malaria share by country-year ---------
d  <- read.csv(file.path(DATA, "U5MR_IHME_GBD_2023.csv"), stringsAsFactors = FALSE)
d  <- d[d$measure == "Deaths" & d$metric == "Number", ]
grab <- function(cause, age, nm) { s <- d[d$cause == cause & d$age == age, c("location", "year", "val")]
  names(s)[3] <- nm; s }
mm <- Reduce(function(a, b) merge(a, b, by = c("location", "year")),
             list(grab("All causes", "<5 years", "ac5"), grab("All causes", "<28 days", "acn"),
                  grab("Malaria",    "<5 years", "ma5"), grab("Malaria",    "<28 days", "man")))
mm$allcause_pn <- mm$ac5 - mm$acn                       # post-neonatal all-cause deaths
mm$malaria_pn  <- mm$ma5 - mm$man                       # post-neonatal malaria deaths
mm$share_pn    <- 100 * mm$malaria_pn / mm$allcause_pn  # % of post-neonatal deaths due to malaria
mm$iso3 <- countrycode::countrycode(mm$location, "country.name", "iso3c", warn = FALSE)

## ---- merge year-specific PfPR2-10 (MAP), GDP p.c., DTP3 ----------------------
pf  <- read.csv(file.path(DATA, "pfpr_by_country_year.csv"))     # iso3, year, pfpr_pct
gdp <- read.csv(file.path(DATA, "wb_gdp_pc.csv"))                 # iso3, year, gdp_pc
dtp <- read.csv(file.path(DATA, "wb_dtp3.csv"))                   # iso3, year, dtp3
m <- merge(mm, pf, by = c("iso3", "year"))
m <- merge(m, gdp, by = c("iso3", "year")); m <- merge(m, dtp, by = c("iso3", "year"))
m$log_gdp <- log(m$gdp_pc); m$year_c <- m$year - 2018            # centre calendar year
m <- m[is.finite(m$share_pn) & is.finite(m$pfpr_pct) & is.finite(m$log_gdp) &
         is.finite(m$dtp3) & m$pfpr_pct > 0.5, ]
m$iso3 <- factor(m$iso3)
cat(sprintf("panel: %d country-years, %d countries, %d-%d\n",
            nrow(m), nlevels(m$iso3), min(m$year), max(m$year)))

## ---- models: full (spline) + prevalence x year interaction ------------------
m_main <- gam(share_pn ~ s(pfpr_pct) + log_gdp + dtp3 + year_c + s(iso3, bs = "re"),
              data = m, method = "REML")
m_int  <- gam(share_pn ~ s(pfpr_pct) + ti(pfpr_pct, year_c, k = c(6, 5)) +
                log_gdp + dtp3 + year_c + s(iso3, bs = "re"),
              data = m, method = "REML")
m_lin  <- gam(share_pn ~ pfpr_pct * year_c + log_gdp + dtp3 + s(iso3, bs = "re"),
              data = m, method = "REML")
sti <- summary(m_int)$s.table["ti(pfpr_pct,year_c)", ]
cat(sprintf("INTERACTION ti(PfPR, year): edf=%.2f, p=%.3f | dAIC(main - int)=%+.1f\n",
            sti["edf"], sti["p-value"], AIC(m_main) - AIC(m_int)))
pl <- summary(m_lin)$p.table["pfpr_pct:year_c", ]
cat(sprintf("Linear PfPR x year coef=%.4f pp/(pt*yr) (p=%.3f): change in the per-point prevalence slope per calendar year\n",
            pl["Estimate"], pl["Pr(>|t|)"]))
cat(sprintf("  prevalence slope at 2018 = %.3f pp/PfPR-pt; year main effect = %.3f pp/yr\n",
            summary(m_lin)$p.table["pfpr_pct", "Estimate"], summary(m_lin)$p.table["year_c", "Estimate"]))

## ---- plot the LINEAR-interaction model (m_lin): straight lines per year ------
yrsL <- c(2013, 2018, 2023)
gl <- do.call(rbind, lapply(yrsL, function(y) data.frame(
  pfpr_pct = seq(min(m$pfpr_pct), max(m$pfpr_pct), length.out = 100),
  year = y, year_c = y - 2018, log_gdp = mean(m$log_gdp), dtp3 = mean(m$dtp3), iso3 = m$iso3[1])))
gl$fit <- as.numeric(predict(m_lin, gl, exclude = "s(iso3)")); gl$Year <- factor(gl$year)
pL <- ggplot() +
  geom_point(data = m, aes(pfpr_pct, share_pn), colour = "grey70", alpha = 0.45, size = 1.6) +
  geom_line(data = gl, aes(pfpr_pct, fit, colour = Year), linewidth = 1.3) +
  geom_hline(yintercept = 0, colour = "grey85") +
  scale_colour_manual(values = c("2013" = "#2166ac", "2018" = "#7a7a7a", "2023" = "#b2182b"), name = "Year") +
  labs(x = expression("National age-standardised "*italic(Pf)*"PR"[2-10]*" (%)"),
       y = "Malaria's share of post-neonatal\nchild deaths, IHME/GBD (%)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
        axis.title = element_text(size = 14), axis.text = element_text(size = 13))
ggsave(file.path(RESULTS, "fig_ihme_share_mlin.png"), pL, width = 8, height = 5.6, dpi = 320)
cat("saved: results/fig_ihme_share_mlin.png (linear PfPR x year model)\n")

## ---- try te(): a single full 2-D tensor smooth over (log PfPR, year) --------
## te() carries BOTH main effects and the interaction, so no separate s()/year_c.
m$log_pfpr <- log(m$pfpr_pct)
m_te <- gam(share_pn ~ te(log_pfpr, year_c, k = c(6, 5)) + log_gdp + dtp3 + s(iso3, bs = "re"),
            data = m, method = "REML")
ste <- summary(m_te)$s.table
te_row <- ste[grep("^te\\(", rownames(ste)), ]
cat(sprintf("te(log PfPR, year): edf=%.2f, p=%.3f; dev.expl=%.0f%%; AIC=%.1f (m_int AIC=%.1f)\n",
            te_row["edf"], te_row["p-value"], 100 * summary(m_te)$dev.expl, AIC(m_te), AIC(m_int)))
xhi2 <- as.numeric(quantile(m$pfpr_pct, 0.95))         # draw over the well-supported range
gte <- do.call(rbind, lapply(c(2013, 2018, 2023), function(y) data.frame(
  pfpr_pct = exp(seq(log(min(m$pfpr_pct)), log(xhi2), length.out = 150)),
  year = y, year_c = y - 2018, log_gdp = mean(m$log_gdp), dtp3 = mean(m$dtp3), iso3 = m$iso3[1])))
gte$log_pfpr <- log(gte$pfpr_pct)
gte$fit <- as.numeric(predict(m_te, gte, exclude = "s(iso3)")); gte$Year <- factor(gte$year)
pte <- ggplot() +
  geom_point(data = m, aes(pfpr_pct, share_pn), colour = "grey70", alpha = 0.45, size = 1.6) +
  geom_line(data = gte, aes(pfpr_pct, fit, colour = Year), linewidth = 1.3) +
  geom_hline(yintercept = 0, colour = "grey85") +
  scale_x_log10(breaks = c(1, 3, 10, 30)) +
  scale_colour_manual(values = c("2013" = "#2166ac", "2018" = "#7a7a7a", "2023" = "#b2182b"), name = "Year") +
  labs(x = expression("National age-standardised "*italic(Pf)*"PR"[2-10]*" (%), log scale"),
       y = "Malaria's share of post-neonatal\nchild deaths, IHME/GBD (%)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
        axis.title = element_text(size = 14), axis.text = element_text(size = 13))
ggsave(file.path(RESULTS, "fig_ihme_share_te.png"), pte, width = 8, height = 5.6, dpi = 320)
cat("saved: results/fig_ihme_share_te.png (te full tensor, log PfPR x year)\n")

## ---- FIGURE: fit to the most recent year (2023) only ------------------------
## Single-year cross-section: no calendar-year term, no country random effect.
## Spline on LOG PfPR (matches the log axis) so the fit tracks the data across
## the prevalence range, including at low prevalence.
d23 <- m[m$year == 2023, ]
m23 <- gam(share_pn ~ s(log(pfpr_pct), k = 5) + log_gdp + dtp3, data = d23, method = "REML")
cat(sprintf("2023-only fit: n=%d countries; s(log PfPR) edf=%.2f; deviance explained=%.0f%%\n",
            nrow(d23), summary(m23)$s.table[1, "edf"], 100 * summary(m23)$dev.expl))
g <- data.frame(pfpr_pct = exp(seq(log(min(d23$pfpr_pct)), log(max(d23$pfpr_pct)), length.out = 200)),
                log_gdp = mean(d23$log_gdp), dtp3 = mean(d23$dtp3))
pr <- predict(m23, g, se.fit = TRUE)
g$fit <- pr$fit; g$lo <- pr$fit - 1.96 * pr$se.fit; g$hi <- pr$fit + 1.96 * pr$se.fit

p <- ggplot() +
  geom_ribbon(data = g, aes(pfpr_pct, ymin = lo, ymax = hi), fill = "grey75", alpha = 0.4) +
  geom_line(data = g, aes(pfpr_pct, fit), colour = "black", linewidth = 1.1) +
  geom_point(data = d23, aes(pfpr_pct, share_pn, colour = location), size = 2.8) +
  geom_hline(yintercept = 0, colour = "grey85") +
  scale_x_log10(breaks = c(1, 3, 10, 30)) +
  guides(colour = guide_legend(title = "Country", ncol = 2, override.aes = list(size = 2.8))) +
  labs(x = expression("National age-standardised "*italic(Pf)*"PR"[2-10]*" (%), log scale"),
       y = "Malaria's share of post-neonatal\nchild deaths, IHME/GBD 2023 (%)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.text = element_text(size = 7), legend.key.size = unit(3.6, "mm"),
        legend.title = element_text(size = 9),
        axis.title = element_text(size = 14), axis.text = element_text(size = 12))

OUT <- file.path(RESULTS, "fig_ihme_share_2023.png")
ggsave(OUT, p, width = 11, height = 6, dpi = 320)
cat(sprintf("saved: %s\n", OUT))
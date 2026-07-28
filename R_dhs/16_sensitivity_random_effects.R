# =============================================================================
# 16_sensitivity_random_effects.R — burden sensitivity to country random effects.
#
# The national burden extrapolation (script 10 / figure 4) uses the population-
# average model with country random effects set to zero. Here we quantify how
# the sub-Saharan Africa total changes if the fitted country-specific random
# effects are INCLUDED instead.
#
# In the attributable-fraction contrast (national PfPR vs a 1% counterfactual,
# same country and year) the random INTERCEPT cancels, so only the country
# random PfPR slope matters. Random effects exist only for the in-sample
# countries; all other SSA countries stay at the population-average value.
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2", "countrycode"))

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
model <- bundle$primary_fits$postneonatal$model
year_center <- unique(bundle$year_center)[1]
ref <- AF_REFERENCE

nts <- merge(
  read.csv(file.path(DATA_DIR, "wb_mortality_timeseries.csv"),
           stringsAsFactors = FALSE)[, c("iso3", "year", "allcause_1mo5y", "births")],
  read.csv(file.path(DATA_DIR, "pfpr_by_country_year.csv"), stringsAsFactors = FALSE),
  by = c("iso3", "year"))
nts$region <- countrycode::countrycode(nts$iso3, "iso3c", "region", warn = FALSE)
nts <- nts[nts$region == "Sub-Saharan Africa" & is.finite(nts$allcause_1mo5y) &
             is.finite(nts$pfpr_pct) & is.finite(nts$births) &
             nts$year >= 2000 & nts$year <= 2024, , drop = FALSE]

lvl <- levels(model$model$country)
insample <- nts$iso3 %in% lvl
cat(sprintf("Random-effects burden sensitivity: %d country-years; REs available for %d of %d SSA countries.\n",
            nrow(nts), length(unique(nts$iso3[insample])), length(unique(nts$iso3))))

## ---- attributable fraction with / without country random effects -----------
newdata <- function(pf10, iso) {
  data.frame(pfpr10 = pf10, year_c = nts$year - year_center, exposure = 1,
             country = factor(iso, levels = lvl),
             G = I(matrix(0, nrow(nts), ncol(model$model$G),
                          dimnames = list(NULL, colnames(model$model$G)))))
}
re_cols <- grep("^s\\(country\\)|^s\\(country,pfpr10\\)",
                colnames(predict(model, newdata(nts$pfpr_pct / 10, lvl[1]), type = "lpmatrix")))

af_contrast <- function(iso, zero_re) {
  Xh <- predict(model, newdata(nts$pfpr_pct / 10, iso), type = "lpmatrix")
  Xl <- predict(model, newdata(rep(ref / 10, nrow(nts)), iso), type = "lpmatrix")
  if (length(zero_re)) { Xh[zero_re, re_cols] <- 0; Xl[zero_re, re_cols] <- 0 }
  pmax(1 - exp(-as.numeric((Xh - Xl) %*% coef(model))), 0)
}
# population-average: zero REs for every row
af_pop <- af_contrast(rep(lvl[1], nrow(nts)), zero_re = seq_len(nrow(nts)))
# REs included: in-sample countries use their own levels; others stay population
af_re <- af_contrast(ifelse(insample, nts$iso3, lvl[1]), zero_re = which(!insample))

deaths_pop <- af_pop * nts$allcause_1mo5y
deaths_re  <- af_re  * nts$allcause_1mo5y

## ---- per-year totals --------------------------------------------------------
pop_by_year <- tapply(deaths_pop, nts$year, sum)   # named by year, ascending
re_by_year  <- tapply(deaths_re,  nts$year, sum)
summary_tab <- data.frame(
  year = as.integer(names(pop_by_year)),
  deaths_re0_population = as.numeric(pop_by_year),
  deaths_re_included    = as.numeric(re_by_year))
summary_tab$diff_pct <- 100 * (summary_tab$deaths_re_included / summary_tab$deaths_re0_population - 1)
write.csv(summary_tab, file.path(RESULTS_DIR, "sensitivity_random_effects_summary.csv"), row.names = FALSE)

decl <- function(v) 100 * (1 - v[summary_tab$year == 2024] / v[summary_tab$year == 2000])
cat(sprintf("Population (RE=0):  2000=%.0fk  2024=%.0fk  decline=%.0f%%\n",
            summary_tab$deaths_re0_population[summary_tab$year == 2000] / 1000,
            summary_tab$deaths_re0_population[summary_tab$year == 2024] / 1000,
            decl(summary_tab$deaths_re0_population)))
cat(sprintf("REs included:       2000=%.0fk  2024=%.0fk  decline=%.0f%%\n",
            summary_tab$deaths_re_included[summary_tab$year == 2000] / 1000,
            summary_tab$deaths_re_included[summary_tab$year == 2024] / 1000,
            decl(summary_tab$deaths_re_included)))
cat(sprintf("Max absolute yearly difference: %.1f%%; cumulative 2000-2024 %.2fM vs %.2fM (%+.1f%%)\n",
            max(abs(summary_tab$diff_pct)), sum(deaths_pop) / 1e6, sum(deaths_re) / 1e6,
            100 * (sum(deaths_re) / sum(deaths_pop) - 1)))

## ---- 2024 country shifts (in-sample only) -----------------------------------
last <- nts$year == 2024 & insample
country_shift <- data.frame(
  iso3 = nts$iso3[last],
  deaths_re0_population = deaths_pop[last],
  deaths_re_included = deaths_re[last])
country_shift$difference <- country_shift$deaths_re_included - country_shift$deaths_re0_population
country_shift <- country_shift[order(-abs(country_shift$difference)), ]
write.csv(country_shift, file.path(RESULTS_DIR, "sensitivity_random_effects_country_2024.csv"), row.names = FALSE)
cat("Largest 2024 country shifts from including REs (thousands):\n")
print(head(transform(country_shift,
                     pop_k = round(deaths_re0_population / 1000, 1),
                     re_k = round(deaths_re_included / 1000, 1),
                     delta_k = round(difference / 1000, 1))[, c("iso3", "pop_k", "re_k", "delta_k")], 6),
      row.names = FALSE)

## ---- overlay figure ---------------------------------------------------------
plot_df <- rbind(
  data.frame(year = summary_tab$year, deaths = summary_tab$deaths_re0_population,
             series = "Country REs set to 0 (reported)"),
  data.frame(year = summary_tab$year, deaths = summary_tab$deaths_re_included,
             series = "Country REs included"))
p <- ggplot2::ggplot(plot_df, ggplot2::aes(year, deaths / 1000, colour = series)) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::scale_colour_manual(values = c("Country REs set to 0 (reported)" = "#08519c",
                                          "Country REs included" = "#d95f0e"), name = NULL) +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::labs(x = NULL, y = "Malaria-attributable under-5 post-neonatal\ndeaths in sub-Saharan Africa (thousands)",
                title = "Burden sensitivity to country random effects",
                subtitle = "Population-average vs country-specific random PfPR slopes (differs by under 3% in any year)") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = c(0.99, 0.98), legend.justification = c(1, 1),
                 legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.75), colour = NA))
ggplot2::ggsave(file.path(RESULTS_DIR, "sensitivity_random_effects_curves.png"),
                p, width = 8, height = 5, dpi = 320, bg = "white")
cat("saved: sensitivity_random_effects_summary.csv + _country_2024.csv + _curves.png\n")

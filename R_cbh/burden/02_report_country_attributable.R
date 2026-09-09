#!/usr/bin/env Rscript
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/burden/settings.R")
year <- cbh_burden_year()
year_text <- function(x) gsub("{year}", as.character(year), x, fixed = TRUE)
library(ggplot2)
library(patchwork)
spec <- cbh_trial_spec()
out <- file.path("results/cbh", spec$id, year_text("country_burden_{year}"))
r <- cbh_read_csv(file.path(out, year_text("country_age_attributable_{year}.csv")))
t <- cbh_read_csv(file.path(out, year_text("country_totals_{year}.csv")))
src <- cbh_read_csv(file.path(out, "ihme_disjoint_age_inputs.csv"))
ages <- cbh_config()$age_bands$age_band
z <- r[r$iso3 == "COD", ]; z <- z[match(ages, z$age_band), ]
stopifnot(nrow(z) == 7, all(is.finite(z$attributable_rate_per100000)))
d <- src[src$iso3 == "COD", ]
value <- function(a) d$rate_per100000[match(a, d$source_age)] / 1e5
# A period life table under piecewise constant hazards. First boundary is
# IHME's 28 days, next is 6 months. This figure explicitly marks the small
# neonatal mapping difference from the DHS model's completed-month boundary.
edges <- c(0, 28 / 365.25 * 12, 6, 12, 24, 36, 48, 60)
H <- z$ihme_rate_per100000 / 1e5 * diff(edges) / 12
H[1] <- value("0 to 6 days (early neonatal)") * 7/365.25 +
        value("7 to 27 days (late neonatal)") * 21/365.25
life <- data.frame(age_band = ages, start_month = head(edges, -1), end_month = tail(edges, -1),
  cumulative_hazard = H, death_probability = -expm1(-H),
  counterfactual_death_probability = -expm1(-H * z$hr_zero_vs_current),
  survival_end = exp(-cumsum(H)), counterfactual_survival_end = exp(-cumsum(H * z$hr_zero_vs_current)))
stopifnot(max(abs(cumprod(1-life$death_probability)-life$survival_end)) < 1e-12,
          max(abs(cumprod(1-life$counterfactual_death_probability)-life$counterfactual_survival_end)) < 1e-12)
cbh_atomic_csv(life, file.path(out, year_text("drc_period_life_table_{year}.csv")))
labels <- c("IHME baseline", "Zero-PfPR counterfactual")
cols <- setNames(c("#263C52", "#CA6130"), labels)
surv <- rbind(data.frame(age = edges, survival = c(1, life$survival_end), scenario = labels[1]),
              data.frame(age = edges, survival = c(1, life$counterfactual_survival_end), scenario = labels[2]))
prob <- rbind(data.frame(life[c("start_month", "end_month")], probability = life$death_probability, scenario = labels[1]),
              data.frame(life[c("start_month", "end_month")], probability = life$counterfactual_death_probability, scenario = labels[2]))
theme_set(theme_minimal(base_size = 12))
sty <- theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
             plot.caption = element_text(hjust = 0), legend.position = "bottom")
a <- ggplot(surv, aes(age, survival, colour = scenario)) + geom_step(linewidth = .8) +
  geom_point(size = 1.4) + scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(0,12,24,36,48,60)) +
  labs(title = "A  Survival to the end of each age band", x = "Age (months)", y = "Survival probability", colour = NULL) + sty
a <- a + guides(colour = "none")
b <- ggplot(prob, aes(x = start_month, xend = end_month, y = probability * 1000,
                       yend = probability * 1000, colour = scenario)) +
  geom_segment(linewidth = 1.2, lineend = "butt") + scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(0,12,24,36,48,60)) +
  labs(title = "B  Conditional death probability", x = "Age (months)",
       y = "Deaths per 1,000 children entering the band", colour = NULL) + sty
z$label <- factor(c("0-27 days", "1-5 months", "6-11 months", "12-23 months",
                    "24-35 months", "36-47 months", "48-59 months"),
                  levels = c("0-27 days", "1-5 months", "6-11 months", "12-23 months", "24-35 months", "36-47 months", "48-59 months"))
c <- ggplot(z, aes(label, attributable_rate_per100000)) +
  geom_hline(yintercept = 0, colour = "grey50") + geom_col(fill = "#147C80", width = .65) +
  geom_errorbar(aes(ymin = attributable_rate_per100000_lower_95, ymax = attributable_rate_per100000_upper_95), width = .15) +
  labs(title = "C  Malaria-attributable mortality rate", x = NULL,
       y = "Deaths per 100,000 person-years") + sty
figure <- ((a | b) / c) + plot_layout(guides = "collect", heights = c(1.15, 1)) +
  plot_annotation(title = year_text("DR Congo: malaria contribution to under-five mortality, {year}"),
    subtitle = sprintf("National population-weighted PfPR2-10: %.2f%% -> 0%% | IHME all-cause rates | Shared-calendar-year mortality model", z$pfpr_pct[1]),
    caption = paste(year_text("A-B: synthetic cohort under {year} hazards; point estimates. C: signed rate difference; intervals cover model and HIV-imputation uncertainty only."),
      "Ages 2-4 share one IHME baseline rate; deaths allocated equally. Neonatal: 0-27 days mapped to the <1-month model effect.",
      year_text("MAP {year} weighted using GPW 2020 population geography. IHME/MAP uncertainty and survey-design uncertainty are not included."), sep = "\n"),
    theme = theme(plot.title = element_text(face = "bold", size = 17), plot.caption = element_text(hjust = 0, size = 10)))
figure <- figure & theme(legend.position = "bottom")
ggsave(file.path(out, year_text("drc_malaria_contribution_{year}.png")), figure, device = ragg::agg_png,
       width = 14, height = 10, dpi = 160, bg = "white")
dr <- t[t$iso3 == "COD", ]
table_rows <- vapply(seq_len(nrow(z)), function(i) sprintf("| %s | %.1f | %.1f | %.1f (%.1f to %.1f) | %.0f |",
  as.character(z$label[i]), z$ihme_rate_per100000[i], z$counterfactual_rate_per100000[i],
  z$attributable_rate_per100000[i], z$attributable_rate_per100000_lower_95[i],
  z$attributable_rate_per100000_upper_95[i], z$attributable_deaths[i]), "")
writeLines(c(year_text("# National malaria-attributable mortality, {year}"), "",
  sprintf("Estimates are available for **%d of %d countries** in the new IHME export. Countries without usable MAP coverage remain missing: %s.",
    sum(t$status == "estimated"), nrow(t), paste(t$country[t$status != "estimated"], collapse = ", ")), "",
  sprintf("**Exposure coverage:** %d estimated countries have MAP covering less than 95%% of population weight within the raster footprint. Their estimates apply the covered-area mean to the whole country and should be treated as provisional. Swaziland/Eswatini has only %.1f%% coverage; its value is especially poorly representative. Missing pixels are not assumed malaria-free. DR Congo has %.1f%% coverage. The 95%% threshold is a reporting flag, not an exclusion rule.",
    sum(t$map_coverage_below_95pct & t$status == "estimated", na.rm = TRUE),
    100*t$map_population_coverage_within_raster[t$iso3 == "SWZ"],
    100*t$map_population_coverage_within_raster[t$iso3 == "COD"]), "",
  year_text("Each country's {year} national PfPR2-10 is evaluated on each fitted age-specific spline. The zero-PfPR hazard ratio is exp[f_g(0)-f_g(P_country)]. Counterfactual mortality rate = IHME rate x HR; attributable rate = IHME rate x (1-HR). Death counts use the same fraction and fixed annual person-time. No adjustment is made to other covariates, calendar year, or country random effects in this contrast: those terms cancel in the current additive model."), "",
  "The seven PfPR effects come from the latest shared-calendar-year mortality model and its ten HIV-incidence uncertainty refits. Contrasts are pooled on the log-HR scale using within-fit covariance plus between-imputation variance and finite-imputation t intervals. Reported point HRs exponentiate pooled mean log HRs. Bounds condition on IHME and MAP point estimates, fitted smoothing parameters, and the age allocation. They do not include source-estimate, survey-design or residual-clustering uncertainty. National total death counts are point estimates; marginal age-band interval endpoints are not summed into total intervals.", "",
  "## Source handling and assumptions", "",
  year_text("- The selected IHME export is the file dated 2026-09-09 10-58-22. Only All causes / Deaths / Both sexes / {year} is used. Six disjoint source age groups are checked against Under 1 and Under 5 totals; those aggregates are never counted twice. Original rate/count lower and upper bounds are retained in the source output."),
  "- Early and late neonatal deaths are added. Their combined annual rate is total deaths divided by the sum of their implied person-years, not the sum or simple mean of rates. The IHME 0-27-day group uses the model's <1-completed-month PfPR effect, an explicit boundary approximation.",
  "- As authorized, ages 24-35, 36-47 and 48-59 months each use the IHME 2-4-year mortality rate and one-third of its deaths/person-time. Each then receives its own fitted PfPR contrast. These three baseline age estimates are assumed, not separately observed in IHME.",
  year_text("- National exposure is calculated from local {year} MAP rasters and country polygons, using GPW 2020 density x grid-cell area x polygon overlap as weights. The earlier extraction used density alone. Both values are saved for comparison. This uses fixed all-age 2020 population geography, not a {year} age-specific population surface. Coverage is conditional on the available raster footprint; missing MAP is never set to zero."),
  "- The national-mean PfPR scenario is the requested approximation; the nonlinear response evaluated at a country mean does not equal a subnational burden aggregation. The fitted PfPR curves are transported to countries outside the mortality fitting sample.",
  "- Negative attributable effects are retained. They mean the fitted association predicts higher mortality at zero PfPR, not protective malaria established by evidence. Zero-PfPR support flags and central-range flags are exported; the model's extrapolation, causal transport and smoothing-stability limitations remain relevant.", "",
  "## DR Congo", "",
  sprintf("National PfPR: **%.2f%%**. IHME under-five deaths: **%s**. Signed model-attributable deaths: **%s** (**%.1f%%**), holding the annual population exposure fixed.",
    dr$pfpr_pct, format(round(dr$ihme_under5_deaths), big.mark = ","),
    format(round(dr$attributable_under5_deaths), big.mark = ","), 100 * dr$attributable_fraction), "",
  "All rates below are per 100,000 person-years. Age-band intervals reflect model/HIV-imputation uncertainty only.", "",
  "| Age | IHME rate | Zero-PfPR rate | Attributable rate (95% interval) | Attributable annual deaths |",
  "|---|---:|---:|---:|---:|", table_rows, "",
  year_text("![DRC malaria contribution](drc_malaria_contribution_{year}.png)"), "",
  year_text("The survival and probability panels describe a synthetic cohort exposed to the {year} mortality schedule. They use q=1-exp(-H), with IHME neonatal hazards integrated separately over 7 and 21 days, followed by ages 28 days to 6 months, 6-12 months, and annual age intervals. This explicitly handles the source neonatal boundary. The survival product is checked numerically. Neither these probabilities nor their differences are used as annual death-count fractions."), "",
  "## Outputs", "",
  year_text("- [All country/age estimates](country_age_attributable_{year}.csv)"),
  year_text("- [Country totals](country_totals_{year}.csv)"),
  year_text("- [National prevalence and coverage](national_pfpr_{year}.csv)"),
  year_text("- [IHME original {year} rows](ihme_source_{year}.csv) and [disjoint age inputs](ihme_disjoint_age_inputs.csv)"),
  "- [Individual-fit contrasts](individual_log_hazard_contrasts.csv), [PfPR support](model_pfpr_support.csv)",
  year_text("- [DRC period life table](drc_period_life_table_{year}.csv), [provenance](provenance.csv)")), file.path(out, "REPORT.md"))
message("Wrote DRC figure and country burden report: ", out)

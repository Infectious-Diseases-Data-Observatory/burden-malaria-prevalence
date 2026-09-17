#!/usr/bin/env Rscript
# Figure 5: matched countries and one population denominator across sources.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
library(ggplot2)
root <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional"))$out
out <- file.path(root, "annual_comparison")
x <- cbh_read_csv(file.path(out, "annual_totals_2004_2024.csv"))
stopifnot(identical(x$year, 2004:2024), all(x$countries == 42))
model_label <- cbh_paper_model_label()
sources <- c(model = model_label, ihme_malaria = "IHME", who_cacode = "UN IGME")
long <- do.call(rbind, lapply(names(sources), function(s) data.frame(
  year = x$year, source = sources[[s]], deaths = x[[paste0(s, "_deaths")]],
  rate_per100000 = x[[paste0(s, "_rate_per100000")]],
  under5_person_years = x$under5_person_years, countries = x$countries)))
stopifnot(all(is.finite(long$rate_per100000)), all(long$rate_per100000 > 0),
  max(abs(long$rate_per100000 - long$deaths / long$under5_person_years * 1e5)) < 1e-8)
long$source <- factor(long$source, levels = unname(sources))
cbh_atomic_csv(long, file.path(out, "figure5_data.csv"))
p <- ggplot(long, aes(year, rate_per100000, colour = source, linetype = source)) +
  geom_line(linewidth = 1.15) +
  scale_colour_manual(values = setNames(c("#16747C", "#253746", "#CC6A30"), unname(sources))) +
  scale_linetype_manual(values = setNames(c("solid", "longdash", "dotdash"), unname(sources))) +
  scale_x_continuous(breaks = c(2004, 2009, 2014, 2019, 2024), limits = c(2004, 2024),
    expand = expansion(mult = c(.015, .025))) +
  scale_y_continuous(limits = c(0, NA), labels = scales::label_comma(),
    expand = expansion(mult = c(0, .05))) +
  labs(x = "Year", y = "Under-five malaria deaths\nper 100,000 child-years", colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 18) +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16),
    legend.text = element_text(size = 16), legend.position = "bottom",
    legend.key.width = grid::unit(1.2, "cm"), panel.grid.minor = element_blank(),
    plot.margin = margin(14, 20, 14, 14))
ggsave(file.path(out, "fig5_annual_malaria_mortality.png"), p,
  device = ragg::agg_png, width = 12, height = 7.5, dpi = 300, bg = "white")
caption <- paste(
  "Annual malaria mortality before age five, 2004–2024, pooled across the same 42 countries covered by the national burden analysis.",
  "For each source, the plotted rate is 100,000 times the sum of national under-five malaria deaths divided by the sum of annual under-five person-years implied by the primary IHME all-cause death counts and rates.",
  "Rates therefore use a common population denominator and are not averages of country rates or probabilities per live birth.",
  paste("The", model_label, "applies the seven separate primary MAP gamma=2 PfPR effects to annual age-specific IHME all-cause deaths, using the zero-PfPR counterfactual."),
  "National annual MAP PfPR[2–10] is weighted by population counts using a fixed GPW 2020 spatial distribution; these exposure weights are separate from the annual mortality denominator.",
  "IHME ages 2–4 share a mortality rate and divide deaths/person-time equally across the model's three annual bands. Signed age-specific contributions are retained.",
  "The UN IGME comparator is the CA-CODE 2026 series from the UN IGME portal, retrieved from UNICEF's CME_CAUSE_OF_DEATH dataflow on 17 September 2026 (under five, both sexes, malaria, deaths); it is not the World Malaria Report all-age series or a fixed proportion of it.",
  "All curves show point estimates. Joint uncertainty across countries and age bands is not available; marginal interval endpoints are not summed.",
  paste("Model-attributable all-cause mortality reductions and the comparator cause-specific estimates are different estimands. The", model_label, "shares IHME all-cause inputs and this comparison is not independent validation."),
  "Zero exposure requires extrapolation. Cape Verde, Lesotho and São Tomé and Príncipe are excluded because usable MAP prevalence is unavailable; the other countries absent from the original national input set are not added.")
writeLines(c("# Figure 5 caption", "", caption), file.path(out, "CAPTION.md"))
audit <- cbh_read_csv(file.path(out, "denominator_audit.csv"))
check <- cbh_read_csv(file.path(out, "agreement_with_existing_primary.csv"))
fmt <- function(z) format(round(z), big.mark = ",", trim = TRUE, scientific = FALSE)
rows <- vapply(seq_len(nrow(x)), function(i) sprintf("| %d | %s | %.1f | %s | %.1f | %s | %.1f |",
  x$year[i], fmt(x$model_deaths[i]), x$model_rate_per100000[i], fmt(x$ihme_malaria_deaths[i]),
  x$ihme_malaria_rate_per100000[i], fmt(x$who_cacode_deaths[i]), x$who_cacode_rate_per100000[i]), "")
writeLines(c("# Figure 5: annual under-five malaria mortality", "",
  "All required inputs are available for 2004–2024. The comparison contains the same 42 countries in all 21 years (882 country-years per source). This is the national burden coverage, not only the 34 countries contributing DHS surveys to model fitting.", "",
  "## Input audit", "",
  "- Primary model: seven verified saved MAP gamma=2 fits; no refitting required.",
  "- MAP: all 21 rasters readable; national means recalculated using population density × cell area, not the legacy density-only country CSV.",
  "- IHME all-cause: all required disjoint age counts and rates, 2004–2024, export dated 9 September 2026.",
  "- IHME malaria: direct under-five, both-sex national counts and rates, export dated 3 September 2026, filtered to 2004–2024.",
  "- UN IGME: direct under-five malaria deaths for every included country-year in the CA-CODE 2026 portal release. The previously cached WHO all-age African Region series multiplied by 0.75 is not used.", "",
  "## Calculation and checks", "",
  "For country c, year t and age band g: D_model(c,t,g) = D_allcause(c,t,g) × {1 − exp[f_g(0) − f_g(P_c,t)]}. Sum over ages and countries. For each source s: rate_s(t) = 100,000 × sum_c D_s(c,t) / sum_c PY_U5(c,t). Here PY_U5 = D_allcause,U5 / (rate_allcause,U5 / 100,000). Keep this same denominator for all sources.", "",
  sprintf("Disjoint age deaths reproduce the all-cause under-five totals (maximum relative error %.3g). The sum of age-implied person-years differs from the directly reported under-five implied denominator by at most %.4f%%; we use the direct under-five denominator for the pooled rate.",
    max(abs(audit$death_sum_relative_error)), 100 * max(abs(audit$person_year_sum_relative_error))), "",
  sprintf("IHME malaria and all-cause count/rate pairs imply denominators differing by at most %.4f%% where the malaria rate is positive. Native IHME malaria rates and the differences are retained in denominator_audit.csv; the Figure 5 IHME rates are recomputed using the common all-cause denominator. Zero malaria count/rate pairs cannot identify population and are not used to infer it.",
    100 * max(abs(audit$malaria_denominator_relative_difference), na.rm = TRUE)), "",
  sprintf("All 126 existing country estimates for 2005, 2015 and 2024 are reproduced; maximum absolute death-count difference %.3g.", max(abs(check$death_difference))), "",
  "The exports identify the data suite as GBD but do not state an unambiguous release or estimation/forecast version. Precise release alignment remains unverified; date and denominator differences are documented rather than assumed away. UN IGME uses a separate cause-of-death model and mortality envelope. National exposure is transported beyond the DHS fitting countries, and MAP population coverage and extrapolation flags are retained in the country/age CSVs. Fixed 2020 spatial weights do not describe changes in within-country population geography over time.", "",
  "Point estimates only: no cross-country or cross-age joint uncertainty has been supplied. Do not sum marginal bounds. Source uncertainty, HIV imputation uncertainty, exposure error and survey design are not propagated.", "",
  "## Annual totals", "", "Rates are deaths per 100,000 under-five child-years.", "",
  paste0("| Year | ", model_label, " deaths | Rate | IHME deaths | Rate | UN IGME deaths | Rate |"),
  "|---|---:|---:|---:|---:|---:|---:|", rows, "",
  "## Sources and reproducibility", "",
  "[WHO indicator metadata](https://www.who.int/data/gho/data/indicators/indicator-details/GHO/number-of-deaths), [UN IGME cause-of-death portal](https://childmortality.org/causes-of-death/data), [CA-CODE 2000–2024 study](https://doi.org/10.1136/bmj-2025-088686). Exact UNICEF API URL, retrieval time and SHA256 are in [source metadata](who_source.json).", "",
  "Run `Rscript R_cbh/burden/04_annual_comparison.R --audit-only` first, then `Rscript R_cbh/burden/04_annual_comparison.R` and `Rscript R_cbh/reporting/07_annual_mortality_comparison.R`. The scripts use local inputs and do not download data or refit models.", "",
  "[Included countries](included_countries.csv) · [Annual totals](annual_totals_2004_2024.csv) · [Country estimates](country_estimates_2004_2024.csv) · [Country-year input audit](country_year_input_audit.csv) · [Denominator audit](denominator_audit.csv) · [Input hashes](input_provenance.csv)", "",
  "![Figure 5](fig5_annual_malaria_mortality.png)", "", caption), file.path(out, "README.md"))
file.copy("data/external/figure5/source.json", file.path(out, "who_source.json"), overwrite = TRUE)
paths <- c(file.path(out, "annual_totals_2004_2024.csv"),
  "R_cbh/reporting/07_annual_mortality_comparison.R", "R_cbh/reporting/labels.R")
cbh_atomic_csv(data.frame(file = paths, md5 = vapply(paths, cbh_file_hash, "")),
  file.path(out, "figure_provenance.csv"))

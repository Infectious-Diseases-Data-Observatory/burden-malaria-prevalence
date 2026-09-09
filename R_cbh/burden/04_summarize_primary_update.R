#!/usr/bin/env Rscript
# Compare primary updates across the three requested years using aggregate files.
source("R_cbh/load_pipeline.R")
library(ggplot2)
years <- c(2005L, 2015L, 2024L)
primary <- "results/cbh/age_band_separate_v1"
historical <- "results/cbh/age_band_hiv_incidence_shared_time_v3"
out <- file.path(primary, "primary_update_2005_2015_2024")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
inputs <- character()
rows <- lapply(years, function(year) {
  new_dir <- file.path(primary, paste0("country_burden_", year))
  old_dir <- file.path(historical, paste0("country_burden_", year))
  new_path <- file.path(new_dir, paste0("model_vs_ihme_malaria_", year, ".csv"))
  old_path <- file.path(old_dir, paste0("country_totals_", year, ".csv"))
  n <- cbh_read_csv(new_path); o <- cbh_read_csv(old_path)
  cbh_unique(n, "iso3", "New country estimates"); cbh_unique(o, "iso3", "Previous country estimates")
  stopifnot(setequal(n$iso3, o$iso3), all(n$year == year), all(n$model_id == "age_band_separate_v1"))
  o <- o[match(n$iso3, o$iso3), ]
  # Hold the source rates/counts, exposure and missing-country policy fixed.
  for (v in c("pfpr_pct", "ihme_under5_deaths", "map_population_coverage_within_raster")) {
    stopifnot(identical(is.na(n[[v]]), is.na(o[[v]])),
      max(abs(n[[v]] - o[[v]]), na.rm = TRUE) < 1e-7)
  }
  stopifnot(identical(n$status, o$status))
  d <- n[c("year", "iso3", "country", "status", "comparison_status", "pfpr_pct",
    "map_population_coverage_within_raster", "map_coverage_below_95pct", "ihme_under5_deaths",
    "attributable_under5_deaths", "counterfactual_under5_deaths", "attributable_fraction",
    "ihme_malaria_deaths", "model_to_ihme_ratio", "model_id", "hiv_treatment")]
  d$previous_joint_attributable_deaths <- o$attributable_under5_deaths
  d$previous_joint_hiv_treatment <- "ten_posterior_draws"
  d$change_from_previous_joint_deaths <- d$attributable_under5_deaths - d$previous_joint_attributable_deaths
  d$percent_change_from_previous_joint <- 100 * d$change_from_previous_joint_deaths / d$previous_joint_attributable_deaths
  inputs <<- c(inputs, new_path, old_path)
  d
})
d <- do.call(rbind, rows)
cbh_unique(d, c("year", "iso3"), "Combined country-year estimates")
cbh_atomic_csv(d, file.path(out, "country_estimates_2005_2015_2024.csv"))
m <- d[d$comparison_status == "matched", ]
stopifnot(all(table(m$year) == 42), all(is.finite(m$attributable_under5_deaths)),
  all(m$attributable_under5_deaths > 0), all(m$ihme_malaria_deaths > 0))
s <- do.call(rbind, lapply(years, function(y) {
  x <- m[m$year == y, ]
  data.frame(year = y, countries = nrow(x),
    separate_age_attributable_deaths = sum(x$attributable_under5_deaths),
    previous_joint_attributable_deaths = sum(x$previous_joint_attributable_deaths),
    ihme_malaria_deaths = sum(x$ihme_malaria_deaths),
    model_to_ihme_ratio = sum(x$attributable_under5_deaths)/sum(x$ihme_malaria_deaths),
    percent_change_from_previous_joint = 100*(sum(x$attributable_under5_deaths)/sum(x$previous_joint_attributable_deaths)-1))
}))
cbh_atomic_csv(s, file.path(out, "year_summary.csv"))
m$coverage <- factor(ifelse(m$map_coverage_below_95pct, "MAP coverage <95%", "MAP coverage >=95%"),
  levels = c("MAP coverage >=95%", "MAP coverage <95%"))
span <- range(m$attributable_under5_deaths, m$ihme_malaria_deaths)
limits <- span * 10^c(-.3, .3)
ticks <- 10^seq(floor(log10(span[1])), ceiling(log10(span[2])))
labels <- function(x) vapply(x, function(v) {
  if (is.na(v)) return(NA_character_)
  if (v >= 1e6) paste0(format(v/1e6, trim = TRUE), "M") else
    if (v >= 1e3) paste0(format(v/1e3, trim = TRUE), "k") else format(v, scientific = FALSE, trim = TRUE)
}, "")
p <- ggplot(m, aes(ihme_malaria_deaths, attributable_under5_deaths)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = coverage), size = 2.2, alpha = .85) +
  geom_text(data = m[m$iso3 %in% c("NGA", "COD"), ], aes(label = iso3),
    nudge_y = .16, size = 3, colour = "#253746") +
  facet_wrap(~year, nrow = 1) + coord_fixed() +
  scale_x_log10(limits = limits, breaks = ticks, labels = labels) +
  scale_y_log10(limits = limits, breaks = ticks, labels = labels) +
  scale_colour_manual(values = c("MAP coverage >=95%" = "#16747C", "MAP coverage <95%" = "#CC6A30")) +
  labs(title = "Updated under-five malaria mortality by country",
    subtitle = "Separate age-band primary models | One fixed median HIV imputation | 42 matched countries per year",
    x = "IHME malaria deaths (log scale)", y = "Our malaria-attributable deaths (log scale)", colour = NULL,
    caption = paste("Dashed line: equal estimates. Our estimate is the all-cause mortality reduction under national PfPR2-10 -> 0%.",
      "Point estimates only. MAP uses fixed 2020 population weights; ages 2-4 retain equal baseline rates and death/person-time shares.",
      "Three countries remain missing because MAP is unavailable. Orange points have limited MAP population coverage.", sep = "\n")) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16), strip.text = element_text(face = "bold"),
    legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 9),
    plot.margin = margin(10, 20, 10, 10))
ggsave(file.path(out, "model_vs_ihme_all_years.png"), p, width = 15, height = 7,
  dpi = 180, device = ragg::agg_png, bg = "white")
fmt <- function(x) format(round(x), big.mark = ",", trim = TRUE)
report <- c("# Country mortality estimates after selecting separate age-band models", "",
  "Updated estimates for 2005, 2015 and 2024 use the seven existing separate age-band models and the fixed posterior-median child HIV-incidence imputation. The earlier joint-model results remain available in their original directory.", "",
  "| Year | Matched countries | Separate-age attributable deaths | Previous joint estimate | IHME malaria deaths | Separate / IHME |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(s)), function(i) sprintf("| %d | %d | %s | %s | %s | %.2f |", s$year[i], s$countries[i],
    fmt(s$separate_age_attributable_deaths[i]), fmt(s$previous_joint_attributable_deaths[i]), fmt(s$ihme_malaria_deaths[i]), s$model_to_ihme_ratio[i]), ""), "",
  "Both estimates use the same IHME all-cause inputs, national PfPR values, population weights and age allocation; these were checked numerically for every country and year. The previous joint estimates also pooled ten HIV-imputation fits, whereas the new estimates use one fixed median imputation. The change column therefore does not isolate model structure from HIV-imputation treatment.", "",
  "Our estimate is the signed reduction in all-cause deaths under zero national PfPR; IHME reports cause-specific malaria deaths. These are different estimands. Zero-PfPR extrapolation, transport to countries outside the fitting sample and incomplete MAP coverage remain limitations. Missing values are never replaced with zero, and negative age-band contributions are retained.", "",
  "Country totals are point estimates. Age-band intervals in the year-specific tables use within-fit covariance conditional on the fixed HIV imputation and fitted smoothing parameters. They exclude HIV-imputation, survey-design, residual-clustering, source and allocation uncertainty. Separate fits do not imply independent sampling errors across age bands; no country-total interval is constructed by summing marginal interval endpoints.", "",
  "![Updated country comparisons](model_vs_ihme_all_years.png)", "",
  "## Files", "",
  "- [All 135 country-year rows, including explicit missing estimates](country_estimates_2005_2015_2024.csv)",
  "- [Year totals](year_summary.csv)",
  vapply(years, function(y) sprintf("- **%d:** [country totals](../country_burden_%d/country_totals_%d.csv), [age-band rates and counts](../country_burden_%d/country_age_attributable_%d.csv), [all-country plot](../country_burden_%d/model_vs_ihme_malaria_by_country_%d.png), [report](../country_burden_%d/REPORT.md).", y,y,y,y,y,y,y,y), ""), "",
  "Reproduce with the three annual burden stages using `--year=2005`, `--year=2015` and `--year=2024` (default `--model=separate`), followed by `Rscript R_cbh/burden/04_summarize_primary_update.R`. No mortality model is refitted by these scripts.")
writeLines(report, file.path(out, "REPORT.md"))
inputs <- c(inputs, "R_cbh/burden/04_summarize_primary_update.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")), file.path(out, "provenance.csv"))
print(s, row.names = FALSE)
message("Combined primary update complete: ", out)

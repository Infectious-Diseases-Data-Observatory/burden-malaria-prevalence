#!/usr/bin/env Rscript
# Compare saved aggregate curves. No model refitting or individual-level data reads.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/sensitivity/sahel/settings.R")
library(ggplot2)
st <- cbh_sahel_settings(); primary <- cbh_primary_settings()
formula_path <- file.path(st$out, "model_formula.txt")
writeLines(trimws(readLines(formula_path), which = "right"), formula_path)
read <- function(x) cbh_read_csv(file.path(st$out, x))
write_csv <- function(d, name) cbh_atomic_csv(d, file.path(st$out, name))
ages <- cbh_config()$age_bands$age_band
diag <- read("fit_diagnostics.csv")
stopifnot(nrow(diag) == 7L, all(diag$converged), all(diag$input_verified),
  all(diag$gamma == 2), all(diag$rank == diag$coefficients),
  all(diag$min_smoothing_hessian_eigenvalue > 0))
sahel <- read("pfpr_curves.csv")
full_path <- file.path(primary$out, "pfpr_curves.csv")
full <- cbh_read_csv(full_path)
pd <- cbh_read_csv(file.path(primary$out, "fit_diagnostics.csv"))
stopifnot(all(pd$gamma == 2), all(pd$converged), all(full$series == "map_full"))
full <- full[names(sahel)]; full$series <- "Full analysis"; sahel$series <- "Sahel (≥12°N)"
d <- rbind(full, sahel)
d$hazard_ratio <- exp(d$log_hazard_ratio)
d$hr_lower_95 <- exp(d$lower_95); d$hr_upper_95 <- exp(d$upper_95)
stopifnot(all(d$hazard_ratio[d$pfpr_pct == 20] == 1),
  all(is.finite(d$hazard_ratio)), all(d$reference_within_observed_support))
write_csv(d, "pfpr_curves_comparison.csv")
write_csv(unique(d[c("series", "age_band", "pfpr_min", "pfpr_p025", "pfpr_p975", "pfpr_max")]),
          "exposure_support.csv")
d$age_label <- factor(d$age_band, levels = ages,
  labels = ifelse(ages == "<1", "<1 month", paste(ages, "months")))
d$series <- factor(d$series, levels = c("Full analysis", "Sahel (≥12°N)"))
colours <- c("Full analysis" = "#215E91", "Sahel (≥12°N)" = "#C34D26")
paper <- theme_minimal(base_size = 20) + theme(panel.grid.minor = element_blank(),
  axis.text = element_text(size = 16), strip.text = element_text(size = 18, face = "bold"),
  legend.position = "bottom", legend.text = element_text(size = 17),
  panel.spacing = grid::unit(24, "pt"), plot.margin = margin(10, 24, 10, 10))
shown <- d[d$within_central_support, ]
p <- ggplot(shown, aes(pfpr_pct, hazard_ratio, colour = series, fill = series)) +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = .4) +
  geom_vline(xintercept = 20, colour = "grey80", linewidth = .4) +
  geom_ribbon(aes(ymin = hr_lower_95, ymax = hr_upper_95), alpha = .14, colour = NA) +
  geom_line(linewidth = 1.05) + facet_wrap(~age_label, ncol = 4) +
  scale_x_continuous(limits = c(0,100), breaks = c(0,20,40,60,80,100)) +
  scale_y_log10(breaks = c(.1,.2,.5,1,2,5,10,20,50), labels = scales::label_number()) +
  scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
  labs(x = "PfPR[2–10] (%)", y = "Hazard ratio relative to PfPR = 20%\n(log scale)",
       colour = NULL, fill = NULL) + paper
ggsave(file.path(st$out, "sahel_vs_full_pfpr_hazard_ratios.png"), p,
       width = 14, height = 8, dpi = 240, device = ragg::agg_png, bg = "white")
p_log <- ggplot(shown, aes(pfpr_pct, log_hazard_ratio, colour = series, fill = series)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = .4) +
  geom_vline(xintercept = 20, colour = "grey80", linewidth = .4) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = .14, colour = NA) +
  geom_line(linewidth = 1.05) + facet_wrap(~age_label, ncol = 4) +
  scale_x_continuous(limits = c(0,100), breaks = c(0,20,40,60,80,100)) +
  scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
  labs(x = "PfPR[2–10] (%)", y = "Log hazard ratio\nrelative to PfPR = 20%", colour = NULL, fill = NULL) + paper
ggsave(file.path(st$out, "sahel_vs_full_pfpr_log_hazard_ratios.png"), p_log,
       width = 14, height = 8, dpi = 240, device = ragg::agg_png, bg = "white")

# Both contrasts are directly identified from the saved curves relative to 20.
# No subtraction of independent variances: the full sample contains the subgroup.
contrasts <- do.call(rbind, lapply(c(0,40), function(x) {
  z <- d[d$pfpr_pct == x, ]
  sign <- if (x == 40) -1 else 1
  loghr <- sign * z$log_hazard_ratio
  data.frame(series = z$series, age_band = z$age_band,
    contrast = if (x == 40) "40% to 20%" else "20% to 0%",
    hazard_ratio = exp(loghr), lower_95 = exp(loghr - 1.96*z$standard_error),
    upper_95 = exp(loghr + 1.96*z$standard_error),
    both_values_within_observed_support = x >= z$pfpr_min & x <= z$pfpr_max &
      z$reference_within_observed_support,
    both_values_within_central_support = pmin(x,20) >= z$pfpr_p025 & pmax(x,20) <= z$pfpr_p975)
}))
write_csv(contrasts, "pfpr_contrasts.csv")
# A geographic audit figure shows every represented survey-region centroid.
regions <- read("region_selection.csv")
regions$selection <- ifelse(regions$selected, "Included", "Excluded")
map <- ggplot(regions, aes(lon, lat, colour = selection)) +
  geom_hline(yintercept = 12, linetype = 2, colour = "grey40") +
  geom_point(alpha = .5, size = 1.8) + coord_quickmap() +
  scale_colour_manual(values = c("Excluded" = "grey70", "Included" = "#C34D26")) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)", colour = NULL) + paper
ggsave(file.path(st$out, "region_selection.png"), map,
       width = 9, height = 8, dpi = 180, device = ragg::agg_png, bg = "white")
sample <- read("sample_summary.csv"); countries <- read("country_sample.csv")
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
z <- sample[sample$series == "sahel_12N", ]
support <- unique(d[d$series == "Sahel (≥12°N)", c("age_band", "pfpr_p025", "pfpr_p975")])
lines <- c("# Sahel-only sensitivity of the primary MAP models", "",
  paste0("Retained ", fmt(z$distinct_children), " children, ", fmt(z$records), " child–age-band records and ",
         fmt(z$deaths), " deaths from ", z$survey_regions, " survey regions in ", z$surveys,
         " surveys and ", z$countries, " countries: ", paste(countries$country, collapse = ", "), "."), "",
  "## Geographic definition", "",
  "Include whole survey regions whose cached boundary centroid is at or north of 12°N and west of 36°E, excluding Ethiopia, Eritrea, Somalia and Djibouti. The western/Horn rule follows the project's existing Sahel convention; the requested 12°N replaces the 11°N cutoff only for this sensitivity. All 1,015 primary survey-region keys matched exactly to survey-specific centroid records; no missing locations were dropped. `region_selection.csv` records inclusion for every primary survey region.", "",
  "This is an approximate geographic proxy for seasonal transmission, not a classification from measured seasonality. Centroids are not child locations or population-weighted locations. Whole regions can cross the cutoff. Nigeria's broad northern-zone centroids are below 12°N, so no Nigerian regions enter; do not interpret this as evidence of absent seasonal transmission in northern Nigeria. The older 11°N monthly mortality analysis is unchanged.", "",
  "## Model and comparison", "",
  "Refit seven separate age-band models with primary MAP prevalence, gamma=2, cr PfPR splines (k=5), cr calendar-year splines (k=6), binomial complementary-log-log likelihood, fREML, discrete bam and the fixed full-band-duration offset. Keep the same confounders, globally computed covariate scaling, fixed median child HIV incidence imputation and unweighted likelihood. Each age/subgroup model has its own time curve, confounder effects and survey/country/region random effects. No outcome, exposure or HIV data are rebuilt.", "",
  "Knots are subgroup quantiles of distinct predictor values, matching the cr construction convention; dimensions are unchanged and knots are recorded in `knots.csv`. This adapts basis locations and re-estimates smoothing parameters for the subgroup. The reference is the saved primary full-sample MAP gamma=2 fit. Prepared-data hashes are checked against its manifest.", "",
  "Figure caption: Adjusted mortality hazard ratios across MAP PfPR₂–₁₀ in the full primary sample (blue) and Sahel subset defined by regional centroids ≥12°N (orange), by completed-month age band. Both curves are relative to 20% PfPR; ribbons are pointwise 95% model-based intervals conditional on fitted smoothing parameters, fixed exposure and the fixed HIV imputation. Each curve is drawn over its own 2.5th–97.5th percentile exposure range. The hazard-ratio figure uses a logarithmic y-axis; a log-hazard-ratio version is also provided. Covariates and fitting choices match the primary models. The subgroup overlaps the full sample, so these curves are dependent; this comparison does not provide an independent-sample interaction test or establish a seasonal causal effect.", "",
  "[Hazard-ratio overlay](sahel_vs_full_pfpr_hazard_ratios.png) · [Log-hazard-ratio overlay](sahel_vs_full_pfpr_log_hazard_ratios.png) · [Region selection](region_selection.png)", "",
  "## Supported 40% to 20% contrasts", "",
  "Hazard ratios below one imply lower mortality at 20% than at 40%. Intervals are conditional pointwise 95% intervals; compare descriptively, without assuming the fits are independent.", "",
  "| Age (months) | Full sample HR (95% interval) | Sahel HR (95% interval) | Sahel central PfPR range (%) |",
  "|---|---:|---:|---:|")
for (age in ages) {
  a <- contrasts[contrasts$age_band == age & contrasts$contrast == "40% to 20%", ]
  b <- support[support$age_band == age, ]
  hr <- function(x) sprintf("%.2f (%.2f–%.2f)", x$hazard_ratio, x$lower_95, x$upper_95)
  lines <- c(lines, sprintf("| %s | %s | %s | %.1f–%.1f |", age,
    hr(a[a$series == "Full analysis", ]), hr(a[a$series == "Sahel (≥12°N)", ]), b$pfpr_p025, b$pfpr_p975))
}
example <- contrasts[contrasts$age_band == "12-23" & contrasts$contrast == "40% to 20%", ]
lines <- c(lines, "", "The Sahel curves show a stronger positive association above 20% PfPR from six months onward, particularly at 12–23 months; the curves below six months are closer. These are descriptive differences in fitted shapes, not a formal test of effect modification by seasonality.", "",
  paste0("For 12–23 months, the 40%→20% hazard ratio is ",
    hr(example[example$series == "Sahel (≥12°N)", ]), " in the Sahel subset versus ",
    hr(example[example$series == "Full analysis", ]), " in the full analysis."), "",
  "The CSV also records 20% to zero contrasts. Zero is outside observed support in both samples; those contrasts extrapolate and are not shown as supported curve segments.", "",
  "## Validation and outputs", "",
  "All seven subgroup models converged, with full rank, finite coefficients/covariances and positive smoothing-Hessian eigenvalues. Every fitted model column was checked against the selected input. Curves were checked against direct link predictions and anchored exactly at HR=1 at 20% PfPR. Input hashes, knots, fit manifests, aggregate sample counts, smoothing EDF and numerical diagnostics are saved alongside this report. Intervals do not account fully for DHS sampling design, HIV/exposure uncertainty or cross-model covariance.", "",
  paste0("Pure fitting time across seven models: ", round(sum(diag$elapsed_seconds), 1), " seconds (excluding reading/writing)."), "",
  "Reproduce from the project root:", "", "```sh", "Rscript R_cbh/sensitivity/sahel/01_fit.R",
  "Rscript R_cbh/sensitivity/sahel/02_report.R", "```", "",
  "Only aggregate tables/figures are published. Individual-level prepared input and fitted model objects remain in the ignored data directory. Primary models, burden estimates and manuscript files are unchanged.")
writeLines(lines, file.path(st$out, "REPORT.md"))
source_files <- c(full_path, file.path(primary$out,"fit_diagnostics.csv"),
  file.path(st$out,"pfpr_curves.csv"), file.path(st$out,"fit_diagnostics.csv"),
  "R_cbh/sensitivity/sahel/02_report.R")
write_csv(data.frame(file = source_files, md5 = vapply(source_files, cbh_file_hash, "")), "report_provenance.csv")
message("Sahel versus primary overlays and report complete")

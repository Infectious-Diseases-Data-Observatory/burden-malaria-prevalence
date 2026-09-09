#!/usr/bin/env Rscript
# Render aggregate sensitivity results without reading research microdata.
source("R_cbh/load_pipeline.R")
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% "--include-south")) stop("Argument: --include-south")
include_south <- "--include-south" %in% args
suffix <- if (include_south) "sensitivity_single_imputation_with_south" else "sensitivity_single_imputation"
out <- file.path("results/cbh/age_band_hiv_incidence_shared_time_v3", suffix)
d <- cbh_read_csv(file.path(out, "pfpr_curves.csv"))
diag <- cbh_read_csv(file.path(out, "fit_diagnostics.csv"))
surveys <- cbh_read_csv(file.path(out, "survey_groups.csv"))
stopifnot(nrow(diag) == 12, all(diag$converged), all(diag$finite_covariance), all(diag$finite_coefficients))
count <- cbh_read_csv(file.path(out, "sample_by_fit_age.csv"))
base_count <- count[count$fit_id == "reference", ]
for (ids in list(paste0("age_", 1:7), c("early", "late"))) {
  check <- aggregate(cbind(rows, deaths) ~ age_band, count[count$fit_id %in% ids, ], sum)
  j <- match(base_count$age_band, check$age_band)
  stopifnot(!anyNA(j), identical(base_count$rows, check$rows[j]), identical(base_count$deaths, check$deaths[j]))
}
stopifnot(length(unique(paste(d$fit_id, d$age_band))) == 42,
  !anyDuplicated(paste(d$fit_id, d$age_band, d$pfpr_pct)),
  all(d$lower_95 <= d$log_hazard_ratio), all(d$upper_95 >= d$log_hazard_ratio))
age <- cbh_config()$age_bands$age_band
cutoff <- median(surveys$survey_year)
early_label <- paste0("Early (", min(surveys$survey_year), "-", cutoff, ")")
late_label <- paste0("Late (", min(surveys$survey_year[surveys$period == "late"]), "-", max(surveys$survey_year), ")")
east_label <- if (include_south) "East, Central and South" else "East and Central"
labels <- c(reference = "Joint model", separate = "Separate age fits", west = "West Africa",
  east_central = east_label, early = early_label, late = late_label)
colors <- setNames(c("#4C566A", "#176B87", "#C35423", "#297E56", "#775AA6", "#D09012"), labels)
d$model <- ifelse(grepl("^age_", d$fit_id), "separate", d$fit_id)
d$series <- factor(unname(labels[d$model]), levels = unname(labels))
d$age_label <- factor(d$age_band, levels = age,
  labels = ifelse(age == "<1", "<1 month", paste0(age, " months")))
stopifnot(!anyNA(d$series))
z <- d[d$pfpr_pct == 20, ]
stopifnot(all(abs(z$log_hazard_ratio) < 1e-10), all(z$standard_error < 1e-10))
central <- d[d$within_central_support, ]
y_limits <- range(c(central$lower_95, central$upper_95))
y_limits <- y_limits + c(-1, 1) * .025 * diff(y_limits)
groups <- list("Separate age-band fits" = c("reference", "separate"),
  "Geographic split" = c("reference", "west", "east_central"),
  "Survey-period split" = c("reference", "early", "late"))
caption <- paste("Single median HIV imputation. Shading: pointwise 95% conditional model intervals.",
  "Each curve is shown over its own central 95% PfPR range; common axes across plots.",
  "Intervals omit HIV-imputation, survey-design and smoothing-parameter uncertainty.", sep = "\n")
style <- function(p) p +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed", linewidth = .35) +
  geom_vline(xintercept = 20, colour = "grey85", linewidth = .3) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = series), alpha = .10, colour = NA) +
  geom_line(linewidth = .8) +
  scale_colour_manual(values = colors, drop = TRUE) +
  scale_fill_manual(values = colors, drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100), expand = expansion(mult = .015)) +
  coord_cartesian(ylim = y_limits) +
  labs(x = "PfPR[2-10] (%)", y = "Log mortality hazard ratio\nrelative to PfPR = 20%", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold"), plot.caption = element_text(hjust = 0, size = 9))
files <- c("separate_age_splines", "geography_splines", "period_splines")
overview <- list()
for (i in seq_along(groups)) {
  x <- central[central$model %in% groups[[i]], ]
  x$comparison <- names(groups)[i]; overview[[i]] <- x
  p <- style(ggplot(x, aes(pfpr_pct, log_hazard_ratio, colour = series, group = series))) +
    facet_wrap(~ age_label, ncol = 4) +
    labs(title = names(groups)[i], subtitle = "Same HIV-imputed data and adjustment set; smoothing re-estimated for each fit", caption = caption)
  ggsave(file.path(out, paste0(files[i], ".png")), p, width = 14, height = 8.5, dpi = 180, device = ragg::agg_png, bg = "white")
}
ov <- do.call(rbind, overview)
ov$comparison <- factor(ov$comparison, levels = names(groups))
p <- style(ggplot(ov, aes(pfpr_pct, log_hazard_ratio, colour = series, group = series))) +
  facet_grid(comparison ~ age_label) +
  labs(title = "PfPR spline sensitivity analyses", subtitle = "Reference: the joint seven-band model with one shared calendar-year spline", caption = caption) +
  theme(strip.text.y = element_text(angle = 0, size = 9))
ggsave(file.path(out, "all_sensitivity_splines.png"), p, width = 23, height = 11, dpi = 180, device = ragg::agg_png, bg = "white")
ggsave(file.path(out, "all_sensitivity_splines.pdf"), p, width = 23, height = 11, device = cairo_pdf, bg = "white")

# A familiar numerical contrast accompanies the curves, with explicit support flags.
hr <- d[d$pfpr_pct == 40, c("fit_id", "model", "series", "age_band", "log_hazard_ratio", "standard_error",
  "pfpr_min", "pfpr_p025", "pfpr_p975", "pfpr_max", "reference_within_observed_support")]
hr$hazard_ratio_40_to_20 <- exp(-hr$log_hazard_ratio)
hr$lower_95 <- exp(-hr$log_hazard_ratio - 1.96 * hr$standard_error)
hr$upper_95 <- exp(-hr$log_hazard_ratio + 1.96 * hr$standard_error)
hr$contrast_within_observed_support <- hr$pfpr_min <= 20 & hr$pfpr_max >= 40
hr$contrast_within_central_support <- hr$pfpr_p025 <= 20 & hr$pfpr_p975 >= 40
cbh_atomic_csv(hr, file.path(out, "pfpr_40_to_20_contrasts.csv"))
reference <- d[d$model == "reference", c("age_band", "pfpr_pct", "log_hazard_ratio", "within_central_support")]
names(reference)[3:4] <- c("reference_log_hr", "reference_central")
delta <- merge(d[d$model != "reference", ], reference, by = c("age_band", "pfpr_pct"))
delta <- delta[delta$within_central_support & delta$reference_central & delta$pfpr_pct %% .5 == 0, ]
delta$difference_from_reference <- delta$log_hazard_ratio - delta$reference_log_hr
dist <- do.call(rbind, lapply(split(delta, paste(delta$fit_id, delta$age_band)), function(x) {
  data.frame(fit_id = x$fit_id[1], age_band = x$age_band[1],
    overlap_min = min(x$pfpr_pct), overlap_max = max(x$pfpr_pct), grid_points = nrow(x),
    rms_log_hr_difference = sqrt(mean(x$difference_from_reference^2)),
    max_absolute_log_hr_difference = max(abs(x$difference_from_reference)))
}))
cbh_atomic_csv(dist, file.path(out, "curve_differences_from_reference.csv"))
support <- unique(d[c("fit_id", "age_band", "pfpr_min", "pfpr_p025", "pfpr_p975", "pfpr_max", "reference_within_observed_support")])
cbh_atomic_csv(support, file.path(out, "pfpr_support.csv"))
md <- c("# PfPR spline sensitivity analyses", "",
  "All comparisons use the same saved posterior-median HIV-incidence imputation and the same complete-case records within each subset. Other confounders, scaling, unweighted likelihood and band-width offset are unchanged. These are separate age, geography and period sensitivities, not a full crossed set of subgroup fits.", "",
  "## Model changes", "",
  "- **Separate age fits:** seven independent models. Each has one PfPR spline (k=5), one time spline (k=6), the existing confounders and survey/country/region random intercepts. Unlike the reference, time functions and all random-effect variances can differ by age. Country replaces country-by-age because each fit contains one age band.",
  "- **Geography:** two joint seven-band fits, each retaining one shared time spline and the reference random-effect structure.",
  paste0("- **Period:** early surveys (year <= ", cutoff, ") versus late surveys (year > ", cutoff, "). The cutoff is the median metadata year across the ", nrow(surveys), " included surveys, with each survey counted once. All ties are assigned to early, and each survey stays intact."),
  "- Smoothing parameters are re-estimated for every new fit. Cubic basis dimensions match the reference; knots adapt to each fitting sample.", "",
  "## Geographic groups", "",
  "Country assignments follow [UN M49](https://unstats.un.org/unsd/methodology/m49/overview/), accessed 9 September 2026. Central Africa means the UN Middle Africa grouping.", "",
  paste0("- West: ", paste(sort(unique(surveys$country[surveys$geography == "west"])), collapse = ", "), "."),
  paste0("- ", east_label, ": ", paste(sort(unique(surveys$country[surveys$geography == "east_central"])), collapse = ", "), "."),
  if (!include_south) "- Namibia, Eswatini and South Africa are excluded only from the geographic comparison; they remain in the reference, age and period fits." else "- Southern Africa is included with East and Central Africa.", "",
  "## Spline comparisons", "",
  "Curves are f_g(P) - f_g(20), the log mortality hazard ratio relative to 20% PfPR. This common reference removes arbitrary spline centering. Curves and ribbons are shown within each fit's central 95% exposure range, on identical axes. Tables include observed and central support flags; predictions outside observed support are not used in the plots. Ribbons are pointwise 95% conditional model intervals, not simultaneous bands.", "",
  "![All sensitivity curves](all_sensitivity_splines.png)", "",
  "Larger individual plots: [age](separate_age_splines.png), [geography](geography_splines.png), [period](period_splines.png). [Vector PDF](all_sensitivity_splines.pdf).", "",
  "## Sample sizes and fitting diagnostics", "",
  "| Fit | Child-band records | Deaths | Countries | Surveys | Fit seconds | Converged |",
  "|---|---:|---:|---:|---:|---:|---|",
  vapply(seq_len(nrow(diag)), function(i) paste0("| ", diag$fit_id[i], " | ", format(diag$rows[i], big.mark = ","), " | ",
    format(diag$deaths[i], big.mark = ","), " | ", diag$countries[i], " | ", diag$surveys[i], " | ",
    round(diag$elapsed_seconds[i], 1), " | ", diag$converged[i], " |"), character(1)), "",
  "## PfPR 40% to 20% hazard ratios", "",
  paste0("| Age (months) | ", paste(labels, collapse = " | "), " |"),
  paste0("|---|", paste(rep("---:", length(labels)), collapse = "|"), "|"),
  vapply(age, function(a) {
    x <- hr[hr$age_band == a, ]; x <- x[match(names(labels), x$model), ]
    paste0("| ", a, " | ", paste(sprintf("%.2f (%.2f–%.2f)%s", x$hazard_ratio_40_to_20, x$lower_95, x$upper_95,
      ifelse(x$contrast_within_central_support, "", "*")), collapse = " | "), " |")
  }, character(1)), "",
  "Values are HR (pointwise 95% interval). * denotes a contrast outside that fit's central 95% exposure range; see CSV for observed-range flags.", "",
  "## Interpretation and limitations", "",
  "- Separate age fits change time adjustment and variance pooling as well as separating the data; differences cannot be attributed only to removing shared PfPR information (the reference already has separate age-specific PfPR curves).",
  "- Geographic and period differences can reflect country/sample composition, exposure support and residual confounding as well as effect heterogeneity. Period is survey year, not the child's exposure year.",
  "- A single HIV imputation is held fixed. Intervals omit imputation, survey-design, residual-clustering, MAP and smoothing-parameter uncertainty. No survey weights are used.",
  "- Comparisons with the reference use overlapping data. Ribbon overlap and the descriptive curve-distance measures are not formal tests of differences.",
  paste0("- All fits report convergence and finite coefficients/covariance. ", sum(diag$min_smoothing_hessian_eigenvalue < -1e-8, na.rm = TRUE),
    " fits (including the reference if applicable) have a negative smoothing-Hessian eigenvalue below -1e-8; see diagnostics. Smoothing stability remains an exploratory limitation."), "",
  "## Reproduction and files", "",
  "Run from the project root:", "", "```sh",
  paste("Rscript R_cbh/sensitivity/01_fit.R", if (include_south) "--include-south" else ""),
  paste("Rscript R_cbh/sensitivity/02_report.R", if (include_south) "--include-south" else ""), "```", "",
  "Private fits remain under the ignored data/derived_cbh/models directory. Outputs here contain aggregate statistics only.", "",
  "- `survey_groups.csv`: exact survey-year and geographic assignments.",
  "- `pfpr_curves.csv`, `pfpr_support.csv`: predictions, intervals and exposure support.",
  "- `pfpr_40_to_20_contrasts.csv`: numerical comparisons.",
  "- `curve_differences_from_reference.csv`: descriptive RMS and maximum absolute log-HR differences, evaluated on the common central exposure range at 0.5 percentage-point spacing.",
  "- `fit_diagnostics.csv`, `sample_by_fit_age.csv`, `*_smooth_summary.csv`: fit checks, samples and smooth EDFs.",
  "- `*_formula.txt`, `specification.txt`, `scaling.csv`, `code_provenance.csv`: reproducible specifications.")
writeLines(md, file.path(out, "REPORT.md"))
message("Saved comparison plots and report: ", out)

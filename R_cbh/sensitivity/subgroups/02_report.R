#!/usr/bin/env Rscript
# Overlay subgroup PfPR curves on the full primary curves; saved aggregates only.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
source("R_cbh/sensitivity/subgroups/settings.R")
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE); stopifnot(all(args %in% "--dhs-only"))
st <- cbh_subgroup_settings(cbh_subgroup_sample(args)); primary <- cbh_primary_settings(st$primary_version)
read <- function(x) cbh_read_csv(file.path(st$out, x))
write_csv <- function(d, name) cbh_atomic_csv(d, file.path(st$out, name))
ages <- cbh_config()$age_bands$age_band
subgroups <- c("sahel", "east_africa", "early", "late")
diag <- read("fit_diagnostics.csv")
stopifnot(nrow(diag) == 28L, all(diag$converged), all(diag$input_verified), all(diag$gamma == 2),
  all(diag$rank == diag$coefficients), all(diag$min_smoothing_hessian_eigenvalue > 0))
sub <- read("pfpr_curves.csv")
full_path <- file.path(primary$out, "pfpr_curves.csv")
full <- cbh_read_csv(full_path)
pd <- cbh_read_csv(file.path(primary$out, "fit_diagnostics.csv"))
stopifnot(all(pd$gamma == 2), all(pd$converged), all(full$series == "map_full"))
full$subgroup <- "map_full"; full <- full[names(sub)]
d <- rbind(full, sub)
d$hazard_ratio <- exp(d$log_hazard_ratio)
stopifnot(all(abs(d$log_hazard_ratio[d$pfpr_pct == 20]) < 1e-10), all(is.finite(d$hazard_ratio)),
  all(d$reference_within_observed_support))
d$series <- factor(st$labels[d$subgroup], levels = unname(st$labels))
write_csv(d, "pfpr_curves_comparison.csv")
write_csv(unique(d[c("subgroup", "series", "age_band", "pfpr_min", "pfpr_p025", "pfpr_p975", "pfpr_max")]),
          "exposure_support.csv")
d$age_label <- factor(d$age_band, levels = ages, labels = ifelse(ages == "<1", "<1 month", paste(ages, "months")))
colours <- setNames(unname(st$colours), unname(st$labels))
paper <- theme_minimal(base_size = 18) + theme(panel.grid.minor = element_blank(),
  axis.text = element_text(size = 14), axis.title = element_text(size = 18),
  strip.text = element_text(size = 16, face = "bold"), legend.position = "bottom",
  legend.text = element_text(size = 15), legend.title = element_blank(),
  panel.spacing = grid::unit(20, "pt"), plot.margin = margin(10, 24, 10, 10))
shown <- d[d$within_central_support, ]
# PfPR axis 0-80%, as in the other supplementary spline figures; every curve's central support lies inside it.
stopifnot(max(shown$pfpr_pct) <= 80)
shown$is_full <- shown$subgroup == "map_full"
p <- ggplot(shown, aes(pfpr_pct, log_hazard_ratio, colour = series)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = .4) +
  geom_vline(xintercept = 20, colour = "grey80", linewidth = .4) +
  geom_ribbon(data = shown[shown$is_full, ], aes(ymin = lower_95, ymax = upper_95),
              fill = "grey55", alpha = .18, colour = NA) +
  geom_line(aes(linewidth = is_full)) +
  scale_linewidth_manual(values = c(`TRUE` = 1.5, `FALSE` = 1), guide = "none") +
  facet_wrap(~age_label, ncol = 4) +
  scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) +
  scale_colour_manual(values = colours) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(linewidth = 1.3))) +
  labs(x = "PfPR[2–10] (%)", y = "Log hazard ratio\nrelative to PfPR = 20%") + paper
stopifnot(is.null(p$labels$title), is.null(p$labels$subtitle))
ggsave(file.path(st$out, "sfig_pfpr_splines_by_subgroup.png"), p, width = 15, height = 8.5,
       dpi = 240, device = ragg::agg_png, bg = "white")

# Contrasts identified directly from each curve; the subgroups are nested in the full sample.
contrasts <- do.call(rbind, lapply(c(0, 40), function(x) {
  z <- d[d$pfpr_pct == x, ]
  sign <- if (x == 40) -1 else 1
  loghr <- sign * z$log_hazard_ratio
  data.frame(subgroup = z$subgroup, series = z$series, age_band = z$age_band,
    contrast = if (x == 40) "40% to 20%" else "20% to 0%",
    hazard_ratio = exp(loghr), lower_95 = exp(loghr - 1.96 * z$standard_error),
    upper_95 = exp(loghr + 1.96 * z$standard_error),
    both_values_within_observed_support = x >= z$pfpr_min & x <= z$pfpr_max & z$reference_within_observed_support,
    both_values_within_central_support = pmin(x, 20) >= z$pfpr_p025 & pmax(x, 20) <= z$pfpr_p975)
}))
write_csv(contrasts, "pfpr_contrasts.csv")
cz <- contrasts[contrasts$contrast == "40% to 20%", ]
cz$age_label <- factor(cz$age_band, levels = rev(ages))
pc <- ggplot(cz, aes(hazard_ratio, age_label, colour = series)) +
  geom_vline(xintercept = 1, colour = "grey60", linetype = 2) +
  geom_pointrange(aes(xmin = lower_95, xmax = upper_95), position = position_dodge(width = .7), size = .5) +
  scale_x_log10(breaks = c(.5, .7, .8, .9, 1, 1.1, 1.25)) + scale_colour_manual(values = colours) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Hazard ratio, PfPR 40% to 20% (95% interval, log scale)", y = "Age (completed months)") + paper
ggsave(file.path(st$out, "pfpr_40_to_20_by_subgroup.png"), pc, width = 10, height = 7, dpi = 240,
       device = ragg::agg_png, bg = "white")

## ---- caption and report --------------------------------------------------------------
sample <- read("sample_summary.csv"); defs <- readLines(file.path(st$out, "subgroup_definitions.csv"))
cutoff <- as.integer(sub("median_survey_year,", "", defs[1]))
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
row <- function(s) sample[sample$subgroup == s, ]
describe <- function(s) sprintf("%s: %s children, %s records, %s deaths, %d survey-regions, %d surveys, %d countries",
  st$labels[[s]], fmt(row(s)$distinct_children), fmt(row(s)$records), fmt(row(s)$deaths),
  row(s)$survey_regions, row(s)$surveys, row(s)$countries)
caption <- paste(
  paste0("Sensitivity of the fitted PfPR–mortality relationship to the analysis subset. Log mortality hazard ratios relative to PfPR[2–10] = 20% by completed-month age band, from the ", cbh_paper_model_label(), " fitted to the full primary sample (black, with its pointwise 95% conditional interval shaded) and refitted separately in four subsets: survey regions with boundary centroids at or north of 12°N and west of 36°E, excluding the Horn of Africa (Sahel); UN M49 Eastern Africa countries; surveys conducted up to the median survey year (", cutoff, "); and surveys conducted after it."),
  "Each subset refit uses the same seven separate age-band models, 17 standardised covariates with the full-sample scaling, fixed HIV incidence imputation, gamma = 2 and basis dimensions as the primary analysis, with knots placed at quantiles of the subset's own predictor values. Each curve is drawn over the central 95% of its own exposure distribution.",
  paste0("Subset sizes. ", paste(vapply(subgroups, describe, ""), collapse = "; "), "."),
  "The subsets overlap the full sample and one another, so the curves are not independent and their differences are descriptive; the comparison does not constitute a test of effect modification. Subset intervals are omitted for legibility and are available in the saved curve table.")
writeLines(c("# Supplementary figure caption: PfPR splines by subgroup", "", caption), file.path(st$out, "CAPTION.md"))
hr <- function(x) sprintf("%.2f (%.2f–%.2f)", x$hazard_ratio, x$lower_95, x$upper_95)
table_lines <- c(paste0("| Age (months) | ", paste(unname(st$labels), collapse = " | "), " |"),
  paste0("|---|", paste(rep("---:", length(st$labels)), collapse = "|"), "|"))
for (age in ages) {
  a <- cz[cz$age_band == age, ]
  table_lines <- c(table_lines, paste0("| ", age, " | ",
    paste(vapply(names(st$labels), function(s) hr(a[a$subgroup == s, ]), ""), collapse = " | "), " |"))
}
lines <- c("# Subgroup sensitivity of the primary PfPR models", "",
  paste0("Seven separate MAP gamma=2 age-band models refitted in four subsets of the fitted 17-variable sample (`", primary$id, "`). Reference: the full primary fits in `", primary$out, "`."), "",
  "## Subsets", "", paste0("- ", vapply(subgroups, describe, "")), "",
  paste0("Definitions: ", paste(defs[-1], collapse = "; "), ". Survey regions enter whole; centroids are survey-specific boundary centroids (DHS: the cached shapefile summary; MICS: the analysis-region polygons, `R_mics/11_region_centroids.R`). The period split counts each survey once."), "",
  "## Figures", "",
  "![Splines by subgroup](sfig_pfpr_splines_by_subgroup.png)", "", caption, "",
  "![40% to 20% contrasts](pfpr_40_to_20_by_subgroup.png)", "",
  "## Hazard ratios, PfPR 40% to 20%", "", "Values below one indicate lower mortality at 20% than at 40%. Intervals are pointwise conditional 95% intervals; the samples are nested, so compare descriptively.", "",
  table_lines, "",
  "## Validation", "",
  sprintf("All 28 subgroup fits converged with full rank, finite coefficients and covariances and positive smoothing-Hessian eigenvalues; %d needed the tighter-tolerance restart (see restarts.csv where present). Fitted model columns were checked against the selected input rows; curves are anchored at zero at 20%% PfPR. Total pure fitting time %.0f seconds.",
    sum(grepl("strict", cbh_read_csv(file.path(st$out, "fit_manifest.csv"))$selected_version)), sum(diag$elapsed_seconds)), "",
  paste0("Reproduce: `Rscript R_cbh/sensitivity/subgroups/01_fit.R` then `Rscript R_cbh/sensitivity/subgroups/02_report.R`", if (st$sample == "dhs") " (both with `--dhs-only`)" else "", ". Only aggregate tables and figures are published; fitted objects stay under the ignored data directory. Primary fits, burden estimates and manuscript files are unchanged."))
writeLines(lines, file.path(st$out, "REPORT.md"))
source_files <- c(full_path, file.path(primary$out, "fit_diagnostics.csv"),
  file.path(st$out, c("pfpr_curves.csv", "fit_diagnostics.csv", "sample_summary.csv")),
  "R_cbh/sensitivity/subgroups/02_report.R", "R_cbh/sensitivity/subgroups/settings.R", "R_cbh/reporting/labels.R")
write_csv(data.frame(file = source_files, md5 = vapply(source_files, cbh_file_hash, "")), "report_provenance.csv")
message("Subgroup overlays and report complete: ", st$out)

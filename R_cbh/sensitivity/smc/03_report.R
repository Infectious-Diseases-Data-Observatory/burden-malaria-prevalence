#!/usr/bin/env Rscript
# Report for the SMC before/after analysis: SMC hazard ratios by age band (main =
# confirmed status; placebo; sensitivities), PfPR curves with and without the SMC term
# on the main sample, sample accounting and interpretation limits. Saved aggregates only.
source("R_cbh/load_pipeline.R")
source("R_cbh/reporting/labels.R")
source("R_cbh/sensitivity/smc/settings.R")
library(ggplot2)
st <- cbh_smc_settings(); read <- function(x) cbh_read_csv(file.path(st$out, x))
ages <- cbh_config()$age_bands$age_band
lab_age <- function(x) ifelse(x == "<1", "<1 month", paste(x, "months"))
diag <- read("fit_diagnostics.csv"); eff <- read("smc_effects.csv"); curves <- read("pfpr_curves.csv")
samp <- read("sample_summary.csv"); ctry <- read("country_summary.csv"); units <- read("admin1_smc_years.csv")
timing <- read("timing_pre_season.csv"); omitted <- read("omitted_smc_countries.csv")
restarts <- if (file.exists(file.path(st$out, "restarts.csv"))) read("restarts.csv") else data.frame(strict_selected = logical())
stopifnot(nrow(diag) == 7L * length(st$variants), all(diag$converged), all(diag$input_verified),
  all(diag$rank == diag$coefficients), all(diag$min_smoothing_hessian_eigenvalue > 0), nrow(eff) == 28L)
paper <- theme_minimal(base_size = 17) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
  legend.title = element_blank(), strip.text = element_text(face = "bold"), plot.margin = margin(10, 20, 10, 10))
shown <- c("main", "placebo", "national", "coverage")
eff$series <- factor(st$labels[eff$variant], levels = unname(st$labels[shown]))
eff$age <- factor(lab_age(eff$age_band), levels = rev(lab_age(ages)))
colours <- setNames(c("#111111", "#9E9E9E", "#C34D26", "#009E73"), unname(st$labels[shown]))
p <- ggplot(eff, aes(hazard_ratio, age, colour = series, shape = series)) + geom_vline(xintercept = 1, colour = "grey55", linetype = 2) +
  geom_pointrange(aes(xmin = lower_95, xmax = upper_95), position = position_dodge(width = .7), size = .45) +
  scale_x_log10(breaks = c(.5, .6, .7, .8, .9, 1, 1.1, 1.25, 1.5, 2)) + scale_colour_manual(values = colours) +
  scale_shape_manual(values = setNames(c(16, 1, 17, 15), unname(st$labels[shown]))) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  labs(x = "Hazard ratio, after versus before SMC switch-on (95% interval, log scale)", y = NULL) + paper
ggsave(file.path(st$out, "smc_hazard_ratios_by_age.png"), p, width = 11.5, height = 7.5, dpi = 220, device = ragg::agg_png, bg = "white")
cv <- curves[curves$variant %in% c("main", "no_smc") & curves$within_central_support, ]
cv$series <- factor(st$labels[cv$variant], levels = unname(st$labels[c("no_smc", "main")]))
cv$age_label <- factor(lab_age(cv$age_band), levels = lab_age(ages))
cc <- setNames(c("#B36B39", "#215E91"), levels(cv$series))
pc <- ggplot(cv, aes(pfpr_pct, log_hazard_ratio, colour = series, fill = series)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = .4) + geom_vline(xintercept = 20, colour = "grey80", linewidth = .4) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = .12, colour = NA) + geom_line(linewidth = 1) +
  facet_wrap(~age_label, ncol = 4) + scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) +
  scale_colour_manual(values = cc) + scale_fill_manual(values = cc) +
  labs(x = "PfPR[2–10] (%)", y = "Log hazard ratio\nrelative to PfPR = 20%") + paper
ggsave(file.path(st$out, "smc_pfpr_curves_with_without_smc.png"), pc, width = 14, height = 8, dpi = 220, device = ragg::agg_png, bg = "white")

## ---- report ---------------------------------------------------------------------------------------------
fmt <- function(x) format(round(x), big.mark = ",", trim = TRUE)
hr <- function(d) sprintf("%.2f (%.2f–%.2f)", d$hazard_ratio, d$lower_95, d$upper_95)
tab <- c(paste0("| Age | ", paste(unname(st$labels[shown]), collapse = " | "), " |"), paste0("|---|", paste(rep("---:", length(shown)), collapse = "|"), "|"))
for (a in ages) tab <- c(tab, paste0("| ", lab_age(a), " | ", paste(vapply(shown, function(v) hr(eff[eff$variant == v & eff$age_band == a, ]), ""), collapse = " | "), " |"))
row <- function(v) samp[samp$variant == v, ]
m <- row("main"); nat <- row("national"); pl <- row("placebo")
sw_band <- eff[eff$variant == "main", ]; sw_band <- sw_band[match(ages, sw_band$age_band), ]
nat_units <- units[units$status_basis == "national", ]
asg <- read("survey_admin1_assignment.csv"); nat_share <- sum(asg$post_national[asg$basis == "national"]) / sum(asg$post_national)
no_post <- ctry$iso3[ctry$main_post_records == 0 & ctry$national_post_records == 0]
re <- diag[diag$variant == "main", ]; re <- re[match(ages, re$age_band), ]
ps <- timing[match(ages, timing$age_band), ]
lines <- c("# SMC introduction and all-cause mortality by age (DHS and MICS)", "",
  paste0("The seven primary age-band models (`primary_map_regional17_dhsmics_gamma2_v5`: 17 covariates with the full-sample scaling, PfPR and calendar-year splines with the primary reference knots, survey, country and survey-region random intercepts, gamma = 2) refitted in the countries that have introduced seasonal malaria chemoprevention according to `data/SMC_rollout/smc_by_dhs_cluster.csv`, with one extra term: an admin-1 indicator equal to 1 when the band entry year is at or after the unit's SMC switch-on year. Everything else is the primary model. SMC countries in the file: ", nat$countries_list, "; ", paste(no_post, collapse = ", "), " have no bands entered after their switch-on (no v5 survey postdates it), so they inform only the covariates, splines and random effects. DRC and Zambia are in the file with no SMC (true zeros) and are excluded."), "",
  "**Main analysis (decided 23 September 2026): confirmed status only.** The codebook records each switch-on year at the finest scope available and warns that national-scope years are not a treatment flag. The main analysis therefore uses only district- or region-confirmed switch-ons (Nigeria, Burkina Faso, Uganda, Cameroon North and Far North); records of units whose switch-on is known only from national-scope records (Mali, Niger, Chad, Côte d'Ivoire, the rest of Cameroon) are kept while no campaign had been recorded in the unit (a known zero) and dropped from the first recorded campaign year onward (for Mali's Kayes, 2013, from district-confirmed clusters). Treating national-scope years as switch-ons is a sensitivity analysis.", "",
  "## Admin-1 switch-on years", "",
  "Admin-1 units are the analysis regions, except in Nigeria, whose analysis regions are the six zones; there the unit is the state, taken from each record's survey cluster (state variable in the DHS recodes; HH7 in MICS 2016). The zone remains the unit of PfPR, covariates and the survey-region random intercept. A unit's switch-on year is the first year by which at least half of its DHS/MIS clusters, pooled over survey rounds, lay in areas with SMC; unmatched (no-GPS) clusters are dropped and clusters in areas never covered count as never treated. DHS rounds and MICS surveys absent from the file take their unit's year by name. See [admin1_smc_years.csv](admin1_smc_years.csv), [admin1_smc_coverage.csv](admin1_smc_coverage.csv), [survey_admin1_assignment.csv](survey_admin1_assignment.csv) and the [coverage maps](smc_coverage_maps_2015_2022.png).", "",
  "![SMC coverage maps](smc_coverage_maps_2015_2022.png)", "",
  sprintf("**Countries not covered.** The file covers ten countries. Other v5 countries with SMC programmes are not in it and are not analysed: %s (%s records and %s deaths in the v5 sample). The codebook states that Mozambique has district-level SMC records for 2020–2025 that were not merged; the others are to be confirmed with the data provider.",
    paste(omitted$country, collapse = ", "), fmt(sum(omitted$records)), fmt(sum(omitted$deaths))), "",
  "## Sample", "",
  "| Variant | Records | Deaths | Post-SMC records | Post-SMC deaths | Surveys | Countries | Admin-1 units (survey-specific) | Units switching within a survey |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
  sprintf("| %s | %s | %s | %s | %s | %d | %d | %d | %d (%s) |", st$labels[samp$variant], fmt(samp$records), fmt(samp$deaths), fmt(samp$post_records),
    fmt(samp$post_deaths), samp$surveys, samp$countries, samp$admin1_units, samp$switching_units, samp$switching_by_country), "",
  "Placebo: post means after the fake switch-on 3 years early, within the pre-switch-on records. Coverage: post means any coverage.", "",
  "| Country | Records | Deaths | Records in main sample | Post-SMC records, main | Post-SMC records, national years | Admin-1 units | Band entry years |", "|---|---:|---:|---:|---:|---:|---:|---|",
  sprintf("| %s | %s | %s | %s | %s | %s | %d | %s |", ctry$iso3, fmt(ctry$records), fmt(ctry$deaths), fmt(ctry$main_records), fmt(ctry$main_post_records),
    fmt(ctry$national_post_records), ctry$admin1_units, ctry$entry_years), "",
  "## SMC hazard ratios by age band", "",
  "Hazard ratio for bands entered after versus before the unit's SMC switch-on (coverage: full versus none), with Wald 95% intervals conditional on the smoothing parameters. The models adjust for PfPR at band entry, so the estimate is the association of SMC with all-cause mortality beyond any change it brings about in MAP prevalence.", "",
  tab, "", "![SMC hazard ratios](smc_hazard_ratios_by_age.png)", "",
  "## PfPR curves with and without the SMC term (main sample)", "", "![PfPR curves](smc_pfpr_curves_with_without_smc.png)", "",
  "## Interpretation limits", "",
  sprintf("- **Identifying variation.** Admin-1 units switching within a survey, by age band (main): %s. Survey-region random intercepts are estimated near zero in some bands (standard deviation by band: %s), so comparisons between units and between surveys also inform the SMC term, not only within-survey switches; survey-level differences can therefore be confounded with SMC timing.",
    paste(sprintf("%s %d", lab_age(ages), sw_band$switching_units), collapse = ", "), paste(sprintf("%.3f", re$re_sd_region), collapse = ", ")),
  sprintf("- **Falsification.** SMC eligibility starts at 3 months, so an effect in the <1 and 1–5 month bands would point to confounding rather than chemoprevention. The placebo (fake switch-on 3 years early, pre-switch-on records only) should show hazard ratios near 1; departures of the size of the main estimates indicate that the design cannot separate SMC from other changes over time."),
  sprintf("- **National-scope years (sensitivity).** In the national-years sensitivity, units whose switch-on is known only from national-scope records supply %.0f%% of post-SMC records; in those countries every unit switches in the same year, so the contrast is a national before/after comparison against a calendar trend shared across countries.", 100 * nat_share),
  sprintf("- **Timing.** A band is post-SMC if entered in or after the switch-on year. In the main sample %s post-SMC records (%s deaths) end before 1 July of the switch-on year and so precede the first campaign season (by band: %s).",
    fmt(sum(timing$pre_season_records)), fmt(sum(timing$pre_season_deaths)), paste(sprintf("%s %s", lab_age(ages), fmt(ps$pre_season_records)), collapse = ", ")),
  "- **Exposure misclassification.** A unit-level switch-on dates the start of campaigns in at least half of the unit's clusters; it does not measure coverage, cycles or eligible ages. The coverage sensitivity uses the share of confirmed clusters covered instead.",
  "- **Confounding.** SMC roll-out coincided with other interventions and with changes in survey programmes; the calendar-year spline is shared across the SMC countries and cannot absorb country-specific trends. The results are descriptive associations, not causal effects.", "",
  sprintf("All %d fits converged with full rank and positive smoothing-Hessian eigenvalues; %d took the tighter-tolerance restart (subgroup-refit acceptance rule: converged, full rank, positive Hessian, fREML no more than 0.001 worse; restarts.csv records the smoothing gradients and fREML values). Reproduce: `Rscript R_cbh/sensitivity/smc/01_assign_smc.R`, `02_fit.R`, `03_report.R`, `04_maps.R`. Cluster-level SMC data and the Nigerian cluster-to-state crosswalk stay in the ignored data folder; only admin-1 aggregates are published.",
    nrow(diag), sum(restarts$strict_selected)))
writeLines(lines, file.path(st$out, "REPORT.md"))
src <- file.path(st$out, c("fit_diagnostics.csv", "smc_effects.csv", "pfpr_curves.csv", "sample_summary.csv", "country_summary.csv",
  "admin1_smc_years.csv", "timing_pre_season.csv", "omitted_smc_countries.csv"))
src <- c(src, "R_cbh/sensitivity/smc/03_report.R", "R_cbh/sensitivity/smc/settings.R")
cbh_atomic_csv(data.frame(file = src, md5 = vapply(src, cbh_file_hash, "")), file.path(st$out, "report_provenance.csv"))
print(eff[, c("variant", "age_band", "hazard_ratio", "lower_95", "upper_95", "p_value", "switching_units")])
message("SMC report written: ", file.path(st$out, "REPORT.md"))

#!/usr/bin/env Rscript
# Supplementary figure (added 24 September 2026): adjusted hazard ratios for the 17 covariates in
# each of the seven primary age-band models, per 1 SD of the covariate on the modelling sample
# (the scaling used in the fit), with 95% intervals from the coefficient covariance conditional on
# the smoothing parameters (Vp), as for the PfPR contrasts. Reads the saved fits; no refitting.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(mgcv)})
st <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION", "regional_mics"))
out <- file.path(st$out, "covariate_effects"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
ages <- cbh_config()$age_bands$age_band
manifest <- fread(file.path(st$out, "fit_manifest.csv"))[series == "map_full"][match(ages, age_band)]
scaling <- fread(file.path(st$out, "covariate_scaling.csv"))
covs <- cbh_primary_regional_spec(st)$covariates
stopifnot(nrow(manifest) == 7L, identical(scaling$variable, covs))

est <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  path <- manifest$model_file[i]
  stopifnot(identical(cbh_file_hash(path), manifest$md5[i]))
  f <- readRDS(path)$fit; terms <- paste0("z_", covs)
  stopifnot(all(terms %in% names(coef(f))), identical(f$family$link, "cloglog"))
  b <- coef(f)[terms]; se <- sqrt(diag(f$Vp)[match(terms, names(coef(f)))])
  z <- data.table(age_band = manifest$age_band[i], variable = covs, log_hr = unname(b), se = unname(se))
  rm(f); gc(FALSE); z
}))
est[, `:=`(hr = exp(log_hr), lower_95 = exp(log_hr - 1.96 * se), upper_95 = exp(log_hr + 1.96 * se))]
est[, excludes_1 := lower_95 > 1 | upper_95 < 1]
est <- merge(est, scaling[, .(variable, sample_mean = mean, sample_sd = sd)], by = "variable")

# Labels give the size of 1 SD in natural units (fold change for log-transformed covariates).
lab <- c(mean_maternal_age_first_birth = "Maternal age at first birth (per %.1f years)",
  mean_maternal_education_years = "Maternal education (per %.1f years)", mean_wealth_quintile = "Wealth-quintile score (per %.2f)",
  urban_pct = "Urban residence (per %.0f pp)", electricity_pct = "Electricity (per %.0f pp)",
  improved_water_pct = "Improved water (per %.0f pp)", improved_sanitation_pct = "Improved sanitation (per %.0f pp)",
  dtp3_pct = "DTP3 coverage (per %.0f pp)", measles_pct = "Measles coverage (per %.0f pp)",
  facility_delivery_pct = "Facility delivery (per %.0f pp)", short_birth_interval_pct = "Birth interval <24 months (per %.0f pp)",
  wasting_pct = "Wasting (per %.1f pp)", stunting_pct = "Stunting (per %.0f pp)",
  log_hiv_incidence = "Child HIV incidence (per %.1f-fold)", log_gdp_pc = "GDP per capita (per %.1f-fold)",
  log_health_expenditure_pc = "Health expenditure per capita (per %.1f-fold)", political_stability = "Political stability (per %.2f units)")
group <- c(mean_maternal_age_first_birth = "Maternal and\nhousehold", mean_maternal_education_years = "Maternal and\nhousehold",
  mean_wealth_quintile = "Maternal and\nhousehold", urban_pct = "Maternal and\nhousehold", electricity_pct = "Maternal and\nhousehold",
  improved_water_pct = "Maternal and\nhousehold", improved_sanitation_pct = "Maternal and\nhousehold",
  dtp3_pct = "Health\nservices", measles_pct = "Health\nservices", facility_delivery_pct = "Health\nservices",
  short_birth_interval_pct = "Health\nservices", wasting_pct = "Nutrition", stunting_pct = "Nutrition",
  log_hiv_incidence = "National", log_gdp_pc = "National", log_health_expenditure_pc = "National", political_stability = "National")
stopifnot(setequal(names(lab), covs), setequal(names(group), covs))
est[, unit := fifelse(startsWith(variable, "log_"), exp(sample_sd), sample_sd)]
est[, label := sprintf(lab[variable], unit)]
est[, group := factor(group[variable], levels = unique(group[covs]))]
est[, label := factor(label, levels = rev(unique(est[order(match(variable, covs))]$label)))]
est[, age_band := factor(age_band, levels = ages, labels = ifelse(ages == "<1", "<1 month", paste(ages, "months")))]
setorder(est, age_band, group, label)
cbh_atomic_csv(as.data.frame(est[, .(age_band, variable, label, group, log_hr, se, hr, lower_95, upper_95, excludes_1, sample_mean, sample_sd)]),
  file.path(out, "covariate_hazard_ratios.csv"))

rng <- range(c(est$lower_95, est$upper_95)); brk <- c(0.7, 0.8, 0.9, 1, 1.1, 1.25, 1.5)
p <- ggplot(est, aes(hr, label)) +
  geom_vline(xintercept = 1, colour = "grey55", linewidth = .4) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), width = 0, orientation = "y", colour = "#1F4E79", linewidth = .45) +
  geom_point(aes(shape = excludes_1), colour = "#1F4E79", fill = "#1F4E79", size = 1.9, stroke = .6) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), labels = c(`TRUE` = "95% interval excludes 1", `FALSE` = "95% interval includes 1"), name = NULL) +
  scale_x_log10(breaks = brk[brk >= rng[1] * .95 & brk <= rng[2] * 1.05], labels = function(x) sub("\\.?0+$", "", sprintf("%.2f", x))) +
  facet_grid(group ~ age_band, scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Adjusted hazard ratio per 1 SD increase (log scale)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"), strip.text.x = element_text(face = "bold"),
    panel.spacing.x = grid::unit(8, "pt"), panel.border = element_rect(fill = NA, colour = "grey85"),
    axis.text.x = element_text(size = 8.5), legend.position = "bottom")
figure <- file.path(out, "sfig_covariate_forest.png")
ggsave(figure, p, width = 16, height = 7.5, dpi = 300, device = ragg::agg_png, bg = "white")

n_ex <- est[, .(n = sum(excludes_1)), by = age_band]
caption <- paste0("Adjusted hazard ratios for the 17 covariates in each of the seven age-band models of the primary analysis (", st$id, "). ",
  "Each hazard ratio is for a 1 SD increase in the covariate on the modelling sample (the SD in natural units is given in the label; for log-transformed national covariates, as a fold change), ",
  "conditional on PfPR[2–10] and calendar year at band entry, the other covariates, and survey, country and survey-region random intercepts. ",
  "Horizontal bars are 95% intervals from the coefficient covariance conditional on the estimated smoothing parameters; filled points mark intervals that exclude 1. ",
  "Regional covariates are survey-region summaries at the time of the survey; national covariates are matched to the calendar year of band entry. ",
  "The covariates were included to adjust the PfPR association for confounding; their own coefficients are mutually adjusted and are not estimates of causal effects.")
writeLines(c("# Supplementary figure: covariate hazard ratios by age band", "", caption), file.path(out, "CAPTION.md"))
inputs <- c(manifest$model_file, file.path(st$out, c("fit_manifest.csv", "covariate_scaling.csv")), "R_cbh/reporting/14_covariate_forest.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")), file.path(out, "provenance.csv"))
print(n_ex); message("Covariate forest plot written: ", figure)

#!/usr/bin/env Rscript
# Figure 3: age-specific fraction of all-cause mortality attributable to PfPR.
# Read the saved primary spline components; no refitting or national IHME inputs.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
library(mgcv)
library(ggplot2)
st <- cbh_primary_settings(); out <- st$out
component_path <- file.path(st$private, "pfpr_components.rds")
components <- readRDS(component_path)
manifest <- cbh_read_csv(file.path(out, "fit_manifest.csv"))
diagnostics <- cbh_read_csv(file.path(out, "fit_diagnostics.csv"))
curves <- cbh_read_csv(file.path(out, "pfpr_curves.csv"))
ages <- cbh_config()$age_bands$age_band
prevalences <- c(10, 20, 30, 40)
stopifnot(length(components) == 7L, nrow(manifest) == 7L,
  all(manifest$series == "map_full"), all(diagnostics$gamma == 2), all(diagnostics$converged))
rows <- lapply(seq_along(ages), function(i) {
  id <- paste0("map_full_age_", i); piece <- components[[id]]
  stopifnot(identical(piece$model_md5, manifest$md5[match(id, manifest$fit_id)]),
            identical(piece$age_band, ages[i]))
  nd <- data.frame(pfpr_pct = prevalences); zero <- nd; zero$pfpr_pct <- 0
  # HR is zero prevalence versus P, not P versus the 20% plotting reference.
  L <- PredictMat(piece$smooth, zero) - PredictMat(piece$smooth, nd)
  loghr <- drop(L %*% piece$coef)
  variance <- rowSums((L %*% piece$covariance) * L)
  stopifnot(all(is.finite(variance)), min(variance) >= -1e-9)
  se <- sqrt(pmax(variance, 0))
  z <- curves[curves$age_band == ages[i], ]
  expected <- z$log_hazard_ratio[z$pfpr_pct == 0] -
    z$log_hazard_ratio[match(prevalences, z$pfpr_pct)]
  stopifnot(max(abs(loghr - expected)) < 1e-9)
  data.frame(age_band = ages[i], pfpr_pct = prevalences,
    log_hr_zero_vs_current = loghr, log_hr_se = se,
    attributable_fraction = -expm1(loghr),
    af_lower_95 = -expm1(loghr + 1.96*se), af_upper_95 = -expm1(loghr - 1.96*se),
    zero_below_observed_support = piece$support[1] > 0,
    current_outside_central95 = prevalences < piece$support[2] | prevalences > piece$support[3])
})
d <- do.call(rbind, rows)
stopifnot(nrow(d) == 28L, all(is.finite(d$attributable_fraction)),
          all(d$attributable_fraction >= 0 & d$attributable_fraction <= 1))
# This is a fraction within each age band; no across-age normalization or clipping.
cbh_atomic_csv(d, file.path(out, "attributable_fraction_by_age.csv"))
d$age_band <- factor(d$age_band, levels = ages)
d$prevalence <- factor(d$pfpr_pct, levels = prevalences, labels = paste0(prevalences, "%"))
p <- ggplot(d, aes(age_band, attributable_fraction, colour = prevalence, group = prevalence)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  scale_colour_manual(values = c("10%" = "#0072B2", "20%" = "#009E73",
                                 "30%" = "#E69F00", "40%" = "#CC79A7")) +
  scale_y_continuous(breaks = seq(0,1,.2), labels = scales::label_percent(accuracy = 1),
                     expand = expansion(mult = 0)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Age (completed months)", y = "Malaria-attributable share\nof all-cause deaths",
       colour = "PfPR[2–10]") +
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    legend.position = "bottom", legend.text = element_text(size = 16),
    legend.title = element_text(size = 18), plot.margin = margin(16,20,12,12))
stopifnot(is.null(p$labels$title), is.null(p$labels$subtitle), is.null(p$labels$caption))
ggsave(file.path(out, "attributable_fraction_by_age.png"), p, width = 11, height = 7,
       dpi = 240, device = ragg::agg_png, bg = "white")
inputs <- c(component_path, file.path(out, c("fit_manifest.csv", "fit_diagnostics.csv", "pfpr_curves.csv")),
            "R_cbh/reporting/06_attributable_fraction_by_age.R", "R_cbh/primary/settings.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")),
               file.path(out, "attributable_fraction_by_age_provenance.csv"))
message("Figure 3 and all 28 age/prevalence estimates saved")

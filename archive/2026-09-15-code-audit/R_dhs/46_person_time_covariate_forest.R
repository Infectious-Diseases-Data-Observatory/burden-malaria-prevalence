# =============================================================================
# 46_person_time_covariate_forest.R — covariate effects under the person-time
# model, by age band.
#
# Each per-band model in script 42 carries the same 18 standardised covariates
# in one ridge-penalised block, with the penalty chosen by fREML. This draws the
# conditional association of each covariate with the band's mortality hazard
# (% change per standard deviation on the transformed scale) for the six bands
# split at four months. Two things to keep in mind when reading it: the
# coefficients are ridge-shrunk, and the survey intercept absorbs anything that
# is constant within a survey, so the national covariates (identical across a
# survey's regions) can only be identified from between-survey variation that
# the intercept has not already taken, and are shrunk towards zero accordingly.
#
# Output
#   results/dhs_rebuild/person_time_covariate_effects.csv
#   results/dhs_rebuild/figure30_person_time_covariate_forest.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages("ggplot2")

bundle <- readRDS(file.path(DERIVED_DIR, "person_time_model_bundle.rds"))
labels <- read.csv(file.path(RESULTS_DIR, "ridge_covariate_effects.csv"),
                   stringsAsFactors = FALSE)[, c("variable", "label", "level")]

rows <- lapply(names(bundle$bands6b_smooth_fits), function(band) {
  fit <- bundle$bands6b_smooth_fits[[band]]
  b <- coef(fit); v <- diag(vcov(fit))
  idx <- grep("^G", names(b))
  data.frame(age_band = band, variable = sub("^G", "", names(b)[idx]),
             estimate = b[idx], se = sqrt(v[idx]),
             ridge_lambda = unname(fit$sp["G"]), stringsAsFactors = FALSE)
})
effects <- do.call(rbind, rows)
effects$pct <- 100 * (exp(effects$estimate) - 1)
effects$lo <- 100 * (exp(effects$estimate - 1.96 * effects$se) - 1)
effects$hi <- 100 * (exp(effects$estimate + 1.96 * effects$se) - 1)
effects <- merge(effects, labels, by = "variable", all.x = TRUE)
effects$label[is.na(effects$label)] <- effects$variable[is.na(effects$label)]
effects$label[effects$variable == "polstab"] <- "Political stability (WGI)"
effects$level[is.na(effects$level)] <- "unknown"
effects$age_band <- factor(effects$age_band, levels = AGE6B)
rownames(effects) <- NULL
write.csv(effects, file.path(RESULTS_DIR, "person_time_covariate_effects.csv"), row.names = FALSE)

# order covariates by their effect in the 4-11 month band, grouped by level
ordering <- effects[effects$age_band == AGE6B[3], ]
ordering <- ordering[order(ordering$level, ordering$pct), ]
effects$label <- factor(effects$label, levels = ordering$label)
effects$level <- factor(effects$level, levels = c("survey-region", "national-exact-year",
                                                  "national-nearest-year", "unknown"),
                        labels = c("Survey-region (DHS)", "National, exact year",
                                   "National, nearest year", "Other"))
effects$level <- droplevels(effects$level)

message("Ridge penalty (smoothing parameter for the covariate block) by band:")
print(unique(effects[, c("age_band", "ridge_lambda")]), row.names = FALSE)
message("\nLargest |effect| per band (% per SD):")
for (band in AGE6B) {
  z <- effects[effects$age_band == band, ]; z <- z[order(-abs(z$pct)), ][1:3, ]
  cat(sprintf("  %-13s %s\n", band, paste(sprintf("%s %+.1f", z$label, z$pct), collapse = "; ")))
}

plot <- ggplot2::ggplot(effects, ggplot2::aes(pct, label, colour = level)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey55") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.3, linewidth = 0.5) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_wrap(~age_band, nrow = 1) +
  ggplot2::scale_colour_manual(values = c("#1D6F8B", "#D95F0E", "#7B3294", "grey40"), name = NULL) +
  ggplot2::labs(x = "Change in mortality hazard per 1 SD of the covariate (%), ridge-shrunk",
                y = NULL,
                title = "Covariate effects under the person-time model, by age band",
                subtitle = paste("One negative-binomial model per band with country and survey intercepts;",
                                 "national covariates are constant within a survey and largely absorbed by its intercept.\n",
                                 "Bars are 95% CIs conditional on the fitted ridge penalty")) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(file.path(RESULTS_DIR, "figure30_person_time_covariate_forest.png"), plot,
                width = 13, height = 5.6, dpi = 200, bg = "white")
message("\nWrote person_time_covariate_effects.csv and figure30_person_time_covariate_forest.png")

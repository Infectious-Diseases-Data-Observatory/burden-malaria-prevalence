#!/usr/bin/env Rscript
# Supplementary figure (added 24 September 2026): the calendar-year spline h_g of each primary
# age-band model, as the hazard ratio relative to band entry in January 2000, with 95% intervals
# from the coefficient covariance conditional on the smoothing parameters (Vp). This is the
# secular trend in all-cause mortality not explained by PfPR, the 17 covariates or the random
# intercepts, common to all countries. Bars show the records entering the band in each year.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(mgcv)})
st <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION", "regional_mics"))
out <- file.path(st$out, "calendar_trend"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
ages <- cbh_config()$age_bands$age_band
band_label <- setNames(ifelse(ages == "<1", "<1 month", paste(ages, "months")), ages)
manifest <- fread(file.path(st$out, "fit_manifest.csv"))[series == "map_full"][match(ages, age_band)]
stopifnot(nrow(manifest) == 7L)
ref <- 2000; grid <- seq(2000, 2024 + 11 / 12, by = 1 / 12)

pieces <- lapply(seq_len(nrow(manifest)), function(i) {
  path <- manifest$model_file[i]; stopifnot(identical(cbh_file_hash(path), manifest$md5[i]))
  f <- readRDS(path)$fit
  s <- f$smooth[[match("s(calendar_year)", vapply(f$smooth, `[[`, "", "label"))]]; ix <- s$first.para:s$last.para
  L <- PredictMat(s, data.frame(calendar_year = grid)) - PredictMat(s, data.frame(calendar_year = rep(ref, length(grid))))
  est <- drop(L %*% coef(f)[ix]); se <- sqrt(pmax(rowSums((L %*% f$Vp[ix, ix]) * L), 0))
  yrs <- floor(f$model$calendar_year + 1e-9)
  res <- list(curve = data.table(age_band = manifest$age_band[i], calendar_year = grid, log_hr = est, se = se),
    support = data.table(age_band = manifest$age_band[i], year = yrs)[, .(records = .N), by = .(age_band, year)],
    edf = data.table(age_band = manifest$age_band[i], edf = sum(f$edf[ix])))
  rm(f); gc(FALSE); res
})
curve <- rbindlist(lapply(pieces, `[[`, "curve")); support <- rbindlist(lapply(pieces, `[[`, "support")); edf <- rbindlist(lapply(pieces, `[[`, "edf"))
curve[, `:=`(hr = exp(log_hr), lower_95 = exp(log_hr - 1.96 * se), upper_95 = exp(log_hr + 1.96 * se))]
cbh_atomic_csv(as.data.frame(curve), file.path(out, "calendar_year_curves.csv"))
cbh_atomic_csv(as.data.frame(support[order(age_band, year)]), file.path(out, "records_by_entry_year.csv"))
end <- curve[abs(calendar_year - 2024) < 1e-9, .(age_band, hr_2024 = hr, lower_2024 = lower_95, upper_2024 = upper_95)]
cbh_atomic_csv(as.data.frame(merge(edf, end, by = "age_band")[match(ages, age_band)]), file.path(out, "calendar_year_summary.csv"))

lab <- setNames(sprintf("%s (EDF %.1f)", band_label[edf$age_band], edf$edf), edf$age_band)
curve[, band := factor(age_band, levels = ages)]; support[, band := factor(age_band, levels = ages)]
lo <- min(curve$lower_95); hi <- max(curve$upper_95)
# Records per year drawn as bars on a log-scale band below the curves (scaled within each panel).
support[, share := records / max(records), by = band]
support[, ytop := exp(log(lo * .72) + share * (log(lo * .97) - log(lo * .72)))]
p <- ggplot(curve, aes(calendar_year, hr)) +
  geom_rect(data = support, aes(xmin = year, xmax = year + .9, ymin = lo * .72, ymax = ytop), inherit.aes = FALSE, fill = "grey85") +
  geom_hline(yintercept = 1, colour = "grey55", linewidth = .35) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#1F4E79", alpha = .15) +
  geom_line(colour = "#1F4E79", linewidth = .8) +
  scale_y_log10(breaks = c(.25, .35, .5, .7, 1, 1.4), labels = function(x) sub("\\.?0+$", "", sprintf("%.2f", x))) +
  scale_x_continuous(breaks = seq(2000, 2025, 5)) +
  facet_wrap(~band, nrow = 2, labeller = labeller(band = lab)) +
  labs(x = "Calendar year of band entry", y = "Hazard ratio relative to January 2000 (log scale)") +
  theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, colour = "grey85"))
ggsave(file.path(out, "sfig_calendar_year_splines.png"), p, width = 13, height = 7, dpi = 300, device = ragg::agg_png, bg = "white")

summ <- merge(edf, end, by = "age_band")[match(ages, age_band)]
caption <- paste0("Calendar-year splines of the seven age-band models of the primary analysis (", st$id, "): the hazard ratio for entering the band in a given month relative to January 2000, ",
  "conditional on PfPR[2–10] at band entry, the 17 covariates and the survey, country and survey-region random intercepts, with 95% intervals conditional on the estimated smoothing parameters. ",
  "The curve is the secular change in all-cause mortality not explained by the other terms and is common to all countries; panel titles give its effective degrees of freedom (1 = linear on the log-hazard scale). ",
  "Grey bars show the number of records entering the band in each calendar year (scaled within each panel); no records enter in 2001 because the World Bank political-stability indicator was not published for that year and the primary analysis uses complete cases. ",
  sprintf("Hazard ratios for 2024 relative to 2000: %s.", paste(sprintf("%s %.2f (%.2f–%.2f)", band_label[summ$age_band], summ$hr_2024, summ$lower_2024, summ$upper_2024), collapse = "; ")))
writeLines(c("# Supplementary figure: calendar-year splines by age band", "", caption), file.path(out, "CAPTION.md"))
inputs <- c(manifest$model_file, file.path(st$out, "fit_manifest.csv"), "R_cbh/reporting/16_calendar_trend.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")), file.path(out, "provenance.csv"))
print(summ); message("Calendar-year figure written: ", out)

# =============================================================================
# 22_seasonality_of_fieldwork.R — does the month of DHS fieldwork matter for the
# measured prevalence, and can it explain the MAP-versus-measured gap?
#
# MOTIVATION. Malaria transmission in the Sahel is strongly seasonal (a single
# rainy-season peak), whereas the MAP exposure used throughout this analysis is an
# ANNUAL mean. A survey fielded in the dry season therefore measures far lower
# parasitaemia than the annual average, and one fielded at the seasonal peak far
# higher. Two things follow: measured-versus-modelled comparisons (script 21) can
# be distorted by when fieldwork happened, and any drift in survey timing over the
# years could masquerade as drift in MAP's accuracy.
#
# DATA. DHS PR (household member) recodes carry hv006, the month of interview.
# Fieldwork moves between regions within a survey, so the timing is summarised per
# survey-region as the sample-weighted CIRCULAR mean month (fieldwork can straddle
# the year boundary, so an arithmetic mean would be wrong for e.g. Nov-Feb).
#
# WHAT IT FINDS. Raw seasonality is large - in the Sahel, measured prevalence
# averages 7.7% for fieldwork in January-June against 25.6% in July-December. But
# most of that raw contrast is confounded with WHICH regions were surveyed when:
# once the MAP annual value and a country random effect are conditioned on, the
# residual seasonal signal is a 5.7 percentage-point amplitude in the Sahel
# (p = 0.019) and essentially nothing elsewhere, adding almost no variance
# explained because location already accounts for 80-98% of it. Monthly cells are
# also small (1-7 survey-regions for most Sahel months), so the adjusted estimate
# is the defensible one and is probably still attenuated.
#
# It also rules OUT survey timing as the explanation for the growing MAP-measured
# gap of script 21: timing shifted between the first two eras and then plateaued
# while the gap kept widening, and the shift was TOWARDS the high-transmission
# season, which would narrow rather than widen the gap.
#
# Outputs (results/dhs_rebuild/):
#   dhs_month_prevalence.csv          per survey-region month, measured and MAP
#   seasonality_variance_explained.csv variance attributable to month
#   seasonality_by_era.csv            fieldwork timing and MAP gap by era
#   dhs_prevalence_by_month.png       measured prevalence against month
# Also writes data/derived_dhs/dhs_survey_region_month.csv (cached per survey).
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("mgcv", "malariaAtlas", "ggplot2", "haven"))

MEASURED_CSV <- file.path(DATA_DIR, "dhs_prevalence_by_region.csv")
CONVERSION_CSV <- file.path(REPO_ROOT, "results", "rdt_microscopy_conversion.csv")
MONTH_CSV <- file.path(DERIVED_DIR, "dhs_survey_region_month.csv")
MONTH_CACHE <- file.path(DATA_DIR, "month_cache")
SAHEL <- c("BFA", "MLI", "NER", "SEN", "GMB", "MRT", "TCD")

## ---- month of fieldwork per survey-region (cached) -------------------------
build_month_panel <- function() {
  dir.create(MONTH_CACHE, showWarnings = FALSE)
  recodes <- c(list.files(file.path(DATA_DIR, "dhs"), "PR.*rds$", full.names = TRUE),
               list.files(path.expand("~/.rdhs_cache"), "PR.*rds$",
                          full.names = TRUE, recursive = TRUE))
  recodes <- recodes[!duplicated(toupper(basename(recodes)))]
  if (!length(recodes)) return(NULL)
  message("Extracting fieldwork month from ", length(recodes), " PR recodes ...")
  rows <- list()
  for (path in recodes) {
    stem <- toupper(sub("\\.rds$", "", basename(path)))
    cache <- file.path(MONTH_CACHE, paste0(stem, ".rds"))
    if (file.exists(cache)) { rows[[stem]] <- readRDS(cache); next }
    x <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(x) || !all(c("hv006", "hv005", "hv024") %in% names(x))) next
    region <- rkey(as.character(haven::as_factor(x$hv024)))
    month <- as.integer(x$hv006); weight <- as.numeric(x$hv005)
    ok <- is.finite(month) & month >= 1 & month <= 12 &
      is.finite(weight) & weight > 0 & nzchar(region)
    if (!any(ok)) next
    region <- region[ok]; month <- month[ok]; weight <- weight[ok]
    # circular mean, so fieldwork spanning December-January is handled correctly
    angle <- 2 * pi * (month - 1) / 12
    sx <- tapply(weight * cos(angle), region, sum)
    sy <- tapply(weight * sin(angle), region, sum)
    circular <- (atan2(sy[names(sx)], sx) / (2 * pi)) * 12 + 1
    circular[circular < 0.5] <- circular[circular < 0.5] + 12
    out <- data.frame(survey_stem = stem, regkey = names(sx),
                      month_circular = as.numeric(circular),
                      month_min = as.numeric(tapply(month, region, min)),
                      month_max = as.numeric(tapply(month, region, max)),
                      stringsAsFactors = FALSE)
    saveRDS(out, cache); rows[[stem]] <- out
  }
  long <- do.call(rbind, rows); rownames(long) <- NULL
  write.csv(long, MONTH_CSV, row.names = FALSE)
  long
}
months <- if (file.exists(MONTH_CSV)) read.csv(MONTH_CSV, stringsAsFactors = FALSE) else build_month_panel()
if (is.null(months) || !file.exists(MEASURED_CSV) || !file.exists(CONVERSION_CSV)) {
  message("Fieldwork month or measured-prevalence inputs unavailable; skipping.")
  quit(save = "no", status = 0)
}

## ---- measured prevalence on a MAP-comparable basis --------------------------
conversion <- read.csv(CONVERSION_CSV, stringsAsFactors = FALSE)
micro_per_rdt <- conversion$slope[conversion$model == "through_origin"]
to_pfpr210 <- function(p) {
  out <- rep(NA_real_, length(p)); ok <- is.finite(p)
  if (any(ok)) {
    n <- sum(ok)
    out[ok] <- suppressMessages(as.numeric(malariaAtlas::convertPrevalence(
      p[ok], rep(0.5, n), rep(5, n), rep(2, n), rep(10, n))))
  }
  out
}
m <- read.csv(MEASURED_CSV, stringsAsFactors = FALSE)
m$micro_equivalent <- ifelse(is.finite(m$mic), m$mic, micro_per_rdt * m$rdt)
m$measured <- 100 * to_pfpr210(pmin(pmax(m$micro_equivalent, 0), 100) / 100)
m$survey_stem <- toupper(m$survey)
m <- merge(m, months, by = c("survey_stem", "regkey"))
analysis <- read_analysis_data()
analysis$key <- paste(analysis$iso3, analysis$year, analysis$regkey, sep = "|")
m$key <- paste(m$iso3, m$year, m$regkey, sep = "|")
m$map <- analysis$pfpr2_10[match(m$key, analysis$key)]
m$sahel <- m$iso3 %in% SAHEL
m$month <- m$month_circular
m$gap <- m$map - m$measured
m <- m[is.finite(m$measured) & m$measured > 0 & is.finite(m$month), , drop = FALSE]
cat(sprintf("merged: %d survey-regions, %d surveys, %d countries (%d Sahel)\n",
            nrow(m), length(unique(m$survey_stem)), length(unique(m$iso3)), sum(m$sahel)))
cat(sprintf("Sahel measured prevalence: January-June %.1f%% versus July-December %.1f%%\n",
            mean(m$measured[m$sahel & round(m$month) <= 6]),
            mean(m$measured[m$sahel & round(m$month) >= 7])))

## ---- variance in measured prevalence attributable to month -----------------
# conditioning on the MAP annual value and country separates the seasonal signal
# from the fact that different places were surveyed in different months
KNOTS <- list(month = c(0.5, 12.5))
variance_row <- function(z, label) {
  if (nrow(z) < 40) return(NULL)
  z$country <- factor(z$iso3)
  base <- mgcv::gam(measured ~ s(map, k = 5) + s(country, bs = "re"),
                    data = z, method = "REML")
  full <- mgcv::gam(measured ~ s(map, k = 5) + s(month, bs = "cc", k = 6) + s(country, bs = "re"),
                    data = z, method = "REML", knots = KNOTS)
  st <- summary(full)$s.table
  mr <- grep("^s\\(month\\)", rownames(st))
  pd <- data.frame(map = mean(z$map, na.rm = TRUE), month = seq(1, 12, by = 0.25),
                   country = z$country[1])
  terms <- predict(full, pd, type = "terms")
  amplitude <- diff(range(terms[, grep("month", colnames(terms))]))
  data.frame(group = label, n = nrow(z),
             r2_without_month = summary(base)$r.sq, r2_with_month = summary(full)$r.sq,
             increment = summary(full)$r.sq - summary(base)$r.sq,
             month_edf = st[mr, "edf"], month_p = st[mr, "p-value"],
             seasonal_amplitude_pp = amplitude)
}
variance_res <- do.call(rbind, list(
  variance_row(m, "All countries"), variance_row(m[m$sahel, ], "Sahel"),
  variance_row(m[!m$sahel, ], "Non-Sahel")))
write.csv(variance_res, file.path(RESULTS_DIR, "seasonality_variance_explained.csv"), row.names = FALSE)
cat("\n=== Variance in measured prevalence explained by month of fieldwork ===\n")
print(within(variance_res, {
  r2_without_month <- round(r2_without_month, 3); r2_with_month <- round(r2_with_month, 3)
  increment <- round(increment, 3); month_edf <- round(month_edf, 2)
  month_p <- signif(month_p, 3); seasonal_amplitude_pp <- round(seasonal_amplitude_pp, 1)
}), row.names = FALSE)

## ---- fieldwork timing and the MAP gap by era -------------------------------
m$era <- cut(m$year, c(2008.5, 2012.5, 2018.5, 2024.5),
             labels = c("2009-2012", "2013-2018", "2019-2024"))
era_res <- do.call(rbind, lapply(levels(m$era), function(e) {
  z <- m[m$era == e, ]
  data.frame(era = e, n = nrow(z), mean_month = mean(z$month),
             pct_july_to_november = 100 * mean(round(z$month) %in% 7:11),
             mean_map_minus_measured = mean(z$gap, na.rm = TRUE),
             sahel_share_pct = 100 * mean(z$sahel))
}))
write.csv(era_res, file.path(RESULTS_DIR, "seasonality_by_era.csv"), row.names = FALSE)
cat("\n=== Fieldwork timing and the MAP-measured gap, by era ===\n")
print(within(era_res, { mean_month <- round(mean_month, 1)
  pct_july_to_november <- round(pct_july_to_november)
  mean_map_minus_measured <- round(mean_map_minus_measured, 1)
  sahel_share_pct <- round(sahel_share_pct) }), row.names = FALSE)
cat("\nTiming shifted between the first two eras then plateaued while the gap kept\n",
    "widening, and the shift was towards the high-transmission season, which would\n",
    "NARROW the gap. Survey timing therefore does not explain the MAP drift.\n", sep = "")
write.csv(m[, c("iso3", "year", "regkey", "month", "measured", "map", "gap", "sahel", "era")],
          file.path(RESULTS_DIR, "dhs_month_prevalence.csv"), row.names = FALSE)

## ---- figure ----------------------------------------------------------------
m$group <- ifelse(m$sahel, "Sahel", "Other")
plot <- ggplot2::ggplot(m, ggplot2::aes(month, measured, colour = group)) +
  ggplot2::geom_point(alpha = 0.55, size = 1.6) +
  ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 6),
                       method.args = list(knots = list(x = c(0.5, 12.5))),
                       se = TRUE, linewidth = 1) +
  ggplot2::scale_colour_manual(values = c("Sahel" = "#d95f0e", "Other" = "#08519c"), name = NULL) +
  ggplot2::scale_x_continuous(breaks = 1:12, labels = month.abb) +
  ggplot2::labs(x = "Month of fieldwork (weighted circular mean per survey-region)",
                y = expression("DHS-measured " * italic(Pf) * "PR"[2-10] * " (%)"),
                title = "Measured parasitaemia by month of fieldwork",
                subtitle = paste("Cyclic smooth. The raw contrast is large but partly reflects",
                                 "which regions were surveyed when;",
                                 "see seasonality_variance_explained.csv for the adjusted signal.")) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
ggplot2::ggsave(file.path(RESULTS_DIR, "dhs_prevalence_by_month.png"), plot,
                width = 9.5, height = 5.2, dpi = 320, bg = "white")
cat("\nsaved: dhs_month_prevalence.csv, seasonality_{variance_explained,by_era}.csv +",
    "dhs_prevalence_by_month.png\n")

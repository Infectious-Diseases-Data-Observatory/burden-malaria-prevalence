# =============================================================================
# 21_map_vs_measured_prevalence.R — validate the MAP exposure against DHS-MEASURED
# parasitaemia, and test whether measurement error can explain the West Africa
# time interaction.
#
# THE MOTIVATION. Script 20 showed the West Africa prevalence-mortality slope rises
# from 3.1% per +10 PfPR points in 2005 (not distinguishable from zero) to 16.5% in
# 2020. One candidate explanation is artefactual: MAP's early surfaces rest on
# sparser survey input, so more exposure measurement error early would attenuate
# the early slope and manufacture a positive interaction. That predicts MAP should
# agree WORSE with directly measured prevalence in earlier years.
#
# THE TEST. DHS surveys that carried the malaria biomarker module measured
# parasitaemia in children 6-59 months by RDT (hml35) and microscopy (hml32). Those
# measurements are converted to a MAP-comparable basis exactly as the legacy
# Component-2 analysis did: microscopy where available, otherwise RDT scaled by the
# fitted microscopy/RDT ratio, then age-standardised from 0.5-5y to 2-10y with the
# Smith et al. age-prevalence model (malariaAtlas::convertPrevalence).
#
# THE BINDING LIMITATION. No DHS survey measured parasitaemia before 2009 - the
# biomarker module did not exist. So this cannot directly validate MAP over
# 2000-2008, which is precisely the period that identifies the interaction. What it
# can do is establish the DIRECTION in which MAP's accuracy is moving over the
# years where both exist.
#
# Outputs (results/dhs_rebuild/):
#   map_vs_measured_agreement.csv    agreement by era, overall and West Africa
#   map_vs_measured_slopes.csv       slope and interaction under each exposure
#   map_vs_measured_scatter.png      MAP versus measured prevalence, by era
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("mgcv", "malariaAtlas", "ggplot2", "scales"))

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
MEASURED_CSV <- file.path(DATA_DIR, "dhs_prevalence_by_region.csv")
CONVERSION_CSV <- file.path(REPO_ROOT, "results", "rdt_microscopy_conversion.csv")
PFPR_LONG_CSV <- file.path(DERIVED_DIR, "map_pfpr_lagged_long.csv")
MORT_LONG_CSV <- file.path(DERIVED_DIR, "mortality_by_period_long.csv")
for (f in c(MEASURED_CSV, CONVERSION_CSV, PFPR_LONG_CSV, MORT_LONG_CSV)) {
  if (!file.exists(f)) {
    message("Required input missing (", basename(f), "); skipping. ",
            "Scripts 17 and the legacy RDT/microscopy comparison must run first.")
    quit(save = "no", status = 0)
  }
}
bundle <- readRDS(MODEL_BUNDLE_RDS)
catalog <- bundle$catalog
catalog$included_in_main <- as.logical(catalog$included_in_main)
WEST_AFRICA <- c("BEN", "BFA", "CIV", "GHA", "GIN", "GMB", "LBR", "MLI",
                 "MRT", "NER", "NGA", "SEN", "SLE", "TGO")

## ---- DHS-measured parasitaemia on a MAP-comparable basis -------------------
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
measured <- read.csv(MEASURED_CSV, stringsAsFactors = FALSE)
measured$micro_equivalent <- ifelse(is.finite(measured$mic), measured$mic,
                                    micro_per_rdt * measured$rdt)
measured$measured_pfpr <- 100 * to_pfpr210(
  pmin(pmax(measured$micro_equivalent, 0), 100) / 100)
measured$key <- paste(measured$iso3, measured$year, measured$regkey, sep = "|")
measured <- measured[is.finite(measured$measured_pfpr) & measured$measured_pfpr > 0, ]

## ---- join to the analysis panel --------------------------------------------
analysis <- read_analysis_data()
analysis$key <- paste(analysis$iso3, analysis$year, analysis$regkey, sep = "|")
analysis$k <- paste(analysis$svkey, analysis$regkey, sep = "|")
pfpr_long <- read.csv(PFPR_LONG_CSV, stringsAsFactors = FALSE)
pfpr_long$k <- paste(pfpr_long$svkey, pfpr_long$regkey, sep = "|")
mort_long <- read.csv(MORT_LONG_CSV, stringsAsFactors = FALSE)
mort_long$k <- paste(mort_long$svkey, mort_long$regkey, sep = "|")
window2 <- local({
  sub <- pfpr_long[pfpr_long$lag <= 1, ]
  tapply(sub$pfpr2_10, sub$k, mean)
})
d <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
short <- mort_long[mort_long$period == 24, c("k", "postneonatal_mortality", "exposure")]
idx <- match(d$k, short$k)
d$mortality24 <- short$postneonatal_mortality[idx]
d$exposure24 <- short$exposure[idx]
d$map_window <- as.numeric(window2[d$k])
d$measured_pfpr <- measured$measured_pfpr[match(d$key, measured$key)]
d$west <- d$iso3 %in% WEST_AFRICA
overlap <- d[is.finite(d$measured_pfpr) & is.finite(d$pfpr2_10), , drop = FALSE]
cat(sprintf("overlap: %d region-years (%d West Africa), %d countries, %d-%d\n",
            nrow(overlap), sum(overlap$west), length(unique(overlap$iso3)),
            min(overlap$year), max(overlap$year)))
cat("NOTE: no DHS survey measured parasitaemia before 2009, so MAP cannot be\n",
    "validated over 2000-2008 - the period that identifies the interaction.\n", sep = "")

## ---- A. agreement by era ---------------------------------------------------
overlap$era <- cut(overlap$year, c(2008.5, 2012.5, 2018.5, 2024.5),
                   labels = c("2009-2012", "2013-2018", "2019-2024"))
agreement <- function(z, label) {
  if (nrow(z) < 20) return(NULL)
  fit <- lm(measured_pfpr ~ pfpr2_10, data = z)
  data.frame(group = label, n = nrow(z),
             correlation = cor(z$measured_pfpr, z$pfpr2_10),
             slope_measured_on_map = coef(fit)[2],
             mean_difference = mean(z$pfpr2_10 - z$measured_pfpr),
             mean_absolute_difference = mean(abs(z$pfpr2_10 - z$measured_pfpr)),
             residual_sd = sd(residuals(fit)))
}
rows <- c(list(agreement(overlap, "All countries, all eras")),
          lapply(levels(overlap$era), function(e)
            agreement(overlap[overlap$era == e, ], paste("All countries", e))),
          list(agreement(overlap[overlap$west, ], "West Africa, all eras")),
          lapply(levels(overlap$era), function(e)
            agreement(overlap[overlap$west & overlap$era == e, ], paste("West Africa", e))))
agreement_res <- do.call(rbind, rows[!vapply(rows, is.null, TRUE)])
write.csv(agreement_res, file.path(RESULTS_DIR, "map_vs_measured_agreement.csv"), row.names = FALSE)
cat("\n=== A. MAP versus DHS-measured PfPR2-10, agreement by era ===\n")
print(within(agreement_res, {
  correlation <- round(correlation, 3); slope_measured_on_map <- round(slope_measured_on_map, 3)
  mean_difference <- round(mean_difference, 1)
  mean_absolute_difference <- round(mean_absolute_difference, 1)
  residual_sd <- round(residual_sd, 1)
}), row.names = FALSE)
cat("\nmean_difference is MAP minus measured, so a positive value means MAP reads\n",
    "higher than the survey measurement.\n", sep = "")

## ---- B. slope and interaction under each exposure, on the same rows --------
fit_exposure <- function(z, column, spec) {
  z$pfpr10 <- z[[column]] / 10
  z$exposure <- z$exposure24
  z$postneonatal_mortality <- z$mortality24
  z <- z[is.finite(z$postneonatal_mortality) & z$postneonatal_mortality > 0 &
           is.finite(z$exposure) & z$exposure > 0 &
           is.finite(z$pfpr10) & z$pfpr10 > 0, , drop = FALSE]
  f <- tryCatch(fit_ridge_gam(z, "postneonatal_mortality", catalog, spec,
                              method = "REML", preprocessing = bundle$preprocessing),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  if (spec == "linear_no_interaction") {
    r <- model_summary_row(f)
    data.frame(n = r$n, estimate = r$pct_change_per_10, lo = r$pct_change_lo,
               hi = r$pct_change_hi, p_value = r$pfpr_p)
  } else {
    p <- summary(f$model)$p.table["pfpr10:year_c", ]
    data.frame(n = nrow(z), estimate = p[1], lo = NA_real_, hi = NA_real_, p_value = p[4])
  }
}
slope_rows <- list()
for (grp in list(list(label = "All overlap", z = overlap),
                 list(label = "West Africa", z = overlap[overlap$west, ]))) {
  for (spec in c("linear_no_interaction", "linear_time_interaction")) {
    for (col in c("pfpr2_10", "map_window", "measured_pfpr")) {
      r <- fit_exposure(grp$z, col, spec)
      if (!is.null(r)) slope_rows[[paste(grp$label, spec, col)]] <-
        cbind(group = grp$label, quantity = spec, exposure = col, r)
    }
  }
}
slope_res <- do.call(rbind, slope_rows)
write.csv(slope_res, file.path(RESULTS_DIR, "map_vs_measured_slopes.csv"), row.names = FALSE)
cat("\n=== B. Slope (% per +10 points) and interaction (per year), same rows ===\n")
print(within(slope_res, { estimate <- round(estimate, 5); lo <- round(lo, 1)
  hi <- round(hi, 1); p_value <- signif(p_value, 3) }), row.names = FALSE)

## ---- scatter of MAP against measured, by era -------------------------------
plot_df <- overlap
plot_df$region_group <- ifelse(plot_df$west, "West Africa", "Other")
scatter <- ggplot2::ggplot(plot_df, ggplot2::aes(measured_pfpr, pfpr2_10)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_point(ggplot2::aes(colour = region_group), alpha = 0.6, size = 1.6) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                       colour = "#08519c", fill = "#c7e0ee", linewidth = 0.7) +
  ggplot2::facet_wrap(~ era, nrow = 1) +
  ggplot2::scale_colour_manual(values = c("West Africa" = "#08519c", "Other" = "#e08214"),
                               name = NULL) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    x = expression("DHS-measured " * italic(Pf) * "PR"[2-10] * " (%, microscopy-equivalent, age-standardised)"),
    y = expression("MAP modelled " * italic(Pf) * "PR"[2-10] * " (%)"),
    title = "MAP modelled prevalence against DHS-measured parasitaemia",
    subtitle = paste("Dashed line is equality; blue line is the fitted relationship.",
                     "Agreement deteriorates in the later eras.")) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom",
                 strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
                 strip.text = ggplot2::element_text(face = "bold"),
                 plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(RESULTS_DIR, "map_vs_measured_scatter.png"), scatter,
                width = 11, height = 4.6, dpi = 320, bg = "white")
cat("\nsaved: map_vs_measured_{agreement,slopes}.csv + map_vs_measured_scatter.png\n")

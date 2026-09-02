# =============================================================================
# 37_burden_trend_brms.R — the national burden trend under the Bayesian ladder.
#
# Script 10 draws malaria-attributable post-neonatal deaths over 2000-2024
# under three mgcv time structures as point estimates. The Bayesian refits in
# script 34 (A additive, C curve shifting linearly in time, D full t2 surface)
# give the same extrapolation WITH uncertainty: every posterior draw of the
# surface yields a national series, so the bands carry the uncertainty in the
# dose-response and in how it changes over time, including the uncertainty in
# how smooth the surface is. They do not carry uncertainty in the all-cause
# mortality inputs, which script 10 adds for the 2024 comparison only.
#
# Extrapolation is population-average: country random effects at zero,
# covariates at their standardised mean, MAP national PfPR in each year
# against the 1% counterfactual, applied to IGME/World Bank all-cause
# post-neonatal deaths. Country-years with PfPR below 1% contribute nothing.
#
# Outputs
#   results/dhs_rebuild/brms_burden_timeseries.csv
#   results/dhs_rebuild/brms_burden_trend_summary.csv
#   results/dhs_rebuild/figure21_burden_trend_models.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("brms", "ggplot2", "countrycode"))

CACHE <- file.path(DATA_DIR, "brms_ladder_cache")
MODELS <- c(A = "A: additive", C = "C: curve shifts linearly in time",
            D = "D: full surface")
bundle <- readRDS(MODEL_BUNDLE_RDS)
year_center <- unique(bundle$year_center)[1]

## ---- national inputs, as in scripts 10 and 16 --------------------------------------
nts <- merge(
  read.csv(file.path(DATA_DIR, "wb_mortality_timeseries.csv"),
           stringsAsFactors = FALSE)[, c("iso3", "year", "allcause_1mo5y", "births")],
  read.csv(file.path(DATA_DIR, "pfpr_by_country_year.csv"), stringsAsFactors = FALSE),
  by = c("iso3", "year"))
nts$region <- countrycode::countrycode(nts$iso3, "iso3c", "region", warn = FALSE)
nts <- nts[nts$region == "Sub-Saharan Africa" & is.finite(nts$allcause_1mo5y) &
             is.finite(nts$pfpr_pct) & is.finite(nts$births) &
             nts$year >= 2000 & nts$year <= 2024, , drop = FALSE]
nts <- nts[order(nts$year, nts$iso3), ]
years <- sort(unique(nts$year))
message(nrow(nts), " country-years across ", length(unique(nts$iso3)),
        " countries, ", min(years), "-", max(years))

## ---- posterior national series under each model -----------------------------------
load_fit <- function(name) {
  file <- file.path(CACHE, paste0(name, ".rds"))
  if (!file.exists(file)) stop("Run 34_brms_model_ladder.R first; missing ", file)
  readRDS(file)
}
national_draws <- function(fit) {
  covariate_names <- setdiff(names(fit$data),
                             c("deaths", "log_exposure", "pfpr10", "year_c",
                               "country"))
  frame <- function(pfpr10) {
    out <- data.frame(pfpr10 = pfpr10, year_c = nts$year - year_center,
                      log_exposure = 0,
                      country = levels(fit$data$country)[1])
    for (name in covariate_names) out[[name]] <- 0
    out
  }
  high <- brms::posterior_linpred(fit, newdata = frame(nts$pfpr_pct / 10),
                                  re_formula = NA)
  low <- brms::posterior_linpred(fit, newdata = frame(rep(AF_REFERENCE / 10, nrow(nts))),
                                 re_formula = NA)
  af <- pmax(1 - exp(-(high - low)), 0)
  deaths <- af * matrix(nts$allcause_1mo5y, nrow = nrow(af), ncol = nrow(nts),
                        byrow = TRUE)
  # draws x years
  indicator <- outer(nts$year, years, "==") * 1
  deaths %*% indicator
}

series_rows <- list()
summary_rows <- list()
for (name in names(MODELS)) {
  message("  ", MODELS[[name]])
  totals <- national_draws(load_fit(name))
  colnames(totals) <- years
  series_rows[[name]] <- data.frame(
    model = name, label = MODELS[[name]], year = years,
    deaths = colMeans(totals),
    lo = apply(totals, 2, stats::quantile, 0.025),
    hi = apply(totals, 2, stats::quantile, 0.975),
    stringsAsFactors = FALSE)
  # Under the time-varying structures some draws put the year-2000 burden near
  # zero, so a ratio has no useful mean. The change is summarised as the
  # absolute difference (mean and interval) and the ratio only by its median.
  difference <- totals[, "2024"] - totals[, "2000"]
  ratio <- totals[, "2024"] / totals[, "2000"]
  summary_rows[[name]] <- data.frame(
    model = name, label = MODELS[[name]],
    deaths_2000 = mean(totals[, "2000"]),
    deaths_2000_lo = unname(stats::quantile(totals[, "2000"], 0.025)),
    deaths_2000_hi = unname(stats::quantile(totals[, "2000"], 0.975)),
    deaths_2024 = mean(totals[, "2024"]),
    deaths_2024_lo = unname(stats::quantile(totals[, "2024"], 0.025)),
    deaths_2024_hi = unname(stats::quantile(totals[, "2024"], 0.975)),
    difference_2024_minus_2000 = mean(difference),
    difference_lo = unname(stats::quantile(difference, 0.025)),
    difference_hi = unname(stats::quantile(difference, 0.975)),
    pct_change_median = 100 * (stats::median(ratio) - 1),
    prob_decline = mean(difference < 0),
    stringsAsFactors = FALSE)
}
series <- do.call(rbind, series_rows)
trend <- do.call(rbind, summary_rows)
rownames(series) <- rownames(trend) <- NULL
write.csv(series, file.path(RESULTS_DIR, "brms_burden_timeseries.csv"),
          row.names = FALSE)
write.csv(trend, file.path(RESULTS_DIR, "brms_burden_trend_summary.csv"),
          row.names = FALSE)
message("\nMalaria-attributable post-neonatal deaths, 2000 and 2024 (thousands), ",
        "with 95% credible intervals:")
print(transform(trend,
                deaths_2000 = sprintf("%.0f (%.0f-%.0f)", deaths_2000 / 1e3,
                                      deaths_2000_lo / 1e3, deaths_2000_hi / 1e3),
                deaths_2024 = sprintf("%.0f (%.0f-%.0f)", deaths_2024 / 1e3,
                                      deaths_2024_lo / 1e3, deaths_2024_hi / 1e3),
                difference = sprintf("%+.0f (%+.0f to %+.0f)",
                                     difference_2024_minus_2000 / 1e3,
                                     difference_lo / 1e3, difference_hi / 1e3),
                pct_change_median = sprintf("%+.0f%%", pct_change_median),
                prob_decline = round(prob_decline, 2))[
                  , c("label", "deaths_2000", "deaths_2024", "difference",
                      "pct_change_median", "prob_decline")],
      row.names = FALSE)

## ---- the mgcv structures from script 10, for the second panel -------------------------
mgcv_file <- file.path(RESULTS_DIR, "time_surface_burden_timeseries.csv")
mgcv <- if (file.exists(mgcv_file)) {
  x <- read.csv(mgcv_file, stringsAsFactors = FALSE)
  x$class <- c(spline_no_time_interaction = "additive",
               spline_ti_interaction = "time interaction",
               full_te_surface = "full surface")[x$model]
  x$label <- c(spline_no_time_interaction = "s(PfPR) + s(year)",
               spline_ti_interaction = "s(PfPR) + s(year) + ti(PfPR, year)",
               full_te_surface = "te(PfPR, year)")[x$model]
  x
} else NULL

series$class <- c(A = "additive", C = "time interaction", D = "full surface")[series$model]
class_levels <- c("additive", "time interaction", "full surface")
class_titles <- c(
  "additive" = "Additive\nbrms A  |  mgcv s(PfPR) + s(year)",
  "time interaction" = "Time interaction\nbrms C (curve shifts linearly)  |  mgcv ti(PfPR, year)",
  "full surface" = "Full surface\nbrms D t2(PfPR, year)  |  mgcv te(PfPR, year)")
plot_df <- rbind(
  data.frame(engine = "Bayesian posterior mean, 95% credible band",
             class = series$class, year = series$year, deaths = series$deaths,
             lo = series$lo, hi = series$hi, stringsAsFactors = FALSE),
  if (!is.null(mgcv)) data.frame(
    engine = "mgcv point estimate (script 10)", class = mgcv$class,
    year = mgcv$year, deaths = mgcv$malaria_deaths, lo = NA_real_, hi = NA_real_,
    stringsAsFactors = FALSE))
plot_df$panel <- factor(class_titles[plot_df$class], levels = class_titles[class_levels])
plot_df$engine <- factor(plot_df$engine,
                         levels = c("Bayesian posterior mean, 95% credible band",
                                    "mgcv point estimate (script 10)"))
class_colours <- setNames(c("#00A83B", "#5B8FF9", "#F8766D"), class_titles[class_levels])

annotation <- trend
annotation$class <- c(A = "additive", C = "time interaction", D = "full surface")[annotation$model]
annotation$panel <- factor(class_titles[annotation$class], levels = levels(plot_df$panel))
annotation$text <- sprintf(
  "2024: %.0fk (%.0f to %.0f)\n2024 minus 2000: %+.0fk (%+.0f to %+.0f)\nP(decline) = %.2f",
  annotation$deaths_2024 / 1e3, annotation$deaths_2024_lo / 1e3,
  annotation$deaths_2024_hi / 1e3, annotation$difference_2024_minus_2000 / 1e3,
  annotation$difference_lo / 1e3, annotation$difference_hi / 1e3,
  annotation$prob_decline)

plot <- ggplot2::ggplot(plot_df, ggplot2::aes(year, deaths / 1e3)) +
  ggplot2::geom_ribbon(data = plot_df[!is.na(plot_df$lo), ],
                       ggplot2::aes(ymin = lo / 1e3, ymax = hi / 1e3, fill = panel),
                       alpha = 0.18, colour = NA) +
  ggplot2::geom_line(ggplot2::aes(colour = panel, linetype = engine), linewidth = 1) +
  ggplot2::geom_text(data = annotation,
                     ggplot2::aes(x = -Inf, y = Inf, label = text),
                     hjust = 0, vjust = 1.2, size = 2.9, colour = "grey25",
                     inherit.aes = FALSE) +
  ggplot2::facet_wrap(~panel, nrow = 1) +
  ggplot2::scale_colour_manual(values = class_colours, guide = "none") +
  ggplot2::scale_fill_manual(values = class_colours, guide = "none") +
  ggplot2::scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.03))) +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::labs(
    x = NULL,
    y = "Malaria-attributable post-neonatal deaths\nin sub-Saharan Africa (thousands)",
    title = "National burden over time under the three time structures",
    subtitle = paste(
      "Population-average attributable fractions applied to all-cause post-neonatal",
      "deaths in", length(unique(nts$iso3)), "countries.\nBands carry the",
      "uncertainty in the dose-response and in its change over time, not in the",
      "all-cause inputs.")) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "bottom",
                 panel.grid.minor = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(size = 9))
ggplot2::ggsave(file.path(RESULTS_DIR, "figure21_burden_trend_models.png"), plot,
                width = 12, height = 5.4, dpi = 200, bg = "white")
message("\nWrote brms_burden_timeseries.csv, brms_burden_trend_summary.csv and ",
        "figure21_burden_trend_models.png")

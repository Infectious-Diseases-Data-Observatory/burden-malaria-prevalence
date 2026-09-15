# =============================================================================
# 16_sensitivity_random_effects.R — burden sensitivity to the country slope.
#
# The national burden extrapolation (script 10 / figure 4) uses the population-
# average model with country random effects set to zero.
#
# This script used to ask what happens if the fitted country random effects are
# INCLUDED instead. That question relied on the primary model carrying a country
# random PfPR slope: in the attributable-fraction contrast (national PfPR versus
# a 1% counterfactual, same country and same year) the random INTERCEPT cancels
# exactly, so the slope was the only random term that could move the answer.
# The primary model no longer has that slope (INCLUDE_COUNTRY_PFPR_SLOPE in
# 00_config.R), which makes the old comparison identically zero by construction
# rather than reassuringly small - a fabricated result.
#
# So the comparison is now between two MODELS rather than two predictions from
# one, which is the question that still has content:
#
#   primary_population      the reported model, country REs set to zero
#   with_slope_population   the same specification refitted WITH the country
#                           random PfPR slope, REs also set to zero
#   with_slope_re_included  that refit with each in-sample country's own random
#                           effects applied; out-of-sample countries stay at the
#                           population-average value
#
# The first pair says what re-admitting the slope does to the reported burden.
# The second pair is the original question, asked of the only model in which it
# is still meaningful.
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2", "countrycode"))

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
primary_model <- bundle$primary_fits$postneonatal$model
year_center <- unique(bundle$year_center)[1]
ref <- AF_REFERENCE

## ---- refit the same specification with the country slope re-admitted -------
analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
catalog$included_in_main <- as.logical(catalog$included_in_main)

with_slope_fit <- fit_ridge_gam(
  analysis, "postneonatal_mortality", catalog,
  specification = bundle$selected_specification,
  method = "REML", preprocessing = bundle$preprocessing,
  include_country_slope = TRUE
)
with_slope_model <- with_slope_fit$model
message("Refitted ", bundle$selected_specification,
        " with the country PfPR slope: AIC ",
        sprintf("%.1f", stats::AIC(with_slope_model)),
        " against ", sprintf("%.1f", stats::AIC(primary_model)),
        " for the primary model.")

nts <- merge(
  read.csv(file.path(DATA_DIR, "wb_mortality_timeseries.csv"),
           stringsAsFactors = FALSE)[, c("iso3", "year", "allcause_1mo5y", "births")],
  read.csv(file.path(DATA_DIR, "pfpr_by_country_year.csv"), stringsAsFactors = FALSE),
  by = c("iso3", "year"))
nts$region <- countrycode::countrycode(nts$iso3, "iso3c", "region", warn = FALSE)
nts <- nts[nts$region == "Sub-Saharan Africa" & is.finite(nts$allcause_1mo5y) &
             is.finite(nts$pfpr_pct) & is.finite(nts$births) &
             nts$year >= 2000 & nts$year <= 2024, , drop = FALSE]

lvl <- levels(primary_model$model$country)
insample <- nts$iso3 %in% lvl
cat(sprintf(
  "Country-slope burden sensitivity: %d country-years; REs available for %d of %d SSA countries.\n",
  nrow(nts), length(unique(nts$iso3[insample])), length(unique(nts$iso3))))

## ---- attributable fraction under a given model ------------------------------
newdata <- function(model, pf10, iso) {
  data.frame(pfpr10 = pf10, year_c = nts$year - year_center, exposure = 1,
             country = factor(iso, levels = lvl),
             G = I(matrix(0, nrow(nts), ncol(model$model$G),
                          dimnames = list(NULL, colnames(model$model$G)))))
}

# Both the intercept and, where present, the slope columns are zeroed for a
# population-average prediction.
random_columns <- function(model) {
  grep("^s\\(country\\)|^s\\(country,pfpr10\\)",
       colnames(predict(model, newdata(model, nts$pfpr_pct / 10, lvl[1]),
                        type = "lpmatrix")))
}

af_contrast <- function(model, iso, zero_re) {
  columns <- random_columns(model)
  Xh <- predict(model, newdata(model, nts$pfpr_pct / 10, iso), type = "lpmatrix")
  Xl <- predict(model, newdata(model, rep(ref / 10, nrow(nts)), iso),
                type = "lpmatrix")
  if (length(zero_re)) {
    Xh[zero_re, columns] <- 0
    Xl[zero_re, columns] <- 0
  }
  pmax(1 - exp(-as.numeric((Xh - Xl) %*% coef(model))), 0)
}

all_rows <- seq_len(nrow(nts))
af_primary <- af_contrast(primary_model, rep(lvl[1], nrow(nts)), all_rows)
af_with_population <- af_contrast(with_slope_model, rep(lvl[1], nrow(nts)),
                                  all_rows)
af_with_re <- af_contrast(with_slope_model,
                          ifelse(insample, nts$iso3, lvl[1]),
                          zero_re = which(!insample))

series <- list(
  primary_population = af_primary * nts$allcause_1mo5y,
  with_slope_population = af_with_population * nts$allcause_1mo5y,
  with_slope_re_included = af_with_re * nts$allcause_1mo5y
)

## ---- per-year totals --------------------------------------------------------
totals <- lapply(series, function(v) tapply(v, nts$year, sum))
summary_tab <- data.frame(year = as.integer(names(totals[[1]])))
for (name in names(totals)) summary_tab[[name]] <- as.numeric(totals[[name]])
summary_tab$pct_slope_vs_primary <- 100 *
  (summary_tab$with_slope_population / summary_tab$primary_population - 1)
summary_tab$pct_re_vs_population <- 100 *
  (summary_tab$with_slope_re_included / summary_tab$with_slope_population - 1)
write.csv(summary_tab,
          file.path(RESULTS_DIR, "sensitivity_random_effects_summary.csv"),
          row.names = FALSE)

decline <- function(v) {
  100 * (1 - v[summary_tab$year == 2024] / v[summary_tab$year == 2000])
}
for (name in names(totals)) {
  cat(sprintf("%-24s 2000=%.0fk  2024=%.0fk  decline=%.0f%%\n", name,
              summary_tab[[name]][summary_tab$year == 2000] / 1000,
              summary_tab[[name]][summary_tab$year == 2024] / 1000,
              decline(summary_tab[[name]])))
}
cat(sprintf(
  "Re-admitting the slope moves the reported burden by at most %.1f%% in a year (%+.1f%% cumulative).\n",
  max(abs(summary_tab$pct_slope_vs_primary)),
  100 * (sum(series$with_slope_population) / sum(series$primary_population) - 1)))
cat(sprintf(
  "Within the with-slope model, including country REs moves it by at most %.1f%% in a year.\n",
  max(abs(summary_tab$pct_re_vs_population))))

## ---- 2024 country shifts (in-sample only) -----------------------------------
last <- nts$year == 2024 & insample
country_shift <- data.frame(
  iso3 = nts$iso3[last],
  primary_population = series$primary_population[last],
  with_slope_population = series$with_slope_population[last],
  with_slope_re_included = series$with_slope_re_included[last])
country_shift$difference <- country_shift$with_slope_re_included -
  country_shift$primary_population
country_shift <- country_shift[order(-abs(country_shift$difference)), ]
write.csv(country_shift,
          file.path(RESULTS_DIR, "sensitivity_random_effects_country_2024.csv"),
          row.names = FALSE)
cat("Largest 2024 country shifts, primary against with-slope plus REs (thousands):\n")
print(head(transform(country_shift,
                     primary_k = round(primary_population / 1000, 1),
                     with_re_k = round(with_slope_re_included / 1000, 1),
                     delta_k = round(difference / 1000, 1))[
                       , c("iso3", "primary_k", "with_re_k", "delta_k")], 6),
      row.names = FALSE)

## ---- overlay figure ---------------------------------------------------------
labels <- c(
  primary_population = "Primary: no country slope, REs = 0 (reported)",
  with_slope_population = "Country slope re-admitted, REs = 0",
  with_slope_re_included = "Country slope re-admitted, REs included"
)
plot_df <- do.call(rbind, lapply(names(labels), function(name) {
  data.frame(year = summary_tab$year, deaths = summary_tab[[name]],
             series = labels[[name]], stringsAsFactors = FALSE)
}))
plot <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(year, deaths / 1000, colour = series)) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::scale_colour_manual(values = setNames(
    c("#08519c", "#d95f0e", "#31a354"), unname(labels)), name = NULL) +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::labs(
    x = NULL,
    y = "Malaria-attributable under-5 post-neonatal\ndeaths in sub-Saharan Africa (thousands)",
    title = "Burden sensitivity to the country-specific PfPR slope",
    subtitle = paste0(
      "The primary model carries a country random intercept only; the slope is ",
      "shown here as a structural sensitivity")) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = c(0.99, 0.98),
                 legend.justification = c(1, 1),
                 legend.background = ggplot2::element_rect(
                   fill = scales::alpha("white", 0.75), colour = NA))
ggplot2::ggsave(file.path(RESULTS_DIR, "sensitivity_random_effects_curves.png"),
                plot, width = 8, height = 5, dpi = 320, bg = "white")
cat("saved: sensitivity_random_effects_summary.csv + _country_2024.csv + _curves.png\n")

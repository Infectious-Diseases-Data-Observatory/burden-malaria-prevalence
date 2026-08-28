# =============================================================================
# 29_subgroup_fits.R — the dose-response fitted separately by era and by region.
#
# Three splits, and the formal interaction test that goes with each:
#
#   1. before 2010 versus 2010 onwards
#   2. West Africa versus Central and East/Southern Africa
#   3. the 2 x 2 of those
#
# Two things are reported for every subset, because they answer different
# questions and the first alone is unstable:
#
#   * the AIC-selected specification and the attributable fractions under it -
#     what one would report if the subset were analysed on its own;
#   * the same quantities under ONE fixed structure held common to every
#     subset. Model selection turns on margins of under two AIC points in this
#     data (see 27_horizon_lag_selection.R), so comparing subsets under whatever
#     each happens to select confounds a real subgroup difference with a
#     coin-flip between near-tied structures.
#
# A caveat on the regional split that does not apply to the era split: region is
# a country-level attribute, so the prevalence-by-region interaction is
# identified only from BETWEEN-country contrasts and inherits the confounding
# that goes with them. Era varies within country and does not have this problem.
#
# Outputs
#   results/dhs_rebuild/subgroup_fits_summary.csv
#   results/dhs_rebuild/subgroup_interaction_tests.csv
#   results/dhs_rebuild/figure11_subgroup_curves.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2"))
set.seed(20260828)

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
catalog$included_in_main <- as.logical(catalog$included_in_main)

analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]

SUBREGION_MAP <- c(
  BEN = "West", BFA = "West", CIV = "West", CPV = "West", GHA = "West",
  GIN = "West", GMB = "West", GNB = "West", LBR = "West", MLI = "West",
  MRT = "West", NER = "West", NGA = "West", SEN = "West", SLE = "West",
  TGO = "West",
  AGO = "Central & East", CAF = "Central & East", CMR = "Central & East",
  COD = "Central & East", COG = "Central & East", GAB = "Central & East",
  GNQ = "Central & East", STP = "Central & East", TCD = "Central & East",
  BDI = "Central & East", BWA = "Central & East", COM = "Central & East",
  DJI = "Central & East", ERI = "Central & East", ETH = "Central & East",
  KEN = "Central & East", LSO = "Central & East", MDG = "Central & East",
  MOZ = "Central & East", MWI = "Central & East", NAM = "Central & East",
  RWA = "Central & East", SDN = "Central & East", SOM = "Central & East",
  SSD = "Central & East", SWZ = "Central & East", TZA = "Central & East",
  UGA = "Central & East", ZAF = "Central & East", ZMB = "Central & East",
  ZWE = "Central & East")

analysis$region_group <- unname(SUBREGION_MAP[analysis$iso3])
analysis$era <- ifelse(analysis$year < 2010, "before 2010", "2010 onwards")
analysis <- analysis[!is.na(analysis$region_group), , drop = FALSE]

COMMON_SPEC <- "spline_no_interaction"
ANCHORS <- c(10, 30, 50)

## ---- one subset ------------------------------------------------------------
fit_subset <- function(data, label, split) {
  if (nrow(data) < 80L || length(unique(data$iso3)) < 5L) {
    message("  ", label, ": only ", nrow(data), " rows / ",
            length(unique(data$iso3)), " countries; skipped")
    return(NULL)
  }
  # A single era spans at most nine calendar years, which cannot support the
  # eight-basis year smooth the full panel uses.
  year_k <- max(3L, min(8L, length(unique(data$year_c)) - 1L))
  fits <- lapply(MODEL_SPECS$specification, function(specification) {
    tryCatch(
      fit_ridge_gam(data, "postneonatal_mortality", catalog,
                    specification = specification, method = "ML",
                    preprocessing = bundle$preprocessing, year_k = year_k),
      error = function(e) NULL
    )
  })
  names(fits) <- MODEL_SPECS$specification
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (!length(fits)) return(NULL)

  aic <- vapply(fits, function(f) stats::AIC(f$model), numeric(1))
  ordered <- sort(aic)
  selected <- names(ordered)[1]

  # Evaluate every attributable fraction at the subset's OWN mean year. Using a
  # common centre would extrapolate a time-interaction fit outside the years the
  # subset actually covers.
  centre <- mean(data$year_c, na.rm = TRUE)
  anchors_of <- function(fit) {
    if (is.null(fit)) return(rep(NA_real_, length(ANCHORS)))
    values <- tryCatch(af_from_model(fit$model, ANCHORS, year_c = centre),
                       error = function(e) NULL)
    if (is.null(values)) rep(NA_real_, length(ANCHORS)) else values$af
  }
  # Mirror the main pipeline: ML chooses the specification, REML fits the model
  # that gets reported. This matters more here than it does on the full panel -
  # in West Africa ML shrinks the prevalence spline to a straight line (edf 1.0)
  # while REML keeps it curved (edf 3.6), which moves the attributable fraction
  # at 10% prevalence from 4.5% to 18.7%. The subgroup does not determine the
  # curve shape, so the smoothing criterion decides it, and the edf is recorded
  # below to keep that visible.
  refit <- function(specification) {
    if (is.null(fits[[specification]])) return(NULL)
    tryCatch(
      fit_ridge_gam(data, "postneonatal_mortality", catalog,
                    specification = specification, method = "REML",
                    preprocessing = bundle$preprocessing, year_k = year_k),
      error = function(e) NULL
    )
  }
  reported_selected <- refit(selected)
  reported_common <- refit(COMMON_SPEC)

  smooth_edf <- function(fit) {
    if (is.null(fit)) return(NA_real_)
    table <- summary(fit$model)$s.table
    row <- grep("pfpr10", rownames(table))
    if (!length(row)) return(NA_real_)
    unname(table[row[1], "edf"])
  }

  af_selected <- anchors_of(reported_selected)
  af_common <- anchors_of(reported_common)

  linear <- refit("linear_no_interaction")
  effect <- rep(NA_real_, 4)
  if (!is.null(linear)) {
    coefficients <- summary(linear$model)$p.table
    row <- grep("^pfpr10$", rownames(coefficients))
    if (length(row)) {
      estimate <- coefficients[row, 1]
      se <- coefficients[row, 2]
      effect <- c(100 * (exp(estimate) - 1),
                  100 * (exp(estimate - 1.96 * se) - 1),
                  100 * (exp(estimate + 1.96 * se) - 1),
                  coefficients[row, 4])
    }
  }

  data.frame(
    split = split, subset = label,
    n = nrow(data), countries = length(unique(data$iso3)),
    surveys = length(unique(data$svkey)),
    median_year = stats::median(data$year),
    selected = selected, year_k = year_k,
    pfpr_edf_ml = smooth_edf(fits[[COMMON_SPEC]]),
    pfpr_edf_reml = smooth_edf(reported_common),
    delta_aic_to_runner_up = if (length(ordered) > 1) {
      unname(ordered[2] - ordered[1])
    } else NA_real_,
    linear_pct_change_per_10 = effect[1],
    linear_lo = effect[2], linear_hi = effect[3], linear_p = effect[4],
    af10_selected = af_selected[1], af30_selected = af_selected[2],
    af50_selected = af_selected[3],
    af10_common = af_common[1], af30_common = af_common[2],
    af50_common = af_common[3],
    stringsAsFactors = FALSE
  )
}

message("Fitting subgroups")
subsets <- list(
  list(label = "all data", split = "reference", rows = rep(TRUE, nrow(analysis))),
  list(label = "before 2010", split = "era", rows = analysis$era == "before 2010"),
  list(label = "2010 onwards", split = "era", rows = analysis$era == "2010 onwards"),
  list(label = "West", split = "region", rows = analysis$region_group == "West"),
  list(label = "Central & East", split = "region",
       rows = analysis$region_group == "Central & East")
)
for (era in c("before 2010", "2010 onwards")) {
  for (region in c("West", "Central & East")) {
    subsets[[length(subsets) + 1L]] <- list(
      label = paste0(region, ", ", era), split = "era x region",
      rows = analysis$era == era & analysis$region_group == region
    )
  }
}

summary_rows <- lapply(subsets, function(s) {
  message("  ", s$label)
  fit_subset(analysis[s$rows, , drop = FALSE], s$label, s$split)
})
summary_table <- do.call(rbind, summary_rows[!vapply(summary_rows, is.null,
                                                     logical(1))])
write.csv(summary_table,
          file.path(RESULTS_DIR, "subgroup_fits_summary.csv"), row.names = FALSE)

## ---- formal interaction tests on the pooled data ---------------------------
# Fitted directly rather than through fit_ridge_gam so the prevalence term can
# carry an interaction with a grouping variable. The rest of the structure -
# ridge block, smooth year, country random intercept, offset - is the primary
# model's.
interaction_test <- function(data, term, label) {
  data <- data[is.finite(data$postneonatal_mortality) &
                 data$postneonatal_mortality > 0 &
                 is.finite(data$exposure) & data$exposure > 0 &
                 is.finite(data$pfpr10) & is.finite(data$year_c), , drop = FALSE]
  data$country <- factor(data$iso3)
  data$deaths <- round(data$postneonatal_mortality / 1000 * data$exposure)
  ridge <- make_ridge_matrix(data, catalog, bundle$preprocessing)
  data$G <- ridge$matrix
  penalty <- list(G = list(diag(ncol(ridge$matrix))))

  formula <- stats::as.formula(sprintf(
    "deaths ~ pfpr10 * %s + G + s(year_c, k = 8) + s(country, bs = 're') + offset(log(exposure))",
    term
  ))
  model <- tryCatch(
    mgcv::gam(formula, family = mgcv::nb(), method = "ML", paraPen = penalty,
              data = data),
    error = function(e) NULL
  )
  if (is.null(model)) return(NULL)
  coefficients <- summary(model)$p.table
  rows <- grep("^pfpr10:", rownames(coefficients))
  if (!length(rows)) return(NULL)
  data.frame(
    test = label, term = rownames(coefficients)[rows],
    estimate = coefficients[rows, 1], se = coefficients[rows, 2],
    p_value = coefficients[rows, 4],
    pct_change_in_slope = 100 * (exp(coefficients[rows, 1]) - 1),
    n = nrow(data), stringsAsFactors = FALSE
  )
}

message("Interaction tests")
analysis$era_region <- factor(
  paste0(analysis$region_group, ", ", analysis$era),
  levels = c("Central & East, before 2010", "Central & East, 2010 onwards",
             "West, before 2010", "West, 2010 onwards")
)
tests <- rbind(
  interaction_test(analysis, "era", "prevalence x era"),
  interaction_test(analysis, "region_group", "prevalence x region"),
  interaction_test(analysis, "era_region", "prevalence x era-region cell")
)
rownames(tests) <- NULL
write.csv(tests, file.path(RESULTS_DIR, "subgroup_interaction_tests.csv"),
          row.names = FALSE)

message("\nSubgroup dose-response (linear no-interaction summary):")
print(summary_table[, c("split", "subset", "n", "countries", "selected",
                        "linear_pct_change_per_10", "linear_lo", "linear_hi",
                        "linear_p")], row.names = FALSE, digits = 3)
message("\nAttributable fractions under the common ", COMMON_SPEC, " structure:")
report <- summary_table[, c("split", "subset", "af10_common", "af30_common",
                            "af50_common")]
for (column in grep("^af", names(report), value = TRUE)) {
  report[[column]] <- round(100 * report[[column]], 1)
}
print(report, row.names = FALSE)
message("\nInteraction tests:")
print(tests[, c("test", "term", "pct_change_in_slope", "p_value")],
      row.names = FALSE, digits = 3)

## ---- figure: the 2 x 2 dose-response curves --------------------------------
curve_for <- function(data, label) {
  year_k <- max(3L, min(8L, length(unique(data$year_c)) - 1L))
  fit <- tryCatch(
    fit_ridge_gam(data, "postneonatal_mortality", catalog,
                  specification = COMMON_SPEC, method = "REML",
                  preprocessing = bundle$preprocessing, year_k = year_k),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  grid <- seq(1, min(60, max(data$pfpr10 * 10, na.rm = TRUE)), by = 1)
  values <- tryCatch(
    af_from_model(fit$model, grid, year_c = mean(data$year_c, na.rm = TRUE)),
    error = function(e) NULL
  )
  if (is.null(values)) return(NULL)
  values$subset <- label
  values
}

curves <- do.call(rbind, lapply(subsets[-1], function(s) {
  d <- analysis[s$rows, , drop = FALSE]
  if (nrow(d) < 80L) return(NULL)
  out <- curve_for(d, s$label)
  if (!is.null(out)) out$split <- s$split
  out
}))

plot <- ggplot2::ggplot(curves,
                        ggplot2::aes(prevalence, 100 * af, colour = subset,
                                     fill = subset)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * lo, ymax = 100 * hi),
                       alpha = 0.15, colour = NA) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::facet_wrap(~split, nrow = 1) +
  ggplot2::labs(
    x = "PfPR2-10 (%)",
    y = "Malaria-attributable fraction of\npost-neonatal mortality (%)",
    title = "Dose-response by era and region",
    subtitle = paste0("All panels use the same ", COMMON_SPEC,
                      " structure fitted by REML, evaluated at each subset's own mean year"),
    colour = NULL, fill = NULL
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(RESULTS_DIR, "figure11_subgroup_curves.png"), plot,
                width = 11, height = 4.2, dpi = 200)
message("\nWrote subgroup_fits_summary.csv, subgroup_interaction_tests.csv and ",
        "figure11_subgroup_curves.png")

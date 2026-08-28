# =============================================================================
# 07_sensitivity_model_structure.R
# Structural sensitivities on the identical primary sample. All alternatives
# retain ridge penalisation for whatever covariate block they include.
#
# Comparisons:
#   * selected primary structure
#   * selected structure without a country random PfPR slope
#   * survey-region covariates only
#   * national covariates only
#   * lower/higher spline basis dimension when the selected model is nonlinear
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2"))

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
catalog <- bundle$catalog
catalog$included_in_main <- as.logical(catalog$included_in_main)

subset_preprocessing <- function(preprocessing, variables) {
  keep <- preprocessing$variables %in% variables
  variables <- preprocessing$variables[keep]
  list(
    variables = variables,
    means = preprocessing$means[variables],
    sds = preprocessing$sds[variables],
    types = preprocessing$types[variables]
  )
}

regional_catalog <- catalog[
  catalog$included_in_main & catalog$level == "survey-region",
  ,
  drop = FALSE
]
national_catalog <- catalog[
  catalog$included_in_main & grepl("^national-", catalog$level),
  ,
  drop = FALSE
]

configurations <- list(
  primary_structure = list(
    catalog = catalog,
    preprocessing = bundle$preprocessing,
    include_country_slope = INCLUDE_COUNTRY_PFPR_SLOPE,
    spline_k = 6
  ),
  # The primary model no longer lets the prevalence slope vary by country, so
  # the informative contrast is now the other way round: what re-admitting it
  # would do.
  with_country_pfpr_slope = list(
    catalog = catalog,
    preprocessing = bundle$preprocessing,
    include_country_slope = TRUE,
    spline_k = 6
  ),
  regional_covariates_only = list(
    catalog = regional_catalog,
    preprocessing = subset_preprocessing(
      bundle$preprocessing,
      regional_catalog$variable
    ),
    include_country_slope = INCLUDE_COUNTRY_PFPR_SLOPE,
    spline_k = 6
  ),
  national_covariates_only = list(
    catalog = national_catalog,
    preprocessing = subset_preprocessing(
      bundle$preprocessing,
      national_catalog$variable
    ),
    include_country_slope = INCLUDE_COUNTRY_PFPR_SLOPE,
    spline_k = 6
  )
)

if (grepl("spline", bundle$selected_specification)) {
  configurations$spline_k4 <- list(
    catalog = catalog,
    preprocessing = bundle$preprocessing,
    include_country_slope = INCLUDE_COUNTRY_PFPR_SLOPE,
    spline_k = 4
  )
  configurations$spline_k8 <- list(
    catalog = catalog,
    preprocessing = bundle$preprocessing,
    include_country_slope = INCLUDE_COUNTRY_PFPR_SLOPE,
    spline_k = 8
  )
}

fits <- lapply(names(configurations), function(name) {
  config <- configurations[[name]]
  message("Fitting structural sensitivity: ", name)
  fit_ridge_gam(
    analysis,
    outcome = "postneonatal_mortality",
    catalog = config$catalog,
    specification = bundle$selected_specification,
    method = "ML",
    preprocessing = config$preprocessing,
    include_country_slope = config$include_country_slope,
    spline_k = config$spline_k
  )
})
names(fits) <- names(configurations)

summary_rows <- do.call(rbind, lapply(names(fits), function(name) {
  fit <- fits[[name]]
  row <- model_summary_row(fit)
  anchors <- af_from_model(
    fit$model,
    prevalence = c(10, 30, 50),
    year_c = 0
  )
  row$sensitivity <- name
  row$ridge_covariates <- ncol(fit$model$model$G)
  row$country_random_slope <- fit$include_country_slope
  row$af10 <- anchors$af[anchors$prevalence == 10]
  row$af30 <- anchors$af[anchors$prevalence == 30]
  row$af50 <- anchors$af[anchors$prevalence == 50]
  row
}))
summary_rows$dAIC <- summary_rows$AIC - min(summary_rows$AIC)
summary_rows <- summary_rows[, c(
  "sensitivity", "specification", "n", "countries",
  "ridge_covariates", "country_random_slope", "AIC", "dAIC",
  "pfpr_smooth_edf", "pct_change_per_10", "pfpr_p",
  "af10", "af30", "af50"
)]

write.csv(
  summary_rows,
  file.path(RESULTS_DIR, "sensitivity_structure_summary.csv"),
  row.names = FALSE
)
saveRDS(fits, file.path(DERIVED_DIR, "sensitivity_structure_models.rds"))

summary_rows$sensitivity <- factor(
  summary_rows$sensitivity,
  levels = rev(summary_rows$sensitivity[order(summary_rows$dAIC)])
)
plot <- ggplot2::ggplot(
  summary_rows,
  ggplot2::aes(dAIC, sensitivity)
) +
  ggplot2::geom_col(fill = "#41ab5d", width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f", dAIC)),
    hjust = -0.15
  ) +
  ggplot2::expand_limits(x = max(summary_rows$dAIC) * 1.12 + 1) +
  ggplot2::labs(
    x = expression(Delta * "AIC within the structural sensitivity set"),
    y = NULL
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(
  file.path(RESULTS_DIR, "sensitivity_structure_comparison.png"),
  plot,
  width = 8.5,
  height = 5.2,
  dpi = 300
)

message("Structural sensitivity results:")
print(summary_rows, row.names = FALSE)

# =============================================================================
# 08_sensitivity_likelihood.R
# Test sensitivity to the negative-binomial count likelihood by modelling the
# log mortality rate with Gaussian errors. The alternative retains the primary
# sample, ridge-penalised covariate block, PfPR/time structure, smooth calendar
# year, and country random intercept.
#
# This is the rebuilt equivalent of the manuscript's log-normal mixed-model
# sensitivity while changing only the outcome likelihood/scale.
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2"))

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
catalog <- bundle$catalog
catalog$included_in_main <- as.logical(catalog$included_in_main)

fit_log_gaussian_gam <- function(
    data,
    outcome,
    specification,
    preprocessing) {
  keep <- is.finite(data[[outcome]]) & data[[outcome]] > 0 &
    is.finite(data$pfpr10) & is.finite(data$year_c) &
    !is.na(data$iso3)
  dd <- data[keep, , drop = FALSE]
  dd$country <- factor(dd$iso3)
  dd$log_rate <- log(dd[[outcome]])

  ridge <- make_ridge_matrix(dd, catalog, preprocessing)
  dd$G <- ridge$matrix
  penalty <- list(G = list(diag(ncol(ridge$matrix))))

  formula <- model_formula(specification)
  formula <- update(formula, log_rate ~ . - offset(log(exposure)))
  model <- mgcv::gam(
    formula,
    family = gaussian(),
    method = "REML",
    paraPen = penalty,
    data = dd
  )
  list(
    model = model,
    data = dd,
    preprocessing = ridge$preprocessing,
    specification = specification,
    outcome = outcome,
    include_country_slope = INCLUDE_COUNTRY_PFPR_SLOPE,
    spline_k = 6,
    year_k = 8
  )
}

outcomes <- c(
  postneonatal = "postneonatal_mortality",
  all_under_5 = "u5mr",
  neonatal = "nnmr"
)

selected_fits <- setNames(
  lapply(outcomes, function(outcome) {
    fit_log_gaussian_gam(
      analysis,
      outcome,
      bundle$selected_specification,
      bundle$preprocessing
    )
  }),
  names(outcomes)
)
linear_fits <- setNames(
  lapply(outcomes, function(outcome) {
    fit_log_gaussian_gam(
      analysis,
      outcome,
      "linear_no_interaction",
      bundle$preprocessing
    )
  }),
  names(outcomes)
)

summary_rows <- do.call(rbind, lapply(names(outcomes), function(name) {
  selected <- model_summary_row(selected_fits[[name]])
  linear <- model_summary_row(linear_fits[[name]])
  anchors <- af_from_model(
    selected_fits[[name]]$model,
    prevalence = c(10, 30, 50),
    year_c = 0
  )
  data.frame(
    outcome = name,
    n = selected$n,
    countries = selected$countries,
    selected_specification = bundle$selected_specification,
    selected_pfpr_smooth_edf = selected$pfpr_smooth_edf,
    selected_pfpr_p = selected$pfpr_p,
    selected_time_interaction_edf = selected$time_interaction_edf,
    selected_time_interaction_p = selected$time_interaction_p,
    linear_pct_change_per_10 = linear$pct_change_per_10,
    linear_pct_change_lo = linear$pct_change_lo,
    linear_pct_change_hi = linear$pct_change_hi,
    linear_pfpr_p = linear$pfpr_p,
    af10 = anchors$af[anchors$prevalence == 10],
    af30 = anchors$af[anchors$prevalence == 30],
    af50 = anchors$af[anchors$prevalence == 50]
  )
}))
write.csv(
  summary_rows,
  file.path(RESULTS_DIR, "sensitivity_likelihood_summary.csv"),
  row.names = FALSE
)
saveRDS(
  list(selected = selected_fits, linear = linear_fits),
  file.path(DERIVED_DIR, "sensitivity_likelihood_models.rds")
)

prevalence <- exp(seq(
  log(min(analysis$pfpr2_10, na.rm = TRUE)),
  log(max(analysis$pfpr2_10, na.rm = TRUE)),
  length.out = 220
))
nb_model <- bundle$primary_fits$postneonatal$model
gaussian_model <- selected_fits$postneonatal$model

curves <- do.call(rbind, lapply(
  list(
    "Negative binomial count GAM" = nb_model,
    "Log-Gaussian rate GAM" = gaussian_model
  ),
  function(model) {
    af <- af_from_model(model, prevalence, year_c = 0)
    data.frame(
      prevalence = prevalence,
      attributable_fraction = 100 * af$af,
      lo = 100 * af$lo,
      hi = 100 * af$hi
    )
  }
))
curves$model <- rep(
  c("Negative binomial count GAM", "Log-Gaussian rate GAM"),
  each = length(prevalence)
)

plot <- ggplot2::ggplot(
  curves,
  ggplot2::aes(prevalence, attributable_fraction, colour = model, fill = model)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lo, ymax = hi),
    alpha = 0.12,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_x_log10(
    breaks = c(1, 2, 5, 10, 20, 50, 80)
  ) +
  ggplot2::labs(
    x = expression(italic(Pf) * "PR"[2-10] * " (%), log scale"),
    y = "Post-neonatal attributable fraction (%, versus 1% PfPR)",
    colour = NULL,
    fill = NULL,
    subtitle = paste(
      "Same sample, ridge covariates, PfPR/time terms and country effects;",
      "likelihood differs"
    )
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(RESULTS_DIR, "sensitivity_likelihood_curves.png"),
  plot,
  width = 9,
  height = 5.8,
  dpi = 300
)

message("Likelihood sensitivity results:")
print(summary_rows, row.names = FALSE)

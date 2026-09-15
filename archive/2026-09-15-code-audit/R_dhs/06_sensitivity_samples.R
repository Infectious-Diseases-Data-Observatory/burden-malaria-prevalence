# =============================================================================
# 06_sensitivity_samples.R
# Refit the selected ridge-GAM specification under alternative sample rules:
#   * primary regional PfPR >=1% sample
#   * include sub-1% regions from otherwise eligible countries
#   * restrict to PfPR 5-40%
#   * complete cases for all covariates eligible for the primary model
#
# These comparisons vary the sample only; model structure and ridge
# preprocessing are held fixed.
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2"))

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
analysis <- read_analysis_data()
catalog <- bundle$catalog
catalog$included_in_main <- as.logical(catalog$included_in_main)

positive_outcomes <- is.finite(analysis$postneonatal_mortality) &
  analysis$postneonatal_mortality > 0 &
  is.finite(analysis$nnmr) & analysis$nnmr > 0 &
  is.finite(analysis$u5mr) & analysis$u5mr > 0 &
  is.finite(analysis$exposure) & analysis$exposure > 0

sample_flags <- list(
  primary_pfpr_ge_1 = as.logical(analysis$main_sample),
  include_sub_1_percent = analysis$country_mean_pfpr_gt_1 & positive_outcomes,
  restricted_pfpr_5_to_40 = as.logical(analysis$main_sample) &
    analysis$pfpr2_10 >= 5 & analysis$pfpr2_10 <= 40,
  complete_case_covariates = as.logical(analysis$main_sample) &
    as.logical(analysis$complete_case_eligible)
)

sensitivity_fits <- lapply(names(sample_flags), function(sample_name) {
  dd <- analysis[sample_flags[[sample_name]], , drop = FALSE]
  message("Fitting sample sensitivity ", sample_name, " (n=", nrow(dd), ").")
  fit_ridge_gam(
    dd,
    outcome = "postneonatal_mortality",
    catalog = catalog,
    specification = bundle$selected_specification,
    method = "REML",
    preprocessing = bundle$preprocessing
  )
})
names(sensitivity_fits) <- names(sample_flags)

summary_rows <- do.call(rbind, lapply(
  names(sensitivity_fits),
  function(sample_name) {
    fit <- sensitivity_fits[[sample_name]]
    row <- model_summary_row(fit)
    anchors <- af_from_model(
      fit$model,
      prevalence = c(10, 30, 50),
      year_c = 0
    )
    row$sample <- sample_name
    row$af10 <- anchors$af[anchors$prevalence == 10]
    row$af30 <- anchors$af[anchors$prevalence == 30]
    row$af50 <- anchors$af[anchors$prevalence == 50]
    row
  }
))
summary_rows <- summary_rows[, c(
  "sample", "specification", "n", "countries", "AIC",
  "pfpr_smooth_edf", "pct_change_per_10", "pct_change_lo",
  "pct_change_hi", "pfpr_p", "af10", "af30", "af50"
)]
write.csv(
  summary_rows,
  file.path(RESULTS_DIR, "sensitivity_sample_summary.csv"),
  row.names = FALSE
)
saveRDS(
  sensitivity_fits,
  file.path(DERIVED_DIR, "sensitivity_sample_models.rds")
)

prevalence_grid <- exp(seq(
  log(max(0.1, min(analysis$pfpr2_10, na.rm = TRUE))),
  log(max(analysis$pfpr2_10, na.rm = TRUE)),
  length.out = 220
))
curves <- do.call(rbind, lapply(
  names(sensitivity_fits),
  function(sample_name) {
    model <- sensitivity_fits[[sample_name]]$model
    high <- newdata_at_mean(model, prevalence_grid / 10, year_c = 0)
    low <- newdata_at_mean(
      model,
      rep(AF_REFERENCE / 10, length(prevalence_grid)),
      year_c = 0
    )
    delta <- link_difference(model, high, low)
    data.frame(
      sample = sample_name,
      prevalence = prevalence_grid,
      relative_change = 100 * (exp(delta$fit) - 1)
    )
  }
))

plot <- ggplot2::ggplot(
  curves,
  ggplot2::aes(prevalence, relative_change, colour = sample)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dotted",
    colour = "grey55"
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_x_log10(
    breaks = c(0.2, 0.5, 1, 2, 5, 10, 20, 50, 80)
  ) +
  ggplot2::labs(
    x = expression(italic(Pf) * "PR"[2-10] * " (%), log scale"),
    y = "Post-neonatal mortality change versus 1% PfPR (%)",
    colour = "Analysis sample"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(RESULTS_DIR, "sensitivity_sample_curves.png"),
  plot,
  width = 9,
  height = 5.8,
  dpi = 300
)

message("Sample sensitivity results:")
print(summary_rows, row.names = FALSE)

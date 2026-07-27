# =============================================================================
# 05_make_main_plots.R
# Create figures only from the saved model bundle and analysis dataset.
# No models are refitted here and no files are copied outside the repository.
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2", "patchwork"))

if (!file.exists(MODEL_BUNDLE_RDS)) {
  stop("Main model bundle missing. Run script 04.")
}
bundle <- readRDS(MODEL_BUNDLE_RDS)
analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]

post_model <- bundle$primary_fits$postneonatal$model
grid <- data.frame(
  prevalence = exp(seq(
    log(max(PFPR_FLOOR, min(analysis$pfpr2_10, na.rm = TRUE))),
    log(max(analysis$pfpr2_10, na.rm = TRUE)),
    length.out = 240
  ))
)

post_prediction <- link_prediction(
  post_model,
  newdata_at_mean(post_model, grid$prevalence / 10, year_c = 0)
)
grid$rate <- 1000 * exp(post_prediction$fit)
grid$rate_lo <- 1000 * exp(post_prediction$fit - 1.96 * post_prediction$se)
grid$rate_hi <- 1000 * exp(post_prediction$fit + 1.96 * post_prediction$se)

af <- af_from_model(post_model, grid$prevalence, year_c = 0)
grid$af <- 100 * af$af
grid$af_lo <- 100 * af$lo
grid$af_hi <- 100 * af$hi

minor_grid <- ggplot2::element_line(colour = "grey92", linewidth = 0.3)
major_grid <- ggplot2::element_line(colour = "grey85", linewidth = 0.4)
x_scale <- ggplot2::scale_x_log10(
  breaks = c(1, 2, 5, 10, 20, 50, 80),
  expand = ggplot2::expansion(mult = c(0.02, 0.03))
)

panel_a <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = analysis,
    ggplot2::aes(pfpr2_10, postneonatal_mortality, size = exposure),
    colour = "grey45",
    alpha = 0.35
  ) +
  ggplot2::geom_ribbon(
    data = grid,
    ggplot2::aes(prevalence, ymin = rate_lo, ymax = rate_hi),
    fill = "#08519c",
    alpha = 0.18
  ) +
  ggplot2::geom_line(
    data = grid,
    ggplot2::aes(prevalence, rate),
    colour = "#08519c",
    linewidth = 1.1
  ) +
  ggplot2::scale_size_area(
    max_size = 6,
    name = "DHS weighted\nbirth exposure",
    breaks = c(500, 2000, 5000),
    labels = scales::comma
  ) +
  x_scale +
  ggplot2::scale_y_log10(breaks = c(10, 20, 50, 100, 200)) +
  ggplot2::coord_cartesian(ylim = c(10, 200)) +
  ggplot2::labs(
    x = expression(italic(Pf) * "PR"[2-10] * " (%), MAP, log scale"),
    y = "All-cause post-neonatal mortality\n(per 1,000 live births, log scale)"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = minor_grid,
    panel.grid.major = major_grid,
    legend.position = c(0.99, 0.02),
    legend.justification = c(1, 0),
    legend.background = ggplot2::element_rect(
      fill = scales::alpha("white", 0.7),
      colour = NA
    )
  )

panel_b <- ggplot2::ggplot(
  grid,
  ggplot2::aes(prevalence, af)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = af_lo, ymax = af_hi),
    fill = "#d73027",
    alpha = 0.16
  ) +
  ggplot2::geom_line(colour = "#d73027", linewidth = 1.1) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dotted",
    colour = "grey55"
  ) +
  x_scale +
  ggplot2::labs(
    x = expression(italic(Pf) * "PR"[2-10] * " (%), MAP, log scale"),
    y = "Malaria-attributable fraction of post-neonatal\ndeaths (%, versus 1% PfPR)"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = minor_grid,
    panel.grid.major = major_grid
  )

figure_1 <- panel_a + panel_b +
  patchwork::plot_annotation(
    tag_levels = "A",
    subtitle = sprintf(
      "Selected ridge-GAM: %s; population curves shown at reference year %d",
      bundle$selected_specification,
      unique(bundle$year_center)[1]
    )
  )
ggplot2::ggsave(
  file.path(RESULTS_DIR, "figure1_pfpr_mortality_and_af.png"),
  figure_1,
  width = 13,
  height = 5.8,
  dpi = 320
)

comparison <- bundle$comparison
comparison$specification <- factor(
  comparison$specification,
  levels = rev(comparison$specification[order(comparison$dAIC)])
)
figure_comparison <- ggplot2::ggplot(
  comparison,
  ggplot2::aes(dAIC, specification)
) +
  ggplot2::geom_col(fill = "#2c7fb8", width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f", dAIC)),
    hjust = -0.15,
    size = 3.5
  ) +
  ggplot2::expand_limits(x = max(comparison$dAIC) * 1.12 + 1) +
  ggplot2::labs(
    x = expression(Delta * "AIC relative to the best model"),
    y = NULL,
    title = "Ridge-penalised GAM model comparison"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(
  file.path(RESULTS_DIR, "figure2_model_comparison.png"),
  figure_comparison,
  width = 8,
  height = 4.8,
  dpi = 300
)

outcome_labels <- c(
  postneonatal = "Post-neonatal (1 month–5 years)",
  all_under_5 = "All under-5",
  neonatal = "Neonatal negative control"
)
relative_curves <- do.call(rbind, lapply(
  names(bundle$primary_fits),
  function(outcome) {
    model <- bundle$primary_fits[[outcome]]$model
    high <- newdata_at_mean(model, grid$prevalence / 10, year_c = 0)
    reference <- newdata_at_mean(
      model,
      rep(AF_REFERENCE / 10, nrow(grid)),
      year_c = 0
    )
    delta <- link_difference(model, high, reference)
    data.frame(
      prevalence = grid$prevalence,
      outcome = outcome_labels[[outcome]],
      relative_change = 100 * (exp(delta$fit) - 1),
      lo = 100 * (exp(delta$fit - 1.96 * delta$se) - 1),
      hi = 100 * (exp(delta$fit + 1.96 * delta$se) - 1)
    )
  }
))

figure_negative_control <- ggplot2::ggplot(
  relative_curves,
  ggplot2::aes(prevalence, relative_change, colour = outcome, fill = outcome)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dotted",
    colour = "grey55"
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lo, ymax = hi),
    alpha = 0.12,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 1.05) +
  x_scale +
  ggplot2::scale_colour_manual(
    values = c(
      "Post-neonatal (1 month–5 years)" = "#08519c",
      "All under-5" = "#41ab5d",
      "Neonatal negative control" = "#d73027"
    ),
    name = NULL
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Post-neonatal (1 month–5 years)" = "#08519c",
      "All under-5" = "#41ab5d",
      "Neonatal negative control" = "#d73027"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    x = expression(italic(Pf) * "PR"[2-10] * " (%), MAP, log scale"),
    y = "Modelled change in mortality versus 1% PfPR (%)",
    subtitle = sprintf(
      "Population curves at reference year %d",
      unique(bundle$year_center)[1]
    )
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(RESULTS_DIR, "figure3_negative_control_curves.png"),
  figure_negative_control,
  width = 9,
  height = 5.8,
  dpi = 300
)

linear_effects <- read.csv(
  file.path(RESULTS_DIR, "linear_effects_by_outcome.csv"),
  stringsAsFactors = FALSE
)
outcome_column_labels <- c(
  postneonatal_mortality = "Post-neonatal (1 month–5 years)",
  u5mr = "All under-5",
  nnmr = "Neonatal negative control"
)
linear_effects$label <- outcome_column_labels[linear_effects$outcome]
linear_effects$label <- factor(
  linear_effects$label,
  levels = rev(unname(outcome_column_labels))
)
figure_linear <- ggplot2::ggplot(
  linear_effects,
  ggplot2::aes(pct_change_per_10, label)
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey55"
  ) +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = pct_change_lo, xmax = pct_change_hi)
  ) +
  ggplot2::geom_point(size = 3, colour = "#08519c") +
  ggplot2::labs(
    x = "% change in mortality per +10 PfPR points\n(ridge-linear comparison; 95% CI)",
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(
  file.path(RESULTS_DIR, "figure4_linear_outcome_comparison.png"),
  figure_linear,
  width = 9,
  height = 4.2,
  dpi = 300
)

write.csv(
  grid,
  file.path(RESULTS_DIR, "main_postneonatal_prediction_curve.csv"),
  row.names = FALSE
)
write.csv(
  relative_curves,
  file.path(RESULTS_DIR, "main_outcome_relative_curves.csv"),
  row.names = FALSE
)

message("Saved four main figures and their prediction tables under ", RESULTS_DIR)

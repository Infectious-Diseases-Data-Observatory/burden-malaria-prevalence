# =============================================================================
# 14_manuscript_figures.R — main-text manuscript figures for main.tex.
#
#   Figure 2 (4 panels): PfPR2-10 vs mortality + attributable-fraction share,
#     for post-neonatal (A, C) and neonatal negative control (B, D).
#   Figure 4: national malaria-attributable post-neonatal deaths 2000-2024 under
#     the selected model (with a simulated 95% band) versus IHME/GBD and WHO.
#
# Figures use only the saved model bundle and the national burden inputs; no
# models are refitted. Run after scripts 04 and 10.
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2", "patchwork", "countrycode", "MASS", "scales"))
set.seed(20260727)

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
post_model <- bundle$primary_fits$postneonatal$model
neo_model  <- bundle$primary_fits$neonatal$model
year_center <- unique(bundle$year_center)[1]

## ============================================================================
## Figure 2: association (A/B scatter + fit) and attributable share (C/D)
## ============================================================================
prev_grid <- exp(seq(
  log(max(PFPR_FLOOR, min(analysis$pfpr2_10, na.rm = TRUE))),
  log(max(analysis$pfpr2_10, na.rm = TRUE)),
  length.out = 240
))
x_scale <- ggplot2::scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 80))

rate_curve <- function(model) {
  pr <- link_prediction(model, newdata_at_mean(model, prev_grid / 10, year_c = 0))
  data.frame(prevalence = prev_grid,
             rate = 1000 * exp(pr$fit),
             lo = 1000 * exp(pr$fit - 1.96 * pr$se),
             hi = 1000 * exp(pr$fit + 1.96 * pr$se))
}
scatter_panel <- function(model, yvar, ylab, family_col, y_lower) {
  g <- rate_curve(model)
  y_upper <- max(c(analysis[[yvar]], g$hi), na.rm = TRUE) * 1.02
  ggplot2::ggplot() +
    ggplot2::geom_point(data = analysis,
      ggplot2::aes(pfpr2_10, .data[[yvar]], size = exposure),
      colour = "grey45", alpha = 0.35) +
    ggplot2::geom_ribbon(data = g, ggplot2::aes(prevalence, ymin = lo, ymax = hi),
      fill = family_col, alpha = 0.18) +
    ggplot2::geom_line(data = g, ggplot2::aes(prevalence, rate),
      colour = family_col, linewidth = 1.1) +
    ggplot2::scale_size_area(max_size = 5.5, name = "Births", breaks = c(500, 2000, 5000),
      labels = scales::comma) +
    x_scale + ggplot2::scale_y_log10() +
    ggplot2::coord_cartesian(ylim = c(y_lower, y_upper)) +
    ggplot2::labs(x = expression(italic(Pf) * "PR"[2-10] * " (%), log scale"), y = ylab) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
      legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.7), colour = NA))
}
af_panel <- function(model, ylab, family_col) {
  af <- af_from_model(model, prev_grid, year_c = 0)
  ggplot2::ggplot(data.frame(prevalence = prev_grid, af = 100 * af$af,
                             lo = 100 * af$lo, hi = 100 * af$hi),
    ggplot2::aes(prevalence, af)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = family_col, alpha = 0.16) +
    ggplot2::geom_line(colour = family_col, linewidth = 1.1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
    x_scale +
    ggplot2::labs(x = expression(italic(Pf) * "PR"[2-10] * " (%), log scale"), y = ylab) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
POST <- "#08519c"; NEO <- "#d73027"
pA <- scatter_panel(post_model, "postneonatal_mortality",
  "All-cause post-neonatal mortality\n(per 1000 live births, log scale)", POST, 10)
pB <- scatter_panel(neo_model, "nnmr",
  "All-cause neonatal mortality\n(per 1000 live births, log scale)", NEO, 5)
pC <- af_panel(post_model, "Malaria-attributable share of\npost-neonatal deaths (%, vs 1% PfPR)", POST)
pD <- af_panel(neo_model, "Malaria-attributable share of\nneonatal deaths (%, vs 1% PfPR)", NEO)
fig2 <- (pA + pB) / (pC + pD) + patchwork::plot_annotation(tag_levels = "A")
ggplot2::ggsave(file.path(RESULTS_DIR, "figure2_association_fourpanel.png"),
  fig2, width = 11, height = 9, dpi = 320)
cat("saved: figure2_association_fourpanel.png\n")

## ============================================================================
## Figure 4: national burden 2000-2024, selected model (+95% band) vs IHME/WHO
## ============================================================================
nts <- merge(
  read.csv(file.path(DATA_DIR, "wb_mortality_timeseries.csv"),
           stringsAsFactors = FALSE)[, c("iso3", "year", "allcause_1mo5y")],
  read.csv(file.path(DATA_DIR, "pfpr_by_country_year.csv"), stringsAsFactors = FALSE),
  by = c("iso3", "year"))
nts$region <- countrycode::countrycode(nts$iso3, "iso3c", "region", warn = FALSE)
nts <- nts[nts$region == "Sub-Saharan Africa" &
             is.finite(nts$allcause_1mo5y) & is.finite(nts$pfpr_pct) &
             nts$year >= 2000 & nts$year <= 2024, , drop = FALSE]

# per-year total malaria-attributable post-neonatal deaths + simulated 95% band
high <- newdata_at_mean(post_model, nts$pfpr_pct / 10, year_c = nts$year - year_center)
low  <- newdata_at_mean(post_model, rep(AF_REFERENCE / 10, nrow(nts)),
                        year_c = nts$year - year_center)
dX <- population_lpmatrix(post_model, high) - population_lpmatrix(post_model, low)
af_point <- pmax(1 - exp(-as.numeric(dX %*% coef(post_model))), 0)
draws <- MASS::mvrnorm(4000, mu = coef(post_model), Sigma = vcov(post_model))
af_sims <- pmax(1 - exp(-(dX %*% t(draws))), 0)          # nrow x sims
deaths_point <- af_point * nts$allcause_1mo5y
deaths_sims  <- af_sims * nts$allcause_1mo5y
ours <- do.call(rbind, lapply(sort(unique(nts$year)), function(y) {
  idx <- nts$year == y
  totals <- colSums(deaths_sims[idx, , drop = FALSE])
  data.frame(year = y, source = "Our model (prevalence)",
             deaths = sum(deaths_point[idx]),
             lo = unname(quantile(totals, 0.025)), hi = unname(quantile(totals, 0.975)))
}))

ihme <- read.csv(file.path(REPO_ROOT, "results", "ihme_u5_deaths_ssa_timeseries.csv"),
                 stringsAsFactors = FALSE)
ihme <- data.frame(year = ihme$year, source = "IHME/GBD (under-5)",
                   deaths = ihme$point, lo = ihme$lo, hi = ihme$hi)
ihme <- ihme[ihme$year >= 2000 & ihme$year <= 2024, ]
who <- read.csv(file.path(REPO_ROOT, "results", "who_wmr2025_africa_deaths.csv"),
                stringsAsFactors = FALSE)
who <- data.frame(year = who$year, source = "WHO (under-5 proxy, 0.75×all-age)",
                  deaths = 0.75 * who$point, lo = 0.75 * who$lo, hi = 0.75 * who$hi)
who <- who[who$year >= 2000 & who$year <= 2024, ]
burden <- rbind(ours, ihme, who)
write.csv(burden, file.path(RESULTS_DIR, "figure4_burden_timeseries.csv"), row.names = FALSE)

cols <- c("Our model (prevalence)" = "#08519c", "IHME/GBD (under-5)" = "#e6550d",
          "WHO (under-5 proxy, 0.75×all-age)" = "#31a354")
fig4 <- ggplot2::ggplot(burden, ggplot2::aes(year, deaths / 1000, colour = source, fill = source)) +
  ggplot2::geom_ribbon(data = burden[burden$source == "Our model (prevalence)", ],
    ggplot2::aes(ymin = lo / 1000, ymax = hi / 1000), alpha = 0.18, colour = NA) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::scale_colour_manual(values = cols, name = NULL) +
  ggplot2::scale_fill_manual(values = cols, name = NULL, guide = "none") +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::labs(x = NULL,
    y = "Malaria-attributable under-5 post-neonatal deaths\nin sub-Saharan Africa (thousands)") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
    legend.position = c(0.99, 0.98), legend.justification = c(1, 1),
    legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.75), colour = NA))
ggplot2::ggsave(file.path(RESULTS_DIR, "figure4_burden_vs_ihme_who.png"),
  fig4, width = 9, height = 6, dpi = 320)
cat("saved: figure4_burden_vs_ihme_who.png\n")

decl <- function(s) { z <- burden[burden$source == s, ]; z <- z[order(z$year), ]
  100 * (1 - z$deaths[z$year == 2024] / z$deaths[z$year == 2000]) }
cat(sprintf("2000->2024 decline: ours %.0f%% (%.0fk->%.0fk); IHME %.0f%%; WHO %.0f%%\n",
  decl("Our model (prevalence)"),
  ours$deaths[ours$year == 2000] / 1000, ours$deaths[ours$year == 2024] / 1000,
  decl("IHME/GBD (under-5)"), decl("WHO (under-5 proxy, 0.75×all-age)")))

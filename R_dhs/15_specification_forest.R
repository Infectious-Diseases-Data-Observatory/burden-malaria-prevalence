# =============================================================================
# 15_specification_forest.R — robustness of the PfPR2-10 association.
#
# Supplementary figure supporting the statement that the post-neonatal
# PfPR2-10 association is essentially unchanged whether we adjust for the
# survey-region covariates alone or the full covariate set, and is robust
# across alternative samples, model structure and a log-Gaussian likelihood.
#
# For every specification we refit the LINEAR ridge summary (pfpr10 as a linear
# term, i.e. the headline "% change per +10 PfPR points" estimate) so all rows
# are directly comparable on one scale, and plot a forest with the primary
# (full-covariate) estimate as a reference band.
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2"))

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
analysis <- read_analysis_data()
catalog <- bundle$catalog
catalog$included_in_main <- as.logical(catalog$included_in_main)

OUT <- "postneonatal_mortality"
LIN <- "linear_no_interaction"

## ---- covariate-block catalogs + preprocessing (as in script 07) ------------
subset_preprocessing <- function(pp, vars) {
  keep <- pp$variables %in% vars
  v <- pp$variables[keep]
  list(variables = v, means = pp$means[v], sds = pp$sds[v], types = pp$types[v])
}
regional_catalog <- catalog[catalog$included_in_main & catalog$level == "survey-region", , drop = FALSE]

## ---- sample flags (as in script 06) ----------------------------------------
positive_outcomes <- is.finite(analysis$postneonatal_mortality) & analysis$postneonatal_mortality > 0 &
  is.finite(analysis$nnmr) & analysis$nnmr > 0 & is.finite(analysis$u5mr) & analysis$u5mr > 0 &
  is.finite(analysis$exposure) & analysis$exposure > 0
main <- as.logical(analysis$main_sample)
sample_flags <- list(
  primary = main,
  include_sub_1_percent = analysis$country_mean_pfpr_gt_1 & positive_outcomes,
  restricted_pfpr_5_to_40 = main & analysis$pfpr2_10 >= 5 & analysis$pfpr2_10 <= 40,
  complete_case = main & as.logical(analysis$complete_case_eligible)
)

## ---- log-Gaussian linear fit (copied from script 08) -----------------------
fit_log_gaussian_gam <- function(data, outcome, specification, preprocessing) {
  keep <- is.finite(data[[outcome]]) & data[[outcome]] > 0 &
    is.finite(data$pfpr10) & is.finite(data$year_c) & !is.na(data$iso3)
  dd <- data[keep, , drop = FALSE]
  dd$country <- factor(dd$iso3)
  dd$log_rate <- log(dd[[outcome]])
  ridge <- make_ridge_matrix(dd, catalog, preprocessing)
  dd$G <- ridge$matrix
  penalty <- list(G = list(diag(ncol(ridge$matrix))))
  formula <- update(model_formula(specification), log_rate ~ . - offset(log(exposure)))
  model <- mgcv::gam(formula, family = gaussian(), method = "REML", paraPen = penalty, data = dd)
  list(model = model, data = dd, preprocessing = ridge$preprocessing,
       specification = specification, outcome = outcome,
       include_country_slope = TRUE, spline_k = 6, year_k = 8)
}

## ---- refit the linear PfPR summary under each specification ----------------
lin <- function(fit) {
  r <- model_summary_row(fit)
  data.frame(est = r$pct_change_per_10, lo = r$pct_change_lo, hi = r$pct_change_hi,
             n = r$n, p = r$pfpr_p)
}
nb_fit <- function(flag, cat, pp, slope = TRUE) {
  fit_ridge_gam(analysis[flag, , drop = FALSE], OUT, cat, LIN,
                method = "REML", preprocessing = pp, include_country_slope = slope)
}

rows <- list()
add <- function(label, group, r, ref = FALSE)
  rows[[length(rows) + 1]] <<- cbind(data.frame(label = label, group = group, ref = ref), r)

add("Full covariate set (primary)", "Primary model",
    lin(nb_fit(sample_flags$primary, catalog, bundle$preprocessing)), ref = TRUE)
add("Survey-region covariates only", "Covariate adjustment",
    lin(nb_fit(sample_flags$primary, regional_catalog,
               subset_preprocessing(bundle$preprocessing, regional_catalog$variable))))
add("Include PfPR₂₋₁₀ < 1% regions", "Analysis sample",
    lin(nb_fit(sample_flags$include_sub_1_percent, catalog, bundle$preprocessing)))
add("Restrict to PfPR₂₋₁₀ 5–40%", "Analysis sample",
    lin(nb_fit(sample_flags$restricted_pfpr_5_to_40, catalog, bundle$preprocessing)))
add("Complete cases (all covariates)", "Analysis sample",
    lin(nb_fit(sample_flags$complete_case, catalog, bundle$preprocessing)))
add("No country PfPR random slope", "Model structure",
    lin(nb_fit(sample_flags$primary, catalog, bundle$preprocessing, slope = FALSE)))
add("Log-Gaussian likelihood (log-rate)", "Likelihood",
    lin(fit_log_gaussian_gam(analysis[sample_flags$primary, , drop = FALSE], OUT, LIN, bundle$preprocessing)))

forest <- do.call(rbind, rows)
write.csv(forest, file.path(RESULTS_DIR, "specification_forest_summary.csv"), row.names = FALSE)
cat("Specification forest (% change in post-neonatal mortality per +10 PfPR2-10 points):\n")
print(within(forest, {est <- round(est, 1); lo <- round(lo, 1); hi <- round(hi, 1)})[, c("label", "n", "est", "lo", "hi")], row.names = FALSE)

## ---- forest plot -----------------------------------------------------------
primary <- forest[forest$ref, ]
forest$group <- factor(forest$group,
  levels = c("Primary model", "Covariate adjustment", "Analysis sample", "Model structure", "Likelihood"))
forest$label <- factor(forest$label, levels = rev(forest$label))   # top-to-bottom as listed
BLU <- "#08519c"

p <- ggplot2::ggplot(forest, ggplot2::aes(est, label)) +
  ggplot2::annotate("rect", xmin = primary$lo, xmax = primary$hi, ymin = -Inf, ymax = Inf,
                    fill = BLU, alpha = 0.10) +
  ggplot2::geom_vline(xintercept = primary$est, linetype = "dashed", colour = BLU) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi, colour = ref),
                         orientation = "y", width = 0.28, linewidth = 0.7) +
  ggplot2::geom_point(ggplot2::aes(colour = ref), size = 2.8) +
  ggplot2::scale_colour_manual(values = c("FALSE" = "grey30", "TRUE" = BLU), guide = "none") +
  ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  ggplot2::labs(
    x = "Change in post-neonatal mortality per +10  PfPR₂₋₁₀ points (%, 95% CI)",
    y = NULL,
    title = "Robustness of the prevalence–mortality association",
    subtitle = "Linear summary refit per specification; band = primary 95% CI") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
    strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 10),
    plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(file.path(RESULTS_DIR, "figure6_specification_forest.png"),
  p, width = 10, height = 5.6, dpi = 320, bg = "white")
cat("saved: figure6_specification_forest.png + specification_forest_summary.csv\n")

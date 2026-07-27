# =============================================================================
# 04_fit_main_models.R
# Compare four ridge-penalised negative-binomial GAM specifications on one
# shared post-neonatal analysis sample:
#   1. linear PfPR, no PfPR-by-time interaction
#   2. spline PfPR, no interaction
#   3. linear PfPR with linear PfPR-by-time interaction
#   4. spline PfPR with tensor PfPR-by-time interaction
#
# Every model includes:
#   * ridge-penalised transformed/standardised eligible covariates
#   * smooth calendar year
#   * country random intercept
#   * country random linear PfPR slope
#   * log-exposure offset
#
# Models are compared using ML on an identical sample. The lowest-AIC
# specification is then refitted using REML for post-neonatal, all-U5 and
# neonatal outcomes. A linear no-interaction model is also refitted for each
# outcome to provide a directly comparable per-10-point negative-control table.
# =============================================================================

source("R_dhs/00_config.R")
required_packages("mgcv")
set.seed(20260727)

analysis <- read_analysis_data()
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
catalog$included_in_main <- as.logical(catalog$included_in_main)
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]

if (!nrow(analysis)) stop("The main-sample flag selects zero rows.")
message(
  "Main model sample: ", nrow(analysis), " region-years; ",
  length(unique(analysis$iso3)), " countries; ",
  sum(catalog$included_in_main), " ridge covariates."
)

comparison_fits <- setNames(
  lapply(MODEL_SPECS$specification, function(specification) {
    message("Fitting comparison model: ", specification)
    fit_ridge_gam(
      analysis,
      outcome = "postneonatal_mortality",
      catalog = catalog,
      specification = specification,
      method = "ML"
    )
  }),
  MODEL_SPECS$specification
)

comparison <- do.call(rbind, lapply(comparison_fits, model_summary_row))
comparison$dAIC <- comparison$AIC - min(comparison$AIC)
comparison <- comparison[order(comparison$AIC), , drop = FALSE]
selected_specification <- comparison$specification[1]

message("Selected post-neonatal specification by AIC: ", selected_specification)
print(
  comparison[, c(
    "specification", "n", "countries", "AIC", "dAIC",
    "pfpr_smooth_edf", "pct_change_per_10", "pfpr_p",
    "time_interaction_beta", "time_interaction_edf", "time_interaction_p"
  )],
  row.names = FALSE
)

selected_preprocessing <- comparison_fits[[selected_specification]]$preprocessing

OUTCOMES <- c(
  postneonatal = "postneonatal_mortality",
  all_under_5 = "u5mr",
  neonatal = "nnmr"
)

primary_fits <- setNames(
  lapply(OUTCOMES, function(outcome) {
    fit_ridge_gam(
      analysis,
      outcome = outcome,
      catalog = catalog,
      specification = selected_specification,
      method = "REML",
      preprocessing = selected_preprocessing
    )
  }),
  names(OUTCOMES)
)

linear_fits <- setNames(
  lapply(OUTCOMES, function(outcome) {
    fit_ridge_gam(
      analysis,
      outcome = outcome,
      catalog = catalog,
      specification = "linear_no_interaction",
      method = "REML",
      preprocessing = selected_preprocessing
    )
  }),
  names(OUTCOMES)
)

primary_summary <- do.call(rbind, lapply(primary_fits, model_summary_row))
primary_summary$role <- "selected specification, REML"
linear_summary <- do.call(rbind, lapply(linear_fits, model_summary_row))
linear_summary$role <- "linear comparison effect, REML"

linear_effects <- linear_summary[, c(
  "outcome", "n", "countries", "pct_change_per_10",
  "pct_change_lo", "pct_change_hi", "pfpr_p"
)]

reference_year_c <- 0
af_anchors <- do.call(rbind, lapply(names(primary_fits), function(outcome) {
  values <- af_from_model(
    primary_fits[[outcome]]$model,
    prevalence = c(10, 30, 50),
    year_c = reference_year_c
  )
  values$outcome <- outcome
  values
}))

comparison$formula <- vapply(
  comparison$specification,
  function(specification) {
    paste(deparse(model_formula(specification)), collapse = " ")
  },
  character(1)
)

write.csv(
  comparison,
  file.path(RESULTS_DIR, "main_model_comparison.csv"),
  row.names = FALSE
)
write.csv(
  rbind(primary_summary, linear_summary),
  file.path(RESULTS_DIR, "main_model_summaries.csv"),
  row.names = FALSE
)
write.csv(
  linear_effects,
  file.path(RESULTS_DIR, "linear_effects_by_outcome.csv"),
  row.names = FALSE
)
write.csv(
  af_anchors,
  file.path(RESULTS_DIR, "attributable_fraction_anchors.csv"),
  row.names = FALSE
)

bundle <- list(
  version = "dhs-rebuild-1",
  selected_specification = selected_specification,
  comparison = comparison,
  comparison_fits = comparison_fits,
  primary_fits = primary_fits,
  linear_fits = linear_fits,
  preprocessing = selected_preprocessing,
  catalog = catalog,
  year_center = unique(analysis$year - analysis$year_c),
  sample_ids = analysis[, c("svkey", "regkey", "iso3", "year")],
  af_anchors = af_anchors
)
saveRDS(bundle, MODEL_BUNDLE_RDS)

message("Saved main model bundle: ", MODEL_BUNDLE_RDS)
message("Linear outcome comparison (% change per +10 PfPR points):")
print(linear_effects, row.names = FALSE)
message("Selected-model AF anchors relative to 1% PfPR:")
print(af_anchors, row.names = FALSE)

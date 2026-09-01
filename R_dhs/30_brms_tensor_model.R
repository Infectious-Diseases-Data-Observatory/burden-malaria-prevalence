# =============================================================================
# 30_brms_tensor_model.R — the prevalence-by-time surface, fitted in Stan.
#
# Why bother, given mgcv already fits this. Because the smoothing parameter is
# the problem. In mgcv, lambda is chosen by optimising a criterion and then
# treated as known, so every interval downstream is conditional on that one
# value. Where the likelihood is flat in lambda that choice is doing real work:
# in West Africa, ML drives the prevalence smooth to a straight line (edf 1.00)
# while REML keeps it curved (edf 3.63), moving the attributable fraction at 10%
# prevalence from 4.5% to 18.7% (see 29_subgroup_fits.R).
#
# Under a Bayesian fit the smoothing parameters are just more parameters. They
# get a prior, they get a posterior, and the uncertainty in HOW SMOOTH the
# surface is propagates into every attributable fraction rather than being
# fixed by an optimiser. That is the whole point of doing it this way.
#
# Two translation notes that matter:
#
#   * brms cannot use mgcv's te(). te() builds a tensor product whose penalties
#     do not decompose into the simple independent quadratic forms a mixed-model
#     or Bayesian representation needs. t2() is mgcv's alternative construction
#     designed exactly for that representation, so t2(pfpr10, year_c) is the
#     faithful counterpart of te(pfpr10, year_c) here, not an approximation of
#     convenience.
#
#   * the ridge block becomes a prior. mgcv puts the 18 standardised covariates
#     in one L2-penalised block and ESTIMATES the penalty. Here they enter as
#     population-level effects with a normal(0, 1) prior, which is a ridge with
#     the penalty FIXED at 1 on the standardised scale rather than estimated.
#     RIDGE_PRIOR_SCALE below controls it; see the note at the foot of the file
#     on adaptive alternatives.
#
# Outputs
#   results/dhs_rebuild/brms_tensor_fit.rds
#   results/dhs_rebuild/brms_tensor_diagnostics.csv
#   results/dhs_rebuild/brms_af_anchors.csv
#   results/dhs_rebuild/figure12_brms_af_surface.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("brms", "ggplot2", "mgcv"))

CHAINS <- 4L
ITERATIONS <- 2000L
WARMUP <- 1000L
ADAPT_DELTA <- 0.95
RIDGE_PRIOR_SCALE <- 1
ANCHORS <- c(10, 30, 50)
FIT_RDS <- file.path(RESULTS_DIR, "brms_tensor_fit.rds")

bundle <- if (file.exists(MODEL_BUNDLE_RDS)) readRDS(MODEL_BUNDLE_RDS) else NULL
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
catalog$included_in_main <- as.logical(catalog$included_in_main)

analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
analysis <- analysis[is.finite(analysis$postneonatal_mortality) &
                       analysis$postneonatal_mortality > 0 &
                       is.finite(analysis$exposure) & analysis$exposure > 0 &
                       is.finite(analysis$pfpr10) & is.finite(analysis$year_c),
                     , drop = FALSE]

# Reuse the pipeline's own standardisation so the covariates are on exactly the
# scale the mgcv fits use and the two are comparable.
ridge <- make_ridge_matrix(analysis, catalog,
                           if (is.null(bundle)) NULL else bundle$preprocessing)
covariates <- as.data.frame(ridge$matrix)
names(covariates) <- make.names(colnames(ridge$matrix))

model_data <- data.frame(
  deaths = round(analysis$postneonatal_mortality / 1000 * analysis$exposure),
  log_exposure = log(analysis$exposure),
  pfpr10 = analysis$pfpr10,
  year_c = analysis$year_c,
  country = factor(analysis$iso3),
  stringsAsFactors = FALSE
)
model_data <- cbind(model_data, covariates)
message("brms tensor model: ", nrow(model_data), " region-years, ",
        nlevels(model_data$country), " countries, ", ncol(covariates),
        " covariates.")

## ---- the model --------------------------------------------------------------
formula <- brms::bf(paste(
  "deaths ~ t2(pfpr10, year_c, k = c(6, 6)) +",
  paste(names(covariates), collapse = " + "),
  "+ (1 | country) + offset(log_exposure)"
))

priors <- c(
  brms::prior_string(sprintf("normal(0, %g)", RIDGE_PRIOR_SCALE), class = "b"),
  # Half-t on the smooth's standard deviations: this IS the penalisation, and
  # unlike a fitted lambda it carries uncertainty into everything downstream.
  brms::prior(student_t(3, 0, 2.5), class = "sds"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(gamma(0.01, 0.01), class = "shape")
)

# The cached fit is only reusable while the panel it was fitted to is unchanged.
# Sampling costs a couple of minutes, so a stale cache is not worth the risk of
# silently reporting a fit of superseded data.
FIT_META_RDS <- sub("\\.rds$", "_meta.rds", FIT_RDS)
fingerprint <- list(rows = nrow(model_data),
                    countries = nlevels(model_data$country),
                    covariates = names(covariates),
                    deaths = sum(model_data$deaths))
cached <- file.exists(FIT_RDS) && file.exists(FIT_META_RDS) &&
  isTRUE(all.equal(readRDS(FIT_META_RDS), fingerprint))
if (file.exists(FIT_RDS) && !cached) {
  message("Cached fit does not match the current panel; refitting.")
}

if (cached) {
  message("Reusing the cached fit at ", basename(FIT_RDS))
  fit <- readRDS(FIT_RDS)
} else {
  message("Sampling: ", CHAINS, " chains x ", ITERATIONS, " iterations")
  fit <- brms::brm(
    formula = formula,
    data = model_data,
    family = brms::negbinomial(),
    prior = priors,
    chains = CHAINS, iter = ITERATIONS, warmup = WARMUP,
    cores = min(CHAINS, max(1L, parallel::detectCores() - 1L)),
    control = list(adapt_delta = ADAPT_DELTA, max_treedepth = 12),
    seed = 20260828, refresh = 100
  )
  saveRDS(fit, FIT_RDS)
  saveRDS(fingerprint, FIT_META_RDS)
}

## ---- did it sample properly -------------------------------------------------
draws <- posterior::as_draws_df(fit)
summary_table <- posterior::summarise_draws(draws)
key <- summary_table[grepl("^b_|^sds_|^sd_|^shape", summary_table$variable), ]
diagnostics <- data.frame(
  divergent_transitions = sum(brms::nuts_params(fit)$Parameter == "divergent__" &
                                brms::nuts_params(fit)$Value > 0),
  max_rhat = max(key$rhat, na.rm = TRUE),
  min_bulk_ess = min(key$ess_bulk, na.rm = TRUE),
  min_tail_ess = min(key$ess_tail, na.rm = TRUE)
)
print(diagnostics, row.names = FALSE)
if (diagnostics$divergent_transitions > 0) {
  warning(diagnostics$divergent_transitions, " divergent transitions; raise ",
          "ADAPT_DELTA before trusting the tails.")
}
if (diagnostics$max_rhat > 1.01) {
  warning("Max R-hat ", signif(diagnostics$max_rhat, 4), " exceeds 1.01.")
}
write.csv(cbind(diagnostics,
                key[key$variable %in% grep("^sds_|^sd_", key$variable,
                                           value = TRUE),
                    c("variable", "mean", "q5", "q95", "rhat")]),
          file.path(RESULTS_DIR, "brms_tensor_diagnostics.csv"),
          row.names = FALSE)

message("\nSmoothness parameters (the penalisation, with its uncertainty):")
print(as.data.frame(key[grepl("^sds_|^sd_", key$variable),
                        c("variable", "mean", "q5", "q95", "rhat")]),
      row.names = FALSE, digits = 3)

## ---- attributable fraction, with smoothing uncertainty included -------------
# Population-average: country random effects excluded, covariates at their means
# (zero, because they are standardised) and a unit offset.
af_draws <- function(prevalence, year) {
  frame <- function(p) {
    out <- data.frame(pfpr10 = p / 10, year_c = year, log_exposure = 0,
                      country = levels(model_data$country)[1])
    for (name in names(covariates)) out[[name]] <- 0
    out
  }
  high <- brms::posterior_linpred(fit, newdata = frame(prevalence),
                                  re_formula = NA)
  low <- brms::posterior_linpred(fit, newdata = frame(rep(AF_REFERENCE,
                                                          length(prevalence))),
                                 re_formula = NA)
  pmax(1 - exp(-(high - low)), 0)
}

years <- sort(unique(round(model_data$year_c)))
anchor_rows <- list()
for (year in years) {
  values <- af_draws(ANCHORS, year)
  for (j in seq_along(ANCHORS)) {
    column <- values[, j]
    anchor_rows[[length(anchor_rows) + 1L]] <- data.frame(
      year = year + unique(bundle$year_center)[1], prevalence = ANCHORS[j],
      af = mean(column), lo = unname(stats::quantile(column, 0.025)),
      hi = unname(stats::quantile(column, 0.975)), stringsAsFactors = FALSE
    )
  }
}
anchors <- do.call(rbind, anchor_rows)
write.csv(anchors, file.path(RESULTS_DIR, "brms_af_anchors.csv"),
          row.names = FALSE)

centre_year <- round(stats::median(years))
at_centre <- anchors[anchors$year == centre_year + unique(bundle$year_center)[1], ]
message("\nAttributable fraction at the central year, with 95% credible intervals:")
print(transform(at_centre, af = round(100 * af, 1), lo = round(100 * lo, 1),
                hi = round(100 * hi, 1))[, c("prevalence", "af", "lo", "hi")],
      row.names = FALSE)

## ---- the surface -------------------------------------------------------------
grid <- expand.grid(prevalence = seq(1, 60, by = 2), year_c = years)
surface_rows <- lapply(split(grid, grid$year_c), function(g) {
  values <- af_draws(g$prevalence, g$year_c[1])
  data.frame(prevalence = g$prevalence,
             year = g$year_c[1] + unique(bundle$year_center)[1],
             af = colMeans(values),
             lo = apply(values, 2, stats::quantile, 0.025),
             hi = apply(values, 2, stats::quantile, 0.975),
             stringsAsFactors = FALSE)
})
surface <- do.call(rbind, surface_rows)

shown <- surface[surface$year %in% c(2000, 2008, 2016, 2024), ]
plot <- ggplot2::ggplot(shown, ggplot2::aes(prevalence, 100 * af,
                                            colour = factor(year),
                                            fill = factor(year))) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * lo, ymax = 100 * hi),
                       alpha = 0.15, colour = NA) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::labs(
    x = "PfPR2-10 (%)",
    y = "Malaria-attributable fraction of\npost-neonatal mortality (%)",
    title = "Prevalence-by-time surface, fitted in Stan",
    subtitle = paste("t2(pfpr10, year_c) with the smoothness estimated rather",
                     "than plugged in;\nbands are 95% credible intervals and",
                     "include uncertainty in how smooth the surface is"),
    colour = "Year", fill = "Year"
  ) +
  ggplot2::theme_minimal(base_size = 10)
ggplot2::ggsave(file.path(RESULTS_DIR, "figure12_brms_af_surface.png"), plot,
                width = 7.5, height = 5, dpi = 200)

message("\nWrote brms_tensor_fit.rds, brms_tensor_diagnostics.csv, ",
        "brms_af_anchors.csv and figure12_brms_af_surface.png")

# A note on the ridge block. normal(0, RIDGE_PRIOR_SCALE) fixes the shrinkage
# rather than estimating it, which is the one place this fit is LESS adaptive
# than the mgcv one. brms offers adaptive alternatives - horseshoe() and R2D2()
# - but both are sparsity-inducing rather than ridge-like, so they change the
# estimand for a covariate block that is meant to be shrunk jointly rather than
# selected. Fixing the scale at 1 on the standardised scale is the conservative
# choice; refit with a larger scale to check the covariates are not driving the
# prevalence surface.

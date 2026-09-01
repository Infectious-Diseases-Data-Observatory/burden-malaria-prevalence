# =============================================================================
# 31_brms_subgroup_fits.R — the subgroup fits, redone in Stan.
#
# 29_subgroup_fits.R found that in West Africa the prevalence smooth is not
# determined by the data: ML drives it to a straight line (edf 1.00) while REML
# keeps it curved (edf 3.63), and the attributable fraction at 10% prevalence
# moves from 4.5% to 18.7% depending only on which criterion picks lambda. Every
# West African subset shows the same thing. Central and East Africa does not
# (edf 2.93 against 3.07).
#
# That is a question about how much the data constrain the smoothness, which is
# exactly what a posterior answers. Here each subgroup is refitted with the
# smoothing parameters given priors rather than optimised, so the ambiguity is
# carried into the interval instead of being resolved by the criterion. The
# check worth watching: does the credible interval for West Africa SPAN the ML
# and REML answers? If it does, the two mgcv fits were both inside the range the
# data actually support, and neither was ever the answer on its own.
#
# The structure is the Bayesian counterpart of script 29's common
# specification - s(pfpr10) with a smooth year term, ridge-style prior on the
# covariate block, country random intercept - and is held IDENTICAL across
# subgroups, including the year basis, so the comparison is like for like and
# one compiled Stan program serves every fit.
#
# Outputs
#   results/dhs_rebuild/brms_subgroup_summary.csv
#   results/dhs_rebuild/brms_subgroup_smoothness.csv
#   results/dhs_rebuild/figure13_brms_subgroup_af.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("brms", "ggplot2", "mgcv"))

CHAINS <- 4L
ITERATIONS <- 2000L
WARMUP <- 1000L
# 0.95 left 0.05-0.85% divergent draws in every subgroup. They sit in the
# near-zero-smoothness funnel - precisely the region that decides whether the
# posterior reaches the straight-line answer ML gives - so the tails there have
# to be trusted, not waved through.
ADAPT_DELTA <- 0.999
RIDGE_PRIOR_SCALE <- 1
ANCHORS <- c(10, 30, 50)
# Fixed for every subgroup: the narrowest era spans only nine calendar years,
# and holding it constant also lets one compiled program serve all the fits.
YEAR_K <- 6L
CACHE <- file.path(DATA_DIR, "brms_subgroup_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

bundle <- readRDS(MODEL_BUNDLE_RDS)
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
catalog$included_in_main <- as.logical(catalog$included_in_main)

analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
analysis <- analysis[is.finite(analysis$postneonatal_mortality) &
                       analysis$postneonatal_mortality > 0 &
                       is.finite(analysis$exposure) & analysis$exposure > 0 &
                       is.finite(analysis$pfpr10) & is.finite(analysis$year_c),
                     , drop = FALSE]

WEST <- c("BEN", "BFA", "CIV", "CPV", "GHA", "GIN", "GMB", "GNB", "LBR", "MLI",
          "MRT", "NER", "NGA", "SEN", "SLE", "TGO")
analysis$region_group <- ifelse(analysis$iso3 %in% WEST, "West",
                                "Central & East")
analysis$era <- ifelse(analysis$year < 2010, "before 2010", "2010 onwards")

ridge <- make_ridge_matrix(analysis, catalog, bundle$preprocessing)
covariates <- as.data.frame(ridge$matrix)
names(covariates) <- make.names(colnames(ridge$matrix))

model_data <- data.frame(
  deaths = round(analysis$postneonatal_mortality / 1000 * analysis$exposure),
  log_exposure = log(analysis$exposure),
  pfpr10 = analysis$pfpr10,
  year_c = analysis$year_c,
  country = factor(analysis$iso3),
  region_group = analysis$region_group,
  era = analysis$era,
  stringsAsFactors = FALSE
)
model_data <- cbind(model_data, covariates)

formula <- brms::bf(paste(
  "deaths ~ s(pfpr10, k = 6) + s(year_c, k =", YEAR_K, ") +",
  paste(names(covariates), collapse = " + "),
  "+ (1 | country) + offset(log_exposure)"
))
priors <- c(
  brms::prior_string(sprintf("normal(0, %g)", RIDGE_PRIOR_SCALE), class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "sds"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(gamma(0.01, 0.01), class = "shape")
)

subsets <- list(
  list(label = "all data", split = "reference",
       rows = rep(TRUE, nrow(model_data))),
  list(label = "before 2010", split = "era",
       rows = model_data$era == "before 2010"),
  list(label = "2010 onwards", split = "era",
       rows = model_data$era == "2010 onwards"),
  list(label = "West", split = "region",
       rows = model_data$region_group == "West"),
  list(label = "Central & East", split = "region",
       rows = model_data$region_group == "Central & East")
)
for (era in c("before 2010", "2010 onwards")) {
  for (region in c("West", "Central & East")) {
    subsets[[length(subsets) + 1L]] <- list(
      label = paste0(region, ", ", era), split = "era x region",
      rows = model_data$era == era & model_data$region_group == region
    )
  }
}

## ---- attributable fraction from a fitted subgroup ---------------------------
af_posterior <- function(fit, data, prevalence) {
  frame <- function(p) {
    out <- data.frame(pfpr10 = p / 10,
                      year_c = mean(data$year_c, na.rm = TRUE),
                      log_exposure = 0,
                      country = levels(droplevels(data$country))[1],
                      stringsAsFactors = FALSE)
    for (name in names(covariates)) out[[name]] <- 0
    out[rep(1, length(p)), , drop = FALSE]
  }
  high <- frame(prevalence)
  high$pfpr10 <- prevalence / 10
  low <- frame(prevalence)
  low$pfpr10 <- AF_REFERENCE / 10
  linear_high <- brms::posterior_linpred(fit, newdata = high, re_formula = NA)
  linear_low <- brms::posterior_linpred(fit, newdata = low, re_formula = NA)
  pmax(1 - exp(-(linear_high - linear_low)), 0)
}

## ---- fit every subgroup on one compiled program -----------------------------
template <- NULL
rows_out <- list()
smoothness_out <- list()

for (s in subsets) {
  data <- model_data[s$rows, , drop = FALSE]
  data$country <- droplevels(data$country)
  if (nrow(data) < 80L || nlevels(data$country) < 5L) {
    message("  ", s$label, ": too small; skipped")
    next
  }
  fingerprint <- list(rows = nrow(data), countries = nlevels(data$country),
                      deaths = sum(data$deaths), year_k = YEAR_K,
                      adapt_delta = ADAPT_DELTA, iterations = ITERATIONS,
                      chains = CHAINS)
  cache_file <- file.path(CACHE, paste0(make.names(s$label), ".rds"))
  meta_file <- sub("\\.rds$", "_meta.rds", cache_file)
  reusable <- file.exists(cache_file) && file.exists(meta_file) &&
    isTRUE(all.equal(readRDS(meta_file), fingerprint))

  if (reusable) {
    message("  ", s$label, ": reusing cached fit")
    fit <- readRDS(cache_file)
  } else {
    message("  ", s$label, ": sampling (", nrow(data), " rows, ",
            nlevels(data$country), " countries)")
    fit <- if (is.null(template)) {
      brms::brm(formula = formula, data = data,
                family = brms::negbinomial(), prior = priors,
                chains = CHAINS, iter = ITERATIONS, warmup = WARMUP,
                cores = min(CHAINS, max(1L, parallel::detectCores() - 1L)),
                control = list(adapt_delta = ADAPT_DELTA, max_treedepth = 14),
                seed = 20260828, refresh = 0)
    } else {
      # Same program, new data: skip the compile.
      tryCatch(
        stats::update(template, newdata = data, recompile = FALSE,
                      chains = CHAINS, iter = ITERATIONS, warmup = WARMUP,
                      cores = min(CHAINS, max(1L, parallel::detectCores() - 1L)),
                      control = list(adapt_delta = ADAPT_DELTA,
                                     max_treedepth = 14),
                      seed = 20260828, refresh = 0),
        error = function(e) {
          message("    update() failed (", conditionMessage(e),
                  "); compiling separately")
          brms::brm(formula = formula, data = data,
                    family = brms::negbinomial(), prior = priors,
                    chains = CHAINS, iter = ITERATIONS, warmup = WARMUP,
                    cores = min(CHAINS,
                                max(1L, parallel::detectCores() - 1L)),
                    control = list(adapt_delta = ADAPT_DELTA,
                                   max_treedepth = 14),
                    seed = 20260828, refresh = 0)
        }
      )
    }
    saveRDS(fit, cache_file)
    saveRDS(fingerprint, meta_file)
  }
  if (is.null(template)) template <- fit

  parameters <- posterior::summarise_draws(posterior::as_draws_df(fit))
  key <- parameters[grepl("^b_|^sds_|^sd_|^shape", parameters$variable), ]
  divergent <- sum(brms::nuts_params(fit)$Parameter == "divergent__" &
                     brms::nuts_params(fit)$Value > 0)

  values <- af_posterior(fit, data, ANCHORS)
  for (j in seq_along(ANCHORS)) {
    column <- values[, j]
    rows_out[[length(rows_out) + 1L]] <- data.frame(
      split = s$split, subset = s$label, n = nrow(data),
      countries = nlevels(data$country), prevalence = ANCHORS[j],
      af = mean(column),
      lo = unname(stats::quantile(column, 0.025)),
      hi = unname(stats::quantile(column, 0.975)),
      divergent_transitions = divergent,
      max_rhat = max(key$rhat, na.rm = TRUE),
      min_ess = min(key$ess_bulk, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  smooth_row <- key[grepl("^sds_spfpr10", key$variable), ]
  if (nrow(smooth_row)) {
    smoothness_out[[length(smoothness_out) + 1L]] <- data.frame(
      subset = s$label, n = nrow(data),
      sds_mean = smooth_row$mean[1], sds_q5 = smooth_row$q5[1],
      sds_q95 = smooth_row$q95[1], stringsAsFactors = FALSE
    )
  }
}

summary_table <- do.call(rbind, rows_out)
write.csv(summary_table,
          file.path(RESULTS_DIR, "brms_subgroup_summary.csv"), row.names = FALSE)
smoothness <- do.call(rbind, smoothness_out)
write.csv(smoothness,
          file.path(RESULTS_DIR, "brms_subgroup_smoothness.csv"),
          row.names = FALSE)

## ---- does the posterior span what ML and REML disagreed about? -------------
mgcv_reference <- read.csv(file.path(RESULTS_DIR, "subgroup_fits_summary.csv"),
                           stringsAsFactors = FALSE)
comparison <- merge(
  summary_table[, c("subset", "prevalence", "af", "lo", "hi")],
  data.frame(
    subset = rep(mgcv_reference$subset, 3),
    prevalence = rep(ANCHORS, each = nrow(mgcv_reference)),
    mgcv_reml = c(mgcv_reference$af10_common, mgcv_reference$af30_common,
                  mgcv_reference$af50_common),
    stringsAsFactors = FALSE
  ),
  by = c("subset", "prevalence")
)
comparison$reml_inside_credible <- comparison$mgcv_reml >= comparison$lo &
  comparison$mgcv_reml <= comparison$hi

message("\nDiagnostics: ", sum(summary_table$divergent_transitions > 0) / 3,
        " subgroups with divergent transitions; worst R-hat ",
        signif(max(summary_table$max_rhat), 4))
message("\nAttributable fraction by subgroup, 95% credible intervals:")
report <- summary_table[, c("subset", "n", "prevalence", "af", "lo", "hi")]
report[c("af", "lo", "hi")] <- lapply(report[c("af", "lo", "hi")],
                                      function(x) round(100 * x, 1))
print(report, row.names = FALSE)
message("\nHow tightly is the smoothness itself determined?")
print(transform(smoothness, sds_mean = round(sds_mean, 2),
                sds_q5 = round(sds_q5, 2), sds_q95 = round(sds_q95, 2)),
      row.names = FALSE)

## ---- figure ------------------------------------------------------------------
summary_table$label <- factor(summary_table$subset,
                              levels = rev(unique(summary_table$subset)))
plot <- ggplot2::ggplot(summary_table,
                        ggplot2::aes(100 * af, label, colour = split)) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = 100 * lo, xmax = 100 * hi),
                          height = 0.25, linewidth = 0.6) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~prevalence, nrow = 1,
                      labeller = ggplot2::labeller(
                        prevalence = function(x) paste0("PfPR ", x, "%"))) +
  ggplot2::labs(
    x = "Malaria-attributable fraction of post-neonatal mortality (%)",
    y = NULL, colour = NULL,
    title = "Subgroup attributable fractions, fitted in Stan",
    subtitle = paste("Identical structure in every subgroup; bands are 95%",
                     "credible intervals\nand include uncertainty in the",
                     "smoothness of the prevalence curve")
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(RESULTS_DIR, "figure13_brms_subgroup_af.png"), plot,
                width = 11, height = 4.6, dpi = 200)

message("\nWrote brms_subgroup_summary.csv, brms_subgroup_smoothness.csv and ",
        "figure13_brms_subgroup_af.png")

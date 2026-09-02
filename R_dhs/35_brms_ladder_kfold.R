# =============================================================================
# 35_brms_ladder_kfold.R — survey-grouped 10-fold cross-validation of the ladder.
#
# Script 34 could not separate the additive (A), linear-in-time (C) and full-
# surface (D) models by PSIS-LOO. But pointwise LOO holds out one region-year at
# a time while the other regions of the same survey stay in training: same
# country, same year, same fieldwork. That asks "can the model predict one more
# region of a survey it already has", which flatters flexible time structure.
# Holding out whole surveys asks whether the prevalence-by-time structure
# generalises to a survey the model has never seen, which is the question the
# interaction is really about.
#
# Folds: 10, formed by loo::kfold_split_grouped on the survey key with a fixed
# seed, identical for all three models so the comparison is paired. Each fold
# refits the model with the priors and sampler settings of the original
# (brms::kfold calls update()). A held-out survey whose country has no other
# survey is predicted with a country intercept drawn from the fitted
# distribution of country intercepts, the same for every model.
#
# Usage
#   LADDER_MODEL=A Rscript R_dhs/35_brms_ladder_kfold.R   run and cache one model
#   Rscript R_dhs/35_brms_ladder_kfold.R                  compare (runs any missing)
#
# Outputs
#   results/dhs_rebuild/brms_ladder_kfold.csv
#   results/dhs_rebuild/brms_ladder_kfold_by_era.csv
#   results/dhs_rebuild/brms_ladder_kfold_folds.csv
#   results/dhs_rebuild/figure18_brms_ladder_kfold.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("brms", "loo", "ggplot2"))

K <- 10L
FOLD_SEED <- 20260828L
CACHE <- file.path(DATA_DIR, "brms_ladder_cache")
MODELS <- c("A", "C", "D")
LABELS <- c(A = "A: additive", C = "C: curve changes linearly in time",
            D = "D: full surface")
FIT_CORES <- min(4L, max(1L, parallel::detectCores() - 1L))
options(mc.cores = FIT_CORES)

load_fit <- function(name) {
  file <- file.path(CACHE, paste0(name, ".rds"))
  if (!file.exists(file)) {
    stop("Run 34_brms_model_ladder.R first; missing ", basename(file))
  }
  readRDS(file)
}

## ---- recover the survey key for each row of the fitted panel ------------------
# brms keeps only formula variables in fit$data, so the panel is rebuilt exactly
# as script 34 built it and checked row by row against the fit.
bundle <- readRDS(MODEL_BUNDLE_RDS)
analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
analysis <- analysis[is.finite(analysis$postneonatal_mortality) &
                       analysis$postneonatal_mortality > 0 &
                       is.finite(analysis$exposure) & analysis$exposure > 0 &
                       is.finite(analysis$pfpr10) & is.finite(analysis$year_c),
                     , drop = FALSE]
reference_fit <- load_fit("A")
stopifnot(
  nrow(analysis) == nrow(reference_fit$data),
  isTRUE(all.equal(
    round(analysis$postneonatal_mortality / 1000 * analysis$exposure),
    reference_fit$data$deaths)),
  isTRUE(all.equal(analysis$pfpr10, reference_fit$data$pfpr10)),
  isTRUE(all.equal(analysis$year_c, reference_fit$data$year_c))
)
panel <- data.frame(svkey = analysis$svkey, iso3 = analysis$iso3,
                    year = analysis$year, pfpr = 10 * analysis$pfpr10,
                    stringsAsFactors = FALSE)
panel$era <- ifelse(panel$year < ERA_CUT, ERA_EARLY, ERA_LATE)

set.seed(FOLD_SEED)
folds <- loo::kfold_split_grouped(K = K, x = panel$svkey)
fold_summary <- do.call(rbind, lapply(split(panel, folds), function(d) {
  data.frame(surveys = length(unique(d$svkey)), regions = nrow(d),
             countries = length(unique(d$iso3)),
             stringsAsFactors = FALSE)
}))
fold_summary$fold <- as.integer(rownames(fold_summary))
write.csv(fold_summary[, c("fold", "surveys", "countries", "regions")],
          file.path(RESULTS_DIR, "brms_ladder_kfold_folds.csv"),
          row.names = FALSE)
message(K, " folds over ", length(unique(panel$svkey)), " surveys: ",
        paste(fold_summary$surveys, collapse = "/"), " surveys per fold")

## ---- one model's folds, cached ------------------------------------------------------
run_kfold <- function(name) {
  cache_file <- file.path(CACHE, paste0("kfold_", name, ".rds"))
  meta_file <- sub("\\.rds$", "_meta.rds", cache_file)
  fingerprint <- list(K = K, folds = folds, seed = FOLD_SEED,
                      fit_meta = readRDS(file.path(CACHE,
                                                   paste0(name, "_meta.rds"))))
  if (file.exists(cache_file) && file.exists(meta_file) &&
      isTRUE(all.equal(readRDS(meta_file), fingerprint))) {
    message("  ", LABELS[[name]], ": reusing cached k-fold")
    return(readRDS(cache_file))
  }
  fit <- load_fit(name)
  message("  ", LABELS[[name]], ": refitting ", K, " times")
  started <- Sys.time()
  kf <- brms::kfold(fit, folds = folds, save_fits = FALSE, compare = FALSE,
                    recompile = FALSE, cores = FIT_CORES)
  message("  ", LABELS[[name]], ": done in ",
          round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
          " minutes")
  saveRDS(kf, cache_file)
  saveRDS(fingerprint, meta_file)
  kf
}

only <- Sys.getenv("LADDER_MODEL", "")
if (nzchar(only)) {
  if (!only %in% MODELS) stop("LADDER_MODEL must be one of A, C, D")
  invisible(run_kfold(only))
  quit(save = "no", status = 0)
}

kfolds <- lapply(MODELS, run_kfold)
names(kfolds) <- MODELS

## ---- paired comparison ---------------------------------------------------------------
comparison <- loo::loo_compare(kfolds)
print(comparison, simplify = FALSE)
pointwise <- sapply(kfolds, function(k) k$pointwise[, "elpd_kfold"])
weights <- loo::stacking_weights(pointwise)
message("\nStacking weights on held-out surveys:")
print(round(weights, 3))

loo_reference <- file.path(RESULTS_DIR, "brms_ladder_loo.csv")
psis <- if (file.exists(loo_reference)) {
  read.csv(loo_reference, stringsAsFactors = FALSE)
} else NULL
kfold_table <- data.frame(
  model = rownames(comparison),
  label = LABELS[rownames(comparison)],
  elpd_kfold = comparison[, "elpd_kfold"],
  se_elpd_kfold = comparison[, "se_elpd_kfold"],
  elpd_diff = comparison[, "elpd_diff"],
  se_diff = comparison[, "se_diff"],
  stacking_weight = as.numeric(weights)[match(rownames(comparison), MODELS)],
  stringsAsFactors = FALSE
)
if (!is.null(psis)) {
  kfold_table$elpd_loo_pointwise <- psis$elpd_loo[match(kfold_table$model,
                                                        psis$model)]
}
# loo_compare's standard error treats region-years as independent. In a
# survey-grouped design the exchangeable unit is the survey, so the paired
# difference is also summarised at that level: sum within survey, then the
# standard error across the 118 surveys. This is the stricter of the two.
survey_level <- function(m, ref = "A") {
  d <- tapply(pointwise[, m] - pointwise[, ref], panel$svkey, sum)
  c(total = sum(d), se_survey = stats::sd(d) / sqrt(length(d)) * length(d),
    surveys_favouring = sum(d > 0), surveys = length(d),
    total_without_top3 = sum(d) - sum(sort(d, decreasing = TRUE)[1:3]))
}
versus_A <- t(sapply(kfold_table$model, survey_level))
kfold_table$elpd_diff_vs_A <- versus_A[, "total"]
kfold_table$se_vs_A_survey_level <- versus_A[, "se_survey"]
kfold_table$surveys_favouring_over_A <- as.integer(versus_A[, "surveys_favouring"])
kfold_table$elpd_diff_vs_A_without_top3 <- versus_A[, "total_without_top3"]
rownames(kfold_table) <- NULL
message("\nAgainst the additive model, survey-level: ")
print(transform(kfold_table[, c("model", "elpd_diff_vs_A", "se_vs_A_survey_level",
                                "surveys_favouring_over_A",
                                "elpd_diff_vs_A_without_top3")],
                elpd_diff_vs_A = round(elpd_diff_vs_A, 1),
                se_vs_A_survey_level = round(se_vs_A_survey_level, 1),
                elpd_diff_vs_A_without_top3 = round(elpd_diff_vs_A_without_top3, 1)),
      row.names = FALSE)
write.csv(kfold_table, file.path(RESULTS_DIR, "brms_ladder_kfold.csv"),
          row.names = FALSE)

## ---- where does any difference come from? ---------------------------------------------
# Paired differences against the additive model, summed within era, with the
# standard error from the pointwise differences.
by_group <- function(group, label) {
  do.call(rbind, lapply(split(seq_len(nrow(panel)), group), function(rows) {
    do.call(rbind, lapply(c("C", "D"), function(m) {
      d <- pointwise[rows, m] - pointwise[rows, "A"]
      data.frame(split = label, group = as.character(group[rows[1]]),
                 model = m, regions = length(rows),
                 surveys = length(unique(panel$svkey[rows])),
                 elpd_diff_vs_A = sum(d),
                 se_diff = sqrt(length(d)) * stats::sd(d),
                 stringsAsFactors = FALSE)
    }))
  }))
}
WEST <- c("BEN", "BFA", "CIV", "CPV", "GHA", "GIN", "GMB", "GNB", "LBR", "MLI",
          "MRT", "NER", "NGA", "SEN", "SLE", "TGO")
panel$region_group <- ifelse(panel$iso3 %in% WEST, "West", "Central & East")
panel$prevalence_band <- cut(panel$pfpr, c(-Inf, 10, 30, Inf),
                             labels = c("PfPR < 10%", "10-30%", "30%+"))
by_era <- rbind(by_group(panel$era, "era"),
                by_group(panel$region_group, "region"),
                by_group(panel$prevalence_band, "prevalence band"))
write.csv(by_era, file.path(RESULTS_DIR, "brms_ladder_kfold_by_era.csv"),
          row.names = FALSE)
message("\nHeld-out elpd relative to the additive model, by subset:")
print(transform(by_era, elpd_diff_vs_A = round(elpd_diff_vs_A, 1),
                se_diff = round(se_diff, 1)), row.names = FALSE)

## ---- per-survey figure ---------------------------------------------------------------------
per_survey <- do.call(rbind, lapply(split(seq_len(nrow(panel)), panel$svkey),
                                    function(rows) {
  do.call(rbind, lapply(c("C", "D"), function(m) {
    data.frame(svkey = panel$svkey[rows[1]], iso3 = panel$iso3[rows[1]],
               year = panel$year[rows[1]], regions = length(rows),
               region_group = panel$region_group[rows[1]],
               model = LABELS[[m]],
               elpd_diff_vs_A = sum(pointwise[rows, m] - pointwise[rows, "A"]),
               stringsAsFactors = FALSE)
  }))
}))
plot <- ggplot2::ggplot(per_survey,
                        ggplot2::aes(year, elpd_diff_vs_A, size = regions,
                                     shape = region_group)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  ggplot2::geom_vline(xintercept = ERA_CUT - 0.5, linetype = "dotted",
                      colour = "grey55") +
  ggplot2::geom_point(alpha = 0.7, colour = "#1D6F8B") +
  ggplot2::facet_wrap(~model, nrow = 1) +
  ggplot2::scale_shape_manual(values = c("Central & East" = 16, "West" = 17),
                              name = "Region") +
  ggplot2::scale_size_area(max_size = 5, name = "Survey regions") +
  ggplot2::labs(
    x = "Survey year",
    y = "Held-out log predictive density,\ndifference from the additive model",
    title = "Survey-grouped 10-fold cross-validation: does time structure help predict an unseen survey?",
    subtitle = paste0(
      "One point per held-out survey (", length(unique(panel$svkey)),
      " surveys). Total: ",
      paste(sprintf("%s %+.1f (se %.1f)", kfold_table$model,
                    kfold_table$elpd_diff, kfold_table$se_diff),
            collapse = "; "),
      " against the best; stacking weights ",
      paste(sprintf("%s %.2f", kfold_table$model, kfold_table$stacking_weight),
            collapse = ", "))) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(RESULTS_DIR, "figure18_brms_ladder_kfold.png"), plot,
                width = 11, height = 4.8, dpi = 200)
message("\nWrote brms_ladder_kfold.csv, brms_ladder_kfold_by_era.csv, ",
        "brms_ladder_kfold_folds.csv and figure18_brms_ladder_kfold.png")

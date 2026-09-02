# =============================================================================
# 36_neonatal_as_covariate.R — neonatal mortality as a covariate for post-neonatal.
#
# Neonatal and post-neonatal mortality come from the same birth histories and
# share everything that shapes child survival in a region apart from malaria
# exposure after the first month: care-seeking, nutrition, reporting quality and
# the survey's own sampling noise. The neonatal outcome is null with respect to
# prevalence (the negative control), so conditioning on it should absorb outcome
# variability without absorbing the malaria effect. If the malaria estimate
# survives, that is further evidence against confounding by general
# child-survival conditions; if it collapses, the two outcomes were sharing the
# association.
#
# Two versions of the covariate are used. The 12-month neonatal rate is the one
# the primary sample carries, but at that window it rests on very few deaths:
# the jackknife errors imply that about three quarters of its between-region
# variance is sampling noise, which both weakens its observed correlation with
# the outcome and limits how much variance it can absorb. The 60-month neonatal
# rate for the same survey regions (from 27_horizon_lag_selection.R) is a far
# less noisy measure of the same background conditions, and is the covariate
# that actually tests the idea. Both enter on the log scale after standardising.
#
# One further caveat: if malaria in pregnancy raises neonatal mortality,
# adjusting for it removes part of the malaria effect. The null negative control
# says that channel is small here.
#
# Parts
#   1  scatter of neonatal against post-neonatal mortality, 12- and 60-month
#      windows, with correlations corrected for sampling error
#   2  mgcv: selected specification and linear summary, with and without each term
#   3  brms: the full t2 surface with and without each term, compared by PSIS-LOO
#
# Outputs
#   results/dhs_rebuild/figure19_neonatal_vs_postneonatal.png
#   results/dhs_rebuild/neonatal_vs_postneonatal_correlation.csv
#   results/dhs_rebuild/neonatal_covariate_mgcv.csv
#   results/dhs_rebuild/neonatal_covariate_brms.csv
#   results/dhs_rebuild/figure20_neonatal_covariate_curves.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "brms", "loo", "ggplot2"))

CHAINS <- 4L
ITERATIONS <- 2000L
WARMUP <- 1000L
ADAPT_DELTA <- 0.99
RIDGE_PRIOR_SCALE <- 1
ANCHORS <- c(10, 30, 50)
CONTRAST_YEARS <- c(2005L, 2020L)
CACHE <- file.path(DATA_DIR, "brms_ladder_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
HORIZON_FILE <- file.path(RESULTS_DIR, "horizon_mortality_estimates.csv")

bundle <- readRDS(MODEL_BUNDLE_RDS)
year_center <- unique(bundle$year_center)[1]
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
catalog$included_in_main <- as.logical(catalog$included_in_main)

# The two covariate versions: column name, label, cache stem.
VERSIONS <- list(
  none = list(column = NULL, label = "Without neonatal mortality", stem = NULL),
  m12 = list(column = "log_nnmr_z", label = "With 12-month neonatal rate",
             stem = "D_nnmr"),
  m60 = list(column = "log_nnmr60_z", label = "With 60-month neonatal rate",
             stem = "D_nnmr60")
)

## ---- Part 1: how correlated are the two outcomes? ------------------------------
horizon <- read.csv(HORIZON_FILE, stringsAsFactors = FALSE)
pairs <- horizon[horizon$horizon_months %in% c(12L, 60L), , drop = FALSE]
pairs$window <- factor(sprintf("%d-month window", pairs$horizon_months),
                       levels = c("12-month window", "60-month window"))

analysis <- read_analysis_data()
main <- analysis[as.logical(analysis$main_sample), , drop = FALSE]

# Each regional rate is estimated from a handful of deaths, so the observed
# correlation is attenuated by sampling error in BOTH variables. With the
# jackknife standard errors the reliability of each variable is
# 1 - mean(se^2) / var(estimate), and dividing the observed correlation by the
# square root of the product of the two reliabilities recovers the correlation
# between the underlying regional rates (classical disattenuation).
correlate <- function(x, y, label, x_se = NULL, y_se = NULL) {
  ok <- is.finite(x) & is.finite(y)
  pos <- ok & x > 0 & y > 0
  r <- stats::cor(x[ok], y[ok])
  disattenuated <- reliability_x <- reliability_y <- NA_real_
  if (!is.null(x_se) && !is.null(y_se)) {
    reliability_x <- 1 - mean(x_se[ok]^2, na.rm = TRUE) / stats::var(x[ok])
    reliability_y <- 1 - mean(y_se[ok]^2, na.rm = TRUE) / stats::var(y[ok])
    disattenuated <- r / sqrt(reliability_x * reliability_y)
  }
  data.frame(
    sample = label, n = sum(ok),
    pearson = r,
    spearman = stats::cor(x[ok], y[ok], method = "spearman"),
    pearson_log = stats::cor(log(x[pos]), log(y[pos])),
    reliability_neonatal = reliability_x,
    reliability_postneonatal = reliability_y,
    pearson_disattenuated = disattenuated,
    ols_slope = unname(stats::coef(stats::lm(y[ok] ~ x[ok]))[2]),
    stringsAsFactors = FALSE)
}
correlation <- rbind(
  do.call(rbind, lapply(split(pairs, pairs$window), function(d)
    correlate(d$nnmr, d$postneonatal, paste("all regions,", d$window[1]),
              x_se = d$nnmr_se, y_se = d$postneonatal_se))),
  correlate(main$nnmr, main$postneonatal_mortality,
            "main sample, 12-month window")
)
rownames(correlation) <- NULL
write.csv(correlation,
          file.path(RESULTS_DIR, "neonatal_vs_postneonatal_correlation.csv"),
          row.names = FALSE)
message("Neonatal against post-neonatal mortality:")
print(transform(correlation, pearson = round(pearson, 2),
                spearman = round(spearman, 2), pearson_log = round(pearson_log, 2),
                reliability_neonatal = round(reliability_neonatal, 2),
                reliability_postneonatal = round(reliability_postneonatal, 2),
                pearson_disattenuated = round(pearson_disattenuated, 2),
                ols_slope = round(ols_slope, 2)), row.names = FALSE)

labels <- correlation[grepl("^all regions", correlation$sample), ]
labels$window <- factor(sub("all regions, ", "", labels$sample),
                        levels = levels(pairs$window))
labels$text <- sprintf(
  "n = %d regions\nPearson r = %.2f\nSpearman = %.2f\ncorrected for sampling error: %.2f",
  labels$n, labels$pearson, labels$spearman, labels$pearson_disattenuated)
scatter <- ggplot2::ggplot(pairs, ggplot2::aes(nnmr, postneonatal)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = postneonatal_lo, ymax = postneonatal_hi),
                         colour = "grey75", linewidth = 0.2, alpha = 0.3, width = 0) +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = nnmr_lo, xmax = nnmr_hi),
                         colour = "grey75", linewidth = 0.2, alpha = 0.3, width = 0,
                         orientation = "y") +
  ggplot2::geom_point(size = 0.8, alpha = 0.6, colour = "#1D6F8B") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                       colour = "#B2182B", linewidth = 0.7) +
  # hjust = 0 keeps every line of the label left-aligned; the padding comes from
  # the axis expansion rather than a negative hjust, which shifts long lines.
  ggplot2::geom_text(data = labels, ggplot2::aes(x = -Inf, y = Inf, label = text),
                     hjust = 0, vjust = 1.2, size = 3, colour = "grey25",
                     inherit.aes = FALSE) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.04))) +
  ggplot2::facet_wrap(~window, nrow = 1, scales = "free") +
  ggplot2::labs(
    x = "Neonatal mortality (deaths per 1000 live births)",
    y = "Post-neonatal under-5 mortality\n(deaths per 1000)",
    title = "Neonatal against post-neonatal mortality, one point per survey region",
    subtitle = paste0("Bars are 95% delete-one-cluster jackknife intervals; ",
                      "red line is the least-squares fit.\n",
                      "The corrected correlation divides r by the square root of ",
                      "the two reliabilities implied by the jackknife errors")) +
  ggplot2::theme_minimal(base_size = 10)
ggplot2::ggsave(file.path(RESULTS_DIR, "figure19_neonatal_vs_postneonatal.png"),
                scatter, width = 9.5, height = 4.4, dpi = 200)

## ---- the two covariates ------------------------------------------------------------
main <- main[is.finite(main$postneonatal_mortality) & main$postneonatal_mortality > 0 &
               is.finite(main$nnmr) & main$nnmr > 0 &
               is.finite(main$exposure) & main$exposure > 0, , drop = FALSE]
main$log_nnmr_z <- as.numeric(scale(log(main$nnmr)))
sixty <- horizon[horizon$horizon_months == 60L, , drop = FALSE]
main$nnmr60 <- sixty$nnmr[match(paste(main$svkey, main$regkey),
                                paste(sixty$svkey, sixty$regkey))]
if (any(!is.finite(main$nnmr60) | main$nnmr60 <= 0)) {
  stop("60-month neonatal rate missing or zero for ",
       sum(!is.finite(main$nnmr60) | main$nnmr60 <= 0), " main-sample rows")
}
main$log_nnmr60_z <- as.numeric(scale(log(main$nnmr60)))
message("\nCovariates: 12- and 60-month log neonatal rates correlate at ",
        round(stats::cor(main$log_nnmr_z, main$log_nnmr60_z), 2),
        " across the ", nrow(main), " region-years")

## ---- Part 2: mgcv, with and without each term ---------------------------------------
mgcv_rows <- list()
for (spec in unique(c(bundle$selected_specification, "linear_no_interaction"))) {
  for (v in names(VERSIONS)) {
    column <- VERSIONS[[v]]$column
    fit <- fit_ridge_gam(
      main, "postneonatal_mortality", catalog, specification = spec,
      method = "REML", preprocessing = bundle$preprocessing,
      extra_terms = column)
    row <- model_summary_row(fit)
    af <- af_from_model(fit$model, ANCHORS)
    ptab <- summary(fit$model)$p.table
    nn <- if (is.null(column)) rep(NA_real_, 4) else ptab[column, ]
    mgcv_rows[[length(mgcv_rows) + 1L]] <- data.frame(
      specification = spec, neonatal_covariate = v, n = row$n,
      AIC = stats::AIC(fit$model),
      deviance_explained = summary(fit$model)$dev.expl,
      nb_theta = fit$model$family$getTheta(TRUE),
      pct_change_per_10 = row$pct_change_per_10,
      pct_change_lo = row$pct_change_lo, pct_change_hi = row$pct_change_hi,
      pfpr_p = row$pfpr_p,
      time_interaction_beta = row$time_interaction_beta,
      time_interaction_p = row$time_interaction_p,
      af10 = af$af[1], af30 = af$af[2], af50 = af$af[3],
      af30_lo = af$lo[2], af30_hi = af$hi[2],
      neonatal_pct_per_sd = 100 * (exp(nn[1]) - 1),
      neonatal_p = nn[4],
      stringsAsFactors = FALSE)
  }
}
mgcv_table <- do.call(rbind, mgcv_rows)
rownames(mgcv_table) <- NULL
write.csv(mgcv_table, file.path(RESULTS_DIR, "neonatal_covariate_mgcv.csv"),
          row.names = FALSE)
message("\nmgcv, ", nrow(main), " region-years:")
print(transform(mgcv_table, AIC = round(AIC, 1),
                deviance_explained = round(deviance_explained, 3),
                nb_theta = round(nb_theta, 1),
                pct_change_per_10 = round(pct_change_per_10, 2),
                pct_change_lo = round(pct_change_lo, 2),
                pct_change_hi = round(pct_change_hi, 2),
                af30 = round(100 * af30, 1),
                af30_lo = round(100 * af30_lo, 1), af30_hi = round(100 * af30_hi, 1),
                neonatal_pct_per_sd = round(neonatal_pct_per_sd, 1),
                neonatal_p = signif(neonatal_p, 2))[
                  , c("specification", "neonatal_covariate", "AIC",
                      "deviance_explained", "nb_theta", "pct_change_per_10",
                      "pct_change_lo", "pct_change_hi", "af30", "af30_lo",
                      "af30_hi", "neonatal_pct_per_sd", "neonatal_p")],
      row.names = FALSE)

## ---- Part 3: brms full surface, with and without each term ------------------------
ridge <- make_ridge_matrix(main, catalog, bundle$preprocessing)
covariates <- as.data.frame(ridge$matrix)
names(covariates) <- make.names(colnames(ridge$matrix))
model_data <- data.frame(
  deaths = round(main$postneonatal_mortality / 1000 * main$exposure),
  log_exposure = log(main$exposure), pfpr10 = main$pfpr10, year_c = main$year_c,
  country = factor(main$iso3), log_nnmr_z = main$log_nnmr_z,
  log_nnmr60_z = main$log_nnmr60_z, stringsAsFactors = FALSE)
model_data <- cbind(model_data, covariates)

base_file <- file.path(CACHE, "D.rds")
if (!file.exists(base_file)) stop("Run 34_brms_model_ladder.R first (needs D.rds)")
base <- readRDS(base_file)
stopifnot(nrow(base$data) == nrow(model_data),
          isTRUE(all.equal(base$data$deaths, model_data$deaths)))

priors <- c(
  brms::prior_string(sprintf("normal(0, %g)", RIDGE_PRIOR_SCALE), class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "sds"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(gamma(0.01, 0.01), class = "shape"))

fit_adjusted <- function(column, stem) {
  formula_adj <- paste(
    "deaths ~ t2(pfpr10, year_c, k = c(6, 6)) +",
    paste(names(covariates), collapse = " + "),
    paste0("+ ", column, " + (1 | country) + offset(log_exposure)"))
  fingerprint <- list(rows = nrow(model_data), deaths = sum(model_data$deaths),
                      formula = formula_adj, chains = CHAINS,
                      iterations = ITERATIONS, adapt_delta = ADAPT_DELTA)
  cache_file <- file.path(CACHE, paste0(stem, ".rds"))
  meta_file <- file.path(CACHE, paste0(stem, "_meta.rds"))
  if (file.exists(cache_file) && file.exists(meta_file) &&
      isTRUE(all.equal(readRDS(meta_file), fingerprint))) {
    message("  ", column, ": reusing cached fit")
    return(readRDS(cache_file))
  }
  message("  ", column, ": sampling")
  fit <- brms::brm(
    formula = brms::bf(formula_adj), data = model_data,
    family = brms::negbinomial(), prior = priors,
    chains = CHAINS, iter = ITERATIONS, warmup = WARMUP,
    cores = min(CHAINS, max(1L, parallel::detectCores() - 1L)),
    control = list(adapt_delta = ADAPT_DELTA, max_treedepth = 12),
    save_pars = brms::save_pars(all = TRUE), seed = 20260828, refresh = 100)
  saveRDS(fit, cache_file)
  saveRDS(fingerprint, meta_file)
  fit
}
message("\nbrms full surface")
fits <- list(none = base)
for (v in c("m12", "m60")) {
  fits[[v]] <- fit_adjusted(VERSIONS[[v]]$column, VERSIONS[[v]]$stem)
}

summarise_fit <- function(fit, column) {
  s <- posterior::summarise_draws(posterior::as_draws_df(fit))
  pick <- function(v) {
    r <- s[s$variable == v, ]
    if (!nrow(r)) return(c(NA_real_, NA_real_, NA_real_))
    c(r$mean, r$q5, r$q95)
  }
  nuts <- brms::nuts_params(fit)
  key <- s[grepl("^b_|^bs_|^sds_|^sd_|^shape", s$variable), ]
  list(sd_country = pick("sd_country__Intercept"), shape = pick("shape"),
       b_nnmr = if (is.null(column)) c(NA_real_, NA_real_, NA_real_) else
         pick(paste0("b_", column)),
       divergent = sum(nuts$Parameter == "divergent__" & nuts$Value > 0),
       max_rhat = max(key$rhat, na.rm = TRUE))
}
af_draws <- function(fit, prevalence, year) {
  frame <- function(p) {
    out <- data.frame(pfpr10 = p / 10, year_c = year - year_center,
                      log_exposure = 0, log_nnmr_z = 0, log_nnmr60_z = 0,
                      country = levels(model_data$country)[1])
    for (name in names(covariates)) out[[name]] <- 0
    out
  }
  high <- brms::posterior_linpred(fit, newdata = frame(prevalence), re_formula = NA)
  low <- brms::posterior_linpred(fit, newdata = frame(rep(AF_REFERENCE, length(prevalence))),
                                 re_formula = NA)
  pmax(1 - exp(-(high - low)), 0)
}

message("\nPSIS-LOO")
loos <- lapply(fits, function(f) {
  l <- brms::loo(f)
  if (sum(l$diagnostics$pareto_k > 0.7) > 0) l <- brms::loo(f, moment_match = TRUE)
  l
})
comparison <- loo::loo_compare(loos)
print(comparison, simplify = FALSE)

rows <- list()
curve_rows <- list()
for (v in names(fits)) {
  fit <- fits[[v]]
  st <- summarise_fit(fit, VERSIONS[[v]]$column)
  anchors <- af_draws(fit, ANCHORS, year_center)
  early <- af_draws(fit, ANCHORS, CONTRAST_YEARS[1])
  late <- af_draws(fit, ANCHORS, CONTRAST_YEARS[2])
  delta <- late - early
  rows[[v]] <- data.frame(
    model = v, label = VERSIONS[[v]]$label,
    elpd_loo = loos[[v]]$estimates["elpd_loo", "Estimate"],
    se_elpd_loo = loos[[v]]$estimates["elpd_loo", "SE"],
    p_loo = loos[[v]]$estimates["p_loo", "Estimate"],
    elpd_diff = comparison[v, "elpd_diff"], se_diff = comparison[v, "se_diff"],
    sd_country = st$sd_country[1], sd_country_q5 = st$sd_country[2],
    sd_country_q95 = st$sd_country[3],
    nb_shape = st$shape[1], nb_shape_q5 = st$shape[2], nb_shape_q95 = st$shape[3],
    neonatal_pct_per_sd = 100 * (exp(st$b_nnmr[1]) - 1),
    neonatal_pct_q5 = 100 * (exp(st$b_nnmr[2]) - 1),
    neonatal_pct_q95 = 100 * (exp(st$b_nnmr[3]) - 1),
    af10 = mean(anchors[, 1]), af30 = mean(anchors[, 2]), af50 = mean(anchors[, 3]),
    af30_lo = unname(stats::quantile(anchors[, 2], 0.025)),
    af30_hi = unname(stats::quantile(anchors[, 2], 0.975)),
    af50_2020_minus_2005 = mean(delta[, 3]),
    af50_diff_lo = unname(stats::quantile(delta[, 3], 0.025)),
    af50_diff_hi = unname(stats::quantile(delta[, 3], 0.975)),
    divergent = st$divergent, max_rhat = st$max_rhat,
    stringsAsFactors = FALSE)
  for (year in c(CONTRAST_YEARS[1], year_center, CONTRAST_YEARS[2])) {
    grid <- seq(1, 60, by = 1)
    values <- af_draws(fit, grid, year)
    curve_rows[[paste(v, year)]] <- data.frame(
      model = v, label = VERSIONS[[v]]$label, year = year, prevalence = grid,
      af = colMeans(values),
      lo = apply(values, 2, stats::quantile, 0.025),
      hi = apply(values, 2, stats::quantile, 0.975), stringsAsFactors = FALSE)
  }
}
brms_table <- do.call(rbind, rows)
rownames(brms_table) <- NULL
write.csv(brms_table, file.path(RESULTS_DIR, "neonatal_covariate_brms.csv"),
          row.names = FALSE)
message("\nFull surface with and without the neonatal terms:")
print(transform(brms_table, elpd_loo = round(elpd_loo, 1), p_loo = round(p_loo, 1),
                elpd_diff = round(elpd_diff, 1), se_diff = round(se_diff, 1),
                sd_country = round(sd_country, 3), nb_shape = round(nb_shape, 1),
                neonatal_pct_per_sd = round(neonatal_pct_per_sd, 1),
                neonatal_pct_q5 = round(neonatal_pct_q5, 1),
                neonatal_pct_q95 = round(neonatal_pct_q95, 1),
                af10 = round(100 * af10, 1), af30 = round(100 * af30, 1),
                af50 = round(100 * af50, 1), af30_lo = round(100 * af30_lo, 1),
                af30_hi = round(100 * af30_hi, 1),
                af50_2020_minus_2005 = round(100 * af50_2020_minus_2005, 1),
                af50_diff_lo = round(100 * af50_diff_lo, 1),
                af50_diff_hi = round(100 * af50_diff_hi, 1))[
                  , c("model", "elpd_loo", "elpd_diff", "se_diff", "p_loo",
                      "sd_country", "nb_shape", "neonatal_pct_per_sd",
                      "neonatal_pct_q5", "neonatal_pct_q95", "af10", "af30",
                      "af30_lo", "af30_hi", "af50", "af50_2020_minus_2005",
                      "af50_diff_lo", "af50_diff_hi", "divergent", "max_rhat")],
      row.names = FALSE)

curves <- do.call(rbind, curve_rows)
curves$label <- factor(curves$label, levels = vapply(VERSIONS, `[[`, "", "label"))
curves$year <- factor(curves$year)
plot <- ggplot2::ggplot(curves, ggplot2::aes(prevalence, 100 * af, colour = year,
                                             fill = year)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * lo, ymax = 100 * hi), alpha = 0.15,
                       colour = NA) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::facet_wrap(~label, nrow = 1) +
  ggplot2::scale_colour_manual(values = c("#8FBEDD", "grey45", "#123F63"), name = "Year") +
  ggplot2::scale_fill_manual(values = c("#8FBEDD", "grey45", "#123F63"), guide = "none") +
  ggplot2::labs(x = "MAP PfPR2-10 (%)",
                y = "Malaria-attributable fraction of\npost-neonatal mortality (%)",
                title = "Full prevalence-by-time surface, with and without neonatal mortality as a covariate",
                subtitle = "Bands are 95% credible intervals; the covariate is held at its mean") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(RESULTS_DIR, "figure20_neonatal_covariate_curves.png"),
                plot, width = 12, height = 4.4, dpi = 200)
message("\nWrote figure19_neonatal_vs_postneonatal.png, ",
        "neonatal_vs_postneonatal_correlation.csv, neonatal_covariate_mgcv.csv, ",
        "neonatal_covariate_brms.csv and figure20_neonatal_covariate_curves.png")

# =============================================================================
# 27_horizon_lag_selection.R
#
# How much do the two timing choices in this analysis actually matter?
#
#   * the mortality reference window (`chmort`'s Period): 12, 24, 36, 48 or 60
#     months before interview;
#   * the lag applied to MAP prevalence: the survey year, or one or two years
#     before it.
#
# Part 1 estimates neonatal and post-neonatal mortality for every survey region
# at all five horizons, WITH confidence intervals.
#
#   DHS.rates::chmort(JK = "Yes") reports jackknife intervals for its own rates,
#   but this analysis defines post-neonatal mortality as U5MR - NNMR, which is
#   not one of them and whose variance cannot be recovered from the two
#   separately (they share the neonatal component and are positively
#   correlated, so adding variances overstates the error).
#
#   So the life table is rebuilt here from per-cluster sums of weighted deaths
#   and exposure. With those cached, deleting a cluster is a subtraction rather
#   than a refit, which makes a delete-one-cluster jackknife cheap enough to run
#   for every region, horizon and outcome - including the difference. The point
#   estimates are validated against chmort itself before anything is reported.
#
# Part 2 plots the horizons against each other with their intervals.
#
# Part 3 refits the five prespecified PfPR/time models at each of the 5 x 3
# horizon-by-lag combinations and records which one AIC selects.
#
# Outputs
#   results/dhs_rebuild/horizon_mortality_estimates.csv
#   results/dhs_rebuild/horizon_estimator_validation.csv
#   results/dhs_rebuild/figure8_horizon_scatter_postneonatal.png
#   results/dhs_rebuild/figure9_horizon_scatter_neonatal.png
#   results/dhs_rebuild/horizon_lag_model_selection.csv
#   results/dhs_rebuild/horizon_lag_selection_table.csv
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("DHS.rates", "ggplot2", "mgcv"))
set.seed(20260828)

HORIZONS <- c(12L, 24L, 36L, 48L, 60L)
LAGS <- 0:2
UNDER5_SEGMENTS <- list(c(0, 1), c(1, 3), c(3, 6), c(6, 12),
                        c(12, 24), c(24, 36), c(36, 48), c(48, 60))
HORIZON_CACHE <- file.path(DATA_DIR, "horizon_mortality_cache")
HORIZON_CSV <- file.path(RESULTS_DIR, "horizon_mortality_estimates.csv")
PFPR_LONG_CSV <- file.path(DERIVED_DIR, "map_pfpr_lagged_long.csv")
dir.create(HORIZON_CACHE, showWarnings = FALSE, recursive = TRUE)

## ---- the life table, resolved to the cluster ---------------------------------
# Weighted deaths and exposure for every (region, age segment, cluster) cell,
# following DHS.rates::chmort's segment definitions exactly.
cluster_sums <- function(data, segments, period, region, cluster) {
  keep <- !is.na(data$v005) & data$v005 != 0
  data <- data[keep, , drop = FALSE]
  region <- region[keep]
  cluster <- cluster[keep]
  weight <- data$v005 / 1e6
  upper_time <- data$v008
  lower_time <- upper_time - period

  regions <- sort(unique(region[!is.na(region) & nzchar(region)]))
  clusters <- sort(unique(cluster))
  dims <- c(length(regions), length(segments), length(clusters))
  deaths_arr <- array(0, dims)
  exposure_arr <- array(0, dims)

  region_index <- match(region, regions)
  cluster_index <- match(cluster, clusters)

  for (s in seq_along(segments)) {
    lower <- segments[[s]][1]
    upper <- segments[[s]][2]
    eligible <- which(data$b7 >= lower | is.na(data$b7))
    if (!length(eligible)) next
    birth <- data$b3[eligible]
    tl <- lower_time[eligible]
    tu <- upper_time[eligible]

    exposure <- rep(NA_real_, length(eligible))
    exposure[birth >= (tl - upper) & birth < (tl - lower)] <- 0.5
    exposure[birth >= (tl - lower) & birth < (tu - upper)] <- 1
    exposure[birth >= (tu - upper) & birth < (tu - lower)] <- 0.5

    died <- data$b7[eligible] >= lower & data$b7[eligible] < upper
    deaths <- rep(NA_real_, length(eligible))
    deaths[birth >= (tl - upper) & birth < (tl - lower) & died] <- 0.5
    deaths[birth >= (tl - lower) & birth < (tu - upper) & died] <- 1
    deaths[birth >= (tu - upper) & birth < (tu - lower) & died] <- 1
    deaths[is.na(data$b7[eligible])] <- 0

    w <- weight[eligible]
    ri <- region_index[eligible]
    ci <- cluster_index[eligible]
    ok <- !is.na(ri) & !is.na(ci)

    add <- function(arr, value) {
      use <- ok & !is.na(value)
      if (!any(use)) return(arr)
      cell <- cbind(ri[use], s, ci[use])
      contribution <- (value * w)[use]
      # rowsum over identical cells, then a single vectorised assignment
      key <- paste(cell[, 1], cell[, 3], sep = "|")
      totals <- rowsum(contribution, key, reorder = FALSE)
      parts <- do.call(rbind, strsplit(rownames(totals), "|", fixed = TRUE))
      idx <- cbind(as.integer(parts[, 1]), s, as.integer(parts[, 2]))
      arr[idx] <- arr[idx] + as.numeric(totals)
      arr
    }
    exposure_arr <- add(exposure_arr, exposure)
    deaths_arr <- add(deaths_arr, deaths)
  }
  list(deaths = deaths_arr, exposure = exposure_arr,
       regions = regions, clusters = clusters)
}

# U5MR, NNMR and post-neonatal (U5MR - NNMR) from one set of segment totals.
rates_from_totals <- function(deaths, exposure) {
  q <- ifelse(exposure > 0, deaths / exposure, NA_real_)
  survive_all <- prod(1 - q, na.rm = TRUE)
  u5mr <- 1000 * abs(1 - survive_all)
  nnmr <- 1000 * ifelse(is.na(q[1]), 0, q[1])
  c(u5mr = u5mr, nnmr = nnmr, postneonatal = u5mr - nnmr)
}

# Delete-one-cluster jackknife within each region.
region_rates_jk <- function(sums) {
  out <- list()
  for (g in seq_along(sums$regions)) {
    deaths <- sums$deaths[g, , , drop = FALSE]
    exposure <- sums$exposure[g, , , drop = FALSE]
    dim(deaths) <- dim(deaths)[-1]
    dim(exposure) <- dim(exposure)[-1]
    total_deaths <- rowSums(deaths)
    total_exposure <- rowSums(exposure)
    if (!is.finite(total_exposure[1]) || total_exposure[1] <= 0) next

    estimate <- rates_from_totals(total_deaths, total_exposure)
    present <- which(colSums(exposure) > 0)
    n <- length(present)
    if (n < 2L) next

    replicates <- vapply(present, function(cl) {
      rates_from_totals(total_deaths - deaths[, cl],
                        total_exposure - exposure[, cl])
    }, numeric(3))
    # pseudo-values: n * full - (n - 1) * leave-one-out
    pseudo <- n * estimate - (n - 1) * replicates
    se <- apply(pseudo, 1, stats::sd) / sqrt(n)

    out[[length(out) + 1L]] <- data.frame(
      regkey = sums$regions[g],
      clusters = n,
      exposure = total_exposure[1],
      u5mr = estimate[["u5mr"]], u5mr_se = se[["u5mr"]],
      nnmr = estimate[["nnmr"]], nnmr_se = se[["nnmr"]],
      postneonatal = estimate[["postneonatal"]],
      postneonatal_se = se[["postneonatal"]],
      stringsAsFactors = FALSE
    )
  }
  if (!length(out)) return(NULL)
  result <- do.call(rbind, out)
  for (outcome in c("u5mr", "nnmr", "postneonatal")) {
    estimate <- result[[outcome]]
    se <- result[[paste0(outcome, "_se")]]
    result[[paste0(outcome, "_lo")]] <- pmax(0, estimate - 1.96 * se)
    result[[paste0(outcome, "_hi")]] <- estimate + 1.96 * se
  }
  result
}

## ---- resolve a survey's regions exactly as the panel does --------------------
survey_region_vector <- function(br, survey_map, survey, registry) {
  boundary_labels <- unique(as.character(survey_map$region))
  region_var <- best_region_var(br, boundary_labels)
  reconciled <- match_region_keys(
    unique(as.character(br[[region_var]])), boundary_labels
  )
  if (length(unique(reconciled$to)) < length(unique(survey_map$regkey))) {
    grouped <- group_regions_to_boundary(br, region_var, boundary_labels,
                                         survey, registry)
    if (!is.null(grouped)) {
      br$region_grouped <- grouped$values
      region_var <- "region_grouped"
      reconciled <- match_region_keys(unique(grouped$values), boundary_labels)
    }
  }
  keys <- rkey(as.character(br[[region_var]]))
  reconciled$to[match(keys, reconciled$from)]
}

## ---- validate against chmort before trusting anything -----------------------
validate_estimator <- function(registry, map, checks = 3L) {
  message("Validating the cluster-resolved life table against DHS.rates::chmort")
  rows <- list()
  candidates <- registry$svkey[registry$svkey %in% map$svkey]
  for (svkey in head(candidates, checks)) {
    survey <- registry[registry$svkey == svkey, , drop = FALSE][1, ]
    br <- tryCatch(readRDS(as.character(survey$local_recode)),
                   error = function(e) NULL)
    if (is.null(br)) next
    survey_map <- map[map$svkey == svkey, , drop = FALSE]
    region <- survey_region_vector(br, survey_map, survey, registry)
    if (all(is.na(region))) { rm(br); gc(FALSE); next }

    br$region_final <- region
    reference <- tryCatch(
      suppressMessages(DHS.rates::chmort(
        br[!is.na(br$region_final), , drop = FALSE],
        Class = "region_final", Period = 60,
        Strata = strata_candidates(br)[1], JK = "Yes"
      )),
      error = function(e) NULL
    )
    if (is.null(reference)) { rm(br); gc(FALSE); next }

    sums <- cluster_sums(br, UNDER5_SEGMENTS, 60, region, br$v021)
    mine <- region_rates_jk(sums)
    if (is.null(mine)) { rm(br); gc(FALSE); next }

    for (label in c("U5MR", "NNMR")) {
      block <- reference[grepl(paste0("^", label), rownames(reference)),
                         c("Class", "R", "SE"), drop = FALSE]
      column <- tolower(label)
      matched <- merge(
        data.frame(Class = block$Class, reference = block$R,
                   reference_se = block$SE, stringsAsFactors = FALSE),
        data.frame(Class = mine$regkey, mine = mine[[column]],
                   mine_se = mine[[paste0(column, "_se")]],
                   exposure = mine$exposure, stringsAsFactors = FALSE),
        by = "Class"
      )
      if (!nrow(matched)) next
      # A clustered sample cannot be more precise than a simple random sample of
      # the same size, so a design effect below 1 marks a broken variance.
      proportion <- matched$mine / 1000
      simple_se <- 1000 * sqrt(proportion * (1 - proportion) / matched$exposure)
      rows[[length(rows) + 1L]] <- data.frame(
        svkey = svkey, outcome = label, regions = nrow(matched),
        max_abs_difference = max(abs(matched$reference - matched$mine)),
        median_deft_chmort = stats::median(matched$reference_se / simple_se),
        median_deft_jackknife = stats::median(matched$mine_se / simple_se),
        regions_chmort_deft_below_1 = sum(matched$reference_se < simple_se),
        regions_jackknife_deft_below_1 = sum(matched$mine_se < simple_se),
        stringsAsFactors = FALSE
      )
    }
    rm(br)
    gc(FALSE)
  }
  if (!length(rows)) return(NULL)
  validation <- do.call(rbind, rows)
  worst <- max(validation$max_abs_difference)
  message(sprintf(
    "  worst disagreement with chmort: %.4f per 1000 across %d comparisons",
    worst, nrow(validation)
  ))
  # chmort rounds its published rates to two decimals, so anything at that
  # scale is rounding rather than a different estimator.
  if (worst > 0.05) {
    warning("The cluster-resolved life table does not reproduce chmort ",
            "(worst difference ", signif(worst, 3), "); Part 1 is unreliable.")
  }
  # Why the intervals here are not chmort's. With Class set, chmort's jackknife
  # returns regional standard errors implying design effects far below 1 - and
  # exactly zero for some regions - which no clustered sample can achieve. It
  # appears to resample the national set of clusters while each regional
  # estimate uses only its own, so most replicates reproduce the full estimate
  # and the variance collapses. The jackknife here resamples within the region.
  message(sprintf(
    "  median design effect: chmort %.2f vs jackknife %.2f (chmort below 1 in %d of %d region-outcomes)",
    stats::median(validation$median_deft_chmort),
    stats::median(validation$median_deft_jackknife),
    sum(validation$regions_chmort_deft_below_1), sum(validation$regions)
  ))
  write.csv(validation,
            file.path(RESULTS_DIR, "horizon_estimator_validation.csv"),
            row.names = FALSE)
  validation
}

## ---- Part 1: mortality at every horizon -------------------------------------
build_horizon_estimates <- function(registry, map) {
  rows <- list()
  for (i in seq_len(nrow(registry))) {
    survey <- registry[i, , drop = FALSE]
    cache <- file.path(HORIZON_CACHE, paste0(survey$svkey, ".rds"))
    if (file.exists(cache)) {
      rows[[survey$svkey]] <- readRDS(cache)
      next
    }
    survey_map <- map[map$svkey == survey$svkey, , drop = FALSE]
    if (!nrow(survey_map)) next
    br <- tryCatch(readRDS(as.character(survey$local_recode)),
                   error = function(e) NULL)
    if (is.null(br)) next
    region <- tryCatch(survey_region_vector(br, survey_map, survey, registry),
                       error = function(e) NULL)
    if (is.null(region) || all(is.na(region))) { rm(br); gc(FALSE); next }

    per_survey <- list()
    for (horizon in HORIZONS) {
      sums <- tryCatch(
        cluster_sums(br, UNDER5_SEGMENTS, horizon, region, br$v021),
        error = function(e) NULL
      )
      if (is.null(sums)) next
      estimates <- tryCatch(region_rates_jk(sums), error = function(e) NULL)
      if (is.null(estimates)) next
      estimates$svkey <- survey$svkey
      estimates$iso3 <- survey$iso3
      estimates$year <- survey$year
      estimates$horizon_months <- horizon
      per_survey[[as.character(horizon)]] <- estimates
    }
    if (length(per_survey)) {
      combined <- do.call(rbind, per_survey)
      rownames(combined) <- NULL
      saveRDS(combined, cache)
      rows[[survey$svkey]] <- combined
    }
    message("  ", survey$svkey, ": ",
            if (length(per_survey)) nrow(per_survey[[1]]) else 0,
            " regions x ", length(per_survey), " horizons")
    rm(br)
    gc(FALSE)
  }
  if (!length(rows)) stop("No horizon estimates could be assembled.")
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
map <- read.csv(MAP_REGION_CSV, stringsAsFactors = FALSE)
registry <- registry[registry$svkey %in% map$svkey, , drop = FALSE]

validate_estimator(registry, map)

message("Estimating mortality at ", length(HORIZONS), " horizons for ",
        nrow(registry), " surveys")
horizon_estimates <- build_horizon_estimates(registry, map)
write.csv(horizon_estimates, HORIZON_CSV, row.names = FALSE)
message("Horizon estimates: ", nrow(horizon_estimates), " region-horizons across ",
        length(unique(horizon_estimates$svkey)), " surveys. Wrote ",
        basename(HORIZON_CSV))

## ---- Part 2: the horizons against each other, with intervals ----------------
# Each shorter window is plotted against the 60-month window the main analysis
# uses, so the question the plot answers is "what would change if we shortened
# the recall period".
horizon_scatter <- function(estimates, outcome, label, file) {
  reference <- estimates[estimates$horizon_months == 60L, , drop = FALSE]
  reference <- reference[, c("svkey", "regkey",
                             outcome,
                             paste0(outcome, "_lo"),
                             paste0(outcome, "_hi"))]
  names(reference) <- c("svkey", "regkey", "ref", "ref_lo", "ref_hi")

  others <- estimates[estimates$horizon_months != 60L, , drop = FALSE]
  others$value <- others[[outcome]]
  others$value_lo <- others[[paste0(outcome, "_lo")]]
  others$value_hi <- others[[paste0(outcome, "_hi")]]
  paired <- merge(
    others[, c("svkey", "regkey", "horizon_months", "value", "value_lo",
               "value_hi")],
    reference, by = c("svkey", "regkey")
  )
  if (!nrow(paired)) return(NULL)
  paired$panel <- factor(
    paired$horizon_months,
    levels = c(12L, 24L, 36L, 48L),
    labels = paste0(c(12L, 24L, 36L, 48L) / 12, "-year window")
  )

  agreement <- do.call(rbind, lapply(split(paired, paired$panel), function(d) {
    data.frame(
      panel = d$panel[1], n = nrow(d),
      correlation = stats::cor(d$value, d$ref, use = "complete.obs"),
      median_ratio = stats::median(d$value / d$ref, na.rm = TRUE),
      median_width_ratio = stats::median(
        (d$value_hi - d$value_lo) / (d$ref_hi - d$ref_lo), na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }))
  print(agreement, row.names = FALSE)

  plot <- ggplot2::ggplot(paired, ggplot2::aes(x = ref, y = value)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey40") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = value_lo, ymax = value_hi),
      colour = "grey70", linewidth = 0.2, alpha = 0.35, width = 0
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = ref_lo, xmax = ref_hi),
      colour = "grey70", linewidth = 0.2, alpha = 0.35, height = 0
    ) +
    ggplot2::geom_point(size = 0.7, alpha = 0.6, colour = "#1D6F8B") +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::labs(
      x = paste0(label, ", 5-year window (per 1000)"),
      y = paste0(label, ", shorter window (per 1000)"),
      title = paste0(label, ": shorter recall windows against the 5-year window"),
      subtitle = paste(
        "One point per survey region; bars are 95% delete-one-cluster",
        "jackknife intervals. Dashed line is equality."
      )
    ) +
    ggplot2::theme_minimal(base_size = 10)
  ggplot2::ggsave(file, plot, width = 11, height = 3.4, dpi = 200)
  message("  wrote ", basename(file))
  agreement
}

message("\nPart 2: plotting the horizons against each other")
postneonatal_agreement <- horizon_scatter(
  horizon_estimates, "postneonatal", "Post-neonatal mortality",
  file.path(RESULTS_DIR, "figure8_horizon_scatter_postneonatal.png")
)
neonatal_agreement <- horizon_scatter(
  horizon_estimates, "nnmr", "Neonatal mortality",
  file.path(RESULTS_DIR, "figure9_horizon_scatter_neonatal.png")
)
agreement <- rbind(
  cbind(outcome = "postneonatal", postneonatal_agreement),
  cbind(outcome = "neonatal", neonatal_agreement)
)
write.csv(agreement, file.path(RESULTS_DIR, "horizon_agreement_summary.csv"),
          row.names = FALSE)

## ---- Part 3: model selection across horizon and lag -------------------------
# The analysis sample is held fixed at the main-sample definition throughout, so
# the only things varying across the 15 cells are the mortality window and the
# prevalence lag - not which regions are in the model.
message("\nPart 3: model selection across ", length(HORIZONS), " horizons x ",
        length(LAGS), " prevalence lags")

analysis <- read_analysis_data()
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
catalog$included_in_main <- as.logical(catalog$included_in_main)
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
analysis$key <- paste(analysis$svkey, analysis$regkey, sep = "|")

if (!file.exists(PFPR_LONG_CSV)) {
  stop("Lagged PfPR panel missing. Run R_dhs/17_sensitivity_timing.R first.")
}
lagged <- read.csv(PFPR_LONG_CSV, stringsAsFactors = FALSE)
lagged$key <- paste(lagged$svkey, lagged$regkey, sep = "|")

horizon_estimates$key <- paste(horizon_estimates$svkey,
                               horizon_estimates$regkey, sep = "|")

cell_results <- list()
for (horizon in HORIZONS) {
  mortality <- horizon_estimates[
    horizon_estimates$horizon_months == horizon, , drop = FALSE
  ]
  for (lag in LAGS) {
    prevalence <- lagged[lagged$lag == lag, c("key", "pfpr2_10")]
    names(prevalence)[2] <- "pfpr_lagged"

    frame <- analysis
    frame$postneonatal_horizon <- mortality$postneonatal[
      match(frame$key, mortality$key)
    ]
    frame$exposure <- mortality$exposure[match(frame$key, mortality$key)]
    frame$pfpr_lagged <- prevalence$pfpr_lagged[
      match(frame$key, prevalence$key)
    ]
    frame$pfpr10 <- frame$pfpr_lagged / 10
    frame <- frame[is.finite(frame$postneonatal_horizon) &
                     frame$postneonatal_horizon > 0 &
                     is.finite(frame$exposure) & frame$exposure > 0 &
                     is.finite(frame$pfpr10), , drop = FALSE]
    if (nrow(frame) < 100L) {
      message("  horizon ", horizon, "m, lag ", lag,
              "y: only ", nrow(frame), " rows; skipped")
      next
    }

    fits <- lapply(MODEL_SPECS$specification, function(specification) {
      tryCatch(
        fit_ridge_gam(frame, "postneonatal_horizon", catalog,
                      specification = specification, method = "ML"),
        error = function(e) NULL
      )
    })
    names(fits) <- MODEL_SPECS$specification
    fits <- fits[!vapply(fits, is.null, logical(1))]
    if (!length(fits)) {
      message("  horizon ", horizon, "m, lag ", lag, "y: no model converged")
      next
    }

    aic <- vapply(fits, function(f) stats::AIC(f$model), numeric(1))
    ordered <- sort(aic)
    linear <- tryCatch(
      fit_ridge_gam(frame, "postneonatal_horizon", catalog,
                    specification = "linear_no_interaction", method = "ML"),
      error = function(e) NULL
    )
    summary_row <- if (is.null(linear)) {
      list(pct = NA_real_, p = NA_real_)
    } else {
      coefficients <- summary(linear$model)$p.table
      row <- grep("^pfpr10$", rownames(coefficients))
      if (!length(row)) list(pct = NA_real_, p = NA_real_) else
        list(pct = 100 * (exp(coefficients[row, 1]) - 1),
             p = coefficients[row, 4])
    }

    cell_results[[length(cell_results) + 1L]] <- data.frame(
      horizon_months = horizon, prevalence_lag_years = lag,
      n = nrow(frame), countries = length(unique(frame$iso3)),
      selected = names(ordered)[1],
      runner_up = if (length(ordered) > 1) names(ordered)[2] else NA_character_,
      delta_aic_to_runner_up = if (length(ordered) > 1)
        unname(ordered[2] - ordered[1]) else NA_real_,
      aic_linear_no_interaction = unname(aic["linear_no_interaction"]),
      aic_selected = unname(ordered[1]),
      linear_pct_change_per_10 = summary_row$pct,
      linear_pfpr_p = summary_row$p,
      stringsAsFactors = FALSE
    )
    message(sprintf("  horizon %2dm, lag %dy: n=%4d -> %s (next %s by %.1f AIC)",
                    horizon, lag, nrow(frame), names(ordered)[1],
                    if (length(ordered) > 1) names(ordered)[2] else "-",
                    if (length(ordered) > 1) ordered[2] - ordered[1] else NA))
  }
}

selection <- do.call(rbind, cell_results)
write.csv(selection, file.path(RESULTS_DIR, "horizon_lag_model_selection.csv"),
          row.names = FALSE)

# The requested table: preferred model as a function of horizon and lag.
table_wide <- reshape(
  selection[, c("horizon_months", "prevalence_lag_years", "selected")],
  idvar = "horizon_months", timevar = "prevalence_lag_years",
  direction = "wide"
)
names(table_wide) <- sub("^selected\\.", "lag_", names(table_wide))
write.csv(table_wide,
          file.path(RESULTS_DIR, "horizon_lag_selection_table.csv"),
          row.names = FALSE)

message("\nPreferred specification by mortality horizon and prevalence lag:")
print(table_wide, row.names = FALSE)
message("\nWrote horizon_lag_model_selection.csv and horizon_lag_selection_table.csv")

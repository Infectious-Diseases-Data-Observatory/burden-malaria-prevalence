# =============================================================================
# 32_brms_three_way.R — one model for prevalence x time x region.
#
# The 2 x 2 subgroup fits split the panel into four cells of 211 to 347
# region-years and fitted each independently. That throws away the fact that the
# covariate effects, the country random effects and the dispersion are shared,
# and it leaves the era-by-region contrast as an eyeball comparison across four
# separate fits rather than a quantity with an interval.
#
# The three-way can instead be written as ONE model. In mgcv's smooth syntax a
# tensor product between two continuous variables can be split by a factor:
#
#   t2(pfpr10, year_c, by = region_group)
#
# gives a separate prevalence-by-time surface for each region, each with its own
# smoothing parameters, while everything else in the model is shared. brms
# passes smooth terms through to mgcv, so this works there too - t2 rather than
# te for the reason set out in 30_brms_tensor_model.R. The factor main effect
# must be included alongside it, because a `by` smooth is centred within each
# level and cannot carry the level's overall height.
#
# Note this is a better-posed version of the three-way than the 2 x 2: era is a
# hard cut at 2013, whereas year_c enters continuously here, so nothing hinges
# on where the cut is placed.
#
# The pooled surface is refitted on identical rows so the two can be compared by
# LOO, which is the question worth asking: does letting the surface differ by
# region buy anything?
#
# Outputs
#   results/dhs_rebuild/brms_three_way_loo.csv
#   results/dhs_rebuild/brms_three_way_surface.csv
#   results/dhs_rebuild/figure16_brms_three_way.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("brms", "ggplot2", "loo"))

CHAINS <- 4L
ITERATIONS <- 2000L
WARMUP <- 1000L
ADAPT_DELTA <- 0.99
RIDGE_PRIOR_SCALE <- 1
ANCHORS <- c(10, 30, 50)
CACHE <- file.path(DATA_DIR, "brms_three_way_cache")
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
ridge <- make_ridge_matrix(analysis, catalog, bundle$preprocessing)
covariates <- as.data.frame(ridge$matrix)
names(covariates) <- make.names(colnames(ridge$matrix))

model_data <- data.frame(
  deaths = round(analysis$postneonatal_mortality / 1000 * analysis$exposure),
  log_exposure = log(analysis$exposure),
  pfpr10 = analysis$pfpr10,
  year_c = analysis$year_c,
  country = factor(analysis$iso3),
  region_group = factor(ifelse(analysis$iso3 %in% WEST, "West",
                               "Central & East"),
                        levels = c("Central & East", "West")),
  stringsAsFactors = FALSE
)
model_data <- cbind(model_data, covariates)
message("Three-way model: ", nrow(model_data), " region-years; ",
        paste(sprintf("%s %d", levels(model_data$region_group),
                      table(model_data$region_group)), collapse = ", "))

covariate_terms <- paste(names(covariates), collapse = " + ")
priors <- c(
  brms::prior_string(sprintf("normal(0, %g)", RIDGE_PRIOR_SCALE), class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "sds"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(gamma(0.01, 0.01), class = "shape")
)

formulas <- list(
  pooled = brms::bf(paste(
    "deaths ~ t2(pfpr10, year_c, k = c(6, 6)) +", covariate_terms,
    "+ (1 | country) + offset(log_exposure)")),
  # A separate surface per region. region_group enters as a main effect because
  # a `by` smooth is centred within each level.
  by_region = brms::bf(paste(
    "deaths ~ t2(pfpr10, year_c, k = c(6, 6), by = region_group) +",
    "region_group +", covariate_terms,
    "+ (1 | country) + offset(log_exposure)"))
)

fit_or_load <- function(name) {
  cache_file <- file.path(CACHE, paste0(name, ".rds"))
  meta_file <- file.path(CACHE, paste0(name, "_meta.rds"))
  fingerprint <- list(rows = nrow(model_data),
                      countries = nlevels(model_data$country),
                      deaths = sum(model_data$deaths),
                      adapt_delta = ADAPT_DELTA, iterations = ITERATIONS)
  if (file.exists(cache_file) && file.exists(meta_file) &&
        isTRUE(all.equal(readRDS(meta_file), fingerprint))) {
    message("  ", name, ": reusing cached fit")
    return(readRDS(cache_file))
  }
  message("  ", name, ": sampling")
  fit <- brms::brm(
    formula = formulas[[name]], data = model_data,
    family = brms::negbinomial(), prior = priors,
    chains = CHAINS, iter = ITERATIONS, warmup = WARMUP,
    cores = min(CHAINS, max(1L, parallel::detectCores() - 1L)),
    control = list(adapt_delta = ADAPT_DELTA, max_treedepth = 14),
    seed = 20260828, refresh = 0
  )
  saveRDS(fit, cache_file)
  saveRDS(fingerprint, meta_file)
  fit
}

fits <- lapply(names(formulas), fit_or_load)
names(fits) <- names(formulas)

## ---- diagnostics -------------------------------------------------------------
for (name in names(fits)) {
  parameters <- posterior::summarise_draws(posterior::as_draws_df(fits[[name]]))
  key <- parameters[grepl("^b_|^sds_|^sd_|^shape", parameters$variable), ]
  divergent <- sum(brms::nuts_params(fits[[name]])$Parameter == "divergent__" &
                     brms::nuts_params(fits[[name]])$Value > 0)
  message(sprintf("  %-10s divergences %3d | max R-hat %.4f | min ESS %.0f",
                  name, divergent, max(key$rhat, na.rm = TRUE),
                  min(key$ess_bulk, na.rm = TRUE)))
}

## ---- does splitting the surface by region buy anything? ---------------------
loo_results <- lapply(fits, function(f) brms::loo(f))
comparison <- loo::loo_compare(loo_results)
print(comparison)
loo_table <- as.data.frame(comparison)
loo_table$model <- rownames(loo_table)
write.csv(loo_table, file.path(RESULTS_DIR, "brms_three_way_loo.csv"),
          row.names = FALSE)

## ---- the two surfaces --------------------------------------------------------
af_surface <- function(fit, prevalence, year, region) {
  frame <- function(p) {
    out <- data.frame(pfpr10 = p / 10, year_c = year, log_exposure = 0,
                      country = levels(model_data$country)[1],
                      region_group = factor(region,
                                            levels = levels(model_data$region_group)),
                      stringsAsFactors = FALSE)
    for (name in names(covariates)) out[[name]] <- 0
    out
  }
  high <- brms::posterior_linpred(fit, newdata = frame(prevalence),
                                  re_formula = NA)
  low <- brms::posterior_linpred(
    fit, newdata = frame(rep(AF_REFERENCE, length(prevalence))),
    re_formula = NA)
  pmax(1 - exp(-(high - low)), 0)
}

year_centre <- unique(bundle$year_center)[1]
shown_years <- c(2004, 2013, 2022)
rows <- list()
support_window <- function(region, year, half_width = 4) {
  in_window <- model_data$region_group == region &
    abs(model_data$year_c + year_centre - year) <= half_width
  model_data$pfpr10[in_window] * 10
}
for (region in levels(model_data$region_group)) {
  for (year in shown_years) {
    # Only where that region has data AT THAT TIME.
    observed <- support_window(region, year)
    if (length(observed) < 20L) next
    grid <- seq(max(AF_REFERENCE, stats::quantile(observed, 0.02)),
                stats::quantile(observed, 0.98), length.out = 35)
    values <- af_surface(fits$by_region, grid, year - year_centre, region)
    rows[[length(rows) + 1L]] <- data.frame(
      region = region, year = year, prevalence = grid,
      af = colMeans(values),
      lo = apply(values, 2, stats::quantile, 0.025),
      hi = apply(values, 2, stats::quantile, 0.975),
      stringsAsFactors = FALSE
    )
  }
}
surface <- do.call(rbind, rows)
write.csv(surface, file.path(RESULTS_DIR, "brms_three_way_surface.csv"),
          row.names = FALSE)

## ---- the contrast, with an interval -----------------------------------------
# The quantity the 2 x 2 could only gesture at: West minus Central & East at
# each anchor and year, as a posterior rather than a comparison of point
# estimates from separate fits.
contrast_rows <- list()
for (year in shown_years) {
  west <- af_surface(fits$by_region, ANCHORS, year - year_centre, "West")
  east <- af_surface(fits$by_region, ANCHORS, year - year_centre,
                     "Central & East")
  difference <- west - east
  for (j in seq_along(ANCHORS)) {
    column <- difference[, j]
    near <- function(region) {
      sum(model_data$region_group == region &
            abs(model_data$year_c + year_centre - year) <= 4 &
            model_data$pfpr10 * 10 >= ANCHORS[j])
    }
    contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
      year = year, prevalence = ANCHORS[j],
      difference = mean(column),
      lo = unname(stats::quantile(column, 0.025)),
      hi = unname(stats::quantile(column, 0.975)),
      probability_west_higher = mean(column > 0),
      west_regions_at_or_above = near("West"),
      east_regions_at_or_above = near("Central & East"),
      stringsAsFactors = FALSE
    )
  }
}
contrast <- do.call(rbind, contrast_rows)
message("\nWest minus Central & East, attributable fraction (percentage points).")
message("The last two columns count regions within four years of the stated ",
        "year at or above that prevalence: where they are near zero the ",
        "contrast is extrapolation, not evidence.")
print(transform(contrast, difference = round(100 * difference, 1),
                lo = round(100 * lo, 1), hi = round(100 * hi, 1),
                probability_west_higher = round(probability_west_higher, 2)),
      row.names = FALSE)
write.csv(contrast, file.path(RESULTS_DIR, "brms_three_way_contrast.csv"),
          row.names = FALSE)

## ---- figure ------------------------------------------------------------------
surface$region <- factor(surface$region, levels = c("Central & East", "West"))
surface$period <- factor(surface$year, levels = shown_years)
period_colours <- setNames(c("#8FBEDD", "#4E86B0", "#123F63"),
                           as.character(shown_years))
plot <- ggplot2::ggplot(surface,
                        ggplot2::aes(prevalence, 100 * af, colour = period,
                                     fill = period, linetype = region,
                                     group = interaction(region, period))) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * lo, ymax = 100 * hi),
                       alpha = 0.12, colour = NA) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::facet_wrap(~region, nrow = 1) +
  ggplot2::scale_colour_manual(values = period_colours, name = "Year") +
  ggplot2::scale_fill_manual(values = period_colours, guide = "none") +
  ggplot2::scale_linetype_manual(
    values = c("Central & East" = "solid", "West" = "dashed"), name = "Region") +
  ggplot2::guides(colour = ggplot2::guide_legend(order = 1),
                  linetype = ggplot2::guide_legend(
                    order = 2, override.aes = list(colour = "grey30"))) +
  ggplot2::labs(
    x = "MAP PfPR2-10 (%)",
    y = "Malaria-attributable fraction of\npost-neonatal mortality (%)",
    title = "Prevalence by time by region, as one model",
    subtitle = paste("t2(pfpr10, year_c, by = region_group): a separate surface",
                     "per region, everything else shared.\nCurves span only the",
                     "prevalence each region observes; bands are 95% credible",
                     "intervals")
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(RESULTS_DIR, "figure16_brms_three_way.png"), plot,
                width = 10, height = 4.8, dpi = 200)

message("\nWrote brms_three_way_loo.csv, brms_three_way_surface.csv, ",
        "brms_three_way_contrast.csv and figure16_brms_three_way.png")

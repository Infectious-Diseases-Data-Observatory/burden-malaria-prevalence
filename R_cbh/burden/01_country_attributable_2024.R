#!/usr/bin/env Rscript
# National 2024 PfPR-to-zero contrasts from the shared-time mortality model.
# No refitting, downloads, or modification of historical model/data outputs.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
for (pkg in c("mgcv", "countrycode", "terra", "sf", "exactextractr"))
  stopifnot(requireNamespace(pkg, quietly = TRUE))
spec <- cbh_trial_spec()
stopifnot(!spec$time_by_age)
out <- file.path("results/cbh", spec$id, "country_burden_2024")
private <- file.path("data/derived_cbh/models", spec$id)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
export <- list.files("data", pattern = "2026-09-09 10-58-22[.]csv$", full.names = TRUE)
stopifnot(length(export) == 1L)
x <- cbh_read_csv(export)
x <- x[x$Year == 2024 & x$Sex == "Both" & x$Condition == "All causes" & x$Measure == "Deaths", ]
x$iso3 <- countrycode::countrycode(x$Location, "country.name", "iso3c", warn = FALSE)
stopifnot(!anyNA(x$iso3), all(is.finite(x$Value)), all(x$Value > 0))
cbh_unique(x, c("iso3", "Age", "Unit"), "2024 IHME export")
countries <- unique(x[c("iso3", "Location")]); countries <- countries[order(countries$iso3), ]
cbh_atomic_csv(x[c("iso3", "Location", "Year", "Age", "Unit", "Value", "Lower", "Upper")],
               file.path(out, "ihme_source_2024.csv"))

# Population counts per raster cell = density x cell area. Fixed GPW 2020
# population geography; exact polygon overlap; only pixels with MAP coverage.
# Historical national extracts weighted by density without cell area.
map_path <- "data/pfpr_2to10_africa_2024.tif"
pop_path <- "data/pop/gpw_v4_population_density_rev11_2020_2.5m.tif"
boundary_path <- "data/africa_admin0.rds"
p <- terra::rast(map_path)
stopifnot(terra::nlyr(p) == 1, terra::is.lonlat(p))
b <- sf::st_make_valid(readRDS(boundary_path))
b <- b[b$iso %in% countries$iso3, ]
density <- terra::resample(terra::crop(terra::rast(pop_path), p), p, method = "bilinear")
weights <- terra::ifel(density >= 0, density, NA) * terra::cellSize(p, unit = "km")
covered <- terra::mask(weights, p)
message("Extracting national population-weighted 2024 PfPR.")
layers <- c(covered, p * covered, weights)
names(layers) <- c("covered", "numerator", "all_pop")
agg <- exactextractr::exact_extract(layers, b, "sum", progress = FALSE)
den <- agg$sum.covered; num <- agg$sum.numerator; all_pop <- agg$sum.all_pop
stopifnot(length(den) == nrow(b), length(num) == nrow(b), length(all_pop) == nrow(b))
pr <- data.frame(iso3 = b$iso, pfpr_pct = 100 * num / den,
  population_weight_on_map = den, population_weight_within_raster = all_pop,
  map_population_coverage_within_raster = den / all_pop)
cbh_unique(pr, "iso3", "National MAP extraction")
pr <- merge(countries, pr, by = "iso3", all.x = TRUE, sort = TRUE)
pr$pfpr_pct[!is.finite(pr$pfpr_pct)] <- NA_real_
stopifnot(all(pr$pfpr_pct[!is.na(pr$pfpr_pct)] >= 0 & pr$pfpr_pct[!is.na(pr$pfpr_pct)] <= 100))
old <- cbh_read_csv("data/pfpr_by_country_year.csv"); old <- old[old$year == 2024, ]
pr$legacy_density_weighted_pfpr_pct <- old$pfpr_pct[match(pr$iso3, old$iso3)]
pr$status <- ifelse(is.na(pr$pfpr_pct), "missing_MAP_no_estimate", "estimated")
pr$map_coverage_below_95pct <- pr$map_population_coverage_within_raster < .95
pr$exposure_year <- 2024L; pr$population_weight_year <- 2020L
cbh_atomic_csv(pr, file.path(out, "national_pfpr_2024.csv"))

# Six disjoint source ages. Early and late neonatal receive the same <1m HR.
leaf <- c("0 to 6 days (early neonatal)", "7 to 27 days (late neonatal)",
          "1 to 5 months", "6-11 months", "12 to 23 months", "2 to 4")
source_rows <- do.call(rbind, lapply(countries$iso3, function(iso) {
  z <- x[x$iso3 == iso, ]
  value <- function(ages, unit, field = "Value") {
    ii <- match(ages, z$Age[z$Unit == unit]); stopifnot(!anyNA(ii))
    z[[field]][which(z$Unit == unit)[ii]]
  }
  deaths <- value(leaf, "Number"); rate <- value(leaf, "Rate (per 100,000)")
  stopifnot(abs(sum(deaths) - value("Under 5", "Number")) < 1e-5 * sum(deaths),
    abs(sum(deaths[1:4]) - value("Under 1", "Number")) < 1e-5 * sum(deaths[1:4]))
  data.frame(iso3 = iso, source_age = leaf, deaths = deaths, rate_per100000 = rate,
    implied_person_years = deaths / (rate / 1e5),
    ihme_deaths_lower = value(leaf, "Number", "Lower"),
    ihme_deaths_upper = value(leaf, "Number", "Upper"),
    ihme_rate_lower = value(leaf, "Rate (per 100,000)", "Lower"),
    ihme_rate_upper = value(leaf, "Rate (per 100,000)", "Upper"))
}))
cbh_atomic_csv(source_rows, file.path(out, "ihme_disjoint_age_inputs.csv"))
ages <- cbh_config()$age_bands$age_band
stopifnot(identical(ages, c("<1", "1-5", "6-11", "12-23", "24-35", "36-47", "48-59")))
bands <- do.call(rbind, lapply(countries$iso3, function(iso) {
  z <- source_rows[source_rows$iso3 == iso, ]; stopifnot(identical(z$source_age, leaf))
  deaths <- c(sum(z$deaths[1:2]), z$deaths[3:5], rep(z$deaths[6] / 3, 3))
  exposure <- c(sum(z$implied_person_years[1:2]), z$implied_person_years[3:5], rep(z$implied_person_years[6] / 3, 3))
  data.frame(iso3 = iso, age_band = ages, ihme_deaths = deaths,
    implied_person_years = exposure, ihme_rate_per100000 = deaths / exposure * 1e5,
    age_mapping = c("0-27 days mapped to model <1 completed month", rep("direct age match", 3),
      rep("IHME 2-4 rate shared; deaths and person-time divided equally", 3)))
}))
bands <- merge(bands, pr, by = "iso3", all.x = TRUE)
bands <- bands[order(bands$iso3, match(bands$age_band, ages)), ]; rownames(bands) <- NULL

# Read each existing fit once, retaining only the PfPR coefficient block and
# spline objects in a compact private cache. No DHS rows are exported.
fd <- cbh_read_csv(file.path("results/cbh", spec$id, "multiple_imputation/fit_diagnostics.csv"))
stopifnot(nrow(fd) == 10L, all(fd$converged))
paths <- file.path(private, sprintf("incidence_draw_%03d.rds", fd$draw))
sig <- cbh_hash(vapply(paths, cbh_file_hash, character(1)))
cache_path <- file.path(private, "burden_pfpr_components.rds")
cache <- if (file.exists(cache_path)) readRDS(cache_path) else NULL
if (is.null(cache) || !identical(cache$signature, sig)) {
  pieces <- lapply(seq_along(paths), function(i) {
    message("Extracting PfPR effects from saved mortality fit ", i, "/", length(paths))
    obj <- readRDS(paths[i]); f <- obj$fit
    stopifnot(isTRUE(f$converged), all(is.finite(coef(f))))
    ss <- Filter(function(s) identical(s$term, "pfpr_pct"), f$smooth)
    stopifnot(length(ss) == 7L, identical(vapply(ss, function(s) s$by.level, ""), ages))
    idx <- unlist(lapply(ss, function(s) s$first.para:s$last.para))
    support <- do.call(rbind, lapply(ages, function(a) {
      pp <- f$model$pfpr_pct[f$model$age_band == a]
      data.frame(age_band = a, minimum = min(pp), maximum = max(pp),
        p025 = unname(quantile(pp, .025)), p975 = unname(quantile(pp, .975)))
    }))
    list(draw = fd$draw[i], coef = coef(f)[idx], covariance = f$Vp[idx, idx], smooths = ss,
      support = support, base_signature = obj$base_input_signature,
      external_signature = obj$external_signature)
  })
  cache <- list(signature = sig, pieces = pieces)
  cbh_atomic_rds(cache, cache_path)
}
stopifnot(length(unique(vapply(cache$pieces, function(z) z$base_signature, ""))) == 1L)
pfpr_matrix <- function(piece, data) {
  do.call(cbind, lapply(piece$smooths, function(s) {
    # PredictMat supplies the centred basis; explicitly apply factor-by mask.
    mm <- mgcv::PredictMat(s, data.frame(pfpr_pct = data$pfpr_pct,
                                        age_band = factor(data$age_band, levels = ages)))
    mm * as.numeric(data$age_band == s$by.level)
  }))
}
valid <- bands[!is.na(bands$pfpr_pct), ]
individual <- lapply(cache$pieces, function(piece) {
  ref <- valid; ref$pfpr_pct <- 0
  L <- pfpr_matrix(piece, ref) - pfpr_matrix(piece, valid)
  variance <- rowSums((L %*% piece$covariance) * L)
  stopifnot(all(variance >= -1e-9))
  # Verify the compact basis against prior full-prediction-matrix contrasts.
  check <- data.frame(age_band = ages, pfpr_pct = 40); to <- check; to$pfpr_pct <- 20
  hr <- exp(drop((pfpr_matrix(piece, to) - pfpr_matrix(piece, check)) %*% piece$coef))
  saved <- cbh_read_csv(file.path("results/cbh", spec$id, "multiple_imputation/individual_pfpr_40_to_20.csv"))
  saved <- saved[saved$draw == piece$draw, ]; saved <- saved[match(ages, saved$age_band), ]
  stopifnot(max(abs(hr - saved$hazard_ratio)) < 1e-7)
  data.frame(iso3 = valid$iso3, age_band = valid$age_band, draw = piece$draw,
    log_hr_zero_vs_current = drop(L %*% piece$coef), variance = pmax(variance, 0))
})
individual <- do.call(rbind, individual)
cbh_atomic_csv(individual, file.path(out, "individual_log_hazard_contrasts.csv"))
pooled <- do.call(rbind, lapply(split(individual, paste(individual$iso3, individual$age_band)), function(z) {
  M <- nrow(z); U <- mean(z$variance); B <- var(z$log_hr_zero_vs_current); T <- U + (1 + 1/M) * B
  df <- if (B > 1e-20) (M-1)*(1+U/((1+1/M)*B))^2 else Inf
  mu <- mean(z$log_hr_zero_vs_current); half <- qt(.975, df) * sqrt(T)
  data.frame(iso3 = z$iso3[1], age_band = z$age_band[1], log_hr_zero_vs_current = mu,
    log_hr_se = sqrt(T), hr_zero_vs_current = exp(mu), hr_lower_95 = exp(mu-half),
    hr_upper_95 = exp(mu+half), monte_carlo_se_log_hr = sqrt(B/M), imputations = M)
}))
res <- merge(bands, pooled, by = c("iso3", "age_band"), all.x = TRUE)
res <- res[order(res$iso3, match(res$age_band, ages)), ]; rownames(res) <- NULL
res$attributable_fraction <- 1 - res$hr_zero_vs_current
res$af_lower_95 <- 1 - res$hr_upper_95; res$af_upper_95 <- 1 - res$hr_lower_95
for (unit in c("rate_per100000", "deaths")) {
  baseline <- res[[paste0("ihme_", unit)]]
  res[[paste0("counterfactual_", unit)]] <- baseline * res$hr_zero_vs_current
  res[[paste0("attributable_", unit)]] <- baseline * res$attributable_fraction
  res[[paste0("attributable_", unit, "_lower_95")]] <- baseline * res$af_lower_95
  res[[paste0("attributable_", unit, "_upper_95")]] <- baseline * res$af_upper_95
  stopifnot(max(abs(baseline - res[[paste0("counterfactual_", unit)]] - res[[paste0("attributable_", unit)]]), na.rm = TRUE) < 1e-7)
}
support <- cache$pieces[[1]]$support
res$zero_below_observed_support <- 0 < support$minimum[match(res$age_band, ages)]
res$current_pfpr_outside_central95 <- res$pfpr_pct < support$p025[match(res$age_band, ages)] |
  res$pfpr_pct > support$p975[match(res$age_band, ages)]
res$negative_attributable_estimate <- res$attributable_fraction < 0
cbh_atomic_csv(res, file.path(out, "country_age_attributable_2024.csv"))
totals <- do.call(rbind, lapply(split(res, res$iso3), function(z) data.frame(
  iso3 = z$iso3[1], country = z$Location[1], pfpr_pct = z$pfpr_pct[1], status = z$status[1],
  map_population_coverage_within_raster = z$map_population_coverage_within_raster[1],
  map_coverage_below_95pct = z$map_coverage_below_95pct[1],
  ihme_under5_deaths = sum(z$ihme_deaths), counterfactual_under5_deaths = sum(z$counterfactual_deaths),
  attributable_under5_deaths = sum(z$attributable_deaths),
  attributable_fraction = sum(z$attributable_deaths) / sum(z$ihme_deaths))))
cbh_atomic_csv(totals, file.path(out, "country_totals_2024.csv"))
cbh_atomic_csv(support, file.path(out, "model_pfpr_support.csv"))
files <- c(export, map_path, pop_path, boundary_path, paths,
           "R_cbh/burden/01_country_attributable_2024.R")
cbh_atomic_csv(data.frame(file = files, md5 = vapply(files, cbh_file_hash, "")), file.path(out, "provenance.csv"))
message("Saved country/age estimates: ", sum(totals$status == "estimated"), "/", nrow(totals), " countries.")
print(totals[totals$iso3 == "COD", ], row.names = FALSE)

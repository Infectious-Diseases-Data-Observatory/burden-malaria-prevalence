# Sensitivity: refit the six age-band person-time models with the improved-water and
# improved-sanitation covariates replaced by DHS StatCompiler values.
setwd("/Users/jameswatson/Documents/Claude Projects/MIS:DHS malaria prevalence")
source("R_dhs/00_config.R")
S <- commandArgs(trailingOnly = TRUE)[1]
BAM_THREADS <- max(1L, min(8L, parallel::detectCores() - 2L))
bundle <- readRDS(MODEL_BUNDLE_RDS)
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE); catalog$included_in_main <- as.logical(catalog$included_in_main)
model_data <- read.csv(file.path(RESULTS_DIR, "person_time_model_data.csv"), stringsAsFactors = FALSE)
api <- read.csv(file.path(S, "wash_api_covariates.csv"), stringsAsFactors = FALSE)
key <- paste(model_data$svkey, model_data$regkey); idx <- match(key, paste(api$svkey, api$regkey))
cat("model rows matched to API covariates:", sum(!is.na(idx)), "of", nrow(model_data), "\n")
prep <- function(md) {
  md$country <- factor(md$iso3); md$survey <- factor(md$svkey); md$window_f <- factor(md$window)
  ridge <- make_ridge_matrix(md, catalog, bundle$preprocessing)
  md$G <- ridge$matrix
  list(data = md, penalty = list(G = list(diag(ncol(ridge$matrix)))))
}
fit_band <- function(p, band) {
  data <- p$data[p$data$age6b == band, ]
  data$segment_f <- droplevels(factor(data$segment))
  parts <- c("pfpr10", if (nlevels(data$segment_f) > 1) "segment_f", "window_f", "s(year_c, k = 8)", "G",
             "s(country, bs = 're')", "s(survey, bs = 're')", "offset(log_pm)")
  mgcv::bam(as.formula(paste("deaths_eff ~", paste(parts, collapse = " + "))), family = mgcv::nb(), method = "fREML",
            paraPen = p$penalty, data = data, discrete = TRUE, nthreads = BAM_THREADS)
}
variants <- list(
  original = model_data,
  api_water = within(model_data, imp_water_analysis <- ifelse(is.na(idx), imp_water_analysis, api$imp_water_api[idx])),
  api_water_sanit = within(model_data, { imp_water_analysis <- ifelse(is.na(idx), imp_water_analysis, api$imp_water_api[idx])
                                          imp_sanit_analysis <- ifelse(is.na(idx), imp_sanit_analysis, api$imp_sanit_api[idx]) }),
  drop_water_sanit = model_data)
catalog_drop <- catalog; catalog_drop$included_in_main[catalog_drop$variable %in% c("imp_water", "imp_sanit")] <- FALSE
rows <- list()
for (v in names(variants)) {
  cat_use <- if (v == "drop_water_sanit") catalog_drop else catalog
  md <- variants[[v]]; md$country <- factor(md$iso3); md$survey <- factor(md$svkey); md$window_f <- factor(md$window)
  ridge <- make_ridge_matrix(md, cat_use, if (v == "drop_water_sanit") NULL else bundle$preprocessing)
  md$G <- ridge$matrix; p <- list(data = md, penalty = list(G = list(diag(ncol(ridge$matrix)))))
  for (band in AGE6B) {
    fit <- fit_band(p, band); b <- coef(fit)["pfpr10"]; se <- sqrt(vcov(fit)["pfpr10", "pfpr10"])
    gi <- grep("^G", names(coef(fit))); gn <- sub("^G", "", names(coef(fit))[gi])
    w <- if ("imp_water" %in% gn) coef(fit)[gi][gn == "imp_water"] else NA; s <- if ("imp_sanit" %in% gn) coef(fit)[gi][gn == "imp_sanit"] else NA
    rows[[paste(v, band)]] <- data.frame(variant = v, band = band, pct_per_10 = 100 * (exp(b) - 1), lo = 100 * (exp(b - 1.96 * se) - 1), hi = 100 * (exp(b + 1.96 * se) - 1),
                                         beta_water = 100 * (exp(w) - 1), beta_sanit = 100 * (exp(s) - 1), stringsAsFactors = FALSE)
  }
  cat("done", v, "\n")
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(S, "wash_sensitivity.csv"), row.names = FALSE)
res$band <- factor(res$band, levels = AGE6B)
cat("\nPfPR effect per +10 points (%), by band and covariate variant:\n")
print(reshape(res[, c("variant", "band", "pct_per_10")], idvar = "band", timevar = "variant", direction = "wide"), row.names = FALSE, digits = 3)
cat("\nImproved-water coefficient (% per SD of logit) by variant:\n")
print(reshape(res[, c("variant", "band", "beta_water")], idvar = "band", timevar = "variant", direction = "wide"), row.names = FALSE, digits = 3)
cat("\nImproved-sanitation coefficient (% per SD of logit) by variant:\n")
print(reshape(res[, c("variant", "band", "beta_sanit")], idvar = "band", timevar = "variant", direction = "wide"), row.names = FALSE, digits = 3)

# =============================================================================
# 09_validate_reproduction.R
# Compare the rebuilt DHS/MIS results with the saved legacy headline results.
#
# This is a migration check, not an assertion of bit-for-bit identity. The
# rebuilt analysis applies a coherent missingness rule and formally selects
# among four PfPR/time specifications. It should nevertheless preserve:
#   * the primary analysis sample size and country count;
#   * the direction and significance classification of the three outcomes;
#   * closely matching ridge-linear effects and nonlinear AF anchors.
#
# The script uses only aggregate summaries and saved fitted models. It never
# reads individual DHS records.
# =============================================================================

source("R_dhs/00_config.R")
required_packages("mgcv")

LEGACY_MODEL_RDS <- file.path(REPO_ROOT, "results", "penalized_models.rds")
LEGACY_LINEAR_CSV <- file.path(REPO_ROOT, "results", "paper_fig2_negcontrol.csv")
NEW_LINEAR_CSV <- file.path(RESULTS_DIR, "linear_effects_by_outcome.csv")
NEW_AF_CSV <- file.path(RESULTS_DIR, "attributable_fraction_anchors.csv")

inputs <- c(
  LEGACY_MODEL_RDS,
  LEGACY_LINEAR_CSV,
  MODEL_BUNDLE_RDS,
  NEW_LINEAR_CSV,
  NEW_AF_CSV
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop(
    "Reproduction inputs are missing:\n",
    paste(" -", missing_inputs, collapse = "\n")
  )
}

legacy_models <- readRDS(LEGACY_MODEL_RDS)
new_bundle <- readRDS(MODEL_BUNDLE_RDS)
legacy_linear <- read.csv(LEGACY_LINEAR_CSV, stringsAsFactors = FALSE)
new_linear <- read.csv(NEW_LINEAR_CSV, stringsAsFactors = FALSE)
new_af <- read.csv(NEW_AF_CSV, stringsAsFactors = FALSE)

legacy_linear$outcome_key <- ifelse(
  grepl("^All under-5", legacy_linear$outcome),
  "all_under_5",
  ifelse(
    grepl("^1mo-5y", legacy_linear$outcome),
    "postneonatal",
    "neonatal"
  )
)
new_linear$outcome_key <- c(
  postneonatal_mortality = "postneonatal",
  u5mr = "all_under_5",
  nnmr = "neonatal"
)[new_linear$outcome]

legacy_population_af <- function(model, prevalence, reference = AF_REFERENCE) {
  n <- length(prevalence)
  country_level <- levels(model$model$country)[1]
  g_names <- colnames(model$model$G)
  make_newdata <- function(pfpr10) {
    out <- data.frame(
      pfpr10 = pfpr10,
      year_c = rep(0, length.out = n),
      exposure = 1,
      country = factor(
        rep(country_level, n),
        levels = levels(model$model$country)
      )
    )
    out$G <- matrix(0, nrow = n, ncol = length(g_names))
    colnames(out$G) <- g_names
    out
  }
  high <- make_newdata(prevalence / 10)
  low <- make_newdata(rep(reference / 10, n))
  Xh <- predict(model, high, type = "lpmatrix")
  Xl <- predict(model, low, type = "lpmatrix")
  random_columns <- grep(
    "^s\\(country\\)|^s\\(country,pfpr10\\)",
    colnames(Xh)
  )
  if (length(random_columns)) {
    Xh[, random_columns] <- 0
    Xl[, random_columns] <- 0
  }
  dX <- Xh - Xl
  delta <- as.numeric(dX %*% coef(model))
  100 * pmax(1 - exp(-delta), 0)
}

anchors <- c(10, 30, 50)
legacy_af <- legacy_population_af(legacy_models$pn, anchors)
rebuilt_af <- new_af$af[
  match(
    paste(anchors, "postneonatal"),
    paste(new_af$prevalence, new_af$outcome)
  )
] * 100

add_check <- function(
    category,
    metric,
    legacy_value,
    rebuilt_value,
    tolerance,
    units) {
  difference <- abs(rebuilt_value - legacy_value)
  data.frame(
    category = category,
    metric = metric,
    legacy_value = legacy_value,
    rebuilt_value = rebuilt_value,
    absolute_difference = difference,
    tolerance = tolerance,
    units = units,
    pass = is.finite(difference) && difference <= tolerance
  )
}

new_primary <- new_bundle$primary_fits$postneonatal$model
checks <- rbind(
  add_check(
    "sample",
    "analysis rows",
    nrow(legacy_models$pn$model),
    nrow(new_primary$model),
    0,
    "survey-region-years"
  ),
  add_check(
    "sample",
    "countries",
    nlevels(legacy_models$pn$model$country),
    nlevels(new_primary$model$country),
    0,
    "countries"
  )
)

for (outcome in c("postneonatal", "all_under_5", "neonatal")) {
  old <- legacy_linear[legacy_linear$outcome_key == outcome, , drop = FALSE]
  new <- new_linear[new_linear$outcome_key == outcome, , drop = FALSE]
  checks <- rbind(
    checks,
    add_check(
      "ridge-linear effect",
      paste0(outcome, ": percent change per 10 PfPR points"),
      old$pct,
      new$pct_change_per_10,
      1,
      "percentage points"
    ),
    add_check(
      "inference",
      paste0(outcome, ": p < 0.05"),
      as.numeric(old$p < 0.05),
      as.numeric(new$pfpr_p < 0.05),
      0,
      "indicator"
    )
  )
}

for (i in seq_along(anchors)) {
  checks <- rbind(
    checks,
    add_check(
      "nonlinear attributable fraction",
      paste0("postneonatal AF at ", anchors[i], "% PfPR"),
      legacy_af[i],
      rebuilt_af[i],
      3,
      "percentage points"
    )
  )
}

write.csv(
  checks,
  file.path(RESULTS_DIR, "reproduction_check.csv"),
  row.names = FALSE
)

all_pass <- all(checks$pass)
summary_lines <- c(
  "DHS/MIS rebuild reproduction check",
  paste0("Status: ", if (all_pass) "PASS" else "REVIEW REQUIRED"),
  paste0("Checks passed: ", sum(checks$pass), "/", nrow(checks)),
  paste0(
    "Legacy headline specification: ridge spline without PfPR-time interaction"
  ),
  paste0(
    "Rebuilt AIC-selected specification: ",
    new_bundle$selected_specification
  ),
  paste0(
    "Rebuilt reference year for AF estimates: ",
    new_bundle$year_center
  ),
  "",
  "Interpretation:",
  paste0(
    "- This confirms reproduction of the main findings, not identical ",
    "coefficients from identical code."
  ),
  paste0(
    "- The rebuilt pipeline harmonises covariate eligibility/imputation and ",
    "compares all four prespecified PfPR/time models on one sample."
  ),
  paste0(
    "- Legacy R files should be archived only after the generated data, plots, ",
    "and validation table have been reviewed."
  )
)
writeLines(
  summary_lines,
  file.path(RESULTS_DIR, "reproduction_summary.txt")
)

print(checks, row.names = FALSE)
message(
  "\nReproduction check: ",
  if (all_pass) "PASS" else "REVIEW REQUIRED",
  " (", sum(checks$pass), "/", nrow(checks), " checks)."
)
if (!all_pass) {
  stop("One or more reproduction checks exceeded the prespecified tolerance.")
}

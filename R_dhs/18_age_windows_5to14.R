# =============================================================================
# 18_age_windows_5to14.R — does PfPR2-10 predict mortality at ages 5-14?
#
# DHS.rates::chmort stops at 60 months: its age segments end at (48,60), so a
# death above age 5 is recorded in the birth history but contributes NOTHING to
# the life table (verified: erasing every 5y+ death leaves U5MR unchanged to five
# decimal places). To reach older ages this script implements chmort_ext(), which
# reproduces chmort's synthetic-cohort period algorithm for ARBITRARY age
# segments: the same component death probabilities, the same 0.5/1/0.5 edge
# weights on exposure and deaths at the window boundaries, the same v005/1e6
# sample weighting, combined as 1000 * (1 - prod(1 - q_i)).
#
# It is validated against chmort itself (NNMR, IMR, CMR, U5MR must agree to
# chmort's two-decimal printed precision) before being used above 60 months.
#
# WHY THIS MATTERS. Malaria biology predicts a specific AGE PATTERN: near-zero in
# neonates (maternal antibodies, fetal haemoglobin), strong post-neonatally, and
# weak at 5-14 once disease-controlling immunity is acquired. A 5-14 gradient as
# steep as the post-neonatal one would instead indicate socioeconomic confounding
# that the covariate block has failed to remove, since poverty raises all-cause
# mortality at every age. So this is both an age-specificity check on the main
# result and a bound on how much malaria burden the under-5 ceiling omits.
#
# TWO LIMITATIONS, both of which bias older-age mortality DOWNWARD and must be
# stated with any 5-14 estimate:
#   * maternal-age truncation — DHS interviews women 15-49, so children born to
#     women who have since aged out of the sample are missing, and that loss
#     grows with the child's age;
#   * maternal-survival selection — children whose mother has died are absent
#     entirely, and orphans have elevated mortality.
# Deaths are also sparse at these ages (a sizeable minority of regions have fewer
# than five observed 5-14 deaths), so the spline is unstable and the LINEAR
# gradient is the defensible summary. Note also that 10q5 is a ten-year
# cumulative probability whereas the post-neonatal outcome spans five years; the
# comparison below is of RELATIVE gradients, not of identical estimands.
#
# Outputs (results/dhs_rebuild/):
#   sensitivity_age_windows_summary.csv     gradient by age window (+ robustness)
#   sensitivity_age_windows.png             forest of the gradients
#   data/derived_dhs/mortality_5to14_by_region.csv   the 5-14 life-table output
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2", "scales"))

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
catalog <- bundle$catalog
catalog$included_in_main <- as.logical(catalog$included_in_main)
analysis <- read_analysis_data()
analysis$k <- paste(analysis$svkey, analysis$regkey, sep = "|")
main <- as.logical(analysis$main_sample)

OLDER_CSV <- file.path(DERIVED_DIR, "mortality_5to14_by_region.csv")
OLDER_CACHE <- file.path(DATA_DIR, "mort514_cache")
UNDER5_SEGMENTS <- list(c(0, 1), c(1, 3), c(3, 6), c(6, 12),
                        c(12, 24), c(24, 36), c(36, 48), c(48, 60))
annual_segments <- function(from, to, by = 12) Map(c, seq(from, to - by, by), seq(from + by, to, by))

## ---- the extended life table -----------------------------------------------
# Component death probabilities for arbitrary age segments (months), following
# DHS.rates::chmort exactly. `group` optionally splits by region.
chmort_ext <- function(data, ageseg, period = 60, group = NULL) {
  data <- data[!is.na(data$v005) & data$v005 != 0, , drop = FALSE]
  data$rweight <- data$v005 / 1e6
  data$tu <- data$v008                 # reference period ends at interview
  data$tl <- data$tu - period
  if (is.null(group)) data$group <- "all" else data$group <- group
  groups <- sort(unique(data$group))
  q <- matrix(NA_real_, length(groups), length(ageseg), dimnames = list(groups, NULL))
  at_risk <- observed <- setNames(numeric(length(groups)), groups)
  for (i in seq_along(ageseg)) {
    lower <- ageseg[[i]][1]; upper <- ageseg[[i]][2]
    # children who either survived (b7 missing) or died at or after this age
    seg <- data[which(data$b7 >= lower | is.na(data$b7)), , drop = FALSE]
    exposure <- rep(NA_real_, nrow(seg))
    exposure[seg$b3 >= (seg$tl - upper) & seg$b3 < (seg$tl - lower)] <- 0.5
    exposure[seg$b3 >= (seg$tl - lower) & seg$b3 < (seg$tu - upper)] <- 1
    exposure[seg$b3 >= (seg$tu - upper) & seg$b3 < (seg$tu - lower)] <- 0.5
    deaths <- rep(NA_real_, nrow(seg))
    in_segment <- seg$b7 >= lower & seg$b7 < upper
    deaths[seg$b3 >= (seg$tl - upper) & seg$b3 < (seg$tl - lower) & in_segment] <- 0.5
    deaths[seg$b3 >= (seg$tl - lower) & seg$b3 < (seg$tu - upper) & in_segment] <- 1
    deaths[seg$b3 >= (seg$tu - upper) & seg$b3 < (seg$tu - lower) & in_segment] <- 1
    deaths[is.na(seg$b7)] <- 0
    weighted_exposure <- tapply(exposure * seg$rweight, seg$group, sum, na.rm = TRUE)
    weighted_deaths <- tapply(deaths * seg$rweight, seg$group, sum, na.rm = TRUE)
    counted <- tapply(as.numeric(!is.na(deaths) & deaths > 0), seg$group, sum, na.rm = TRUE)
    present <- names(weighted_exposure)
    q[present, i] <- as.numeric(weighted_deaths[present]) / as.numeric(weighted_exposure[present])
    if (i == 1) at_risk[present] <- as.numeric(weighted_exposure[present])
    observed[present] <- observed[present] +
      ifelse(is.na(counted[present]), 0, as.numeric(counted[present]))
  }
  list(q = q, at_risk = at_risk, observed = observed, groups = groups)
}
combine_q <- function(q) 1000 * abs(1 - prod(1 - q, na.rm = TRUE))

## ---- locate the raw Births Recodes -----------------------------------------
recode_paths <- function() {
  c(list.files(file.path(DATA_DIR, "dhs"), pattern = "rds$", full.names = TRUE),
    list.files(path.expand("~/.rdhs_cache"), pattern = "rds$",
               full.names = TRUE, recursive = TRUE))
}

## ---- validate chmort_ext against chmort itself ------------------------------
validate_against_chmort <- function() {
  if (!requireNamespace("DHS.rates", quietly = TRUE)) return(invisible(NULL))
  paths <- recode_paths()
  if (!length(paths)) return(invisible(NULL))
  br <- readRDS(paths[1])
  if (!all(c("v005", "v008", "b3", "b7") %in% names(br))) return(invisible(NULL))
  fit <- chmort_ext(br[, c("v005", "v008", "b3", "b7")], UNDER5_SEGMENTS, period = 60)
  q <- fit$q[1, ]
  mine <- c(NNMR = combine_q(q[1]), IMR = combine_q(q[1:4]),
            CMR = combine_q(q[5:8]), U5MR = combine_q(q[1:8]))
  reference <- suppressMessages(DHS.rates::chmort(br, Period = 60))
  cat("validation of chmort_ext against DHS.rates::chmort (under-5 segments):\n")
  for (nm in names(mine)) {
    cat(sprintf("  %-5s extended %8.4f   chmort %8.4f   difference %.4f\n",
                nm, mine[[nm]], reference[nm, "R"], mine[[nm]] - reference[nm, "R"]))
  }
  worst <- max(abs(mine - reference[names(mine), "R"]))
  if (worst > 0.01) {
    warning("chmort_ext does not reproduce chmort (worst difference ", signif(worst, 3),
            "); the 5-14 estimates below are NOT trustworthy.", call. = FALSE, immediate. = TRUE)
  } else {
    cat("  agreement is within chmort's two-decimal printed precision.\n")
  }
  invisible(NULL)
}

## ---- build the per-region 5-14 life table -----------------------------------
build_older_mortality <- function() {
  required_packages(c("DHS.rates", "haven"))
  dir.create(OLDER_CACHE, showWarnings = FALSE)
  registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
  registry <- registry[registry$svkey %in% unique(analysis$svkey[main]), , drop = FALSE]
  paths <- recode_paths()
  stems <- toupper(sub("\\.rds$", "", basename(paths)))
  locate <- function(no_extension) {
    hit <- which(stems == toupper(no_extension))
    if (length(hit)) paths[hit[1]] else NA_character_
  }
  segments_5to14 <- annual_segments(60, 180)
  segments_5to9  <- annual_segments(60, 120)
  segments_10to14 <- annual_segments(120, 180)

  extract <- function(svrow) {
    cache <- file.path(OLDER_CACHE, paste0(svrow$svkey, ".rds"))
    if (file.exists(cache)) return(readRDS(cache))
    path <- locate(svrow$no_extension)
    if (is.na(path)) stop("recode unavailable")
    br <- readRDS(path)
    if (!all(c("v024", "v005", "v008", "b3", "b7") %in% names(br))) stop("recode lacks required columns")
    d <- data.frame(v005 = br$v005, v008 = br$v008, b3 = br$b3, b7 = br$b7)
    regions <- rkey(as.character(haven::as_factor(br$v024)))
    keep <- nzchar(regions) & !is.na(d$v005) & d$v005 != 0
    d <- d[keep, , drop = FALSE]; regions <- regions[keep]
    all14 <- chmort_ext(d, segments_5to14, 60, group = regions)
    a59  <- chmort_ext(d, segments_5to9,   60, group = regions)
    a1014 <- chmort_ext(d, segments_10to14, 60, group = regions)
    out <- data.frame(
      svkey = svrow$svkey, regkey = all14$groups,
      q514 = apply(all14$q, 1, combine_q),
      q59 = apply(a59$q, 1, combine_q)[all14$groups],
      q1014 = apply(a1014$q, 1, combine_q)[all14$groups],
      exposure514 = as.numeric(all14$at_risk[all14$groups]),
      deaths514 = as.numeric(all14$observed[all14$groups]),
      stringsAsFactors = FALSE)
    out <- out[is.finite(out$q514) & out$q514 > 0 &
                 is.finite(out$exposure514) & out$exposure514 > 0, , drop = FALSE]
    saveRDS(out, cache)
    out
  }
  message("Building 5-14 life tables for ", nrow(registry), " surveys ...")
  rows <- list()
  for (i in seq_len(nrow(registry))) {
    svrow <- registry[i, , drop = FALSE]
    result <- tryCatch(extract(svrow), error = function(e) {
      message("  skip ", svrow$svkey, ": ", conditionMessage(e)); NULL })
    if (!is.null(result) && nrow(result)) rows[[svrow$svkey]] <- result
  }
  long <- do.call(rbind, rows)
  rownames(long) <- NULL
  write.csv(long, OLDER_CSV, row.names = FALSE)
  long
}

validate_against_chmort()
older <- if (file.exists(OLDER_CSV)) {
  read.csv(OLDER_CSV, stringsAsFactors = FALSE)
} else {
  tryCatch(build_older_mortality(), error = function(e) {
    message("5-14 life tables unavailable (", conditionMessage(e),
            "); raw DHS Births Recodes are required."); NULL })
}
if (is.null(older)) {
  message("Skipping the age-window comparison.")
  quit(save = "no", status = 0)
}
older$k <- paste(older$svkey, older$regkey, sep = "|")
idx <- match(analysis$k, older$k)
analysis$q514 <- older$q514[idx]
analysis$exposure514 <- older$exposure514[idx]
analysis$deaths514 <- older$deaths514[idx]

## ---- age-window comparison on a common sample -------------------------------
# every row must support all three outcomes, so differences are about age and not
# about which region-years survive each definition
common <- main &
  is.finite(analysis$nnmr) & analysis$nnmr > 0 &
  is.finite(analysis$postneonatal_mortality) & analysis$postneonatal_mortality > 0 &
  is.finite(analysis$q514) & analysis$q514 > 0 &
  is.finite(analysis$exposure514) & analysis$exposure514 > 0 &
  is.finite(analysis$exposure) & analysis$exposure > 0 &
  is.finite(analysis$pfpr10)
cat(sprintf("\ncommon sample across all three age windows: %d region-years\n", sum(common)))
cat(sprintf("5-14: median 10q5 %.1f per 1000; median observed deaths per region %.0f; %.0f%% of regions have fewer than 5\n",
            median(analysis$q514[common]), median(analysis$deaths514[common]),
            100 * mean(analysis$deaths514[common] < 5)))

gradient <- function(rows, outcome, exposure_col = NULL) {
  dd <- analysis[rows, , drop = FALSE]
  if (!is.null(exposure_col)) dd$exposure <- dd[[exposure_col]]
  dd <- dd[is.finite(dd[[outcome]]) & dd[[outcome]] > 0 &
             is.finite(dd$exposure) & dd$exposure > 0, , drop = FALSE]
  linear <- model_summary_row(fit_ridge_gam(
    dd, outcome, catalog, "linear_no_interaction",
    method = "REML", preprocessing = bundle$preprocessing))
  data.frame(n = nrow(dd), pct_change_per_10 = linear$pct_change_per_10,
             lo = linear$pct_change_lo, hi = linear$pct_change_hi, p_value = linear$pfpr_p)
}

rows <- list(
  cbind(label = "Neonatal (<1 month)", group = "Age window",
        gradient(common, "nnmr")),
  cbind(label = "Post-neonatal (1-59 months)", group = "Age window",
        gradient(common, "postneonatal_mortality")),
  cbind(label = "Ages 5-14 years", group = "Age window",
        gradient(common, "q514", "exposure514")),
  cbind(label = "5-14: restrict PfPR 5-40%", group = "5-14 robustness",
        gradient(common & analysis$pfpr2_10 >= 5 & analysis$pfpr2_10 <= 40, "q514", "exposure514")),
  cbind(label = "5-14: exclude PfPR < 5%", group = "5-14 robustness",
        gradient(common & analysis$pfpr2_10 >= 5, "q514", "exposure514")),
  cbind(label = "5-14: regions with >=10 deaths", group = "5-14 robustness",
        gradient(common & analysis$deaths514 >= 10, "q514", "exposure514"))
)
res <- do.call(rbind, rows)
write.csv(res, file.path(RESULTS_DIR, "sensitivity_age_windows_summary.csv"), row.names = FALSE)
cat("\n=== PfPR2-10 gradient by age window (identical model specification) ===\n")
print(within(res, { pct_change_per_10 <- round(pct_change_per_10, 1); lo <- round(lo, 1)
  hi <- round(hi, 1); p_value <- signif(p_value, 2) }), row.names = FALSE)
cat("\nThe spline is NOT used for the 5-14 window: deaths are sparse there and the\n",
    "fitted smooth is non-monotone, so the linear gradient is the defensible summary.\n", sep = "")

## ---- forest figure ----------------------------------------------------------
res$label <- factor(res$label, levels = rev(res$label))
res$group <- factor(res$group, levels = c("Age window", "5-14 robustness"))
POST <- "#08519c"
plot <- ggplot2::ggplot(res, ggplot2::aes(pct_change_per_10, label)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi, colour = group),
                         orientation = "y", width = 0.25, linewidth = 0.7) +
  ggplot2::geom_point(ggplot2::aes(colour = group), size = 3) +
  ggplot2::scale_colour_manual(values = c("Age window" = POST, "5-14 robustness" = "grey45"),
                               guide = "none") +
  ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  ggplot2::labs(
    x = expression("Change in mortality per +10 " * italic(Pf) * "PR"[2-10] * " points (%, 95% CI)"),
    y = NULL,
    title = "Age-specificity of the prevalence-mortality association",
    subtitle = "Same ridge NB-GAM and covariates; common sample of survey-region-years") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_blank(),
                 strip.placement = "outside",
                 strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
                 strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 9),
                 plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(RESULTS_DIR, "sensitivity_age_windows.png"), plot,
                width = 9, height = 4.8, dpi = 320, bg = "white")
cat("saved: sensitivity_age_windows_summary.csv + sensitivity_age_windows.png\n")

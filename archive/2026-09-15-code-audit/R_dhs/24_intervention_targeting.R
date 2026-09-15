# =============================================================================
# 24_intervention_targeting.R — did malaria interventions reach high-transmission
# regions first, and did that differ between West and East/Southern Africa?
#
# MOTIVATION. Scripts 20-23 narrowed the explanations for the West Africa
# steepening / East & Southern flattening: it is not HIV, not the meningitis belt,
# not MAP measurement error, not country composition, and not a GROWING deprivation
# entanglement (script 23 shows that entanglement is high but flat in West Africa).
# The remaining candidate from the literature is differential TARGETING - whether
# nets and treatment reached the highest-burden places first or last.
#
# INDICATORS, from the recodes rather than any published aggregate:
#   ITN  children aged 0-4 who slept under an insecticide-treated net last night
#        (PR recode hml12 = "only treated" or "both treated and untreated"),
#        household-weighted by hv005.
#   ACT  among children with fever in the preceding two weeks (BR recode h22),
#        the share given an artemisinin combination (ml13e), weighted by v005.
# Both exist only in surveys carrying the malaria module, so coverage is largely
# 2005 onwards and roughly half the analysis panel carries them.
#
# METHOD, as in script 23: both PfPR and coverage are residualised on survey fixed
# effects, so the correlation asks whether, WITHIN one country at one moment, the
# higher-transmission regions received more or less coverage. A positive value
# means high-burden regions were favoured.
#
# FINDING, which splits by intervention:
#   ITN targeting is POSITIVE in West Africa (+0.25 overall, rising to +0.48 by
#   2018-2024) and around zero or negative in East & Southern Africa (-0.05
#   overall, -0.30 in 2012-2017). This CONTRADICTS the suggestion that nets went
#   high-burden-first in the east and last in the west - the signs are reversed.
#   ACT targeting shows the opposite and matches the mortality divergence: strongly
#   positive in East & Southern Africa (+0.47 overall, +0.68 and +0.70 in the last
#   two eras) but weak in West Africa (+0.15), where coverage also stayed very low
#   (10.6% of febrile children against 12.4%, and only 3.9% in Central Africa).
#
# INTERPRETATION. Treatment, not prevention, is where the regions diverge, and in
# the direction that fits: East & Southern Africa steered ACTs towards its
# high-transmission regions, which would narrow the mortality gap and flatten the
# gradient (as observed), while West Africa did not, leaving high-transmission
# regions with the same poor treatment access throughout (gradient steepens, as
# observed). The improving West African ITN targeting works AGAINST its steepening,
# so the account is not clean.
#
# Outputs (results/dhs_rebuild/):
#   intervention_targeting.csv       within-survey correlations by region and era
# Also writes data/derived_dhs/dhs_itn_act_by_region.csv (cached per survey).
# =============================================================================
source("R_dhs/00_config.R")
required_packages("haven")

ITN_ACT_CSV <- file.path(DERIVED_DIR, "dhs_itn_act_by_region.csv")
ITN_ACT_CACHE <- file.path(DATA_DIR, "itn_act_cache")
WEST <- c("BEN","BFA","CIV","GHA","GIN","GMB","LBR","MLI","MRT","NER","NGA","SEN","SLE","TGO")
CENTRAL <- c("AGO","CMR","COD","COG","GAB","TCD")
EAST <- c("BDI","COM","ETH","KEN","MDG","MOZ","MWI","NAM","RWA","SWZ","UGA","ZMB","ZWE")

analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]

build_itn_act <- function() {
  dir.create(ITN_ACT_CACHE, showWarnings = FALSE)
  registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
  registry <- registry[registry$svkey %in% unique(analysis$svkey), , drop = FALSE]
  find_files <- function(pattern) {
    f <- c(list.files(file.path(DATA_DIR, "dhs"), pattern, full.names = TRUE),
           list.files(path.expand("~/.rdhs_cache"), pattern, full.names = TRUE, recursive = TRUE))
    f[!duplicated(toupper(basename(f)))]
  }
  prs <- find_files("PR.*rds$"); brs <- find_files("BR.*rds$")
  stem_of <- function(p) toupper(sub("\\.rds$", "", basename(p)))
  # a survey's PR and BR files share the country prefix and phase characters
  key6 <- function(s) paste0(substr(s, 1, 2), substr(s, 5, 6))
  weighted_share <- function(flag, weight, group) {
    num <- tapply(as.numeric(flag) * weight, group, sum, na.rm = TRUE)
    den <- tapply(weight, group, sum, na.rm = TRUE)
    100 * num[names(den)] / den
  }
  message("Extracting ITN and ACT coverage for ", nrow(registry), " surveys ...")
  rows <- list()
  for (i in seq_len(nrow(registry))) {
    svkey <- registry$svkey[i]; nx <- toupper(registry$no_extension[i])
    cache <- file.path(ITN_ACT_CACHE, paste0(svkey, ".rds"))
    if (file.exists(cache)) { rows[[svkey]] <- readRDS(cache); next }
    out <- NULL
    pr <- prs[key6(stem_of(prs)) == key6(nx)]
    if (length(pr)) {
      x <- tryCatch(readRDS(pr[1]), error = function(e) NULL)
      if (!is.null(x) && all(c("hml12", "hml16", "hv005") %in% names(x))) {
        age <- suppressWarnings(as.integer(as.character(haven::as_factor(x$hml16))))
        net <- as.character(haven::as_factor(x$hml12))
        weight <- as.numeric(x$hv005) / 1e6
        for (v in intersect(c("hv024", "shstate", "sstate"), names(x))) {
          region <- rkey(as.character(haven::as_factor(x[[v]])))
          ok <- is.finite(age) & age <= 4 & is.finite(weight) & weight > 0 &
            nzchar(region) & !is.na(net)
          if (!any(ok)) next
          share <- weighted_share(grepl("only treated|both treated", net[ok]),
                                  weight[ok], region[ok])
          out <- rbind(out, data.frame(svkey = svkey, regkey = names(share),
            itn_pct = as.numeric(share),
            itn_n = as.numeric(tapply(rep(1, sum(ok)), region[ok], sum)),
            stringsAsFactors = FALSE))
        }
      }
    }
    br <- brs[stem_of(brs) == nx]
    if (length(br)) {
      b <- tryCatch(readRDS(br[1]), error = function(e) NULL)
      if (!is.null(b) && all(c("h22", "v005") %in% names(b))) {
        fever <- as.character(haven::as_factor(b$h22)) == "yes"
        act <- if ("ml13e" %in% names(b)) as.character(haven::as_factor(b$ml13e)) == "yes" else NA
        weight <- as.numeric(b$v005) / 1e6
        for (v in intersect(c("v024", "shstate", "sstate"), names(b))) {
          region <- rkey(as.character(haven::as_factor(b[[v]])))
          ok <- !is.na(fever) & fever & is.finite(weight) & weight > 0 & nzchar(region)
          if (!any(ok) || all(is.na(act))) next
          got <- ifelse(is.na(act[ok]), FALSE, act[ok])
          share <- weighted_share(got, weight[ok], region[ok])
          add <- data.frame(svkey = svkey, regkey = names(share),
            act_pct = as.numeric(share),
            fever_n = as.numeric(tapply(rep(1, sum(ok)), region[ok], sum)),
            stringsAsFactors = FALSE)
          out <- if (is.null(out)) add else merge(out, add, by = c("svkey", "regkey"), all = TRUE)
        }
      }
    }
    if (!is.null(out) && nrow(out)) {
      out <- out[!duplicated(out$regkey), , drop = FALSE]
      saveRDS(out, cache); rows[[svkey]] <- out
    }
  }
  long <- do.call(rbind, lapply(rows, function(z) {
    for (cc in c("itn_pct", "itn_n", "act_pct", "fever_n")) if (!cc %in% names(z)) z[[cc]] <- NA_real_
    z[, c("svkey", "regkey", "itn_pct", "itn_n", "act_pct", "fever_n")]
  }))
  rownames(long) <- NULL
  write.csv(long, ITN_ACT_CSV, row.names = FALSE)
  long
}
coverage <- if (file.exists(ITN_ACT_CSV)) read.csv(ITN_ACT_CSV, stringsAsFactors = FALSE) else
  tryCatch(build_itn_act(), error = function(e) {
    message("ITN/ACT extraction unavailable (", conditionMessage(e), ")."); NULL })
if (is.null(coverage)) { message("Skipping."); quit(save = "no", status = 0) }

analysis$k <- paste(analysis$svkey, analysis$regkey, sep = "|")
coverage$k <- paste(coverage$svkey, coverage$regkey, sep = "|")
idx <- match(analysis$k, coverage$k)
analysis$itn <- coverage$itn_pct[idx]
analysis$act <- coverage$act_pct[idx]
analysis$fever_n <- coverage$fever_n[idx]
analysis$region_group <- ifelse(analysis$iso3 %in% WEST, "West Africa",
  ifelse(analysis$iso3 %in% CENTRAL, "Central Africa",
    ifelse(analysis$iso3 %in% EAST, "East & Southern", NA)))
analysis <- analysis[!is.na(analysis$region_group), , drop = FALSE]
analysis$era <- cut(analysis$year, c(1999.5, 2011.5, 2017.5, 2024.5),
                    labels = c("2005-2011", "2012-2017", "2018-2024"))
cat(sprintf("matched: ITN %d, ACT %d of %d analysis region-years\n",
            sum(is.finite(analysis$itn)), sum(is.finite(analysis$act)), nrow(analysis)))

residualise <- function(v, survey) {
  out <- rep(NA_real_, length(v)); ok <- is.finite(v)
  if (sum(ok) < 10) return(out)
  out[ok] <- residuals(lm(v[ok] ~ factor(survey[ok]))); out
}
within_cor <- function(z, column) {
  a <- residualise(z$pfpr2_10, z$svkey); b <- residualise(z[[column]], z$svkey)
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 20) return(c(NA_real_, NA_real_))
  c(cor(a[ok], b[ok]), sum(ok))
}
res <- list()
for (column in c("itn", "act")) {
  for (g in c("West Africa", "Central Africa", "East & Southern")) {
    z <- analysis[analysis$region_group == g & is.finite(analysis[[column]]), , drop = FALSE]
    if (nrow(z) < 20) next
    for (e in c("ALL", levels(analysis$era))) {
      y <- if (e == "ALL") z else z[z$era == e, , drop = FALSE]
      if (nrow(y) < 25) next
      r <- within_cor(y, column)
      if (is.na(r[1])) next
      res[[paste(column, g, e)]] <- data.frame(indicator = toupper(column),
        region_group = g, era = e, n = r[2],
        mean_coverage_pct = mean(y[[column]]), correlation = r[1])
    }
  }
}
targeting <- do.call(rbind, res)
write.csv(targeting, file.path(RESULTS_DIR, "intervention_targeting.csv"), row.names = FALSE)
cat("\n=== Within-survey correlation of PfPR with intervention coverage ===\n")
cat("    positive = higher-transmission regions received MORE coverage\n\n")
print(within(targeting, { mean_coverage_pct <- round(mean_coverage_pct, 1)
  correlation <- round(correlation, 3) }), row.names = FALSE)
cat("\nmedian regional fever sample (the ACT denominator):\n")
print(round(tapply(analysis$fever_n, analysis$region_group, median, na.rm = TRUE)))

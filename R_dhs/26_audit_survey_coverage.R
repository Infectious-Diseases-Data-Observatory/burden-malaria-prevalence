# =============================================================================
# 26_audit_survey_coverage.R
#
# Audit the survey panel actually analysed against the full registry of
# eligible DHS/MIS Births Recodes.
#
# The analysis panel is built with `03_build_analysis_dataset.R
# --from-legacy-aggregate`, which inherits its survey list from the legacy
# aggregate `results/component2_region_data_full.csv`. That aggregate is smaller
# than the registry, so some eligible surveys never reach the models. This
# script quantifies the gap and, where possible, attributes each omission to a
# cause:
#
#   no_map_extraction  the per-survey MAP regional extraction was never cached,
#                      so the survey has no regional prevalence to merge on
#   region_key_mismatch the recode's v024 labels and the boundary's DHSREGEN
#                      do not reconcile under rkey(), so the regional join
#                      drops rows (fully, when the overlap is zero)
#   boundary_granularity the boundary carries a different number of units from
#                      the recode (e.g. zones versus regions)
#   no_map_series      MAP publishes no national prevalence series for the
#                      country at all (malaria-free; a legitimate exclusion)
#
# Outputs
#   results/survey_coverage_audit.csv        one row per eligible survey
#   results/survey_coverage_audit_by_country.csv
# =============================================================================

source("R_dhs/00_config.R")
suppressPackageStartupMessages({
  library(countrycode)
  library(haven)
})

registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
aggregate_panel <- read.csv(
  file.path(REPO_ROOT, "results", "component2_region_data_full.csv"),
  stringsAsFactors = FALSE
)
analysis <- read_analysis_data()

registry$country <- countrycode(registry$iso3, "iso3c", "country.name",
                                warn = FALSE)
registry$in_aggregate <- registry$svkey %in% aggregate_panel$svkey
registry$in_analysis <- registry$svkey %in% analysis$svkey
registry$in_main_sample <- registry$svkey %in%
  analysis$svkey[as.logical(analysis$main_sample)]

# --- inputs each omitted survey would need -----------------------------------
map_cached <- sub("\\.rds$", "", list.files(MAP_REGION_DIR, pattern = "rds$"))
registry$map_extraction_cached <- registry$svkey %in% map_cached
registry$recode_on_disk <- file.exists(as.character(registry$local_recode))
registry$boundary_on_disk <- file.exists(
  file.path(BOUNDARY_DIR, paste0(registry$SurveyId, ".rds"))
)

# Countries for which MAP publishes no national prevalence series are outside
# the malaria-endemic frame entirely, so their surveys are excluded by design.
pfpr_national <- read.csv(file.path(DATA_DIR, "pfpr_by_country_year.csv"),
                          stringsAsFactors = FALSE)
registry$has_map_series <- registry$iso3 %in% pfpr_national$iso3

# --- reconcile the regional keys ---------------------------------------------
# Only for surveys that are absent, and only where both inputs are readable:
# reading every Births Recode would be needlessly slow.
registry$n_boundary_units <- NA_integer_
registry$n_recode_units <- NA_integer_
registry$n_region_overlap <- NA_integer_

to_check <- which(!registry$in_aggregate & registry$recode_on_disk &
                    registry$boundary_on_disk)
for (i in to_check) {
  boundary <- tryCatch(
    readRDS(file.path(BOUNDARY_DIR, paste0(registry$SurveyId[i], ".rds"))),
    error = function(e) NULL
  )
  recode <- tryCatch(readRDS(as.character(registry$local_recode[i])),
                     error = function(e) NULL)
  if (is.null(boundary) || is.null(recode)) next
  if (!"DHSREGEN" %in% names(boundary) || !"v024" %in% names(recode)) next

  boundary_keys <- unique(rkey(boundary$DHSREGEN))
  recode_keys <- unique(rkey(as.character(haven::as_factor(recode$v024))))
  registry$n_boundary_units[i] <- length(boundary_keys)
  registry$n_recode_units[i] <- length(recode_keys)
  registry$n_region_overlap[i] <- length(intersect(boundary_keys, recode_keys))
}

# --- attribute a cause --------------------------------------------------------
classify <- function(row) {
  if (row$in_aggregate) return(NA_character_)
  if (!row$has_map_series) return("no_map_series")
  if (!row$recode_on_disk) return("recode_absent")
  if (!row$boundary_on_disk) return("boundary_absent")
  overlap <- row$n_region_overlap
  n_boundary <- row$n_boundary_units
  n_recode <- row$n_recode_units
  if (isTRUE(is.finite(overlap))) {
    if (isTRUE(n_boundary != n_recode)) return("boundary_granularity")
    if (isTRUE(overlap < n_recode)) return("region_key_mismatch")
  }
  if (!row$map_extraction_cached) return("no_map_extraction")
  "unattributed"
}
registry$omission_cause <- vapply(
  seq_len(nrow(registry)), function(i) classify(registry[i, ]), character(1)
)

# --- report -------------------------------------------------------------------
message("=== Survey coverage audit ===")
message(sprintf("Eligible Births Recodes %d-%d: %d surveys, %d countries",
                DHS_START_YEAR, DHS_END_YEAR, nrow(registry),
                length(unique(registry$iso3))))
message(sprintf("  in the legacy aggregate : %3d", sum(registry$in_aggregate)))
message(sprintf("  in the analysis dataset : %3d", sum(registry$in_analysis)))
message(sprintf("  in the main sample      : %3d", sum(registry$in_main_sample)))
message(sprintf("  MISSING                 : %3d (%.0f%%)",
                sum(!registry$in_aggregate),
                100 * mean(!registry$in_aggregate)))

omitted <- registry[!registry$in_aggregate, ]
message("\nOmitted surveys by attributed cause:")
print(table(omitted$omission_cause))

message("\nOmitted surveys:")
print(omitted[order(omitted$country, omitted$year),
              c("country", "year", "SurveyType", "svkey", "omission_cause",
                "n_boundary_units", "n_recode_units", "n_region_overlap")],
      row.names = FALSE)

recoverable <- omitted[omitted$has_map_series, ]
message(sprintf(
  paste("\n%d of the %d omitted surveys are in malaria-endemic countries with a",
        "MAP series,\nso they are recoverable rather than excluded by design.",
        "They span %d countries."),
  nrow(recoverable), nrow(omitted), length(unique(recoverable$iso3))
))

by_country <- data.frame(
  iso3 = registry$iso3, country = registry$country,
  eligible = 1L, analysed = as.integer(registry$in_aggregate)
)
by_country <- aggregate(cbind(eligible, analysed) ~ iso3 + country,
                        data = by_country, FUN = sum)
by_country$omitted <- by_country$eligible - by_country$analysed
by_country <- by_country[order(-by_country$omitted, by_country$country), ]
message("\nBy country (countries with any omission):")
print(by_country[by_country$omitted > 0, ], row.names = FALSE)

write.csv(registry, file.path(RESULTS_DIR, "survey_coverage_audit.csv"),
          row.names = FALSE)
write.csv(by_country,
          file.path(RESULTS_DIR, "survey_coverage_audit_by_country.csv"),
          row.names = FALSE)
message("\nWrote results/survey_coverage_audit.csv and ",
        "results/survey_coverage_audit_by_country.csv")

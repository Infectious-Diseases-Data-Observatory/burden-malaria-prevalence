# =============================================================================
# 02c_fetch_statcompiler_covariates.R
# Pull three household covariates for every survey in the registry from the DHS
# Program API (the StatCompiler indicators), at national level and by region:
#
#   WS_SRCE_H_IMP  households whose main drinking-water source is improved (%)
#   WS_TLET_H_IMP  households with an improved sanitation facility (%)
#   CN_NUTS_C_WH2  children under 5 wasted, weight-for-height below -2 SD (%)
#
# These replace the label-matching household aggregates the pipeline used to
# build from the births recodes (see results/dhs_rebuild/covariate_audit): the
# published indicators apply DHS's own category coding, so they are immune to
# the vocabulary, numeric-code and negation problems the audit found, and they
# exist for surveys whose recodes carry no anthropometry in the births file.
#
# Requires the rdhs configuration used by script 01 (the API needs a login for
# survey-level queries). Output is a long table, one row per survey x indicator
# x (national | region label):
#   data/derived_dhs/statcompiler_covariates.csv
#
# Run once; script 03 reads the cached file. Re-run to refresh after new surveys
# are added to the registry.
# =============================================================================

source("R_dhs/00_config.R")
required_packages("rdhs")

registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
surveys <- unique(registry$SurveyId[!is.na(registry$SurveyId) & nzchar(registry$SurveyId)])
message("Requesting ", length(STATCOMPILER_INDICATORS), " indicators for ", length(surveys), " surveys")

pull <- function(breakdown) {
  # the API accepts long survey lists, but keep requests in batches so one
  # failure does not lose everything
  batches <- split(surveys, ceiling(seq_along(surveys) / 40))
  do.call(rbind, lapply(batches, function(ids) {
    out <- tryCatch(
      rdhs::dhs_data(indicatorIds = unname(STATCOMPILER_INDICATORS), surveyIds = ids,
                     breakdown = breakdown),
      error = function(e) stop("DHS API request failed (", breakdown, "): ", conditionMessage(e)))
    if (is.null(out) || !nrow(out)) return(NULL)
    out$level <- breakdown
    as.data.frame(out, stringsAsFactors = FALSE)
  }))
}
national <- pull("national")
subnational <- pull("subnational")
raw <- rbind(national, subnational)
keep <- c("SurveyId", "IndicatorId", "level", "CharacteristicCategory", "CharacteristicLabel",
          "Value", "DenominatorWeighted", "DenominatorUnweighted")
keep <- intersect(keep, names(raw))
table <- raw[, keep]
table$variable <- names(STATCOMPILER_INDICATORS)[match(table$IndicatorId, STATCOMPILER_INDICATORS)]
table$svkey <- registry$svkey[match(table$SurveyId, registry$SurveyId)]
# subnational rows come at several nesting depths, marked by leading dots, and
# some labels carry qualifiers ("(>2010)", "2009"); keep the raw label and a
# cleaned one for matching
table$label_clean <- statcompiler_clean_label(table$CharacteristicLabel)
table <- table[order(table$svkey, table$variable, table$level, table$CharacteristicLabel), ]
write.csv(table, STATCOMPILER_CSV, row.names = FALSE)

coverage <- aggregate(SurveyId ~ variable + level, data = table,
                      FUN = function(x) length(unique(x)))
names(coverage)[3] <- "surveys"
print(coverage, row.names = FALSE)
message("Wrote ", STATCOMPILER_CSV, ": ", nrow(table), " rows")

suppressMessages(library(rdhs))
S <- commandArgs(trailingOnly = TRUE)[1]
setwd("/Users/jameswatson/Documents/Claude Projects/MIS:DHS malaria prevalence")
source("R_dhs/00_config.R")
reg <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
d <- read.csv(file.path(DERIVED_DIR, "dhs_analysis_dataset.csv"), stringsAsFactors = FALSE)
ind <- rdhs::dhs_indicators()
hits <- ind[grepl("improved", ind$Label, ignore.case = TRUE) & grepl("water|sanit|toilet", ind$Label, ignore.case = TRUE), c("IndicatorId", "Label", "Definition")]
print(hits[, 1:2], row.names = FALSE)
ids <- c("WS_SRCE_H_IMP", "WS_SRCE_P_IMP", "WS_TLET_H_IMP", "WS_TLET_P_IMP")
ids <- ids[ids %in% ind$IndicatorId]
surveys <- reg$SurveyId[reg$svkey %in% unique(d$svkey)]
nat <- rdhs::dhs_data(indicatorIds = ids, surveyIds = surveys, breakdown = "national")
sub <- tryCatch(rdhs::dhs_data(indicatorIds = ids, surveyIds = surveys, breakdown = "subnational"), error = function(e) { message("subnational failed: ", conditionMessage(e)); NULL })
write.csv(nat, file.path(S, "dhs_api_water_national.csv"), row.names = FALSE)
if (!is.null(sub)) write.csv(sub, file.path(S, "dhs_api_water_subnational.csv"), row.names = FALSE)
cat("national rows", nrow(nat), " subnational rows", if (is.null(sub)) NA else nrow(sub), "\n")
cat("surveys with a national household improved-water value:", length(unique(nat$SurveyId[nat$IndicatorId == "WS_SRCE_H_IMP"])), "of", length(surveys), "\n")

# Recorded child deaths by age at death across the processed Births Recodes.
# b5 is the survival status (labelled "yes"/"no" or 1/0); b7 is age at death in months.
reg <- read.csv("data/derived_dhs/survey_registry.csv", stringsAsFactors = FALSE)
man <- read.csv("data/derived_cbh/survey_manifest.csv", stringsAsFactors = FALSE)
reg <- reg[reg$svkey %in% man$survey[man$status == "built"], ]
dead_flag <- function(b5) {
  ch <- tolower(trimws(as.character(b5)))
  out <- ch %in% c("no", "0")
  out[is.na(ch)] <- FALSE
  out
}
tot <- data.frame(survey = reg$svkey, deaths = NA_integer_, under5 = NA_integer_,
                  over5 = NA_integer_, unknown = NA_integer_)
for (i in seq_len(nrow(reg))) {
  d <- readRDS(reg$local_recode[i])
  died <- dead_flag(d[["b5"]])
  b7 <- suppressWarnings(as.numeric(as.character(d[["b7"]])))
  tot$deaths[i]  <- sum(died)
  tot$under5[i]  <- sum(died & is.finite(b7) & b7 < 60)
  tot$over5[i]   <- sum(died & is.finite(b7) & b7 >= 60)
  tot$unknown[i] <- sum(died & !is.finite(b7))
  rm(d); gc(FALSE)
}
stopifnot(all(tot$deaths == tot$under5 + tot$over5 + tot$unknown))
write.csv(tot, "results/cbh/primary_map_regional17_gamma2_v3/study_flow/recorded_death_ages_by_survey.csv", row.names = FALSE)
cat(sprintf("surveys %d | deaths %d | under5 %d (%.1f%%) | 60+ months %d (%.1f%%) | unknown age %d\n",
    nrow(tot), sum(tot$deaths), sum(tot$under5), 100 * sum(tot$under5) / sum(tot$deaths),
    sum(tot$over5), 100 * sum(tot$over5) / sum(tot$deaths), sum(tot$unknown)))

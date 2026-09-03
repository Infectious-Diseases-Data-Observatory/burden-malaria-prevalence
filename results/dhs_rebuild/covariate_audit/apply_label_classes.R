# Apply the judged JMP classes to every survey's label table and compare the
# resulting household shares with the pipeline rule and with DHS StatCompiler.
setwd("/Users/jameswatson/Documents/Claude Projects/MIS:DHS malaria prevalence")
source("R_dhs/00_config.R")
S <- commandArgs(trailingOnly = TRUE)[1]
required_packages("jsonlite")
cls <- jsonlite::fromJSON(file.path(S, "label_classes.json"), simplifyVector = TRUE)
reg <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
nat <- read.csv(file.path(S, "dhs_api_water_national.csv"), stringsAsFactors = FALSE)
nat$svkey <- reg$svkey[match(nat$SurveyId, reg$SurveyId)]
water <- read.csv(file.path(S, "v113_labels_by_survey.csv"), stringsAsFactors = FALSE)
other <- read.csv(file.path(S, "v116_v119_labels_by_survey.csv"), stringsAsFactors = FALSE)
sanit <- other[other$variable == "v116", ]

summarise <- function(tab, classes, indicator, name, delivered_improved = TRUE) {
  tab$jmp_class <- classes$jmp_class[match(tab$label, classes$label)]
  tab$jmp_class[tab$numeric_code] <- "numeric_code"
  improved_set <- c("improved", if (delivered_improved) "improved_delivered_or_packaged")
  keep <- !tab$jmp_class %in% c("exclude_not_a_response", "numeric_code", NA)
  by <- split(tab[keep, ], tab$svkey[keep])
  out <- do.call(rbind, lapply(names(by), function(k) {
    z <- by[[k]]
    data.frame(svkey = k,
               rule = 100 * sum(z$w_hh[z$matched]) / sum(z$w_hh),
               jmp = 100 * sum(z$w_hh[z$jmp_class %in% improved_set]) / sum(z$w_hh),
               ambiguous = 100 * sum(z$w_hh[z$jmp_class == "ambiguous"]) / sum(z$w_hh),
               stringsAsFactors = FALSE)
  }))
  api <- nat[nat$IndicatorId == indicator, c("svkey", "Value")]
  out$dhs_api <- api$Value[match(out$svkey, api$svkey)]
  out$year <- reg$year[match(out$svkey, reg$svkey)]
  out$rule_minus_api <- out$rule - out$dhs_api
  out$jmp_minus_api <- out$jmp - out$dhs_api
  numeric <- unique(tab$svkey[tab$numeric_code & tab$w_hh > 0])
  numeric_share <- tapply(tab$w_hh * tab$numeric_code, tab$svkey, sum) / tapply(tab$w_hh, tab$svkey, sum)
  out$numeric_code_share <- 100 * as.numeric(numeric_share[out$svkey])
  cat(sprintf("\n%s (%d surveys): median |rule - DHS| %.1f, median |JMP classes - DHS| %.1f; surveys > 5 pts off: rule %d, JMP %d; > 10: rule %d, JMP %d\n",
              name, nrow(out), median(abs(out$rule_minus_api), na.rm = TRUE), median(abs(out$jmp_minus_api), na.rm = TRUE),
              sum(abs(out$rule_minus_api) > 5, na.rm = TRUE), sum(abs(out$jmp_minus_api) > 5, na.rm = TRUE),
              sum(abs(out$rule_minus_api) > 10, na.rm = TRUE), sum(abs(out$jmp_minus_api) > 10, na.rm = TRUE)))
  show <- out[abs(out$rule_minus_api) > 5 | abs(out$jmp_minus_api) > 5, ]
  show <- show[order(show$jmp_minus_api), c("svkey", "year", "rule", "jmp", "dhs_api", "numeric_code_share", "ambiguous")]
  print(show, row.names = FALSE, digits = 3)
  write.csv(out, file.path(S, paste0("label_classes_", name, "_by_survey.csv")), row.names = FALSE)
  # labels where the judged class disagrees with the rule, with their weight
  dis <- aggregate(w_hh ~ label + matched + jmp_class, data = tab[!tab$numeric_code, ], FUN = sum)
  dis$rule_improved <- dis$matched
  dis$jmp_improved <- dis$jmp_class %in% improved_set
  dis <- dis[dis$rule_improved != dis$jmp_improved & !dis$jmp_class %in% c("exclude_not_a_response"), ]
  dis$share <- 100 * dis$w_hh / sum(tab$w_hh)
  cat("  labels where the rule and the JMP class disagree (share of all households):\n")
  print(head(dis[order(-dis$w_hh), c("label", "rule_improved", "jmp_class", "share")], 30), row.names = FALSE, digits = 2)
  invisible(out)
}
w <- summarise(water, cls$water, "WS_SRCE_H_IMP", "water")
w2 <- summarise(water, cls$water, "WS_SRCE_H_IMP", "water_delivered_unimproved", delivered_improved = FALSE)
s <- summarise(sanit, cls$sanitation, "WS_TLET_H_IMP", "sanitation")
cat("\nverifier corrections:", nrow(cls$corrections$water), "water,", nrow(cls$corrections$sanitation), "sanitation\n")

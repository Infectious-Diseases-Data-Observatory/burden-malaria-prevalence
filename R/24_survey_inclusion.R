# =============================================================================
# 24_survey_inclusion.R — inclusion/exclusion accounting for the expanded panel.
# Enumerated universe (BR, SSA, 2000+, DHS/MIS/AIS)  ->  panel  ->  fit samples.
# Writes results/survey_inclusion.csv and prints a country x year summary.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(rdhs)}); options(rappdir_permission = TRUE)
CACHE <- path.expand("~/.rdhs_cache/datasets_reformatted")
BDIR  <- file.path(DATA, "dhs_boundaries")

## ---- 1. enumerate the candidate universe (same filter as R/18) --------------
ds <- dhs_datasets(fileFormat = "FL"); brs <- ds[ds$FileType == "Births Recode", ]
brs$iso3 <- countrycode::countrycode(brs$CountryName, "country.name", "iso3c", warn = FALSE)
brs$reg  <- countrycode::countrycode(brs$iso3, "iso3c", "region", warn = FALSE)
brs$yr   <- as.integer(brs$SurveyYear)
sv <- brs[!is.na(brs$reg) & brs$reg == "Sub-Saharan Africa" & brs$yr >= 2000 &
            brs$SurveyType %in% c("DHS","MIS","AIS"), ]
sv$svkey <- paste0(substr(sv$FileName, 1, 2), substr(sv$FileName, 5, 8))
sv$noext <- sub("\\..*$", "", sv$FileName)
sv$country <- countrycode::countrycode(sv$iso3, "iso3c", "country.name", warn = FALSE)
u <- sv[!duplicated(sv$svkey), c("iso3","country","yr","SurveyType","svkey","SurveyId","noext")]
names(u)[names(u) == "yr"] <- "year"

## ---- 2. what reached the panel and the fit samples --------------------------
d  <- read.csv(file.path(RESULTS, "component2_region_data_expanded.csv"), stringsAsFactors = FALSE)
nreg <- table(d$svkey)
u$in_panel   <- u$svkey %in% d$svkey
u$n_regions  <- as.integer(nreg[u$svkey]); u$n_regions[is.na(u$n_regions)] <- 0L

dh <- tryCatch(read.csv(file.path(RESULTS, "component2_region_data_expanded_health.csv"), stringsAsFactors = FALSE),
               error = function(e) NULL)
if (!is.null(dh)) {
  cc <- complete.cases(dh[, c("m1mo5y","pfpr10","dtp3","dtp3_reg","facility","educ_yrs","wealth_q",
                              "log_gdp","pct_urban","year_c")]) & is.finite(dh$exposure) & dh$exposure > 0 & dh$pfpr2_10 >= 1
  u$in_health_fit <- u$svkey %in% dh$svkey[cc]
} else u$in_health_fit <- NA

## ---- 3. best-guess exclusion reason for those not in the panel --------------
reason <- function(r) {
  if (r$in_panel) return("")
  if (!file.exists(file.path(CACHE, paste0(r$noext, ".rds")))) return("microdata not downloaded/no access")
  bf <- file.path(BDIR, paste0(r$SurveyId, ".rds"))
  if (!file.exists(bf)) return("no admin-1 boundary (DHS SDR)")
  "no usable mortality/prevalence overlap"
}
u$reason <- vapply(seq_len(nrow(u)), function(i) reason(u[i, ]), "")
u <- u[order(u$country, u$year), ]
write.csv(u[, c("iso3","country","year","SurveyType","svkey","in_panel","n_regions","in_health_fit","reason")],
          file.path(RESULTS, "survey_inclusion.csv"), row.names = FALSE)

## ---- 4. print summary -------------------------------------------------------
cat(sprintf("Enumerated universe (BR, SSA, 2000+, DHS/MIS/AIS): %d surveys, %d countries, %d-%d\n",
            nrow(u), length(unique(u$iso3)), min(u$year), max(u$year)))
cat(sprintf("  -> in expanded panel:            %d surveys, %d countries, %d region-years\n",
            sum(u$in_panel), length(unique(u$iso3[u$in_panel])), sum(u$n_regions)))
cat(sprintf("  -> in region-health fit sample:  %d surveys, %d countries\n",
            sum(u$in_health_fit, na.rm = TRUE), length(unique(u$iso3[u$in_health_fit]))))
cat(sprintf("  -> EXCLUDED from panel:          %d surveys\n\n", sum(!u$in_panel)))

cat("=== INCLUDED (expanded panel), by country: years [type] ===\n")
for (is in sort(unique(u$iso3[u$in_panel]))) {
  s <- u[u$iso3 == is & u$in_panel, ]
  cat(sprintf("  %-26s %s\n", s$country[1], paste(sprintf("%d[%s]", s$year, s$SurveyType), collapse = ", ")))
}
cat("\n=== EXCLUDED from panel: country year [type] — reason ===\n")
ex <- u[!u$in_panel, ]
for (i in seq_len(nrow(ex))) cat(sprintf("  %-26s %d [%s] — %s\n", ex$country[i], ex$year[i], ex$SurveyType[i], ex$reason[i]))

cat("\n=== also dropped from the region-HEALTH model (in panel, not in health fit) ===\n")
hd <- u[u$in_panel & !u$in_health_fit, ]
if (nrow(hd)) for (i in seq_len(nrow(hd))) cat(sprintf("  %-26s %d [%s]\n", hd$country[i], hd$year[i], hd$SurveyType[i])) else cat("  (none)\n")
cat("\nsaved: results/survey_inclusion.csv\n")

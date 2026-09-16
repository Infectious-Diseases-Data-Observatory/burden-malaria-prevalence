#!/usr/bin/env Rscript
# =============================================================================
# Sahel mortality by calendar month: build the aggregate cells.
#
# Question: within Sahelian survey regions, how does the all-cause death rate of
# children under two years vary with the month of the year?
#
# DHS birth histories record month and year of birth (B3, century-month code)
# and age at death with its unit (B6: days, months or years), not a date of
# death. The month of death is therefore derived as birth month plus age at
# death in completed months. That is exact to the month for deaths reported in
# days (assigned to the birth month when 0-30 days, and to floor(days / 30.44)
# months otherwise) and correct to within about half a month for deaths reported
# in completed months, because the day of birth within the month is unknown. A
# death reported in whole years cannot be placed within the year: children whose
# death is reported as "1 year" are excluded from the whole analysis (their
# exposure between 12 and 23 months is unknown too); deaths reported at two or
# more years fall outside the analysed ages, and those children contribute
# their full exposure at 0-23 months.
#
# Exposure. Every valid child contributes one person-month for each calendar
# month in the 60 months before the interview month during which the child was
# alive at the start of the month and aged 0-23 completed months. The month of
# death counts as exposure and carries the death, as in a piecewise-exponential
# model. Children whose birth month was imputed by DHS (flag B10) are dropped:
# DHS imputes the month at random within the admissible range, which would blur
# any seasonal pattern. Design weights (V005) are carried alongside unweighted
# counts.
#
# Geography and calendars follow the CBH builder exactly: survey rules for the
# region variable and any Ethiopian calendar offset, the deterministic region
# matcher with curated synonyms and overrides, and the same validity checks on
# dates, survival status, weights and death ages. Sahelian regions are those
# whose boundary centroid lies at or above 11 degrees north and west of the Horn
# of Africa (R_cbh/seasonality/settings.R).
#
# Outputs (results/cbh/seasonality_sahel_v1), all aggregate:
#   cells_survey_month_age.csv         survey x calendar month x age in months: deaths and person-months
#   cells_survey_region_month_band.csv survey x region x calendar month x age band
#   births_survey_month.csv            births in the exposure window by calendar month
#   sahel_regions.csv                  the regions used, with centroid and mean MAP prevalence
#   exclusions.csv                     children excluded, by survey and reason
#   build_provenance.csv, session_info.txt
# =============================================================================
source("R_cbh/load_pipeline.R")
source("R_cbh/seasonality/settings.R")

cfg <- cbh_config()
st <- cbh_seasonality_settings()
dir.create(st$out, recursive = TRUE, showWarnings = FALSE)
checked <- cbh_check_inputs(cfg)
registry <- checked$registry

## ---- Sahelian regions ------------------------------------------------------------------
centroids <- cbh_read_csv(st$centroids)
cbh_require(centroids, c("SurveyId", "regkey", "lat", "lon"), "Region centroids")
centroids$svkey <- registry$svkey[match(centroids$SurveyId, registry$SurveyId)]
regions <- checked$boundaries
regions$iso3 <- registry$iso3[match(regions$svkey, registry$svkey)]
regions$survey_year <- registry$year[match(regions$svkey, registry$svkey)]
j <- match(paste(regions$svkey, regions$regkey), paste(centroids$svkey, centroids$regkey))
regions$lat <- centroids$lat[j]; regions$lon <- centroids$lon[j]
if (anyNA(regions$lat)) {
  stop("Regions without a centroid: ", paste(head(paste(regions$svkey, regions$regkey)[is.na(regions$lat)]), collapse = ", "))
}
regions$sahel <- regions$lat >= st$latitude_min & regions$lon < st$longitude_max &
  !regions$iso3 %in% st$exclude_countries
sahel <- regions[regions$sahel, c("svkey", "iso3", "survey_year", "region", "regkey", "lat", "lon")]
surveys <- sort(unique(sahel$svkey))
message(nrow(sahel), " Sahelian regions in ", length(surveys), " surveys (centroid >= ",
        st$latitude_min, "N, west of ", st$longitude_max, "E, excluding ",
        paste(st$exclude_countries, collapse = "/"), ")")

# mean MAP PfPR2-10 over the exposure years, for description only
map <- checked$inputs$map
if (all(c("svkey", "regkey", "year", "pfpr2_10") %in% names(map))) {
  m <- map[map$svkey %in% surveys, ]
  m <- m[m$year >= registry$year[match(m$svkey, registry$svkey)] - 5 &
           m$year <= registry$year[match(m$svkey, registry$svkey)], ]
  mean_pfpr <- tapply(m$pfpr2_10, paste(m$svkey, m$regkey), mean, na.rm = TRUE)
  sahel$mean_pfpr_pct <- as.numeric(mean_pfpr[paste(sahel$svkey, sahel$regkey)])
} else {
  sahel$mean_pfpr_pct <- NA_real_
}

## ---- helpers -------------------------------------------------------------------------------
# B10 tells whether the month of birth was reported or imputed. Labels beginning
# "month" (month and year complete; month, year and day; month and age; month
# only) mean the month was reported; numeric codes 1, 2 and 7 are the same cases.
birth_month_reported <- function(b10) {
  lab <- tolower(trimws(cbh_label(b10)))
  num <- cbh_num(b10)
  has_text <- !is.na(lab) & grepl("[a-z]", lab)
  out <- rep(NA, length(lab))
  out[has_text] <- grepl("^month", lab[has_text])
  out[!has_text & !is.na(num)] <- num[!has_text & !is.na(num)] %in% c(1, 2, 7)
  out
}

build_survey <- function(k) {
  survey <- registry[registry$svkey == k, , drop = FALSE]
  rule <- checked$rules[checked$rules$svkey == k, , drop = FALSE]
  boundary <- checked$boundaries[checked$boundaries$svkey == k, , drop = FALSE]
  br <- cbh_read_recode(survey$local_recode, c(rule$region_var, rule$strata_var))
  cbh_require(br, c("v005", "v008", "b3", "b5", "b7"), "Births Recode")
  geo <- cbh_geography(br, survey, rule, boundary, checked$overrides, registry)
  n <- nrow(br)
  offset <- as.integer(rule$cmc_offset_months)
  born <- cbh_column(br, "b3", TRUE) + offset
  interview <- cbh_column(br, "v008", TRUE) + offset
  alive <- cbh_decode_category(br$b5, c("no", "yes"), c(0, 1))
  weight <- cbh_column(br, "v005", TRUE) / 1e6
  b7 <- cbh_column(br, "b7", TRUE)
  b6 <- if ("b6" %in% names(br)) br$b6 else cbh_column(br, "b6")
  age <- cbh_death_age(b6, b7)
  died <- !is.na(alive) & alive == 0

  # completed months of age at death, where the month can be placed
  death_age <- rep(NA_real_, n)
  days <- died & age$unit %in% "days"
  death_age[days] <- ifelse(age$reported_days[days] <= 30, 0, floor(age$reported_days[days] / 30.4375))
  in_months <- died & age$unit %in% c("months", "imputed_months")
  death_age[in_months] <- age$lo[in_months]
  in_years <- died & age$unit %in% "years"
  # two or more years: the child survived the analysed ages; "1 year" cannot be placed
  death_age[in_years & age$lo >= 24] <- age$lo[in_years & age$lo >= 24]
  unresolved <- in_years & age$lo < 24

  invalid <- rep("", n)
  mark <- function(bad, reason) invalid[which(bad & invalid == "")] <<- reason
  mark(!is.finite(born) | !is.finite(interview) | born != floor(born) |
         interview != floor(interview) | born <= 0 | interview < born, "invalid_dates")
  mark(is.na(alive), "unknown_survival_status")
  mark(!is.finite(weight) | weight <= 0, "invalid_weight")
  mark(died & !is.finite(age$lo), "unknown_age_at_death")
  mark(died & is.finite(death_age) & death_age > interview - born, "death_after_interview")
  mark(!is.na(alive) & alive == 1 & is.finite(b7) & b7 >= 0 & b7 < 600, "alive_with_death_age")
  mark(unresolved, "death_month_unresolved_year_unit")
  reported <- birth_month_reported(cbh_column(br, "b10"))
  mark(!is.na(reported) & !reported, "birth_month_imputed")
  if ("caseid" %in% names(br)) {
    mother_key <- trimws(cbh_column(br, "caseid"))
  } else {
    mother_key <- cbh_key(br, c("v001", "v002", "v003"))
  }
  index_name <- if ("bidx" %in% names(br)) "bidx" else "bord"
  ck <- paste(mother_key, cbh_column(br, index_name, TRUE), sep = "|")
  mark(duplicated(ck) | duplicated(ck, fromLast = TRUE), "duplicate_child_identifier")
  mark(is.na(geo$regkey), "region_unmatched")
  in_sahel <- !is.na(geo$regkey) & paste(k, geo$regkey) %in% paste(sahel$svkey, sahel$regkey)
  mark(invalid == "" & !in_sahel, "outside_sahel")

  valid <- invalid == ""
  exclusions <- as.data.frame(table(ifelse(valid, "used", invalid)), stringsAsFactors = FALSE)
  names(exclusions) <- c("status", "children")
  exclusions$survey <- k
  exclusions$birth_month_flag_unavailable <- sum(is.na(reported))

  # exposure months: alive at the start of the month, aged 0-23 completed months,
  # within the 60 months before the interview month
  death_cmc <- ifelse(died & is.finite(death_age), born + death_age, NA_real_)
  start <- pmax(born, interview - st$lookback_months)
  stop <- pmin(born + st$max_age_months, interview - 1)
  stop <- ifelse(is.finite(death_cmc), pmin(stop, death_cmc), stop)
  months <- stop - start + 1
  use <- which(valid & months > 0)
  if (!length(use)) return(list(cells = NULL, births = NULL, exclusions = exclusions))
  idx <- rep(use, months[use])
  cmc <- start[idx] + sequence(months[use]) - 1
  age_months <- as.integer(cmc - born[idx])
  month <- as.integer((cmc - 1) %% 12 + 1)
  death_here <- as.numeric(is.finite(death_cmc[idx]) & death_cmc[idx] == cmc)
  w <- weight[idx]
  key <- paste(geo$regkey[idx], month, age_months, sep = "|")
  sums <- cbind(person_months_w = w, person_months_n = 1, deaths_w = w * death_here, deaths_n = death_here)
  agg <- rowsum(sums, key)
  parts <- do.call(rbind, strsplit(rownames(agg), "|", fixed = TRUE))
  cells <- data.frame(svkey = k, iso3 = survey$iso3, survey_year = survey$year,
                      regkey = parts[, 1], month = as.integer(parts[, 2]),
                      age_months = as.integer(parts[, 3]), agg, stringsAsFactors = FALSE)
  rownames(cells) <- NULL
  # births in the exposure window, by calendar month of birth
  born_in_window <- valid & born >= interview - st$lookback_months & born <= interview - 1
  bm <- as.integer((born[born_in_window] - 1) %% 12 + 1)
  births <- data.frame(svkey = k, iso3 = survey$iso3, month = 1:12,
                       births_w = as.numeric(tapply(weight[born_in_window], factor(bm, levels = 1:12), sum, default = 0)),
                       births_n = as.numeric(table(factor(bm, levels = 1:12))), stringsAsFactors = FALSE)
  births$births_w[is.na(births$births_w)] <- 0
  list(cells = cells, births = births, exclusions = exclusions)
}

## ---- run -------------------------------------------------------------------------------------
cells <- births <- exclusions <- list()
for (i in seq_along(surveys)) {
  k <- surveys[i]
  result <- build_survey(k)
  cells[[k]] <- result$cells; births[[k]] <- result$births; exclusions[[k]] <- result$exclusions
  used <- result$exclusions$children[result$exclusions$status == "used"]
  message(sprintf("[%d/%d] %s: %s children used, %.0f deaths under 24 months", i, length(surveys), k,
                  if (length(used)) used else 0L, if (is.null(result$cells)) 0 else sum(result$cells$deaths_n)))
}
cells <- cbh_bind(cells); births <- cbh_bind(births); exclusions <- cbh_bind(exclusions)
cells$band <- st$bands$band[findInterval(cells$age_months, st$bands$age_lo)]

by_survey_age <- aggregate(cbind(person_months_w, person_months_n, deaths_w, deaths_n) ~
                             svkey + iso3 + survey_year + band + month + age_months, data = cells, FUN = sum)
by_region_band <- aggregate(cbind(person_months_w, person_months_n, deaths_w, deaths_n) ~
                              svkey + iso3 + survey_year + regkey + band + month, data = cells, FUN = sum)
region_totals <- aggregate(cbind(person_months_w, deaths_w) ~ svkey + regkey, data = cells, FUN = sum)
sahel <- merge(sahel, region_totals, by = c("svkey", "regkey"), all.x = TRUE)
sahel <- sahel[order(sahel$iso3, sahel$svkey, sahel$regkey), ]

cbh_atomic_csv(by_survey_age, file.path(st$out, "cells_survey_month_age.csv"))
cbh_atomic_csv(by_region_band, file.path(st$out, "cells_survey_region_month_band.csv"))
cbh_atomic_csv(births, file.path(st$out, "births_survey_month.csv"))
cbh_atomic_csv(sahel, file.path(st$out, "sahel_regions.csv"))
cbh_atomic_csv(exclusions, file.path(st$out, "exclusions.csv"))
provenance <- data.frame(input = c(cfg$registry, cfg$survey_rules, cfg$region_overrides, cfg$boundary_regions,
                                   st$centroids, "R_cbh/seasonality/01_build_cells.R", "R_cbh/seasonality/settings.R"))
provenance$md5 <- vapply(provenance$input, cbh_file_hash, character(1))
provenance$built <- format(Sys.time(), tz = "UTC", usetz = TRUE)
cbh_atomic_csv(provenance, file.path(st$out, "build_provenance.csv"))
writeLines(capture.output(utils::sessionInfo()), file.path(st$out, "session_info.txt"))

summary_tab <- aggregate(cbind(person_months_w, deaths_w, deaths_n) ~ band, data = cells, FUN = sum)
summary_tab$rate_per_1000_child_years <- 12000 * summary_tab$deaths_w / summary_tab$person_months_w
message("\nSahel, ", length(surveys), " surveys, ", nrow(sahel), " regions:")
print(transform(summary_tab, person_years_w = round(person_months_w / 12), deaths_w = round(deaths_w),
                rate_per_1000_child_years = round(rate_per_1000_child_years, 1))[
                  , c("band", "person_years_w", "deaths_w", "deaths_n", "rate_per_1000_child_years")], row.names = FALSE)
excl <- aggregate(children ~ status, data = exclusions, FUN = sum)
message("\nChildren by status:"); print(excl[order(-excl$children), ], row.names = FALSE)
message("\nWrote ", st$out)

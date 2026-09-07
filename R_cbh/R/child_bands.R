cbh_decode_category <- function(x, labels, codes) {
  text <- tolower(trimws(cbh_label(x)))
  out <- unname(codes[match(text, labels)])
  numeric <- cbh_num(x)
  use <- is.na(out) & numeric %in% codes
  out[use] <- numeric[use]
  out
}

# B6 retains the reported units. B7 is the DHS-imputed age in months.
# A years-reported death is assigned to its whole-year band, not to a
# fictitious exact death month. B7 provides a documented fallback.
cbh_death_age <- function(b6, b7) {
  text <- tolower(trimws(cbh_label(b6)))
  raw <- cbh_num(b6)
  unit <- floor(raw / 100)
  value <- raw %% 100
  unit[!unit %in% 1:3 | value >= 97] <- NA_real_
  for (u in 1:3) {
    name <- c("days", "months", "years")[u]
    hit <- !is.na(text) & grepl(paste0("^", name, "?:? *[0-9]+$"), text)
    unit[hit] <- u
    value[hit] <- cbh_num(sub("^[^0-9]*", "", text[hit]))
  }
  zero <- !is.na(text) & text == "died on day of birth"
  unit[zero] <- 1; value[zero] <- 0
  unit[value >= 97] <- NA_real_
  lo <- hi <- rep(NA_real_, length(raw))
  days <- unit == 1 & !is.na(unit)
  # DHS month-based mortality convention: reported 0–30 days are neonatal.
  lo[days] <- ifelse(value[days] <= 30, 0, 1)
  hi[days] <- ifelse(value[days] <= 30, 1, 6)
  months <- unit == 2 & !is.na(unit)
  lo[months] <- value[months]; hi[months] <- value[months] + 1
  years <- unit == 3 & !is.na(unit)
  lo[years] <- 12 * value[years]; hi[years] <- 12 * (value[years] + 1)
  fallback <- is.na(lo) & is.finite(b7) & b7 >= 0 & b7 < 600 & b7 == floor(b7)
  lo[fallback] <- b7[fallback]; hi[fallback] <- lo[fallback] + 1
  list(lo = lo, hi = hi,
       source = ifelse(fallback, "b7_imputed_months", ifelse(is.na(lo), "unknown", "b6_reported")),
       unit = ifelse(fallback, "imputed_months", c("days", "months", "years")[unit]),
       reported_days = ifelse(days, value, NA_real_))
}

cbh_make_bands <- function(br, survey, rules, geography, cfg) {
  cbh_require(br, c("v005", "v008", "b3", "b5", "b7"), "Births Recode")
  n <- nrow(br)
  if (!n) stop("The recode is empty.")
  offset <- as.integer(rules$cmc_offset_months)
  if (length(offset) != 1L || is.na(offset)) stop("An explicit CMC offset is required.")
  born <- cbh_column(br, "b3", TRUE) + offset
  interview <- cbh_column(br, "v008", TRUE) + offset
  observed_year <- stats::median(cbh_cmc_year(interview), na.rm = TRUE)
  if (!is.finite(observed_year) || abs(observed_year - survey$year) > 2) {
    stop("Interview dates disagree with the survey year; review the explicit calendar rule.")
  }
  alive <- cbh_decode_category(br$b5, c("no", "yes"), c(0, 1))
  weight <- cbh_column(br, "v005", TRUE) / 1e6
  b7 <- cbh_column(br, "b7", TRUE)
  b6 <- cbh_column(br, "b6")
  if ("b6" %in% names(br)) b6 <- br$b6
  age <- cbh_death_age(b6, b7)
  died <- !is.na(alive) & alive == 0
  invalid <- rep("", n)
  mark <- function(bad, reason) invalid[which(bad & invalid == "")] <<- reason
  mark(!is.finite(born) | !is.finite(interview) | born != floor(born) |
         interview != floor(interview) | born <= 0 | interview < born, "invalid_dates")
  mark(is.na(alive), "unknown_survival_status")
  mark(!is.finite(weight) | weight <= 0, "invalid_weight")
  mark(died & !is.finite(age$lo), "unknown_age_at_death")
  mark(died & is.finite(age$lo) & age$lo > interview - born, "death_after_interview")
  mark(!is.na(alive) & alive == 1 & is.finite(b7) & b7 >= 0 & b7 < 600,
       "alive_with_death_age")
  death_band <- findInterval(age$lo, c(cfg$age_bands$age_lo, 60))
  # Reported intervals must not cross a chosen band boundary.
  upper_band <- findInterval(age$hi - 1e-7, c(cfg$age_bands$age_lo, 60))
  mark(died & death_band != upper_band & age$lo < 60, "ambiguous_death_band")
  b7_band <- findInterval(b7, c(cfg$age_bands$age_lo, 60))
  death_band_disagreement <- died & age$source == "b6_reported" &
    is.finite(b7) & b7 >= 0 & b7 < 600 & b7 == floor(b7) & b7_band != death_band
  death_band[!died] <- 8L
  # Stable within an input snapshot, pseudonymous, no original case IDs exported.
  if ("caseid" %in% names(br)) {
    mother_key <- trimws(cbh_column(br, "caseid"))
  } else {
    cbh_require(br, c("v001", "v002", "v003"), "Mother identifiers")
    mother_key <- cbh_key(br, c("v001", "v002", "v003"))
    mother_key[!complete.cases(br[c("v001", "v002", "v003")])] <- NA_character_
  }
  index_name <- if ("bidx" %in% names(br)) "bidx" else "bord"
  cbh_require(br, index_name, "Child identifiers")
  child_index <- cbh_column(br, index_name, TRUE)
  mark(is.na(mother_key) | !nzchar(mother_key) | !is.finite(child_index) |
         child_index <= 0 | child_index != floor(child_index),
       "missing_child_identifier")
  ck <- paste(mother_key, child_index, sep = "|")
  duplicate <- duplicated(ck) | duplicated(ck, fromLast = TRUE)
  mark(duplicate, "duplicate_child_identifier")

  sex <- cbh_decode_category(br$b4 %cbh_or% rep(NA_character_, n), c("male", "female"), c(1, 2))
  urban <- cbh_decode_category(br$v025 %cbh_or% rep(NA_character_, n), c("urban", "rural"), c(1, 2))
  wealth <- cbh_decode_category(br$v190 %cbh_or% rep(NA_character_, n),
                               c("poorest", "poorer", "middle", "richer", "richest"), 1:5)
  education <- cbh_column(br, "v133", TRUE)
  education[!is.finite(education) | education < 0 | education > 30] <- NA_real_
  maternal_age <- (cbh_column(br, "b3", TRUE) - cbh_column(br, "v011", TRUE)) / 12
  maternal_age[!is.finite(maternal_age) | maternal_age < 10 | maternal_age > 55] <- NA_real_
  order <- cbh_column(br, "bord", TRUE)
  order[!is.finite(order) | order < 1 | order > 40] <- NA_real_
  preceding <- cbh_column(br, "b11", TRUE)
  preceding[!is.finite(preceding) | preceding < 0 | preceding >= 900] <- NA_real_
  preceding[!is.na(order) & order == 1] <- NA_real_
  b0 <- tolower(cbh_column(br, "b0"))
  multiple <- ifelse(grepl("multiple|twin|triplet", b0), 1,
                     ifelse(b0 == "single birth", 0, NA_real_))
  b0n <- cbh_column(br, "b0", TRUE)
  numeric_b0 <- is.finite(b0n) & b0n >= 0 & b0n <= 8
  multiple[numeric_b0] <- as.numeric(b0n[numeric_b0] > 0)
  stratum_var <- rules$strata_var
  if (is.na(stratum_var) || !nzchar(stratum_var)) {
    usable <- c("v022", "v023")
    usable <- usable[vapply(usable, function(v) v %in% names(br) && any(!is.na(br[[v]])), logical(1))]
    stratum_var <- if (length(usable)) usable[1] else NA_character_
  }
  stratum <- if (is.na(stratum_var)) rep(NA_character_, n) else cbh_column(br, stratum_var)
  cluster <- cbh_column(br, "v021")
  encode_id <- function(x, prefix) ifelse(is.na(x) | !nzchar(trimws(x)), NA_character_,
                                         paste(survey$svkey, prefix, match(x, unique(x)), sep = ":"))
  children <- data.frame(child_id = paste(survey$svkey, "c", seq_len(n), sep = ":"),
    mother_id = encode_id(mother_key, "m"), psu = encode_id(cluster, "p"),
    stratum = encode_id(stratum, "s"), stratum_variable = stratum_var,
    survey = survey$svkey, country = survey$iso3, survey_year = survey$year,
    regkey = geography$regkey, region = geography$region,
    survey_weight = weight, birth_cmc = born, interview_cmc = interview,
    cmc_offset_months = offset, calendar_conversion = rules$calendar,
    sex = ifelse(sex == 1, "male", ifelse(sex == 2, "female", NA_character_)),
    urban = ifelse(urban == 1, 1, ifelse(urban == 2, 0, NA_real_)),
    maternal_age_birth = maternal_age, maternal_education_years = education,
    wealth_quintile = wealth, birth_order = order, multiple_birth = multiple,
    first_birth = as.numeric(order == 1), preceding_birth_interval_months = preceding,
    birth_date_flag = cbh_column(br, "b10"), death_age_flag = cbh_column(br, "b13"),
    death_age_source = ifelse(died, age$source, NA_character_),
    death_age_unit = ifelse(died, age$unit, NA_character_),
    death_band_b6_b7_disagree = death_band_disagreement,
    neonatal_death_days = ifelse(died, age$reported_days, NA_real_))
  rows <- list(); flow <- list()
  for (a in seq_len(nrow(cfg$age_bands))) {
    band <- cfg$age_bands[a, ]
    entry <- born + band$age_lo
    end <- born + band$age_hi
    reached <- invalid == "" & death_band >= a
    window <- reached & entry >= interview - cfg$entry_lookback_months & entry < interview
    complete <- window & end <= interview
    year <- cbh_cmc_year(entry)
    supported <- complete & year >= cfg$first_entry_year & year <= cfg$last_entry_year
    use <- which(supported %in% TRUE)
    d <- children[use, , drop = FALSE]
    d$age_band <- rep(band$age_band, nrow(d))
    d$age_band_index <- rep(a, nrow(d))
    d$age_lo <- rep(band$age_lo, nrow(d)); d$age_hi <- rep(band$age_hi, nrow(d))
    d$band_years <- rep((band$age_hi - band$age_lo) / 12, nrow(d))
    d$log_band_years <- log(d$band_years)
    d$band_entry_cmc <- entry[use]; d$band_end_cmc <- end[use]
    d$entry_year <- year[use]; d$calendar_year <- cbh_cmc_time(entry[use])
    d$months_before_interview_at_entry <- interview[use] - entry[use]
    d$death <- as.integer(died[use] & death_band[use] == a)
    rows[[a]] <- d
    flow[[a]] <- data.frame(survey = survey$svkey, age_band = band$age_band,
      valid_children_reaching_band = sum(reached, na.rm = TRUE),
      entered_in_lookback = sum(window, na.rm = TRUE),
      incomplete_potential_band = sum(window & !complete, na.rm = TRUE),
      full_band_before_interview = sum(complete, na.rm = TRUE),
      entry_year_outside_range = sum(complete & !supported, na.rm = TRUE),
      eligible_band_rows = length(use), deaths = sum(d$death))
  }
  reasons <- as.data.frame(table(ifelse(invalid == "", "valid", invalid)), stringsAsFactors = FALSE)
  names(reasons) <- c("status", "children")
  reasons$survey <- survey$svkey
  list(data = cbh_bind(rows), flow = cbh_bind(flow), child_checks = reasons)
}

`%cbh_or%` <- function(x, y) if (is.null(x)) y else x

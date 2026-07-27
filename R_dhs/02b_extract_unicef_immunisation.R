# =============================================================================
# 02b_extract_unicef_immunisation.R
# Extract a compact country-year WUENIC panel from UNICEF's global dataflow CSV.
#
# Required local input:
#   data/fusion_GLOBAL_DATAFLOW_UNICEF_1.0_all.csv
#
# Output:
#   data/derived_dhs/unicef_immunisation_country_year.csv
#
# The source is a public aggregate dataset. Only the three required vaccine
# series and their country-year provenance flags are retained.
# =============================================================================

source("R_dhs/00_config.R")

if (!file.exists(UNICEF_GLOBAL_CSV)) {
  stop("UNICEF global dataflow file not found: ", UNICEF_GLOBAL_CSV)
}

indicator_map <- c(
  IM_HIB3 = "hib3_wuenic",
  IM_PCVC = "pcv3_wuenic",
  IM_ROTAC = "rotac_wuenic"
)

extract_indicator_lines <- function(path, indicators, chunk_size = 100000L) {
  connection <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  header <- readLines(connection, n = 1L, warn = FALSE)
  if (length(header) != 1L) stop("UNICEF dataflow file has no header.")

  patterns <- paste0(",", indicators, ",")
  selected <- list()
  chunk_index <- 0L
  repeat {
    lines <- readLines(connection, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    keep <- Reduce(
      `|`,
      lapply(patterns, function(pattern) grepl(pattern, lines, fixed = TRUE))
    )
    if (any(keep)) {
      chunk_index <- chunk_index + 1L
      selected[[chunk_index]] <- lines[keep]
    }
  }
  c(header, unlist(selected, use.names = FALSE))
}

message("Extracting WUENIC vaccine rows from: ", UNICEF_GLOBAL_CSV)
selected_lines <- extract_indicator_lines(
  UNICEF_GLOBAL_CSV,
  names(indicator_map)
)
raw <- read.csv(
  text = paste(selected_lines, collapse = "\n"),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

required <- c(
  "REF_AREA", "INDICATOR", "SEX", "TIME_PERIOD", "OBS_VALUE",
  "UNIT_MULTIPLIER", "UNIT_MEASURE", "OBS_STATUS", "DATA_SOURCE", "AGE"
)
missing_columns <- setdiff(required, names(raw))
if (length(missing_columns)) {
  stop(
    "UNICEF file is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

raw <- raw[
  raw$INDICATOR %in% names(indicator_map) &
    raw$SEX == "_T" &
    raw$UNIT_MEASURE == "PCNT" &
    raw$AGE == "M12T23" &
    raw$OBS_STATUS == "E",
  required,
  drop = FALSE
]
raw$year <- suppressWarnings(as.integer(raw$TIME_PERIOD))
raw$value <- suppressWarnings(as.numeric(raw$OBS_VALUE))
raw$iso3 <- toupper(raw$REF_AREA)
raw <- raw[
  nchar(raw$iso3) == 3L &
    is.finite(raw$year) &
    is.finite(raw$value) &
    raw$value >= 0 & raw$value <= 100,
  ,
  drop = FALSE
]

key <- paste(raw$iso3, raw$INDICATOR, raw$year, sep = "|")
duplicate_keys <- unique(key[duplicated(key) | duplicated(key, fromLast = TRUE)])
if (length(duplicate_keys)) {
  conflicting <- vapply(duplicate_keys, function(duplicate_key) {
    length(unique(raw$value[key == duplicate_key])) > 1L
  }, logical(1))
  if (any(conflicting)) {
    stop(
      "Conflicting duplicate UNICEF estimates for: ",
      paste(duplicate_keys[conflicting], collapse = ", ")
    )
  }
  raw <- raw[!duplicated(key), , drop = FALSE]
}

countries <- sort(unique(raw$iso3))
panel <- expand.grid(
  iso3 = countries,
  year = DHS_START_YEAR:DHS_END_YEAR,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

for (indicator in names(indicator_map)) {
  variable <- unname(indicator_map[[indicator]])
  status_variable <- paste0(variable, "_status")
  panel[[variable]] <- NA_real_
  panel[[status_variable]] <- "missing_after_series_start"
  indicator_rows <- raw[raw$INDICATOR == indicator, , drop = FALSE]

  for (iso3 in countries) {
    panel_rows <- panel$iso3 == iso3
    country_rows <- indicator_rows[
      indicator_rows$iso3 == iso3,
      ,
      drop = FALSE
    ]
    if (!nrow(country_rows)) {
      panel[[variable]][panel_rows] <- 0
      panel[[status_variable]][panel_rows] <-
        "no_series_assumed_not_introduced"
      next
    }

    before_series <- panel_rows & panel$year < min(country_rows$year)
    panel[[variable]][before_series] <- 0
    panel[[status_variable]][before_series] <-
      "pre_series_assumed_not_introduced"

    match_index <- match(
      paste(panel$iso3[panel_rows], panel$year[panel_rows], sep = "|"),
      paste(country_rows$iso3, country_rows$year, sep = "|")
    )
    observed <- is.finite(match_index)
    target_rows <- which(panel_rows)[observed]
    panel[[variable]][target_rows] <- country_rows$value[match_index[observed]]
    panel[[status_variable]][target_rows] <- "wuenic_estimate"
  }
}

source_labels <- unique(raw$DATA_SOURCE[nzchar(raw$DATA_SOURCE)])
panel$unicef_immunisation_source <- if (length(source_labels)) {
  paste(source_labels, collapse = "; ")
} else {
  "UNICEF global dataflow"
}

summary_rows <- lapply(unname(indicator_map), function(variable) {
  status <- panel[[paste0(variable, "_status")]]
  data.frame(
    variable = variable,
    source_indicator = names(indicator_map)[indicator_map == variable],
    countries = length(unique(panel$iso3)),
    country_years = nrow(panel),
    wuenic_estimates = sum(status == "wuenic_estimate"),
    assumed_pre_introduction_zeros = sum(grepl("assumed", status)),
    missing_after_series_start = sum(!is.finite(panel[[variable]])),
    min_coverage = min(panel[[variable]], na.rm = TRUE),
    max_coverage = max(panel[[variable]], na.rm = TRUE)
  )
})
summary <- do.call(rbind, summary_rows)

write.csv(panel, UNICEF_IMMUNISATION_CSV, row.names = FALSE)
write.csv(summary, UNICEF_IMMUNISATION_SUMMARY_CSV, row.names = FALSE)

message(
  "Saved ", nrow(panel), " country-years for ", length(countries),
  " countries: ", UNICEF_IMMUNISATION_CSV
)
print(summary, row.names = FALSE)

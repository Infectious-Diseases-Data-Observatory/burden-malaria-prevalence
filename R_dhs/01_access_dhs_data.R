# =============================================================================
# 01_access_dhs_data.R
# Enumerate eligible DHS/MIS Births Recodes and, only when explicitly asked,
# download missing recodes and survey-region boundary files.
#
# Default (inventory only):
#   Rscript R_dhs/01_access_dhs_data.R
#
# Explicit access/download (requires an authorised rdhs configuration):
#   Rscript R_dhs/01_access_dhs_data.R --download
#
# This script does not print or export individual-level records.
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("rdhs", "countrycode", "sf"))
options(rappdir_permission = TRUE)

args <- commandArgs(trailingOnly = TRUE)
download_requested <- "--download" %in% args

enumerate_surveys <- function() {
  datasets <- rdhs::dhs_datasets(fileFormat = "FL")
  births <- datasets[datasets$FileType == "Births Recode", , drop = FALSE]
  births$iso3 <- countrycode::countrycode(
    births$CountryName, "country.name", "iso3c", warn = FALSE
  )
  births$world_bank_region <- countrycode::countrycode(
    births$iso3, "iso3c", "region", warn = FALSE
  )
  births$year <- as.integer(births$SurveyYear)
  keep <- !is.na(births$world_bank_region) &
    births$world_bank_region == "Sub-Saharan Africa" &
    births$year >= DHS_START_YEAR &
    births$year <= DHS_END_YEAR &
    births$SurveyType %in% c("DHS", "MIS")
  out <- births[keep, , drop = FALSE]
  out$no_extension <- toupper(sub("\\..*$", "", out$FileName))
  out$svkey <- survey_key(out$FileName)
  out <- out[!duplicated(out$svkey), , drop = FALSE]
  out$boundary_file <- file.path(BOUNDARY_DIR, paste0(out$SurveyId, ".rds"))
  out$local_recode <- vapply(
    out$no_extension, local_recode_path, character(1)
  )
  out$recode_available <- !is.na(out$local_recode)
  out$boundary_available <- file.exists(out$boundary_file)
  out[order(out$iso3, out$year), , drop = FALSE]
}

download_missing_recodes <- function(registry) {
  missing <- registry[!registry$recode_available, , drop = FALSE]
  if (!nrow(missing)) {
    message("All eligible Births Recodes are already available locally.")
    return(invisible(registry))
  }
  message(
    "Requesting ", nrow(missing),
    " authorised DHS Births Recodes through rdhs."
  )
  paths <- rdhs::get_datasets(
    dataset_filenames = missing$FileName,
    download_option = "rds",
    reformat = TRUE,
    clear_cache = FALSE
  )
  for (i in seq_len(nrow(missing))) {
    key <- missing$no_extension[i]
    source_path <- paths[[key]]
    if (is.character(source_path) && length(source_path) == 1 &&
        file.exists(source_path)) {
      destination <- file.path(DHS_LOCAL_DIR, paste0(key, ".rds"))
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      if (!file.exists(destination)) file.copy(source_path, destination)
    }
  }
  invisible(registry)
}

download_missing_boundaries <- function(registry) {
  missing <- registry[!registry$boundary_available, , drop = FALSE]
  if (!nrow(missing)) {
    message("All available survey boundaries are already cached.")
    return(invisible(registry))
  }
  message("Requesting up to ", nrow(missing), " DHS survey boundary files.")
  for (i in seq_len(nrow(missing))) {
    boundary <- tryCatch(
      rdhs::download_boundaries(
        surveyId = missing$SurveyId[i],
        method = "sf"
      )[[1]],
      error = function(e) NULL
    )
    if (!is.null(boundary)) saveRDS(boundary, missing$boundary_file[i])
  }
  invisible(registry)
}

registry <- enumerate_surveys()
message(
  "Eligible survey universe: ", nrow(registry), " surveys; ",
  length(unique(registry$iso3)), " countries; ",
  min(registry$year), "-", max(registry$year), "."
)
message(
  "Available before optional download: ",
  sum(registry$recode_available), " recodes; ",
  sum(registry$boundary_available), " boundaries."
)

if (download_requested) {
  download_missing_recodes(registry)
  download_missing_boundaries(registry)
  registry <- enumerate_surveys()
}

write.csv(
  registry[, c(
    "iso3", "CountryName", "year", "SurveyType", "SurveyId", "FileName",
    "no_extension", "svkey", "local_recode", "boundary_file",
    "recode_available", "boundary_available"
  )],
  SURVEY_REGISTRY_CSV,
  row.names = FALSE
)

message(
  "Saved survey registry: ", SURVEY_REGISTRY_CSV,
  if (download_requested) " (post-download inventory)." else " (inventory only)."
)

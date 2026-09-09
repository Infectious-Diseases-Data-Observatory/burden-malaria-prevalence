# Shared CLI validation for annual burden and plotting stages.
cbh_burden_year <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args)) return(2024L)
  if (length(args) != 1L || !grepl("^--year=[0-9]{4}$", args))
    stop("Use --year=YYYY (2000-2024), or no argument for 2024.")
  year <- as.integer(sub("^--year=", "", args))
  if (is.na(year) || year < 2000L || year > 2024L) stop("Year must be between 2000 and 2024.")
  year
}

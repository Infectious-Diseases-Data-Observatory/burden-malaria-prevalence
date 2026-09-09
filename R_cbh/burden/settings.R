# Separate-age primary by default; retain explicit access to the historical joint model.
cbh_burden_options <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (any(!grepl("^--year=[0-9]{4}$|^--model=(separate|joint)$", args)))
    stop("Use --year=YYYY (2000-2024) and/or --model=separate|joint.")
  ya <- args[grepl("^--year=", args)]; ma <- args[grepl("^--model=", args)]
  if (length(ya) > 1L || length(ma) > 1L) stop("Specify each option at most once.")
  year <- if (length(ya)) as.integer(sub("^--year=", "", ya)) else 2024L
  if (is.na(year) || year < 2000L || year > 2024L) stop("Year must be between 2000 and 2024.")
  model <- if (length(ma)) sub("^--model=", "", ma) else "separate"
  source_id <- "age_band_hiv_incidence_shared_time_v3"
  list(year = year, model = model, source_id = source_id,
    result_id = if (model == "separate") "age_band_separate_v1" else source_id,
    label = if (model == "separate") "Separate age-band models; one median HIV imputation" else
      "Joint shared-calendar-year model; ten HIV imputations",
    interval_label = if (model == "separate") "conditional model uncertainty only; HIV imputation held fixed" else
      "model and HIV-imputation uncertainty only")
}
cbh_burden_year <- function(args = commandArgs(trailingOnly = TRUE)) cbh_burden_options(args)$year

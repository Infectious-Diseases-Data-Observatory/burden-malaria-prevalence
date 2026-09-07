# Read only shards named in the completed manifest; never glob stale outputs.
cbh_load_analysis <- function(output_dir, surveys = NULL, required_covariates = character(),
                              model_ready_only = TRUE, allow_failed_surveys = FALSE) {
  meta <- readRDS(file.path(output_dir, "manifest.rds"))
  if (!isTRUE(meta$complete)) stop("The dataset build did not complete.")
  m <- meta$manifest
  if (!is.null(surveys)) {
    if (length(setdiff(surveys, m$survey))) stop("Requested survey is absent from this build manifest.")
    m <- m[m$survey %in% surveys, , drop = FALSE]
  }
  if (any(m$status == "failed") && !allow_failed_surveys) stop("The build has failed surveys; review survey_manifest.csv before explicitly allowing partial output.")
  skipped <- m[!m$status %in% c("built", "cached"), c("survey", "status"), drop = FALSE]
  m <- m[m$status %in% c("built", "cached"), , drop = FALSE]
  if (!nrow(m)) stop("No built surveys are available.")
  rows <- report <- vector("list", nrow(m))
  for (i in seq_len(nrow(m))) {
    obj <- readRDS(file.path(output_dir, m$file[i]))
    if (!identical(obj$signature, m$signature[i])) stop("Shard signature differs from manifest; rebuild before reading.")
    d <- obj$data
    cbh_require(d, required_covariates, "Requested adjustment set")
    core <- if (model_ready_only) d$model_ready else rep(TRUE, nrow(d))
    complete <- rep(TRUE, nrow(d))
    for (v in required_covariates) {
      complete <- complete & !is.na(d[[v]])
      if (is.numeric(d[[v]])) complete <- complete & is.finite(d[[v]])
    }
    keep <- core & complete
    report[[i]] <- data.frame(survey = m$survey[i], eligible_rows = nrow(d),
      missing_core_exposure_or_region = sum(!d$model_ready),
      missing_requested_covariates = sum(core & !complete), selected_rows = sum(keep),
      selected_deaths = sum(d$death[keep]))
    d <- d[keep, , drop = FALSE]
    # Normalize after the final row selection. This is not a survey variance estimator.
    d$analysis_weight <- if (nrow(d)) d$survey_weight / mean(d$survey_weight) else numeric()
    rows[[i]] <- d
  }
  d <- cbh_bind(rows)
  d$age_band <- factor(d$age_band, levels = meta$config$age_bands$age_band, ordered = FALSE)
  for (v in c("survey", "country", "region", "sex", "mother_id", "psu", "stratum")) d[[v]] <- factor(d[[v]])
  d$country_age <- interaction(d$country, d$age_band, drop = TRUE)
  attr(d, "selection_report") <- cbh_bind(report)
  attr(d, "skipped_surveys") <- skipped
  attr(d, "schema_version") <- meta$schema_version
  if (nrow(skipped)) warning(nrow(skipped), " surveys are unavailable; see attr(d, 'skipped_surveys').", call. = FALSE)
  d
}

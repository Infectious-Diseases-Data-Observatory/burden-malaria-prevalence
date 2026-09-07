cbh_check_inputs <- function(cfg, surveys = NULL) {
  registry <- cbh_read_csv(cfg$registry)
  cbh_require(registry, c("svkey", "iso3", "year", "local_recode"), "Survey registry")
  cbh_unique(registry, "svkey", "Survey registry")
  rules <- cbh_read_csv(cfg$survey_rules)
  cbh_require(rules, c("svkey", "region_var", "cmc_offset_months", "calendar",
                       "group_donor", "group_fine", "group_coarse", "strata_var"), "Survey rules")
  cbh_unique(rules, "svkey", "Survey rules")
  boundaries <- unique(cbh_read_csv(cfg$boundary_regions)[c("svkey", "region", "regkey")])
  cbh_unique(boundaries, c("svkey", "regkey"), "Boundary region definitions")
  overrides <- cbh_read_csv(cfg$region_overrides)
  cbh_require(overrides, c("svkey", "region_var", "source_label", "regkey", "region_id"), "Region overrides")
  cbh_unique(overrides, c("svkey", "region_var", "source_label"), "Region overrides")
  inputs <- cbh_external_inputs(cfg)
  selected <- registry
  if (!is.null(surveys)) {
    if (length(setdiff(surveys, registry$svkey))) stop("Unknown survey requested.")
    selected <- registry[registry$svkey %in% surveys, , drop = FALSE]
  }
  if (!nrow(selected)) stop("No surveys selected.")
  if (any(!selected$svkey %in% rules$svkey)) stop("Selected surveys need explicit survey rules.")
  if (any(!grepl("^[A-Za-z0-9_-]+$", selected$svkey))) stop("Unsafe survey key.")
  checks <- data.frame(survey = selected$svkey, country = selected$iso3, survey_year = selected$year,
    recode_available = file.exists(selected$local_recode),
    boundary_available = selected$svkey %in% boundaries$svkey,
    annual_map_available = selected$svkey %in% inputs$map$svkey)
  list(registry = registry, selected = selected, rules = rules, boundaries = boundaries,
       overrides = overrides, inputs = inputs, checks = checks)
}

cbh_validate_bands <- function(d, cfg) {
  cbh_unique(d, c("child_id", "age_band"), "Child-band output")
  if (!nrow(d)) return(invisible(TRUE))
  stopifnot(all(d$death %in% 0:1), all(is.finite(d$band_years) & d$band_years > 0),
    all(d$band_end_cmc <= d$interview_cmc),
    all(d$band_entry_cmc >= d$interview_cmc - cfg$entry_lookback_months),
    all(d$band_entry_cmc < d$interview_cmc),
    all(d$entry_year == cbh_cmc_year(d$band_entry_cmc)),
    all(d$exposure_year == d$entry_year),
    all(d$band_years == (d$age_hi - d$age_lo) / 12),
    all(d$entry_year >= cfg$first_entry_year & d$entry_year <= cfg$last_entry_year),
    all(is.finite(d$survey_weight) & d$survey_weight > 0))
  dead <- d[d$death == 1L, c("child_id", "age_band_index"), drop = FALSE]
  if (anyDuplicated(dead$child_id)) stop("A child has multiple death outcomes.")
  j <- match(d$child_id, dead$child_id)
  if (any(d$age_band_index > dead$age_band_index[j], na.rm = TRUE)) stop("A band follows a child's death.")
  invisible(TRUE)
}

cbh_private_output <- function(cfg) {
  # Resolve existing ancestors, including symlinks, before creating any directory.
  parent <- cfg$output_dir; tail <- character()
  while (!dir.exists(parent)) {
    if (basename(parent) %in% c(".", "..", "")) stop("Use a simple output directory under project data/.")
    tail <- c(basename(parent), tail)
    next_parent <- dirname(parent)
    if (identical(parent, next_parent)) stop("Cannot resolve output directory.")
    parent <- next_parent
  }
  resolved <- do.call(file.path, as.list(c(normalizePath(parent), tail)))
  private_root <- normalizePath(file.path(cfg$root, "data"), mustWork = TRUE)
  if (!startsWith(paste0(resolved, "/"), paste0(private_root, "/"))) {
    stop("Generated child data must remain under the project's ignored data/ directory.")
  }
  resolved
}

cbh_build <- function(cfg, surveys = NULL, force = FALSE) {
  checked <- cbh_check_inputs(cfg, surveys)
  out <- cbh_private_output(cfg)
  dir.create(file.path(out, "child_bands"), recursive = TRUE, showWarnings = FALSE)
  code <- sort(c(list.files(file.path(cfg$root, "R_cbh"), "[.]R$", full.names = TRUE),
                 list.files(file.path(cfg$root, "R_cbh", "R"), "[.]R$", full.names = TRUE)))
  sources <- unique(c(cfg$registry, cfg$survey_rules, cfg$region_overrides, cfg$boundary_regions,
                      cfg$annual_map, cfg$statcompiler,
                      vapply(cfg$annual_panels, `[[`, character(1), "path"), code))
  source_hashes <- data.frame(path = sources, md5 = vapply(sources, cbh_file_hash, character(1)))
  settings <- cfg; settings$output_dir <- NULL
  base_signature <- cbh_hash(list(settings, source_hashes))
  started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  # A reader refuses an interrupted build rather than silently loading an old run.
  cbh_atomic_rds(list(complete = FALSE, started = started), file.path(out, "manifest.rds"))
  manifests <- flows <- checks <- crosswalks <- missingness <- list()
  for (i in seq_len(nrow(checked$selected))) {
    survey <- checked$selected[i, , drop = FALSE]
    rule <- checked$rules[checked$rules$svkey == survey$svkey, , drop = FALSE]
    boundary <- checked$boundaries[checked$boundaries$svkey == survey$svkey, , drop = FALSE]
    path <- file.path(out, "child_bands", paste0(survey$svkey, ".rds"))
    donor_path <- checked$registry$local_recode[checked$registry$svkey %in% rule$group_donor]
    raw_hash <- cbh_file_hash(survey$local_recode)
    signature <- cbh_hash(list(base_signature, survey, rule, raw_hash,
                              vapply(donor_path, cbh_file_hash, character(1))))
    m <- data.frame(survey = survey$svkey, country = survey$iso3, survey_year = survey$year,
      status = "pending", rows = 0L, model_ready_rows = 0L, deaths = 0L,
      model_ready_deaths = 0L, signature = signature, raw_md5 = raw_hash,
      file = file.path("child_bands", paste0(survey$svkey, ".rds")), error = "")
    available <- checked$checks[i, ]
    if (!available$recode_available) {
      m$status <- "missing_recode"
    } else if (!available$boundary_available || !available$annual_map_available) {
      m$status <- "missing_map_geography"
    } else {
      result <- tryCatch({
        object <- if (!force && file.exists(path)) tryCatch(readRDS(path), error = function(e) NULL) else NULL
        cached <- !is.null(object) && identical(object$signature, signature)
        if (!cached) {
          br <- cbh_read_recode(survey$local_recode, c(rule$region_var, rule$strata_var))
          geo <- cbh_geography(br, survey, rule, boundary, checked$overrides, checked$registry)
          bands <- cbh_make_bands(br, survey, rule, geo, cfg)
          bands$data <- cbh_attach_external(bands$data, survey, boundary, checked$inputs, cfg)
          cbh_validate_bands(bands$data, cfg)
          object <- c(bands, list(crosswalk = geo$crosswalk, signature = signature,
                                  schema_version = cfg$schema_version))
          cbh_atomic_rds(object, path)
          rm(br, geo, bands)
        }
        cbh_validate_bands(object$data, cfg)
        list(object = object, cached = cached)
      }, error = function(e) list(error = conditionMessage(e)))
      if (!is.null(result$error)) {
        m$status <- "failed"; m$error <- result$error
      } else {
        d <- result$object$data
        m$status <- if (result$cached) "cached" else "built"
        m$rows <- nrow(d); m$model_ready_rows <- sum(d$model_ready)
        m$deaths <- sum(d$death); m$model_ready_deaths <- sum(d$death[d$model_ready])
        flows[[i]] <- result$object$flow; checks[[i]] <- result$object$child_checks
        crosswalks[[i]] <- result$object$crosswalk
        missingness[[i]] <- data.frame(survey = survey$svkey, variable = names(d), rows = nrow(d),
          missing = vapply(d, function(x) sum(is.na(x)), integer(1)))
        rm(d)
      }
      rm(result)
    }
    manifests[[i]] <- m
    message(sprintf("[%d/%d] %s: %s; %s eligible rows, %s with PfPR", i, nrow(checked$selected),
                    survey$svkey, m$status, m$rows, m$model_ready_rows))
    invisible(gc(FALSE))
  }
  manifest <- cbh_bind(manifests)
  reports <- list(survey_manifest = manifest, eligibility_flow = cbh_bind(flows),
                  child_checks = cbh_bind(checks), region_crosswalk = cbh_bind(crosswalks),
                  variable_missingness = cbh_bind(missingness), source_files = source_hashes,
                  input_checks = checked$checks)
  for (name in names(reports)) cbh_atomic_csv(reports[[name]], file.path(out, paste0(name, ".csv")))
  cbh_atomic_rds(list(complete = TRUE, started = started,
    finished = format(Sys.time(), tz = "UTC", usetz = TRUE),
    schema_version = cfg$schema_version, config = cfg, manifest = manifest,
    source_files = source_hashes, session_info = utils::sessionInfo()), file.path(out, "manifest.rds"))
  invisible(manifest)
}

# All joins are exact in calendar year. Values remain in their source units.
cbh_external_inputs <- function(cfg) {
  map <- cbh_read_csv(cfg$annual_map)
  cbh_require(map, c("svkey", "regkey", "year", "pfpr2_10"), "Annual MAP panel")
  cbh_unique(map, c("svkey", "regkey", "year"), "Annual MAP panel")
  panels <- lapply(names(cfg$annual_panels), function(name) {
    spec <- cfg$annual_panels[[name]]
    d <- cbh_read_csv(spec$path, required = FALSE)
    if (!is.null(d)) {
      cbh_require(d, c("iso3", "year", unname(spec$columns), unname(spec$source_status)), name)
      cbh_unique(d, c("iso3", "year"), name)
    }
    d
  })
  names(panels) <- names(cfg$annual_panels)
  stat <- cbh_read_csv(cfg$statcompiler, required = FALSE)
  if (!is.null(stat)) cbh_require(stat, c("svkey", "variable", "level", "CharacteristicLabel", "Value"), "STATcompiler")
  list(map = map, panels = panels, statcompiler = stat)
}

cbh_attach_external <- function(d, survey, boundaries, inputs, cfg) {
  n <- nrow(d)
  map <- inputs$map
  j <- match(cbh_key(d, c("survey", "regkey", "entry_year")),
             cbh_key(map, c("svkey", "regkey", "year")))
  p <- cbh_num(map$pfpr2_10[j])
  valid <- is.finite(p) & p >= 0 & p <= 100
  d$pfpr_pct <- ifelse(valid, p, NA_real_)
  d$pfpr_status <- ifelse(is.na(d$regkey), "unmatched_region",
                          ifelse(is.na(j), "missing_region_year",
                                 ifelse(valid, "observed", "missing_or_invalid_value")))
  d$exposure_year <- d$entry_year
  for (name in names(cfg$annual_panels)) {
    spec <- cfg$annual_panels[[name]]
    panel <- inputs$panels[[name]]
    j <- if (is.null(panel)) rep(NA_integer_, n) else
      match(cbh_key(d, c("country", "entry_year")), cbh_key(panel, c("iso3", "year")))
    for (out in names(spec$columns)) {
      v <- if (is.null(panel)) rep(NA_real_, n) else cbh_num(panel[[spec$columns[[out]]]][j])
      valid <- is.finite(v)
      if (out %in% spec$percent) valid <- valid & v >= 0 & v <= 100
      if (out %in% c("gdp_pc", "health_expenditure_pc")) valid <- valid & v > 0
      status <- ifelse(is.na(j), if (is.null(panel)) "missing_panel" else "missing_country_year",
                       ifelse(valid, "observed", "missing_or_invalid_value"))
      if (out %in% names(spec$source_status) && !is.null(panel)) {
        source_status <- as.character(panel[[spec$source_status[[out]]]][j])
        allowed <- !is.na(source_status) & source_status %in% spec$allowed_status
        rejected <- !is.na(j) & !allowed
        status[rejected] <- paste0("not_observed:", ifelse(is.na(source_status[rejected]), "unknown", source_status[rejected]))
        valid <- valid & allowed
      }
      d[[out]] <- ifelse(valid, v, NA_real_)
      d[[paste0(out, "_status")]] <- status
    }
  }
  # No arbitrary positive constant is added to a zero prevalence.
  for (pair in list(c("hiv_prev_pct", "log_hiv_prev"), c("gdp_pc", "log_gdp_pc"),
                   c("health_expenditure_pc", "log_health_expenditure_pc"))) {
    v <- d[[pair[1]]]
    d[[pair[2]]] <- log(ifelse(is.finite(v) & v > 0, v, NA_real_))
  }
  # These are measured at survey, not at band entry. They are candidates only.
  # Never silently replace missing regional measurements by national values.
  stat <- inputs$statcompiler
  for (variable in c("imp_water", "imp_sanit", "wasting")) {
    sub <- if (is.null(stat)) NULL else stat[stat$svkey == survey$svkey & stat$variable == variable, , drop = FALSE]
    for (level in c("subnational", "national")) {
      out <- paste0(variable, if (level == "subnational") "_survey_region_pct" else "_survey_national_pct")
      v <- rep(NA_real_, n); status <- rep("missing", n)
      if (!is.null(sub) && nrow(sub)) {
        s <- sub[sub$level == level & !is.na(sub$level), , drop = FALSE]
        if (level == "subnational" && nrow(s)) {
          crosswalk <- cbh_match_regions(s$CharacteristicLabel, boundaries)
          s$key <- crosswalk$regkey[match(s$CharacteristicLabel, crosswalk$source_label)]
          ambiguous <- duplicated(s$key) | duplicated(s$key, fromLast = TRUE)
          j <- match(d$regkey, s$key)
          usable <- !is.na(j) & !is.na(d$regkey) & !ambiguous[j]
          v[usable] <- cbh_num(s$Value[j[usable]])
          status[!is.na(j) & ambiguous[j]] <- "ambiguous_source"
        } else if (level == "national" && nrow(s) == 1L) {
          v <- rep(cbh_num(s$Value), n)
        } else if (level == "national" && nrow(s) > 1L) status[] <- "ambiguous_source"
        valid <- is.finite(v) & v >= 0 & v <= 100
        status[valid] <- "observed_at_survey"
        v[!valid] <- NA_real_
      }
      d[[out]] <- v
      d[[paste0(out, "_status")]] <- status
    }
  }
  d$survey_covariate_year <- rep(survey$year, n)
  d$model_ready <- !is.na(d$region) & is.finite(d$pfpr_pct)
  stopifnot(nrow(d) == n)
  d
}

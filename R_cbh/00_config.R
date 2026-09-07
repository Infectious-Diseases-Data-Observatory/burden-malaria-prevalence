# Configuration only. No reads, downloads, installs or fits when sourced.
cbh_config <- function(root = getwd()) {
  root <- normalizePath(root, mustWork = TRUE)
  list(
    root = root,
    output_dir = file.path(root, "data", "derived_cbh"),
    registry = file.path(root, "data/derived_dhs/survey_registry.csv"),
    boundary_regions = file.path(root, "data/derived_dhs/map_pfpr_by_survey_region.csv"),
    annual_map = file.path(root, "data/derived_dhs/map_pfpr_window_years.csv"),
    survey_rules = file.path(root, "R_cbh/config/survey_rules.csv"),
    region_overrides = file.path(root, "R_cbh/config/region_overrides.csv"),
    entry_lookback_months = 60L,
    first_entry_year = 2000L,
    last_entry_year = 2024L,
    exposure_timing = "band_entry_calendar_year",
    age_bands = data.frame(
      age_band = c("<1", "1-5", "6-11", "12-23", "24-35", "36-47", "48-59"),
      age_lo = c(0L, 1L, 6L, 12L, 24L, 36L, 48L),
      age_hi = c(1L, 6L, 12L, 24L, 36L, 48L, 60L)),
    # Candidate X variables are retained, not selected or imputed here.
    # An unavailable optional panel produces NA and a source-status flag.
    annual_panels = list(
      hiv = list(path = file.path(root, "data/derived_dhs/hiv_prevalence_country_year.csv"),
                 columns = c(hiv_prev_pct = "hiv_prev"), percent = "hiv_prev_pct"),
      immunisation = list(path = file.path(root, "data/derived_dhs/unicef_immunisation_country_year.csv"),
                          columns = c(hib3_pct = "hib3_wuenic", pcv3_pct = "pcv3_wuenic",
                                      rotavirus_pct = "rotac_wuenic"),
                          percent = c("hib3_pct", "pcv3_pct", "rotavirus_pct"),
                          source_status = c(hib3_pct = "hib3_wuenic_status",
                                            pcv3_pct = "pcv3_wuenic_status",
                                            rotavirus_pct = "rotac_wuenic_status"),
                          allowed_status = "wuenic_estimate"),
      gdp = list(path = file.path(root, "data/wb_gdp_pc.csv"),
                 columns = c(gdp_pc = "gdp_pc")),
      health_spending = list(path = file.path(root, "data/wb_hexp_pc.csv"),
                             columns = c(health_expenditure_pc = "hexp_pc")),
      governance = list(path = file.path(root, "data/wb_polstab.csv"),
                        columns = c(political_stability = "polstab"))
    ),
    statcompiler = file.path(root, "data/derived_dhs/statcompiler_covariates.csv"),
    # Output is local, access-controlled microdata. Do not put it in results/.
    schema_version = "cbh-band-entry-v1"
  )
}

# Default revised primary data policy; old no-fallback audit is reproducible.
cbh_covariate_settings <- function(complete_case_only=FALSE,expanded=FALSE) {
  expanded <- expanded || complete_case_only
  id <- if(complete_case_only) "regional_adjustment_v1" else if(expanded) "regional_adjustment_unicef_v2" else "regional_adjustment_reduced_v3"
  base <- "data/derived_cbh/regional_adjustment"
  list(id=id,expanded=expanded,regional_mean=!expanded,unicef_fallback=!complete_case_only,base=base,
    private=if(complete_case_only)base else file.path(base,if(expanded)"unicef_v2" else "reduced_v3"),
    out=file.path("results/cbh",id),panel=file.path(base,"unicef/country_year_estimates.csv"))
}

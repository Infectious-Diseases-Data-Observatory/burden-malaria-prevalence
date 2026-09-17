# Default revised primary data policy; old no-fallback audit is reproducible.
cbh_covariate_settings <- function(complete_case_only=FALSE) {
  id <- if(complete_case_only) "regional_adjustment_v1" else "regional_adjustment_unicef_v2"
  base <- "data/derived_cbh/regional_adjustment"
  list(id=id,unicef_fallback=!complete_case_only,base=base,
    private=if(complete_case_only)base else file.path(base,"unicef_v2"),
    out=file.path("results/cbh",id),panel=file.path(base,"unicef/country_year_estimates.csv"))
}

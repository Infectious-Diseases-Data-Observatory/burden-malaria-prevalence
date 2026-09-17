# Explicit national vaccine substitution, with provenance. Source after load_pipeline.R.
cbh_unicef_indicators <- function() c(dtp3_pct="IM_DTP3",measles_pct="IM_MCV1",
  hib3_pct="IM_HIB3",pcv3_pct="IM_PCVC",rotavirus_pct="IM_ROTAC")

cbh_unicef_fill <- function(d,panel,variables=names(cbh_unicef_indicators())) {
  cbh_require(panel,c("country","year","variable","value","source"),"UNICEF vaccine panel")
  cbh_unique(panel,c("country","year","variable"),"UNICEF vaccine panel")
  for(v in variables) {
    year_variable <- if(v %in% c("dtp3_pct","measles_pct")) "survey_year" else "entry_year"
    cbh_require(d,c("country",year_variable,v),"Vaccine substitution input")
    p <- panel[panel$variable==v,,drop=FALSE]
    j <- match(paste(d$country,d[[year_variable]]),paste(p$country,p$year))
    value <- p$value[j]
    available <- is.finite(value) & value>=0 & value<=100
    missing <- !is.finite(d[[v]])
    fill <- missing & available
    original <- d[[v]]
    d[[v]][fill] <- value[fill]
    d[[paste0(v,"_before_imputation")]] <- original
    d[[paste0(v,"_imputed")]] <- fill
    d[[paste0(v,"_source_year")]] <- d[[year_variable]]
    d[[paste0(v,"_imputation_source")]] <- ifelse(fill,p$source[j],
      ifelse(missing,"unresolved_missing",if(v %in% c("dtp3_pct","measles_pct"))
        "retained_regional_estimate" else "retained_national_estimate"))
    stopifnot(identical(d[[v]][!missing],original[!missing]),
      all(!fill | is.finite(d[[v]])))
  }
  d
}

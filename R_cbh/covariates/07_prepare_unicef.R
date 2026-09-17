#!/usr/bin/env Rscript
# Extract actual estimates from the existing public UNICEF snapshot. No network,
# no fabricated zeros, no nearest-year filling, no modification of legacy panel.
source("R_cbh/load_pipeline.R")
source("R_cbh/covariates/unicef.R")
library(data.table)
path <- "data/fusion_GLOBAL_DATAFLOW_UNICEF_1.0_all.csv"
out <- "data/derived_cbh/regional_adjustment/unicef"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
extract <- function(path) {
  con <- file(path,"r",encoding="UTF-8");on.exit(close(con))
  header <- readLines(con,n=1,warn=FALSE)
  found <- list(); pattern <- paste0(",(",paste(cbh_unicef_indicators(),collapse="|"),"),")
  repeat {
    lines <- readLines(con,n=100000L,warn=FALSE)
    if(!length(lines)) break
    keep <- grepl(pattern,lines)
    if(any(keep)) found[[length(found)+1L]] <- lines[keep]
  }
  if(!length(found))stop("No selected UNICEF vaccine records found")
  read.csv(text=paste(c(header,unlist(found)),collapse="\n"),stringsAsFactors=FALSE,
    na.strings=c("","NA"),check.names=FALSE)
}
message("Reading five vaccine series from local UNICEF snapshot")
raw <- extract(path)
cbh_require(raw,c("REF_AREA","INDICATOR","SEX","TIME_PERIOD","OBS_VALUE","UNIT_MULTIPLIER",
  "UNIT_MEASURE","OBS_STATUS","DATA_SOURCE","AGE"),"UNICEF source")
registry <- cbh_read_csv(cbh_config()$registry)
keep <- with(raw,REF_AREA %in% registry$iso3 & SEX=="_T" & UNIT_MEASURE=="PCNT" &
  AGE=="M12T23" & OBS_STATUS=="E" & grepl("WHO/UNICEF estimates of national immunization coverage",DATA_SOURCE,fixed=TRUE))
raw <- raw[which(keep),,drop=FALSE]
multiplier <- suppressWarnings(as.numeric(raw$UNIT_MULTIPLIER))
stopifnot(all(!is.na(multiplier) & multiplier==0))
p <- data.frame(country=raw$REF_AREA,year=suppressWarnings(as.integer(raw$TIME_PERIOD)),
  variable=names(cbh_unicef_indicators())[match(raw$INDICATOR,cbh_unicef_indicators())],
  value=suppressWarnings(as.numeric(raw$OBS_VALUE)),source="UNICEF_WUENIC_country_year",
  source_indicator=raw$INDICATOR,source_release=raw$DATA_SOURCE)
stopifnot(!anyNA(p),all(is.finite(p$value) & p$value>=0 & p$value<=100))
p <- unique(p)
cbh_unique(p,c("country","year","variable"),"Selected WUENIC estimates")
stopifnot(all(names(cbh_unicef_indicators()) %in% p$variable))
p <- p[order(p$country,p$variable,p$year),]
cbh_atomic_csv(p,file.path(out,"country_year_estimates.csv"))
cbh_atomic_csv(raw,file.path(out,"selected_source_rows.csv"))
cbh_atomic_csv(data.frame(file=c(path,"R_cbh/covariates/unicef.R","R_cbh/covariates/07_prepare_unicef.R"),
  md5=vapply(c(path,"R_cbh/covariates/unicef.R","R_cbh/covariates/07_prepare_unicef.R"),cbh_file_hash,"")),
  file.path(out,"input_manifest.csv"))
print(as.data.table(p)[,.(countries=uniqueN(country),country_years=.N,first_year=min(year),last_year=max(year)),by=variable])

#!/usr/bin/env Rscript
# Explicit public-data refresh; never reads or uploads respondent data.
# Run from project root. The audit itself uses this local snapshot offline.
library(jsonlite)
source("R_cbh/load_pipeline.R")
out <- "data/derived_cbh/regional_adjustment"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
ids <- c("CH_VACC_C_DP3", "CH_VACC_C_MSL", "RH_DELP_C_DHF",
  "CN_IYCB_C_EXB", "CN_BFSS_C_EBF", "HC_ELEC_H_ELC",
  "WS_SRCE_H_IMP", "WS_TLET_H_IMP", "FE_BINT_C_I07", "FE_BINT_C_I18")
registry <- cbh_read_csv(cbh_config()$registry)
surveys <- unique(na.omit(registry$SurveyId))
provenance <- list()
fetch <- function(endpoint, query, name) {
  rows <- list(); page <- 1L
  repeat {
    url <- paste0("https://api.dhsprogram.com/rest/dhs/", endpoint,
      "?f=json&perpage=1000&page=", page, "&", query)
    path <- file.path(out,paste0(name,"_page",page,".json"))
    if(!file.exists(path) || "--refresh" %in% commandArgs(TRUE))
      download.file(url,path,quiet=TRUE,method="libcurl")
    x <- jsonlite::fromJSON(path)
    stopifnot(is.data.frame(x$Data) || length(x$Data)==0,
      x$RecordsReturned==NROW(x$Data))
    if(NROW(x$Data)) rows[[page]] <- x$Data
    provenance[[length(provenance)+1L]] <<- data.frame(url=url,file=path,
      md5=cbh_file_hash(path),retrieved_utc=format(file.info(path)$mtime,tz="UTC",usetz=TRUE))
    if(page >= x$TotalPages) break
    page <- page+1L
  }
  if(!length(rows)) return(data.frame())
  z <- as.data.frame(data.table::rbindlist(rows,fill=TRUE))
  stopifnot(nrow(z)==x$RecordCount)
  z
}
metadata <- fetch("indicators",paste0("indicatorIds=",paste(ids,collapse=",")),"indicators")
cbh_atomic_csv(metadata,file.path(out,"indicator_metadata.csv"))
all <- list()
for(level in c("subnational","national")) {
  batches <- split(surveys,ceiling(seq_along(surveys)/20))
  for(i in seq_along(batches)) {
    message(level, " batch ",i,"/",length(batches))
    x <- fetch("data",paste0("indicatorIds=",paste(ids,collapse=","),
      "&surveyIds=",paste(batches[[i]],collapse=","),"&breakdown=",level),paste0(level,"_",i))
    if(nrow(x)) { x$level <- level; all[[length(all)+1L]] <- x }
  }
}
d <- as.data.frame(data.table::rbindlist(all,fill=TRUE))
# API DataId is not globally unique: retain it, but validate semantic keys.
cbh_unique(d,c("SurveyId","IndicatorId","CharacteristicId","ByVariableId","DataId"),"Published source records")
stopifnot(all(d$SurveyId %in% surveys),all(d$IndicatorId %in% ids))
cbh_atomic_csv(d,file.path(out,"published_indicators.csv"))
cbh_atomic_csv(do.call(rbind,provenance),file.path(out,"source_manifest.csv"))
message("Saved ",nrow(d)," published aggregate records.")

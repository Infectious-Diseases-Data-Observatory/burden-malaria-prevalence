#!/usr/bin/env Rscript
# Rebuild only the regional covariate overlay; leave base child-bands/fits intact.
source("R_cbh/load_pipeline.R")
source("R_cbh/covariates/regional.R")
out <- "data/derived_cbh/regional_adjustment"
dir.create(file.path(out,"surveys"),recursive=TRUE,showWarnings=FALSE)
cfg <- cbh_config()
registry <- cbh_read_csv(cfg$registry)
rules <- cbh_read_csv(cfg$survey_rules)
boundaries <- unique(cbh_read_csv(cfg$boundary_regions)[c("svkey","region","regkey")])
overrides <- cbh_read_csv(cfg$region_overrides)
meta <- readRDS(file.path(cfg$output_dir,"manifest.rds"))
stopifnot(meta$complete)
m <- meta$manifest[meta$manifest$status %in% c("built","cached"),]
outputs <- vector("list",nrow(m))
for(i in seq_len(nrow(m))) {
  s <- registry[registry$svkey==m$survey[i],,drop=FALSE]
  rule <- rules[rules$svkey==s$svkey,,drop=FALSE]
  boundary <- boundaries[boundaries$svkey==s$svkey,,drop=FALSE]
  paths <- c(s$local_recode,cfg$survey_rules,cfg$region_overrides,cfg$boundary_regions,
    registry$local_recode[registry$svkey %in% rule$group_donor],
    "R_cbh/covariates/regional.R","R_cbh/covariates/02_extract_regional.R",
    "R_cbh/R/geography.R","R_cbh/R/utils.R","R_cbh/R/child_bands.R")
  hashes <- data.frame(file=paths,md5=vapply(paths,cbh_file_hash,""))
  sig <- cbh_hash(list(s,hashes))
  path <- file.path(out,"surveys",paste0(s$svkey,".rds"))
  old <- if(file.exists(path)) readRDS(path) else NULL
  if(!is.null(old) && identical(sig,old$signature)) outputs[[i]] <- old else {
    br <- cbh_read_recode(s$local_recode,c(rule$region_var,"b19","b9","m4","m39a",
      "v409","v409a","v410","v410a","v411","v411a","v412","v412a","v412b","v412c",
      paste0("v413",c("","a","b","c","d")),paste0("v414",letters[1:23])))
    geo <- cbh_geography(br,s,rule,boundary,overrides,registry)
    x <- cbh_regional_recode(br,s,rule,geo)
    x$signature <- sig; x$source_hashes <- hashes
    cbh_atomic_rds(x,path); outputs[[i]] <- x
    rm(br,geo,x); gc(FALSE)
  }
  message(i,"/",nrow(m),": ",s$svkey)
}
cbh_atomic_csv(cbh_bind(lapply(outputs,`[[`,"values")),file.path(out,"recode_regional.csv"))
cbh_atomic_csv(cbh_bind(lapply(outputs,`[[`,"feeding")),file.path(out,"feeding_extraction.csv"))
cbh_atomic_csv(unique(cbh_bind(lapply(outputs,`[[`,"source_hashes"))),file.path(out,"recode_source_manifest.csv"))

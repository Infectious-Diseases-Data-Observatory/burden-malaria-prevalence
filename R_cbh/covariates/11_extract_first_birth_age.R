#!/usr/bin/env Rscript
# Availability audit only; preserve current modelling inputs and fits.
source("R_cbh/load_pipeline.R")
cfg <- cbh_config()
out <- file.path(cfg$output_dir,"regional_adjustment/planned17_audit")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
reg <- cbh_read_csv(cfg$registry); rules <- cbh_read_csv(cfg$survey_rules)
bounds <- unique(cbh_read_csv(cfg$boundary_regions)[c("svkey","region","regkey")])
ov <- cbh_read_csv(cfg$region_overrides)
meta <- readRDS(file.path(cfg$output_dir,"manifest.rds"));stopifnot(meta$complete)
m <- meta$manifest[meta$manifest$status %in% c("built","cached"),]
rows <- list(); paths <- c(cfg$registry,cfg$survey_rules,cfg$boundary_regions,cfg$region_overrides,
  "R_cbh/covariates/11_extract_first_birth_age.R","R_cbh/R/utils.R","R_cbh/R/geography.R")
for(i in seq_len(nrow(m))) {
  s <- reg[reg$svkey==m$survey[i],,drop=FALSE]
  rule <- rules[rules$svkey==s$svkey,,drop=FALSE]
  br <- cbh_read_recode(s$local_recode,c(rule$region_var,"v212"))
  geo <- cbh_geography(br,s,rule,bounds[bounds$svkey==s$svkey,],ov,reg)
  age <- cbh_column(br,"v008",TRUE)-cbh_column(br,"b3",TRUE)
  recent <- is.finite(age) & age>=0 & age<60
  mother <- trimws(cbh_column(br,"caseid"))
  if(all(is.na(mother)))mother <- cbh_key(br,c("v001","v002","v003"))
  valid_id <- !is.na(mother) & nzchar(mother)
  first <- recent & valid_id & !duplicated(ifelse(recent & valid_id,mother,NA_character_))
  w <- cbh_column(br,"v005",TRUE)/1e6
  v <- cbh_column(br,"v212",TRUE)
  # DHS age is in completed years. Missing/special codes are not ages.
  v[!is.finite(v) | v<8 | v>49] <- NA_real_
  for(r in unique(na.omit(geo$regkey))) {
    ix <- which(geo$regkey==r & first & is.finite(w) & w>0)
    obs <- ix[is.finite(v[ix])]
    rows[[length(rows)+1L]] <- data.frame(survey=s$svkey,regkey=r,
      mean_maternal_age_first_birth=if(length(obs))weighted.mean(v[obs],w[obs]) else NA_real_,
      eligible_mothers=length(ix),observed_mothers=length(obs),missing_mothers=length(ix)-length(obs),
      v212_available="v212" %in% names(br))
  }
  paths <- c(paths,s$local_recode,reg$local_recode[reg$svkey %in% rule$group_donor])
  rm(br,geo);gc(FALSE)
  if(i%%10==0)message("First-birth age: ",i,"/",nrow(m)," surveys")
}
cbh_atomic_csv(cbh_bind(rows),file.path(out,"first_birth_age.csv"))
paths <- unique(paths)
cbh_atomic_csv(data.frame(file=paths,md5=vapply(paths,cbh_file_hash,"")),file.path(out,"first_birth_age_provenance.csv"))
message("First-birth-age regional extraction complete")

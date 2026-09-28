#!/usr/bin/env Rscript
# Regional covariates with DHS indicator definitions from local recodes, for surveys
# whose published subnational rows are excluded (published_survey_exclusions.csv).
# Run after 02 and before 03/13. Licensed recodes are read in place; only regional
# and national aggregates are written, under data/derived_cbh/regional_adjustment.
source("R_cbh/load_pipeline.R")
library(data.table)
cfg <- cbh_config()
out <- file.path(cfg$output_dir,"regional_adjustment")
reg <- cbh_read_csv(cfg$registry); rules <- cbh_read_csv(cfg$survey_rules)
bounds <- unique(cbh_read_csv(cfg$boundary_regions)[c("svkey","region","regkey")])
overrides <- cbh_read_csv(cfg$region_overrides)
excl <- cbh_read_csv("R_cbh/covariates/published_survey_exclusions.csv")
stopifnot(!anyDuplicated(excl$survey),all(excl$survey %in% reg$svkey))
pub <- cbh_read_csv(file.path(out,"published_indicators.csv"))
nutpub <- cbh_read_csv(file.path(out,"planned17_audit/nutrition_published.csv"))
statc <- cbh_read_csv("data/derived_dhs/statcompiler_covariates.csv",required=FALSE)
# Survey-specific inputs. Nutrition national rows are not in the regional API cache:
# the benchmark is the denominator-weighted mean of an exhaustive published partition.
spec <- list(SL61FL=list(nutrition_partition=c("Eastern","Northern","North Western","Southern")))
stopifnot(setequal(names(spec),excl$survey))

# Classified labels; any other label stops. Missing/don't know stay in the
# denominator as "no", as in the published DHS tables.
vacc_yes <- c("vaccination date on card","reported by mother","vaccination marked on card")
vacc_no <- c("no","don't know","missing")
facility <- c("government hospital","government health center","government health post",
  "other public sector","private hospital/clinic","other private sector")
nonfacility <- c("respondent's home","other home","other","missing")
# JMP 2017 ladder as in current DHS tables: bottled/sachet, tanker and cart are improved.
water_imp <- c("piped into dwelling","piped to yard/plot","public tap/standpipe","tube well or borehole",
  "protected well","protected spring","rainwater","tanker truck","cart with small tank","bottled water or sachets")
water_unimp <- c("unprotected well","unprotected spring","river/dam/lake/ponds/stream/canal/irrigation channel","other","missing")
# Improved facility type regardless of sharing (WS_TLET_H_IMP).
san_imp <- c("flush to piped sewer system","flush to septic tank","flush to pit latrine","flush, don't know where",
  "ventilated improved pit latrine (vip)","pit latrine with slab","composting toilet")
san_unimp <- c("flush to somewhere else","pit latrine without slab/open pit","bucket toilet",
  "hanging toilet/latrine","no facility/bush/field","other","missing")
anthro_invalid <- c("flagged cases","height out of plausible limits","weight out of plausible limits",
  "age in days out of plausible limits","missing")
lab <- function(d,v) tolower(trimws(cbh_column(d,v)))
classify <- function(x,mask,sets,what) {
  bad <- setdiff(unique(x[mask]),unlist(sets))
  if(length(bad)) stop(what,": unclassified labels: ",paste(bad,collapse=", "))
}

# One weighted value per boundary region plus the national value; every
# denominator row must have a reviewed region and a classified response.
summarise <- function(s,variable,population,regkey,mask,value,w,missing_coded_no) {
  ok <- mask & is.finite(w) & w>0
  if(anyNA(regkey[ok])) stop(s$svkey,": ",variable," denominator rows without a reviewed region")
  if(anyNA(value[ok])) stop(s$svkey,": ",variable," has an unclassified denominator row")
  expected <- sort(unique(bounds$regkey[bounds$svkey==s$svkey]))
  if(!identical(sort(unique(regkey[ok])),expected)) stop(s$svkey,": ",variable," does not cover every boundary region")
  regional <- cbh_bind(lapply(expected,function(r) {ix <- which(ok & regkey==r)
    data.frame(survey=s$svkey,country=s$iso3,regkey=r,variable=variable,
      value=weighted.mean(value[ix],w[ix]),eligible_n=length(ix),observed_n=length(ix),missing_n=0L,
      weighted_n=sum(w[ix]),population=population,source="DHS_recode_published_equivalent",
      small_denominator=length(ix)<25,low_precision=length(ix)<50)}))
  ix <- which(ok)
  list(rows=regional,nat=data.frame(survey=s$svkey,variable=variable,recode_national=weighted.mean(value[ix],w[ix]),
    observed_n=length(ix),weighted_n=sum(w[ix]),missing_coded_no=sum(missing_coded_no[ix])))
}

results <- sources <- list()

for(sv in excl$survey) {
  s <- reg[reg$svkey==sv,,drop=FALSE]; rule <- rules[rules$svkey==sv,,drop=FALSE]
  stopifnot(nrow(s)==1L,nrow(rule)==1L,rule$region_var=="v024",is.na(rule$group_donor) | !nzchar(rule$group_donor))
  boundary <- bounds[bounds$svkey==sv,,drop=FALSE]
  pr_path <- file.path(dirname(s$local_recode),sub("^([A-Za-z]{2})BR","\\1PR",basename(s$local_recode)))
  stopifnot(file.exists(s$local_recode),file.exists(pr_path),pr_path!=s$local_recode)
  sources[[sv]] <- data.frame(survey=sv,file=c(s$local_recode,pr_path),md5=c(cbh_file_hash(s$local_recode),cbh_file_hash(pr_path)))

  # Births (BR, women's weight v005): vaccination, facility delivery, birth interval.
  br <- cbh_read_recode(s$local_recode,c(rule$region_var,"b19","h3","h5","h7","h9","m15"))
  geo <- cbh_geography(br,s,rule,boundary,overrides,reg)
  num <- function(v) cbh_column(br,v,TRUE)
  age <- if("b19" %in% names(br) && any(is.finite(num("b19")))) num("b19") else num("v008")-num("b3")
  w <- num("v005")/1e6
  alive <- cbh_decode_category(cbh_column(br,"b5"),c("no","yes"),c(0,1))
  stopifnot(!anyNA(alive))
  child <- is.finite(age) & age>=12 & age<=23 & alive==1
  # DTP3 = at least 3 recorded DTP doses (h3, h5, h7), as in the DHS tabulation
  # code: a third dose without an earlier one does not count (SL61FL national
  # 77.85 against published 77.9; h7 alone gives 78.45).
  x <- lapply(c(h3="h3",h5="h5",h7="h7",h9="h9"),function(v) {y <- lab(br,v); classify(y,child,list(vacc_yes,vacc_no),paste(sv,v)); y})
  doses <- Reduce(`+`,lapply(x[c("h3","h5","h7")],function(y) y %in% vacc_yes))
  unknown <- Reduce(`|`,lapply(x[c("h3","h5","h7")],function(y) y %in% c("don't know","missing")))
  results[[length(results)+1L]] <- summarise(s,"dtp3_pct","living_children_12_23_months",
    geo$regkey,child,100*(doses>=3),w,unknown & doses<3)
  results[[length(results)+1L]] <- summarise(s,"measles_pct","living_children_12_23_months",
    geo$regkey,child,100*(x$h9 %in% vacc_yes),w,x$h9 %in% c("don't know","missing"))
  births <- is.finite(age) & age>=0 & age<60
  x <- lab(br,"m15"); classify(x,births,list(facility,nonfacility),paste(sv,"m15"))
  results[[length(results)+1L]] <- summarise(s,"facility_delivery_pct","live_births_last_60_months",geo$regkey,
    births,100*(x %in% facility),w,x %in% "missing")
  b11 <- num("b11")
  results[[length(results)+1L]] <- summarise(s,"short_birth_interval_pct","non_first_births_last_60_months",geo$regkey,
    births & is.finite(b11),100*(b11>=7 & b11<=23),w,rep(FALSE,nrow(br)))
  rm(br,geo); gc(FALSE)

  # Households and anthropometry (PR, household weight hv005). No HR file: one PR
  # row per household (hv001, hv002) carries the household items.
  pr <- readRDS(pr_path); names(pr) <- tolower(names(pr))
  cbh_require(pr,c("hv001","hv002","hv005","hv024","hv103","hv201","hv205","hv206","hc1","hc70","hc72"),paste(sv,"PR"))
  geo <- cbh_geography(data.frame(v024=cbh_column(pr,"hv024")),s,rule,boundary,overrides,reg)
  w <- cbh_column(pr,"hv005",TRUE)/1e6
  hh <- !duplicated(paste(cbh_column(pr,"hv001",TRUE),cbh_column(pr,"hv002",TRUE)))
  hhpop <- "households_PR_one_row_per_hv001_hv002"
  x <- lab(pr,"hv206"); classify(x,hh,list("yes",c("no","missing")),paste(sv,"hv206"))
  results[[length(results)+1L]] <- summarise(s,"electricity_pct",hhpop,geo$regkey,hh,100*(x=="yes"),w,x=="missing")
  x <- lab(pr,"hv201"); classify(x,hh,list(water_imp,water_unimp),paste(sv,"hv201"))
  results[[length(results)+1L]] <- summarise(s,"improved_water_pct",hhpop,geo$regkey,hh,100*(x %in% water_imp),w,x=="missing")
  x <- lab(pr,"hv205"); classify(x,hh,list(san_imp,san_unimp),paste(sv,"hv205"))
  results[[length(results)+1L]] <- summarise(s,"improved_sanitation_pct",hhpop,geo$regkey,hh,100*(x %in% san_imp),w,x=="missing")
  # DHS nutrition: de facto children 0-59 months with a valid z-score (flagged,
  # implausible and missing measurements are labels here and excluded).
  hc1 <- cbh_column(pr,"hc1",TRUE)
  under5 <- lab(pr,"hv103")=="yes" & is.finite(hc1) & hc1>=0 & hc1<60
  for(v in c("hc70","hc72")) {
    z <- cbh_column(pr,v,TRUE); text <- lab(pr,v); coded <- !is.finite(z) & !is.na(text)
    classify(text,under5 & coded,list(anthro_invalid),paste(sv,v))
    valid <- under5 & is.finite(z) & z<9996
    results[[length(results)+1L]] <- summarise(s,c(hc70="stunting_pct",hc72="wasting_pct")[[v]],
      "de_facto_children_0_59_months_valid_measurement",geo$regkey,valid,100*(z< -200),w,rep(FALSE,nrow(pr)))
  }
  rm(pr,geo); gc(FALSE)
}
rows <- as.data.table(cbh_bind(lapply(results,`[[`,"rows"))); nat <- as.data.table(cbh_bind(lapply(results,`[[`,"nat")))
cbh_unique(as.data.frame(rows),c("survey","regkey","variable"),"Recode replacements")

# Published national benchmarks from the cached DHS API rows; the recode national
# value must lie within 2 points (published regions are unusable, totals are not).
ind <- c(dtp3_pct="CH_VACC_C_DP3",measles_pct="CH_VACC_C_MSL",facility_delivery_pct="RH_DELP_C_DHF",
  improved_water_pct="WS_SRCE_H_IMP",improved_sanitation_pct="WS_TLET_H_IMP",electricity_pct="HC_ELEC_H_ELC")
bench <- list()
for(sv in excl$survey) {
  id <- reg$SurveyId[reg$svkey==sv]
  x <- pub[pub$SurveyId==id & pub$level=="national",,drop=FALSE]
  one <- function(i,by="") {y <- x[x$IndicatorId==i & (if(nzchar(by)) x$ByVariableLabel %in% by else is.na(x$ByVariableLabel) | x$ByVariableLabel==""),]
    stopifnot(nrow(y)==1L); as.numeric(y$Value)}
  for(v in names(ind)) bench[[length(bench)+1L]] <- data.frame(survey=sv,variable=v,published_national=
    one(ind[[v]],if(v=="facility_delivery_pct")"Five years preceding the survey" else ""),benchmark="DHS API national row")
  bench[[length(bench)+1L]] <- data.frame(survey=sv,variable="short_birth_interval_pct",
    published_national=one("FE_BINT_C_I07")+one("FE_BINT_C_I18"),benchmark="DHS API national rows FE_BINT_C_I07 + FE_BINT_C_I18")
  # Exhaustive: the partition's unweighted denominators equal the recode's measured children.
  part <- spec[[sv]]$nutrition_partition
  for(i in c("CN_NUTS_C_HA2","CN_NUTS_C_WH2")) {
    v <- c(CN_NUTS_C_HA2="stunting_pct",CN_NUTS_C_WH2="wasting_pct")[[i]]
    y <- nutpub[nutpub$SurveyId==id & nutpub$IndicatorId==i & nutpub$CharacteristicLabel %in% part,]
    stopifnot(nrow(y)==length(part),sum(y$DenominatorUnweighted)==nat[survey==sv & variable==v,observed_n])
    bench[[length(bench)+1L]] <- data.frame(survey=sv,variable=v,
      published_national=sum(y$Value*y$DenominatorWeighted)/sum(y$DenominatorWeighted),
      benchmark=paste("Denominator-weighted exhaustive published regions",paste(part,collapse="/")))
    # Cross-check against the cached StatCompiler national wasting row (R_dhs 02c).
    z <- if(is.null(statc)) statc else statc[statc$SurveyId==id & statc$IndicatorId==i & statc$level=="national",]
    if(!is.null(z) && nrow(z)) stopifnot(nrow(z)==1L,abs(z$Value-bench[[length(bench)]]$published_national)<0.1)
  }
}
nat <- merge(nat,rbindlist(bench),by=c("survey","variable"),all=TRUE)
nat[,difference_pp:=recode_national-published_national]
stopifnot(nrow(nat)==9L*nrow(excl),!anyNA(nat$difference_pp))
print(nat[,.(survey,variable,recode_national=round(recode_national,1),published_national=round(published_national,1),
  difference_pp=round(difference_pp,2),observed_n)])
if(any(abs(nat$difference_pp)>2)) stop("Recode national value differs from the published national value by more than 2 points")

nutvars <- c("stunting_pct","wasting_pct")
cols <- c("survey","regkey","variable","value","eligible_n","observed_n","missing_n","weighted_n",
  "population","source","small_denominator","low_precision")
cbh_atomic_csv(rows[!variable %in% nutvars],file.path(out,"recode_replacements.csv"))
cbh_atomic_csv(rows[variable %in% nutvars,cols,with=FALSE],file.path(out,"recode_replacements_nutrition.csv"))
cbh_atomic_csv(nat,file.path(out,"recode_replacements_national.csv"))
paths <- c("R_cbh/covariates/17_recode_replacements.R","R_cbh/covariates/published_survey_exclusions.csv",
  cfg$registry,cfg$survey_rules,cfg$region_overrides,cfg$boundary_regions,
  file.path(out,c("published_indicators.csv","planned17_audit/nutrition_published.csv")),
  "data/derived_dhs/statcompiler_covariates.csv","R_cbh/R/geography.R","R_cbh/R/utils.R")
cbh_atomic_csv(rbind(rbindlist(sources),data.frame(survey="",file=paths,md5=vapply(paths,cbh_file_hash,""))),
  file.path(out,"recode_replacements_sources.csv"))
message("Recode replacements: ",nrow(rows)," survey-region values for ",nrow(excl)," survey(s); national checks within 2 points.")

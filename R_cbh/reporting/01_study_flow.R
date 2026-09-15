#!/usr/bin/env Rscript
# Primary CBH flow, adapted from R_dhs/11_study_flow.R. Aggregate inputs only.
source("R_cbh/load_pipeline.R")
library(ggplot2)
source("R_cbh/primary/settings.R")
root <- cbh_primary_settings()$out
out <- file.path(root,"study_flow")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
input_paths <- c(
  survey_manifest="data/derived_cbh/survey_manifest.csv",
  child_checks="data/derived_cbh/child_checks.csv",
  eligibility="data/derived_cbh/eligibility_flow.csv",
  selection="results/cbh/age_band_hiv_incidence_shared_time_v3/selection.csv",
  diagnostics=file.path(root,"fit_diagnostics.csv"),
  fits=file.path(root,"fit_manifest.csv"))
m <- cbh_read_csv(input_paths[["survey_manifest"]])
c <- cbh_read_csv(input_paths[["child_checks"]])
e <- cbh_read_csv(input_paths[["eligibility"]])
s <- cbh_read_csv(input_paths[["selection"]])
d <- cbh_read_csv(input_paths[["diagnostics"]]);d <- d[d$series=="map_full",]
f <- cbh_read_csv(input_paths[["fits"]]);f <- f[f$series=="map_full",]
ages <- cbh_config()$age_bands$age_band
cbh_unique(m,"survey","Survey manifest")
cbh_unique(e,c("survey","age_band"),"Eligibility ledger")
cbh_unique(s,"survey","Model selection")
stopifnot(nrow(d)==7,nrow(f)==7,setequal(d$age_band,ages),all(d$gamma==2),
  all(d$converged),all(d$input_verified),!any(m$status=="failed"))
built <- m[m$status %in% c("built","cached"),]
omitted <- m[!m$status %in% c("built","cached"),]
stopifnot(all(omitted$status=="missing_map_geography"),setequal(built$survey,s$survey),
  setequal(unique(e$survey),built$survey),all(e$entered_in_lookback==e$incomplete_potential_band+e$full_band_before_interview),
  all(e$full_band_before_interview==e$entry_year_outside_range+e$eligible_band_rows))
stopifnot(setequal(unique(c$survey),built$survey),
  sum(c$children[c$status=="valid"])==sum(e$valid_children_reaching_band[e$age_band=="<1"]))
# Reconcile every survey as well as final totals; no disjoint-covariate assumption.
ag <- aggregate(e$eligible_band_rows,list(survey=e$survey),sum)
j <- match(s$survey,ag$survey)
k <- match(s$survey,built$survey)
stopifnot(identical(as.numeric(s$eligible_rows),as.numeric(ag$x[j])),
  all(s$eligible_rows==built$rows[k]),all(s$pfpr_available_rows==built$model_ready_rows[k]),
  all(s$pfpr_available_rows==s$excluded_incomplete+s$complete_case_rows),
  sum(d$rows)==sum(s$complete_case_rows),sum(d$deaths)==sum(s$complete_case_deaths),
  all(d$surveys==sum(s$complete_case_rows>0)),
  all(d$countries==length(unique(s$country[s$complete_case_rows>0]))))
n <- list(registry=nrow(m),registry_countries=length(unique(m$country)),omitted=nrow(omitted),
  built=nrow(built),built_countries=length(unique(built$country)),entered=sum(e$entered_in_lookback),
  incomplete=sum(e$incomplete_potential_band),outside_year=sum(e$entry_year_outside_range),
  eligible=sum(e$eligible_band_rows),missing_map=sum(s$eligible_rows-s$pfpr_available_rows),
  missing_covariates=sum(s$excluded_incomplete),primary=sum(s$complete_case_rows),
  deaths=sum(s$complete_case_deaths),surveys=sum(s$complete_case_rows>0),
  countries=length(unique(s$country[s$complete_case_rows>0])))
n$excluded_eligibility <- n$incomplete+n$outside_year
n$excluded_missing <- n$missing_map+n$missing_covariates
n$recorded_births <- sum(c$children)
n$valid_births <- sum(c$children[c$status=="valid"])
n$recent_valid_births <- sum(e$entered_in_lookback[e$age_band=="<1"])
stopifnot(n$entered-n$excluded_eligibility==n$eligible,n$eligible-n$excluded_missing==n$primary)
fmt <- function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE)
counts <- data.frame(item=names(n),count=unlist(n),row.names=NULL)
cbh_atomic_csv(counts,file.path(out,"flow_counts.csv"))
cbh_atomic_csv(omitted[c("survey","country","survey_year","status")],file.path(out,"excluded_surveys.csv"))
cbh_atomic_csv(s,file.path(out,"survey_selection.csv"))
# Coordinates leave substantial white space; 16-19 pt text at native 11-inch width.
# Main labels contain only stages and counts. Eligibility/covariate details are in the caption.
nodes <- data.frame(
  id=c("registry","entered","eligible","sample","omit_surveys","omit_eligibility","omit_missing"),
  x=c(rep(3.8,4),9.4,9.4,9.4),
  y=c(11.45,9.35,7.25,5.15,10.4,8.3,6.2),
  w=c(rep(6.9,4),3.65,3.65,3.65),
  h=c(1.35,1.35,1.35,1.65,1.2,1.2,1.2),
  role=c(rep("main",4),rep("excluded",3)),
  heading=c("Survey registry","Band entries within 5 years","Eligible child–age bands","Primary analysis sample",
    paste(fmt(n$omitted),"surveys excluded"),paste(fmt(n$excluded_eligibility),"excluded"),paste(fmt(n$excluded_missing),"excluded")),
  detail=c(sprintf("%s surveys · %s countries",fmt(n$registry),fmt(n$registry_countries)),
    paste(fmt(n$entered),"child–age bands"),paste(fmt(n$eligible),"records"),
    sprintf("%s records · %s deaths\n%s surveys · %s countries",fmt(n$primary),fmt(n$deaths),fmt(n$surveys),fmt(n$countries)),
    "No MAP geography","Band or year ineligible","Missing MAP or covariates"))
cbh_atomic_csv(nodes,file.path(out,"figure_labels.csv"))
fill <- c(main="#EAF2F8",excluded="#F6F3EF")
border <- c(main="#346A8C",excluded="#A29588")
p <- ggplot()+theme_void()+coord_cartesian(xlim=c(0,11.5),ylim=c(4.1,12.25),expand=FALSE,clip="off")+
  theme(plot.margin=margin(8,12,8,12))
segment <- function(x,y,xe,ye,arrow=FALSE) {
  if(arrow) annotate("segment",x=x,y=y,xend=xe,yend=ye,colour="#687780",linewidth=.7,
    arrow=grid::arrow(length=grid::unit(2.5,"mm"),type="closed"))
  else annotate("segment",x=x,y=y,xend=xe,yend=ye,colour="#687780",linewidth=.7)
}
for(i in 1:3) p <- p+segment(3.8,nodes$y[i]-nodes$h[i]/2,3.8,nodes$y[i+1]+nodes$h[i+1]/2,TRUE)
# Exclusion branches leave the connecting lines, not the retained sample boxes.
for(y in c(10.4,8.3,6.2)) p <- p+segment(3.8,y,9.4-3.65/2,y,TRUE)
for(i in seq_len(nrow(nodes))) {
  z <- nodes[i,]; small <- z$role=="excluded"
  p <- p+annotate("rect",xmin=z$x-z$w/2,xmax=z$x+z$w/2,ymin=z$y-z$h/2,ymax=z$y+z$h/2,
    fill=fill[z$role],colour=border[z$role],linewidth=.75)
  has_detail <- nzchar(z$detail)
  title_y <- z$y+if(has_detail) if(z$id=="sample") .43 else .22 else 0
  detail_y <- z$y+if(z$id=="sample") -.22 else -.23
  p <- p+annotate("text",x=z$x,y=title_y,label=z$heading,fontface="bold",family="sans",
    size=if(small) 5.5 else 6.15,lineheight=1.1,colour="#183746")
  if(has_detail) p <- p+annotate("text",x=z$x,y=detail_y,label=z$detail,family="sans",
    size=if(small) 5.25 else 5.8,lineheight=1.16,colour="#243D49")
}
ggsave(file.path(out,"study_flow_diagram.png"),p,width=11,height=7.6,dpi=300,device=ragg::agg_png,bg="white")
caption <- c("# Figure caption — primary sample inclusion","",
  sprintf("**Figure. Sample inclusion for the primary MAP analysis.** The survey registry contains %s surveys in %s countries. %s surveys lack MAP geography, leaving %s processed surveys in %s countries. Counts below the registry refer to child–age-band records, not unique children, and begin after birth-history validity checks and confirmation that the child reached the band alive.",fmt(n$registry),fmt(n$registry_countries),fmt(n$omitted),fmt(n$built),fmt(n$built_countries)),"",
  sprintf("The processed surveys contain %s recorded births across their complete birth histories; %s pass history-validity checks, including %s births within 60 months before interview. These birth totals exclude the surveys skipped for missing MAP geography. Children born earlier can still contribute later age-band entries within the five-year window.",fmt(n$recorded_births),fmt(n$valid_births),fmt(n$recent_valid_births)),"",
  sprintf("Among %s band entries within the 60 months before interview, %s have an incomplete potential band and %s have entry years outside 2000–2024. Requiring the full potential band to end by interview for deaths and survivors alike leaves %s eligible records. A further %s lack usable regional MAP/geography, and %s lack at least one required covariate after joining child HIV incidence. Covariate exclusions count records once, even if several values are missing. The final sample contains %s records and %s deaths from %s surveys in %s countries.",fmt(n$entered),fmt(n$incomplete),fmt(n$outside_year),fmt(n$eligible),fmt(n$missing_map),fmt(n$missing_covariates),fmt(n$primary),fmt(n$deaths),fmt(n$surveys),fmt(n$countries)),"",
  "Seven completed-month bands are analyzed separately: <1, 1–5, 6–11, 12–23, 24–35, 36–47 and 48–59. The outcome is death during the band. Annual regional MAP PfPR₂–₁₀ and national covariates are assigned at band-entry year; exposure is not averaged over time until death. There is no prevalence floor or requirement for a death in each region/band. Eligible MAP records after 2015 are retained.","",
  "Adjustment variables are sex, multiple birth, birth order, maternal age, maternal education, wealth quintile, urban residence, log child HIV incidence, log GDP per capita, log health expenditure per capita and political stability. Child HIV incidence uses the fixed posterior-median imputation informed by adolescent incidence; other model covariates require complete cases. Vaccines and maternal death fraction are excluded.","",
  "Counts are read from the saved build/eligibility/selection ledgers and reconciled with the seven selected full-sample MAP gamma=2 fits. See [counts](flow_counts.csv), [survey selection](survey_selection.csv), [provenance](provenance.csv) and the [analysis plan](../../../../docs/ANALYSIS_PLAN.md).")
writeLines(caption,file.path(out,"CAPTION.md"))
provenance <- c(unname(input_paths),"R_cbh/reporting/01_study_flow.R","R_cbh/00_config.R")
cbh_atomic_csv(data.frame(file=provenance,md5=vapply(provenance,cbh_file_hash,"")),file.path(out,"provenance.csv"))
message("Saved current primary study flow and caption: ",out)

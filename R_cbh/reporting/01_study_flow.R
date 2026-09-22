#!/usr/bin/env Rscript
# Primary CBH flow, adapted from R_dhs/11_study_flow.R. Aggregate inputs only.
source("R_cbh/load_pipeline.R")
library(ggplot2)
source("R_cbh/primary/settings.R")
settings <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional"))
root <- settings$out
out <- file.path(root,"study_flow")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
input_paths <- c(
  survey_manifest="data/derived_cbh/survey_manifest.csv",
  child_checks="data/derived_cbh/child_checks.csv",
  eligibility="data/derived_cbh/eligibility_flow.csv",
  selection=if(settings$regional) file.path(root,"prepared_selection_by_survey.csv") else "results/cbh/age_band_hiv_incidence_shared_time_v3/selection.csv",
  sample=file.path(root,"primary_sample.csv"),
  survey_coverage=file.path(root,"survey_coverage.csv"),
  diagnostics=file.path(root,"fit_diagnostics.csv"),
  fits=file.path(root,"fit_manifest.csv"),
  attribution="results/cbh/planned17_covariate_missingness/exclusion_attribution.csv")
m <- cbh_read_csv(input_paths[["survey_manifest"]])
c <- cbh_read_csv(input_paths[["child_checks"]])
e <- cbh_read_csv(input_paths[["eligibility"]])
s <- cbh_read_csv(input_paths[["selection"]])
if(settings$regional) {
  j <- match(s$survey,m$survey);stopifnot(!anyNA(j),all(s$country==m$country[j]))
  s$eligible_rows <- m$rows[j];s$pfpr_available_rows <- m$model_ready_rows[j]
  s$complete_case_rows <- s$records;s$complete_case_deaths <- s$deaths
  s$excluded_incomplete <- s$pfpr_available_rows-s$complete_case_rows
  stopifnot(all(s$excluded_incomplete>=0))
}
d <- cbh_read_csv(input_paths[["diagnostics"]]);d <- d[d$series=="map_full",]
f <- cbh_read_csv(input_paths[["fits"]]);f <- f[f$series=="map_full",]
sample <- cbh_read_csv(input_paths[["sample"]])
coverage <- cbh_read_csv(input_paths[["survey_coverage"]])
ages <- cbh_config()$age_bands$age_band
cbh_unique(m,"survey","Survey manifest")
cbh_unique(e,c("survey","age_band"),"Eligibility ledger")
cbh_unique(s,"survey","Model selection")
stopifnot(nrow(d)==7,nrow(f)==7,setequal(d$age_band,ages),all(d$gamma==2),
  all(d$converged),all(d$input_verified),!any(m$status=="failed"))
# Confirm that these are still the selected primary prepared inputs. Hashing does
# not inspect individual records; all displayed counts come from aggregate ledgers.
prepared_path <- settings$data
stopifnot(all(f$prepared_data_md5==cbh_file_hash(prepared_path)),nrow(sample)==1L,
  sample$records==sum(d$rows),sample$deaths==sum(d$deaths),
  sample$surveys==nrow(coverage),sample$countries==length(unique(coverage$country)),
  setequal(coverage$survey,s$survey[s$complete_case_rows>0]))
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
n$primary_distinct_children <- sample$distinct_children
stopifnot(n$entered-n$excluded_eligibility==n$eligible,n$eligible-n$excluded_missing==n$primary)
fmt <- function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE)
counts <- data.frame(item=names(n),count=unlist(n),row.names=NULL)
cbh_atomic_csv(counts,file.path(out,"flow_counts.csv"))
cbh_atomic_csv(omitted[c("survey","country","survey_year","status")],file.path(out,"excluded_surveys.csv"))
cbh_atomic_csv(s,file.path(out,"survey_selection.csv"))
# Disjoint attribution of covariate exclusions to the first missing covariate
# (R_cbh/covariates/14_exclusion_attribution.R); its total must equal the ledger.
attr <- cbh_read_csv(input_paths[["attribution"]])
stopifnot(sum(attr$records)==n$missing_covariates)
registry <- cbh_read_csv(cbh_config()$registry)
stopifnot(all(s$survey %in% registry$svkey),all(registry$SurveyType[registry$svkey %in% s$survey[s$complete_case_rows>0]]=="DHS"))
attr <- attr[order(-attr$records),]
short <- c(wasting_pct="Wasting",stunting_pct="Stunting",facility_delivery_pct="Facility delivery",
  electricity_pct="Electricity",political_stability="Political stability",log_hiv_incidence="HIV incidence",
  log_health_expenditure_pc="Health expenditure",mean_wealth_quintile="Wealth score",log_gdp_pc="GDP")
items <- sprintf("%s %s",ifelse(attr$reason %in% names(short),short[attr$reason],attr$label),fmt(attr$records))
if(length(items)>6) items <- c(items[1:5],paste("other",fmt(sum(attr$records[-(1:5)]))))
reason_text <- paste(vapply(split(items,ceiling(seq_along(items)/2)),paste,"",collapse=" · "),collapse="\n")
omitted_countries <- paste(sort(unique(omitted$country)),collapse=", ")
# Spell the country out in the exclusion box; the code alone reads poorly in a sentence.
omitted_names <- paste(sort(unique(countrycode::countrycode(omitted$country,"iso3c","country.name",
  warn=FALSE))),collapse=", ")
stopifnot(nzchar(omitted_names),!anyNA(omitted_names))
n$with_map <- n$eligible-n$missing_map
# Five retained stages on the left; one box per exclusion reason on the right,
# each leaving the connecting line between the stages it separates.
nodes <- data.frame(
  id=c("registry","entered","eligible","with_map","sample","omit_surveys","omit_eligibility","omit_map","omit_covariates"),
  x=c(rep(3.7,5),rep(9.75,4)),
  y=c(13.75,11.55,9.4,7.25,5.0,12.69,10.44,8.33,6.2),
  w=c(rep(6.7,5),rep(4.7,4)),
  h=c(1.35,1.5,1.35,1.35,1.65,1.35,1.75,1.35,2.0),
  role=c(rep("main",5),rep("excluded",4)),
  heading=c("Survey registry","Band entries within 5 years of interview","Eligible child–age bands",
    "With regional MAP prevalence","Primary analysis sample",
    paste(fmt(n$omitted),"surveys excluded"),paste(fmt(n$excluded_eligibility),"records excluded"),
    paste(fmt(n$missing_map),"records excluded"),paste(fmt(n$missing_covariates),"records excluded")),
  detail=c(sprintf("%s surveys · %s countries",fmt(n$registry),fmt(n$registry_countries)),
    sprintf("%s child–age bands\n%s surveys · %s countries",fmt(n$entered),fmt(n$built),fmt(n$built_countries)),
    paste(fmt(n$eligible),"records"),paste(fmt(n$with_map),"records"),
    sprintf("%s records · %s deaths\n%s surveys · %s countries",fmt(n$primary),fmt(n$deaths),fmt(n$surveys),fmt(n$countries)),
    sprintf("%s is malaria free, so MAP\npublishes no prevalence surface",omitted_names),
    sprintf("%s: band not complete by interview\n%s: entry year outside 2000–2024",fmt(n$incomplete),fmt(n$outside_year)),
    "No MAP PfPR for the region\nin the band-entry year",
    paste0("Required covariate unavailable:\n",reason_text)))
cbh_atomic_csv(nodes,file.path(out,"figure_labels.csv"))
fill <- c(main="#EAF2F8",excluded="#F6F3EF")
border <- c(main="#346A8C",excluded="#A29588")
p <- ggplot()+theme_void()+coord_cartesian(xlim=c(0,12.3),ylim=c(4.0,14.55),expand=FALSE,clip="off")+
  theme(plot.margin=margin(8,12,8,12))
segment <- function(x,y,xe,ye,arrow=FALSE) {
  if(arrow) annotate("segment",x=x,y=y,xend=xe,yend=ye,colour="#687780",linewidth=.7,
    arrow=grid::arrow(length=grid::unit(2.5,"mm"),type="closed"))
  else annotate("segment",x=x,y=y,xend=xe,yend=ye,colour="#687780",linewidth=.7)
}
for(i in 1:4) p <- p+segment(3.7,nodes$y[i]-nodes$h[i]/2,3.7,nodes$y[i+1]+nodes$h[i+1]/2,TRUE)
for(i in 6:9) p <- p+segment(3.7,nodes$y[i],nodes$x[i]-nodes$w[i]/2,nodes$y[i],TRUE)
for(i in seq_len(nrow(nodes))) {
  z <- nodes[i,]; small <- z$role=="excluded"
  p <- p+annotate("rect",xmin=z$x-z$w/2,xmax=z$x+z$w/2,ymin=z$y-z$h/2,ymax=z$y+z$h/2,
    fill=fill[z$role],colour=border[z$role],linewidth=.75)
  n_lines <- length(strsplit(z$detail,"\n")[[1]])
  line_h <- if(small) .21 else .25
  title_y <- z$y+z$h/2-.36
  detail_y <- title_y-.30-line_h*n_lines/2
  p <- p+annotate("text",x=z$x,y=title_y,label=z$heading,fontface="bold",family="sans",
    size=if(small) 5.2 else 5.9,lineheight=1.1,colour="#183746")
  p <- p+annotate("text",x=z$x,y=detail_y,label=z$detail,family="sans",
    size=if(small) 4.4 else 5.4,lineheight=1.1,colour="#243D49")
}
ggsave(file.path(out,"study_flow_diagram.png"),p,width=11.5,height=9.6,dpi=300,device=ragg::agg_png,bg="white")
caption <- c("# Figure caption — primary sample inclusion","",
  sprintf("**Figure. Sample inclusion for the primary MAP analysis.** The survey registry contains %s surveys in %s countries. %s surveys are excluded because Lesotho is malaria free and the Malaria Atlas Project publishes no prevalence surface for it, leaving %s processed surveys in %s countries. Counts below the registry refer to child–age-band records, not unique children, and begin after birth-history validity checks and confirmation that the child reached the band alive.",fmt(n$registry),fmt(n$registry_countries),fmt(n$omitted),fmt(n$built),fmt(n$built_countries)),"",
  sprintf("The processed surveys contain %s recorded births across their complete birth histories; %s pass history-validity checks, including %s births within 60 months before interview. These birth totals exclude the surveys skipped for missing MAP geography. Children born earlier can still contribute later age-band entries within the five-year window.",fmt(n$recorded_births),fmt(n$valid_births),fmt(n$recent_valid_births)),"",
  sprintf("Among %s band entries within the 60 months before interview, %s have an incomplete potential band and %s have entry years outside 2000–2024. Requiring the full potential band to end by interview for deaths and survivors alike leaves %s eligible records. A further %s lack a regional MAP PfPR value for the band-entry year, leaving %s, and %s lack at least one required covariate after the declared HIV, vaccination and available-region substitutions. Each covariate exclusion is attributed to the first missing covariate in a declared order (national annual series, then regional summaries), so the reasons shown are disjoint: %s. The final sample contains %s records and %s deaths from %s surveys in %s countries.",fmt(n$entered),fmt(n$incomplete),fmt(n$outside_year),fmt(n$eligible),fmt(n$missing_map),fmt(n$with_map),fmt(n$missing_covariates),
    paste(sprintf("%s %s",tolower(attr$label),fmt(attr$records)),collapse="; "),fmt(n$primary),fmt(n$deaths),fmt(n$surveys),fmt(n$countries)),"",
  sprintf("%s processed surveys contribute no records to the final sample because a required covariate is unavailable for all of their regions; this includes all %s MIS surveys in the registry, which do not publish the facility-delivery and anthropometry indicators, so the analysed sample is DHS only.",fmt(sum(s$complete_case_rows==0)),fmt(sum(registry$SurveyType=="MIS"))),"",
  sprintf("These final records represent %s distinct children, including children born more than five years before interview whose later band entries meet the lookback rule. This figure describes the full primary sample, not the Sahel-only sensitivity subset.",fmt(n$primary_distinct_children)),"",
  "Seven completed-month bands are analyzed separately: <1, 1–5, 6–11, 12–23, 24–35, 36–47 and 48–59. The outcome is death during the band. Annual regional MAP PfPR₂–₁₀ and national covariates are assigned at band-entry year; exposure is not averaged over time until death. There is no prevalence floor or requirement for a death in each region/band. Eligible MAP records after 2015 are retained.","",
  if(settings$regional) "The 17 adjustment variables comprise survey-region summaries of maternal age at first birth, education, wealth, urban residence, DTP3/measles coverage, facility delivery, short birth interval, improved water, improved sanitation, electricity, wasting and stunting, plus national annual log child HIV incidence, log GDP per capita, log health expenditure per capita and political stability. Child HIV incidence uses the fixed posterior-median imputation informed by adolescent incidence. Missing regional DTP3/measles use exact country/survey-year UNICEF values; remaining regional gaps use the mean of available regions in the same survey. Whole-survey gaps remain missing. Hib3, PCV, rotavirus and exclusive breastfeeding are excluded." else "Adjustment variables are sex, multiple birth, birth order, maternal age, maternal education, wealth quintile, urban residence, log child HIV incidence, log GDP per capita, log health expenditure per capita and political stability. Child HIV incidence uses the fixed posterior-median imputation informed by adolescent incidence; other model covariates require complete cases. Vaccines and maternal death fraction are excluded.","",
  "Counts are read from the saved build/eligibility/selection ledgers and reconciled with the seven selected full-sample MAP gamma=2 fits. See [counts](flow_counts.csv), [survey selection](survey_selection.csv), [provenance](provenance.csv) and the [analysis plan](../../../../docs/ANALYSIS_PLAN.md).")
writeLines(caption,file.path(out,"CAPTION.md"))
provenance <- c(unname(input_paths),prepared_path,"R_cbh/reporting/01_study_flow.R","R_cbh/primary/settings.R","R_cbh/00_config.R")
cbh_atomic_csv(data.frame(file=provenance,md5=vapply(provenance,cbh_file_hash,"")),file.path(out,"provenance.csv"))
message("Saved current primary study flow and caption: ",out)

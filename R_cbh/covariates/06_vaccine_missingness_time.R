#!/usr/bin/env Rscript
# Aggregate existing audit only; no microdata reads, data refreshes or fits.
source("R_cbh/load_pipeline.R")
library(data.table)
library(ggplot2)
out <- "results/cbh/regional_adjustment_v1"
path <- file.path(out,"vaccine_source_status.csv")
d <- fread(path)
stopifnot(all(d$status %in% c("observed",
  "not_observed:pre_series_assumed_not_introduced",
  "not_observed:no_series_assumed_not_introduced")))
summarise <- function(x,groups) x[,.(records=sum(records),
  missing_records=sum(records[status!="observed"]),
  pre_series_records=sum(records[grepl("pre_series",status)]),
  no_series_records=sum(records[grepl("no_series",status)]),
  countries=uniqueN(country),surveys=uniqueN(survey)),by=groups]
annual <- summarise(d,c("year","variable"))
annual[,missing_pct:=100*missing_records/records]
stopifnot(all(annual[,.(denominators=uniqueN(records)),by=year]$denominators==1L),
  all(annual[,.(n=sum(records)),by=variable]$n==5885022L),
  all(annual$missing_records==annual$pre_series_records+annual$no_series_records))
expected <- fread(file.path(out,"missingness_summary.csv"))[baseline=="previous_primary"]
totals <- annual[,.(missing_records=sum(missing_records)),by=variable]
stopifnot(all(totals$missing_records==expected$missing_records[match(totals$variable,expected$variable)]))
d[,period_start:=5L*(year%/%5L)]
period <- summarise(d,c("period_start","variable"))
period[,period:=paste0(period_start,"–",pmin(period_start+4L,max(d$year)))]
period[,missing_pct:=100*missing_records/records]
country_year <- summarise(d,c("country","year","variable"))
country_year[,missing_pct:=100*missing_records/records]
setorder(annual,year,variable);setorder(period,period_start,variable)
cbh_atomic_csv(annual,file.path(out,"vaccine_missingness_by_entry_year.csv"))
cbh_atomic_csv(period,file.path(out,"vaccine_missingness_by_period.csv"))
cbh_atomic_csv(country_year,file.path(out,"vaccine_missingness_by_country_year.csv"))
# Explicit empty years break lines; an absent year is not 0% missingness.
grid <- CJ(year=2000:2024,variable=unique(d$variable))
plot_data <- merge(grid,annual,by=c("year","variable"),all.x=TRUE)
plot_data[,vaccine:=factor(variable,levels=c("hib3_pct","pcv3_pct","rotavirus_pct"),
  labels=c("Hib3","PCV","Rotavirus"))]
p <- ggplot(plot_data,aes(year,missing_pct,color=vaccine,group=vaccine))+
  geom_line(linewidth=1.1,na.rm=TRUE)+geom_point(size=2,na.rm=TRUE)+
  scale_color_manual(values=c(Hib3="#0072B2",PCV="#D55E00",Rotavirus="#009E73"))+
  scale_x_continuous(breaks=seq(2000,2024,4),limits=c(2000,2024))+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20),labels=function(x)paste0(x,"%"))+
  labs(x="Band-entry year",y="Records missing vaccine coverage",color=NULL)+
  theme_classic(base_size=16)+theme(legend.position="bottom",axis.title=element_text(size=18))
ggsave(file.path(out,"vaccine_missingness_over_time.png"),p,width=10,height=6,dpi=200)
wide <- dcast(period,period+records~variable,value.var="missing_pct")
lines <- c("# Vaccine missingness over time","",
  "Missingness among the 5,885,022 records in the previously fitted primary sample, grouped by **band-entry year**, not survey year. These are unweighted record percentages; each period pools its numerator and denominator rather than averaging yearly percentages. Country/survey composition changes over time.","",
  "| Entry period | Records | Hib3 missing | PCV missing | Rotavirus missing |",
  "|---|---:|---:|---:|---:|",
  sprintf("| %s | %s | %.1f%% | %.1f%% | %.1f%% |",wide$period,
    format(wide$records,big.mark=",",trim=TRUE),wide$hib3_pct,wide$pcv3_pct,wide$rotavirus_pct),"",
  "Hib3 has no missing assigned values from 2012 onward in this sample. PCV is missing for 5.8% of 2015 records and 0.3% of 2021 records. Rotavirus is missing for 31.0% in 2015 and 16.6% in 2021. All three have 0% missingness among the represented 2022–2023 records. There are no 2001 or 2024 band-entry records in this previously fitted sample, so missingness for those years is undefined, not zero.","",
  "The national source contains no PCV/rotavirus series for Comoros, Gabon and Guinea. Their lack of later records in this sample must not be mistaken for resolution of their source gaps. Other missing values are pre-series placeholders labelled assumed not introduced; the pipeline currently rejects them as unverified zeros. Confirming genuine pre-introduction zero coverage could substantially reduce historical missingness.","",
  "![Vaccine missingness over time](vaccine_missingness_over_time.png)","",
  "[Annual counts and source-gap breakdown](vaccine_missingness_by_entry_year.csv) · [Period summary](vaccine_missingness_by_period.csv) · [Country/year counts](vaccine_missingness_by_country_year.csv). The annual CSV includes the number of represented countries and surveys. Missing records equal pre-series plus no-series records, and totals were checked against the preceding overall audit.","",
  "Reproduce with `Rscript R_cbh/covariates/06_vaccine_missingness_time.R`. No datasets, fits, paper figures or TeX files are changed.")
writeLines(lines,file.path(out,"VACCINE_MISSINGNESS_OVER_TIME.md"))
cbh_atomic_csv(data.frame(file=c(path,"R_cbh/covariates/06_vaccine_missingness_time.R"),
  md5=vapply(c(path,"R_cbh/covariates/06_vaccine_missingness_time.R"),cbh_file_hash,"")),
  file.path(out,"vaccine_time_input_manifest.csv"))
print(wide)

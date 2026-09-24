#!/usr/bin/env Rscript
# Figure 4 (added 24 September 2026): probability of dying before age 5 by country, all causes
# (IHME age-band death rates) and the part caused by malaria under the PfPR-ACM model.
# For each burden country and year the seven age-band hazards come from IHME all-cause
# death rates (the neonatal band from the early and late neonatal rates, as in the DRC life
# table of R_cbh/primary/02_effects.R); the zero-PfPR hazards multiply each band by the
# model's hazard ratio for PfPR 0 versus the country's current PfPR. The synthetic-cohort
# probability is 5q0 = 1 - exp(-sum of band cumulative hazards). Intervals draw each band's
# log hazard ratio independently (the bands are separate models). Aggregates only.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
library(data.table); library(ggplot2)
st <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional_mics")); root <- st$out
out <- file.path(root,"under5_probability"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
ages <- cbh_config()$age_bands$age_band
edges <- c(0,28/365.25*12,6,12,24,36,48,60)
paths <- c(estimates=file.path(root,"burden/country_age_estimates.csv"),igme="data/igme_u5mr_by_country.csv",
  ihme_inputs=file.path("results/cbh/age_band_separate_v1",sprintf("country_burden_%d",st$years),"ihme_disjoint_age_inputs.csv"))
stopifnot(all(file.exists(paths)))
est <- fread(paths[["estimates"]])[status=="estimated"]
stopifnot(setequal(unique(est$year),st$years),all(est[,.N,by=.(iso3,year)]$N==7L))
set.seed(st$seed); n_draws <- 4000L
res <- rbindlist(lapply(st$years,function(y) {
  src <- fread(file.path("results/cbh/age_band_separate_v1",sprintf("country_burden_%d",y),"ihme_disjoint_age_inputs.csv"))
  rbindlist(lapply(sort(unique(est[year==y]$iso3)),function(iso) {
    z <- est[year==y & iso3==iso][match(ages,age_band)]
    H <- z$ihme_rate_per100000/1e5*diff(edges)/12
    s <- src[iso3==iso]
    early <- s[source_age=="0 to 6 days (early neonatal)"]$rate_per100000; late <- s[source_age=="7 to 27 days (late neonatal)"]$rate_per100000
    stopifnot(length(early)==1,length(late)==1,all(is.finite(H)),all(is.finite(z$log_hr_zero_vs_current)),all(z$log_hr_se>0))
    H[1] <- early/1e5*7/365.25+late/1e5*21/365.25
    H0 <- H*exp(z$log_hr_zero_vs_current)
    draws <- vapply(seq_len(7),function(g) H[g]*exp(rnorm(n_draws,z$log_hr_zero_vs_current[g],z$log_hr_se[g])),numeric(n_draws))
    q0_draws <- -expm1(-rowSums(draws)); q <- -expm1(-sum(H)); q0 <- -expm1(-sum(H0))
    data.table(year=y,iso3=iso,country=z$Location[1],pfpr_pct=weighted.mean(z$pfpr_pct,z$implied_person_years),
      q5_allcause=q,q5_no_malaria=q0,q5_no_malaria_lower=unname(quantile(q0_draws,.025)),q5_no_malaria_upper=unname(quantile(q0_draws,.975)),
      q5_attributable=q-q0,q5_attributable_lower=q-unname(quantile(q0_draws,.975)),q5_attributable_upper=q-unname(quantile(q0_draws,.025)),
      share_attributable=(q-q0)/q,any_zero_below_support=any(z$zero_below_observed_support))
  }))
}))
igme <- fread(paths[["igme"]])[year==max(st$years),.(iso3,igme_u5mr=u5mr)]
res <- merge(res,igme,by="iso3",all.x=TRUE); setorder(res,year,-q5_allcause)
cbh_atomic_csv(as.data.frame(res),file.path(out,"under5_death_probability.csv"))
cur <- res[year==max(st$years)]
check <- cur[is.finite(igme_u5mr),.(n=.N,median_ratio=median(1000*q5_allcause/igme_u5mr),cor=cor(1000*q5_allcause,igme_u5mr))]

fmt1 <- function(x) sprintf("%.0f",1000*x)
## ---- Figure 4 ----------------------------------------------------------------------------------------------
# Probability of dying before age 5 per 1,000 live births (x axis from 0): all causes from the IHME
# age-band death rates (light bar) and the part caused by malaria under the PfPR-ACM model (dark bar,
# 95% interval), with both values listed on the right.
display <- c(Swaziland="Eswatini",`Cote d'Ivoire`="Côte d'Ivoire",Congo="Republic of the Congo")
cur[,country:=fifelse(country %in% names(display),display[country],country)]
setorder(cur,-q5_allcause)
cur[,label:=factor(country,levels=rev(country))]
long <- rbind(cur[,.(label,value=1000*q5_allcause,series="All causes (IHME)")],
  cur[,.(label,value=1000*q5_attributable,series="Caused by malaria (PfPR-ACM model)")])
long[,series:=factor(series,levels=c("All causes (IHME)","Caused by malaria (PfPR-ACM model)"))]
fills <- c(`All causes (IHME)`="#C9D6E3",`Caused by malaria (PfPR-ACM model)`="#1F4E79")
xmax <- ceiling(max(1000*cur$q5_allcause)/20)*20
p <- ggplot()+
  geom_col(data=long[series=="All causes (IHME)"],aes(value,label,fill=series),width=.75)+
  geom_col(data=long[series!="All causes (IHME)"],aes(value,label,fill=series),width=.45)+
  geom_errorbarh(data=cur,aes(xmin=1000*q5_attributable_lower,xmax=1000*q5_attributable_upper,y=label),height=.3,colour="grey20",linewidth=.4)+
  geom_text(data=cur,aes(x=xmax+3,y=label,label=sprintf("%.0f",1000*q5_allcause)),hjust=0,size=3.5,colour="grey25")+
  geom_text(data=cur,aes(x=xmax+15,y=label,label=sprintf("%.0f",1000*q5_attributable)),hjust=0,size=3.5,colour="#1F4E79")+
  annotate("text",x=c(xmax+3,xmax+15),y=nrow(cur)+1.1,label=c("All\ncauses","Malaria"),hjust=0,vjust=0,size=3.4,lineheight=.9,colour=c("grey25","#1F4E79"))+
  scale_fill_manual(values=fills,name=NULL)+
  scale_x_continuous(limits=c(0,xmax+26),breaks=seq(0,xmax,20),expand=c(0,0))+
  coord_cartesian(clip="off")+
  labs(x=sprintf("Probability of dying before age 5, %d (per 1,000 live births)",max(st$years)),y=NULL)+
  theme_minimal(base_size=15)+theme(panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),legend.position="bottom",
    axis.text.y=element_text(size=11.5),legend.text=element_text(size=13),plot.margin=margin(30,12,8,8))
figure <- file.path(out,"under5_death_probability_2024.png")
ggsave(figure,p,width=10,height=12,dpi=260,device=ragg::agg_png,bg="white")
top <- cur[order(-q5_attributable)][1:3]
caption <- paste0("Probability of dying before age 5 in ",max(st$years)," per 1,000 live births, by country (",nrow(cur)," sub-Saharan African countries in the burden comparison): all causes, from IHME age-specific all-cause death rates (light bars), and the part caused by malaria under the ", cbh_paper_model_label(), " (dark bars, with 95% intervals); both values are listed on the right. ",
  "Probabilities are period life-table (synthetic-cohort) calculations, 5q0 = 1 − exp(−Σ_g H_g), over the model's seven age bands (<1, 1–5, 6–11, 12–23, 24–35, 36–47 and 48–59 completed months), with band cumulative hazards H_g from the IHME death rates. ",
  "The malaria-caused probability is 5q0 minus the probability with malaria transmission removed, for which each band's hazard is multiplied by the model's hazard ratio for PfPR[2–10] = 0 versus the country's current population-weighted MAP prevalence. ",
  sprintf("Intervals draw each band's log hazard ratio independently; they condition on the IHME rates and the fitted smoothing parameters. Zero prevalence lies below the observed exposure support, so the counterfactual is an extrapolation. Against UN IGME under-five mortality for %d, the IHME-based all-cause probabilities have median ratio %.2f and correlation %.2f; IHME and UN IGME differ markedly for some countries. ",max(st$years),check$median_ratio,check$cor),
  sprintf("The largest malaria-caused probabilities are in %s (%s per 1,000), %s (%s) and %s (%s); the median is %s per 1,000, %.0f%% of the all-cause probability.",
    top$country[1],fmt1(top$q5_attributable[1]),top$country[2],fmt1(top$q5_attributable[2]),top$country[3],fmt1(top$q5_attributable[3]),fmt1(median(cur$q5_attributable)),100*median(cur$share_attributable)))
writeLines(c("# Figure 4 caption","",caption),file.path(out,"CAPTION.md"))
inputs <- c(unname(paths),"R_cbh/reporting/13_under5_death_probability.R","R_cbh/primary/settings.R")
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"provenance.csv"))
print(cur[1:10,.(country,pfpr=round(pfpr_pct,1),q5=round(1000*q5_allcause,1),malaria=round(1000*q5_attributable,1),lower=round(1000*q5_attributable_lower,1),upper=round(1000*q5_attributable_upper,1))])
print(check); message("Figure 4 (under-5 death probability) written: ",figure)

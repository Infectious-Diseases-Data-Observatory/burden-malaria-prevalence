#!/usr/bin/env Rscript
# Plotting only: read aggregates, never fit models or prepare exposures.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
library(ggplot2)
settings <- cbh_primary_settings();out <- settings$out
writeLines(trimws(readLines(file.path(out,"model_formula.txt")),which="right"),file.path(out,"model_formula.txt"))
ages <- cbh_config()$age_bands$age_band
age_factor <- function(x) factor(x,levels=ages,labels=paste(ages,"months"))
read <- function(x) cbh_read_csv(file.path(out,x))
save_plot <- function(p,name,width,height) ggsave(file.path(out,name),p,width=width,height=height,
  dpi=240,device=ragg::agg_png,bg="white")
d <- read("pfpr_curves.csv");hr <- read("pfpr_40_to_20_contrasts.csv")
diag <- read("fit_diagnostics.csv");sample <- read("primary_sample.csv")
edf <- read("pfpr_edf.csv");sums <- read("burden/year_summary.csv")
totals <- read("burden/country_totals.csv");bands <- read("burden/country_age_estimates.csv")
comparison <- read("burden/comparison_with_previous_primary.csv")
cc <- read("comparison_with_previous_curves.csv")
stopifnot(nrow(diag)==7,all(diag$series=="map_full"),all(diag$converged),all(diag$input_verified),
  all(diag$gamma==2),all(diag$min_smoothing_hessian_eigenvalue>0))
base <- theme_minimal(base_size=14)+theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face="bold"),strip.text=element_text(face="bold"),
  plot.caption=element_text(hjust=0,size=11),legend.position="bottom")
d$age_label <- age_factor(d$age_band)
p <- ggplot(d[d$within_central_support,],aes(pfpr_pct,log_hazard_ratio))+
  geom_hline(yintercept=0,colour="grey65",linewidth=.4)+
  geom_vline(xintercept=20,colour="grey80",linewidth=.4)+
  geom_ribbon(aes(ymin=lower_95,ymax=upper_95),fill="#215E91",alpha=.17)+
  geom_line(colour="#215E91",linewidth=1)+facet_wrap(~age_label,ncol=4)+
  scale_x_continuous(limits=c(0,100),breaks=c(0,20,40,60,80,100))+
  labs(title="PfPR and mortality by age",subtitle="Primary MAP analysis · separate age-band models · gamma = 2",
    x="PfPR[2–10] (%)",y="Log hazard ratio relative to PfPR = 20%",
    caption="Curves show the central 95% of each band's exposure distribution. Shading: conditional pointwise 95% intervals.")+base
save_plot(p,"pfpr_splines.png",14,8)
hr$age_label <- factor(hr$age_band,levels=rev(ages))
p <- ggplot(hr,aes(hazard_ratio_40_to_20,age_label))+
  geom_vline(xintercept=1,colour="grey60",linetype=2)+
  geom_segment(aes(x=lower_95,xend=upper_95,yend=age_label),colour="#215E91",linewidth=.8)+
  geom_point(colour="#215E91",size=3)+scale_x_log10()+
  labs(title="Mortality change if PfPR falls from 40% to 20%",
    subtitle="Primary MAP analysis · gamma = 2",x="Mortality hazard ratio (95% interval)",y="Age (completed months)",
    caption="Intervals condition on fitted smoothing parameters, fixed HIV imputation and exposure.")+base
save_plot(p,"pfpr_40_to_20.png",9,6)
# Country totals are signed point estimates, with missing countries retained in CSVs.
z <- totals[is.finite(totals$attributable_under5_deaths),]
long <- rbind(data.frame(z[c("iso3","country","year")],deaths=z$attributable_under5_deaths,source="Primary MAP model"),
  data.frame(z[c("iso3","country","year")],deaths=z$ihme_malaria_deaths,source="IHME malaria"))
long$country <- factor(long$country,levels=rev(sort(unique(long$country))))
p <- ggplot(long,aes(deaths,country,colour=source,shape=source))+
  geom_point(position=position_dodge(width=.55),size=2.2)+facet_wrap(~year,nrow=1)+
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=100),breaks=c(0,1000,10000,100000),
    labels=scales::label_number(scale_cut=scales::cut_short_scale()))+
  scale_colour_manual(values=c("Primary MAP model"="#215E91","IHME malaria"="#A15132"))+
  labs(title="Malaria-attributable deaths by country",subtitle="Primary MAP gamma = 2 versus IHME cause-specific malaria deaths",
    x="Deaths before age 5 (pseudo-log scale)",y=NULL,colour=NULL,shape=NULL,
    caption="42 countries with national MAP exposure. Model estimates are all-cause reductions under zero PfPR; IHME is a cause-specific comparator.")+base+
  theme(axis.text.y=element_text(size=11))
save_plot(p,"burden/country_deaths_all_years.png",15,14)
p <- ggplot(z,aes(ihme_malaria_deaths,attributable_under5_deaths))+
  geom_abline(slope=1,intercept=0,colour="grey65",linetype=2)+geom_point(colour="#215E91",size=2.5,alpha=.8)+
  geom_text(data=z[z$iso3 %in% c("COD","NGA","AGO","UGA","TZA"),],aes(label=iso3),
    size=3.5,hjust=-.1,vjust=-.5,check_overlap=TRUE)+facet_wrap(~year,nrow=1)+
  coord_equal(xlim=c(min(0,z$attributable_under5_deaths),1.13*max(z$attributable_under5_deaths,z$ihme_malaria_deaths)),
    ylim=c(min(0,z$attributable_under5_deaths),1.13*max(z$attributable_under5_deaths,z$ihme_malaria_deaths)),expand=FALSE)+
  scale_x_continuous(labels=scales::label_number(scale_cut=scales::cut_short_scale()),expand=expansion(mult=c(.04,.13)))+
  scale_y_continuous(labels=scales::label_number(scale_cut=scales::cut_short_scale()),expand=expansion(mult=c(.04,.13)))+
  labs(title="Country mortality estimates versus IHME",subtitle="Primary MAP analysis · gamma = 2",
    x="IHME malaria deaths before age 5",y="Model-attributable deaths before age 5",
    caption="Each point is a country. Dashed line: equality. Country totals are point estimates.")+base
save_plot(p,"burden/country_vs_ihme.png",15,6)
b <- bands[is.finite(bands$attributable_deaths),]
b <- aggregate(b$attributable_deaths,b[c("year","age_band")],sum);names(b)[3] <- "deaths"
b$age_label <- age_factor(b$age_band)
p <- ggplot(b,aes(age_label,deaths))+geom_col(fill="#215E91",width=.7)+facet_wrap(~year,nrow=1)+
  scale_y_continuous(labels=scales::label_number(scale_cut=scales::cut_short_scale()))+
  labs(title="Malaria-attributable deaths by age",subtitle="Totals across the same 42 countries · primary MAP gamma = 2",
    x="Age (completed months)",y="Attributable deaths",
    caption="Signed point estimates; no cross-age independence assumption or aggregate confidence interval.")+base+
  theme(axis.text.x=element_text(angle=45,hjust=1))
save_plot(p,"burden/deaths_by_age.png",15,6)
life <- read("burden/drc_life_table.csv")
long <- rbind(data.frame(life[c("year","end_month")],survival=life$observed_survival,scenario="IHME all-cause"),
  data.frame(life[c("year","end_month")],survival=life$zero_pfpr_survival,scenario="Zero PfPR"))
start <- unique(long[c("year","scenario")]);start$end_month <- 0;start$survival <- 1
long <- rbind(long,start[names(long)])
p <- ggplot(long,aes(end_month,survival,colour=scenario))+geom_line(linewidth=1)+facet_wrap(~year,nrow=1)+
  scale_colour_manual(values=c("IHME all-cause"="#59636B","Zero PfPR"="#215E91"))+
  scale_y_continuous(labels=scales::label_percent(accuracy=1))+
  labs(title="DRC: survival to age 5 under zero PfPR",subtitle="Synthetic cohort calibrated to national IHME mortality rates",
    x="Age (months)",y="Survival probability",colour=NULL,
    caption="Survival is the product of band survival probabilities. Lines join band endpoints; this is separate from annual attributable death counts.")+base
save_plot(p,"burden/drc_survival.png",15,6)
g <- read("grouped_outcome_checks.csv")
g <- g[g$grouping=="pfpr_pct",]
g$pfpr_midpoint <- as.numeric(g$group)+5;g$age_label <- age_factor(g$age_band)
g$residual_per1000 <- 1000*(g$observed_probability-g$fitted_probability)
p <- ggplot(g,aes(pfpr_midpoint,residual_per1000))+geom_hline(yintercept=0,colour="grey65")+
  geom_line(colour="#215E91")+geom_point(aes(size=records),colour="#215E91",alpha=.8)+
  facet_wrap(~age_label,ncol=4)+scale_size_area(max_size=5,guide="none")+
  labs(title="Observed minus fitted mortality by PfPR",subtitle="In-sample check of band death probabilities · primary MAP gamma = 2",
    x="PfPR[2–10] bin midpoint (%)",y="Observed minus fitted deaths per 1,000 band entries",
    caption="Ten-percentage-point bins; cells with fewer than 100 records are omitted. Point area reflects records. This is not held-out validation.")+base
save_plot(p,"outcome_residuals_by_pfpr.png",14,8)
fmt <- function(x) format(round(x),big.mark=",",trim=TRUE)
md <- c("# Primary MAP analysis rerun","",
  sprintf("Seven age-band models were freshly fitted at gamma=2 using %s children, %s child-band records and %s deaths (%s surveys; %s countries). Dataset preparation, HIV imputation and national exposure/IHME input preparation were not rerun.",
    fmt(sample$distinct_children),fmt(sample$records),fmt(sample$deaths),fmt(sample$surveys),fmt(sample$countries)),"",
  "The unweighted binomial complementary-log-log model uses a fixed log band-width offset; separate cr PfPR (k=5) and calendar-year (k=6) splines, the saved scaled confounders, and survey/country/region random intercepts in each age band. Reference knot locations, sample selection and posterior-median child HIV incidence are unchanged. See [formula](model_formula.txt), [knots](../../../R_cbh/primary/reference_knots.csv) and [analysis plan](../../../docs/ANALYSIS_PLAN.md).","",
  "![PfPR curves](pfpr_splines.png)","",
  "| Completed months | Records | Deaths | PfPR EDF | HR: 40% to 20% (95% interval) | Fit seconds |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(ages,function(a) {
    x <- diag[diag$age_band==a,];h <- hr[hr$age_band==a,];e <- edf$edf[edf$age_band==a]
    sprintf("| %s | %s | %s | %.2f | %.3f (%.3f–%.3f) | %.1f |",a,fmt(x$rows),fmt(x$deaths),e,
      h$hazard_ratio_40_to_20,h$lower_95,h$upper_95,x$elapsed_seconds)
  },""),"",
  "[40% to 20% figure](pfpr_40_to_20.png) · [All curve estimates](pfpr_curves.csv) · [Attributable-fraction anchors](pfpr_attributable_fraction_anchors.csv)","",
  "## National mortality","",
  "For each age band, HR = exp[f_g(0) − f_g(P_country)]. The IHME all-cause rate/count is multiplied by HR for the counterfactual and by (1 − HR) for the attributable contribution. National MAP prevalence and IHME inputs are held at the existing values for each year. Ages 2–4 share the IHME rate and split deaths/person-time equally; early/late neonatal inputs share the <1-month model effect.","",
  "| Year | Countries | Attributable deaths | IHME malaria deaths |","|---|---:|---:|---:|",
  vapply(seq_len(nrow(sums)),function(i)sprintf("| %d | %d | %s | %s |",sums$year[i],sums$countries[i],
    fmt(sums$attributable_under5_deaths[i]),fmt(sums$ihme_malaria_deaths[i])),""),"",
  "![Country estimates](burden/country_vs_ihme.png)","",
  "[Every country's figure](burden/country_deaths_all_years.png) · [Country totals](burden/country_totals.csv) · [Age-band rates and deaths](burden/country_age_estimates.csv) · [Age contributions](burden/deaths_by_age.png) · [DRC survival](burden/drc_survival.png)","",
  "All 45 countries remain in the tables; Cape Verde, Lesotho and São Tomé and Príncipe lack national MAP estimates. Totals use the same 42 estimable countries. Missing estimates are not zero; negative attributable estimates are retained. National prevalence uses the saved GPW 2020 population weights in each scenario year. Evaluating a nonlinear spline at national mean prevalence differs from averaging local attributable effects.","",
  "## Verification and limits","",
  sprintf("All seven fits converged, had full coefficient rank, finite covariance and positive smoothing-Hessian eigenvalues. Every outcome, predictor and offset was checked against the prepared sample. Compact spline estimates and variances were checked against full model predictions. Total selected fitting time: %.1f minutes. This mgcv build does not support OpenMP, so its recorded nthreads warning means computation used one thread.",sum(diag$elapsed_seconds)/60),"",
  sprintf("Against the previous MAP gamma=2 run, the maximum absolute change in a curve log-hazard ratio was %.3g, and the largest country/year change was %.3g deaths. [Curve comparison](comparison_with_previous_curves.csv) · [Country comparison](burden/comparison_with_previous_primary.csv).",max(abs(cc$log_hr_difference)),max(abs(comparison$difference_deaths),na.rm=TRUE)),"",
  "Band intervals condition on smoothing parameters, one HIV imputation, exposure and IHME baselines. Country totals have no confidence intervals because cross-age sampling covariance is not estimated. Design, within-child dependence, HIV/exposure uncertainty, transport assumptions and IHME allocation uncertainty remain unresolved. The inherited regional MAP extraction issues documented in the [code audit](../../../docs/CODE_AUDIT.md) were not altered by this rerun. Country-held-out validation and gamma=2 structure/geography/period sensitivity fits remain separate planned analyses.","",
  "Model-attributable all-cause reductions and IHME cause-specific malaria counts are different estimands; agreement is descriptive. Zero/current exposure support flags accompany all age-band national estimates. The DRC survival curves describe a synthetic cohort, not additional annual death counts.","",
  "[Residual figure](outcome_residuals_by_pfpr.png) · [Fitted outcome checks](fitted_outcome_checks.csv) · [Grouped checks](grouped_outcome_checks.csv) · [Fit diagnostics](fit_diagnostics.csv) · [Fit manifest](fit_manifest.csv) · [Fit inputs/provenance](fit_input_provenance.csv) · [Burden provenance](burden/provenance.csv) · [Software versions](session_info.txt)","",
  "Reproduce the primary analysis, without data setup: `Rscript run_all.R --primary`. This forces fresh fitting. To regenerate aggregates/figures from verified current fits, use `Rscript R_cbh/primary/run.R --report-only`."
)
writeLines(md,file.path(out,"REPORT.md"))
files <- c("pfpr_curves.csv","pfpr_40_to_20_contrasts.csv","fit_diagnostics.csv","primary_sample.csv","pfpr_edf.csv",
  "burden/year_summary.csv","burden/country_totals.csv","burden/country_age_estimates.csv","burden/drc_life_table.csv",
  "burden/comparison_with_previous_primary.csv","comparison_with_previous_curves.csv","grouped_outcome_checks.csv","fitted_outcome_checks.csv")
paths <- c(file.path(out,files),"R_cbh/primary/03_report.R","R_cbh/primary/settings.R")
cbh_atomic_csv(data.frame(file=paths,md5=vapply(paths,cbh_file_hash,"")),file.path(out,"report_provenance.csv"))
message("Primary MAP figures and report complete")

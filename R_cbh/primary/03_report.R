#!/usr/bin/env Rscript
# Plotting only: read aggregates, never fit models or prepare exposures.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/age_band_table.R")
source("R_cbh/reporting/labels.R")
source("R_cbh/reporting/log_axes.R")
library(ggplot2)
settings <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional_mics"));out <- settings$out
model_label <- cbh_paper_model_label()
writeLines(trimws(readLines(file.path(out,"model_formula.txt")),which="right"),file.path(out,"model_formula.txt"))
ages <- cbh_config()$age_bands$age_band
age_factor <- function(x) factor(x,levels=ages,labels=ifelse(ages=="<1","<1 month",paste(ages,"months")))
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
# Larger typography for the three manuscript figures, at unchanged plot widths.
paper <- base+theme(text=element_text(size=20),axis.title=element_text(size=20),
  axis.text=element_text(size=16),strip.text=element_text(size=18,face="bold"),
  legend.title=element_text(size=18),legend.text=element_text(size=16),
  panel.spacing=grid::unit(24,"pt"),plot.margin=margin(10,24,10,10))
d$age_label <- age_factor(d$age_band)
# PfPR axis 0-80%, as in the other supplementary spline figures; every curve's central support lies inside it.
stopifnot(max(d$pfpr_pct[d$within_central_support])<=80)
p <- ggplot(d[d$within_central_support,],aes(pfpr_pct,log_hazard_ratio))+
  geom_hline(yintercept=0,colour="grey65",linewidth=.4)+
  geom_vline(xintercept=20,colour="grey80",linewidth=.4)+
  geom_ribbon(aes(ymin=lower_95,ymax=upper_95),fill="#215E91",alpha=.17)+
  geom_line(colour="#215E91",linewidth=1)+facet_wrap(~age_label,ncol=4)+
  scale_x_continuous(limits=c(0,80),breaks=c(0,20,40,60,80))+
  labs(x="PfPR[2–10] (%)",y="Log hazard ratio\nrelative to PfPR = 20%")+paper
save_plot(p,"pfpr_splines.png",14,8)
hr$age_label <- factor(hr$age_band,levels=rev(ages))
p <- ggplot(hr,aes(hazard_ratio_40_to_20,age_label))+
  geom_vline(xintercept=1,colour="grey60",linetype=2)+
  geom_segment(aes(x=lower_95,xend=upper_95,yend=age_label),colour="#215E91",linewidth=.8)+
  geom_point(colour="#215E91",size=3)+scale_x_log10(minor_breaks=cbh_log10_minor_breaks)+
  labs(title="Mortality change if PfPR falls from 40% to 20%",
    subtitle="Primary MAP analysis · gamma = 2",x="Mortality hazard ratio (95% interval)",y="Age (completed months)",
    caption="Intervals condition on fitted smoothing parameters, fixed HIV imputation and exposure.")+base+cbh_log10_grid_theme()
save_plot(p,"pfpr_40_to_20.png",9,6)
# Country totals are signed point estimates, with missing countries retained in CSVs.
z <- totals[is.finite(totals$attributable_under5_deaths),]
long <- rbind(data.frame(z[c("iso3","country","year")],deaths=z$attributable_under5_deaths,source=model_label),
  data.frame(z[c("iso3","country","year")],deaths=z$ihme_malaria_deaths,source="IHME malaria"))
long$country <- factor(long$country,levels=rev(sort(unique(long$country))))
p <- ggplot(long,aes(deaths,country,colour=source,shape=source))+
  geom_point(position=position_dodge(width=.55),size=2.2)+facet_wrap(~year,nrow=1)+
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=100),breaks=c(0,1000,10000,100000),
    labels=scales::label_number(scale_cut=scales::cut_short_scale()))+
  scale_colour_manual(values=setNames(c("#215E91","#A15132"),c(model_label,"IHME malaria")))+
  labs(title="Malaria-attributable deaths by country",subtitle="Primary MAP gamma = 2 versus IHME cause-specific malaria deaths",
    x="Deaths before age 5 (pseudo-log scale)",y=NULL,colour=NULL,shape=NULL,
    caption="42 countries with national MAP exposure. Model estimates are all-cause reductions under zero PfPR; IHME is a cause-specific comparator.")+base+
  theme(axis.text.y=element_text(size=11))
save_plot(p,"burden/country_deaths_all_years.png",15,14)
# True log10 scales require positive values: keep missing/nonpositive estimates
# in the source tables and explicitly record any rows unavailable for this plot.
plot_rows <- totals[is.finite(totals$ihme_malaria_deaths) & totals$ihme_malaria_deaths>0 &
  is.finite(totals$attributable_under5_deaths) & totals$attributable_under5_deaths>0,]
excluded <- totals[!paste(totals$iso3,totals$year) %in% paste(plot_rows$iso3,plot_rows$year),]
cbh_atomic_csv(excluded,file.path(out,"burden/log10_plot_exclusions.csv"))
stopifnot(nrow(plot_rows)>0)
log_limits <- c(1e3,10^ceiling(log10(max(plot_rows$ihme_malaria_deaths,plot_rows$attributable_under5_deaths))))
below_display <- plot_rows[plot_rows$ihme_malaria_deaths<log_limits[1] |
  plot_rows$attributable_under5_deaths<log_limits[1],]
cbh_atomic_csv(below_display,file.path(out,"burden/log10_below_display_minimum.csv"))
log_breaks <- 10^seq(log10(log_limits[1]),log10(log_limits[2]))
log_labels <- function(x) parse(text=paste0("10^",round(log10(x))))
p <- ggplot(plot_rows,aes(ihme_malaria_deaths,attributable_under5_deaths))+
  geom_abline(slope=1,intercept=0,colour="grey65",linetype=2)+
  geom_point(colour="#215E91",size=2.5,alpha=.8)+
  ggrepel::geom_text_repel(data=plot_rows[plot_rows$iso3 %in% c("COD","NGA","AGO","UGA","TZA"),],
    aes(label=iso3),size=5,seed=20260915,max.overlaps=Inf,min.segment.length=0)+
  facet_wrap(~year,nrow=1)+coord_equal(xlim=log_limits,ylim=log_limits,expand=FALSE)+
  scale_x_log10(breaks=log_breaks,labels=log_labels,minor_breaks=cbh_log10_minor_breaks)+
  scale_y_log10(breaks=log_breaks,labels=log_labels,minor_breaks=cbh_log10_minor_breaks)+
  labs(x="IHME malaria deaths before age 5",y=paste0(model_label," deaths\nbefore age 5"))+paper+cbh_log10_grid_theme()+
  theme(panel.spacing.x=grid::unit(40,"pt"),plot.margin=margin(10,26,10,10))
# Both axes have exactly the same base-10 transform and range; no pseudocount.
stopifnot(identical(p$scales$get_scales("x")$trans$name,"log-10"),
  identical(p$scales$get_scales("y")$trans$name,"log-10"))
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
age_table <- cbh_primary_age_table(out,ages)
md <- c("# Primary MAP analysis rerun","",
  sprintf("Seven age-band models were freshly fitted at gamma=2 using %s children, %s child-band records and %s deaths (%s surveys; %s countries). Dataset preparation, HIV imputation and national exposure/IHME input preparation were not rerun.",
    fmt(sample$distinct_children),fmt(sample$records),fmt(sample$deaths),fmt(sample$surveys),fmt(sample$countries)),"",
  "The unweighted binomial complementary-log-log model uses a fixed log band-width offset; separate cr PfPR (k=5) and calendar-year (k=6) splines, the saved scaled confounders, and survey/country/region random intercepts in each age band. Reference knot locations, sample selection and posterior-median child HIV incidence are unchanged. See [formula](model_formula.txt), [knots](../../../R_cbh/primary/reference_knots.csv) and [analysis plan](../../../docs/ANALYSIS_PLAN.md).","",
  "![PfPR curves](pfpr_splines.png)","",
  age_table,"",
  "[LaTeX table fragment](tables/age_band_results.latex.txt) · [Table data](tables/age_band_results.csv)","",
  "[40% to 20% figure](pfpr_40_to_20.png) · [All curve estimates](pfpr_curves.csv) · [Attributable-fraction anchors](pfpr_attributable_fraction_anchors.csv)","",
  "## National mortality","",
  "For each age band, HR = exp[f_g(0) − f_g(P_country)]. The IHME all-cause rate/count is multiplied by HR for the counterfactual and by (1 − HR) for the attributable contribution. National MAP prevalence and IHME inputs are held at the existing values for each year. Ages 2–4 share the IHME rate and split deaths/person-time equally; early/late neonatal inputs share the <1-month model effect.","",
  "| Year | Countries | Attributable deaths | IHME malaria deaths |","|---|---:|---:|---:|",
  vapply(seq_len(nrow(sums)),function(i)sprintf("| %d | %d | %s | %s |",sums$year[i],sums$countries[i],
    fmt(sums$attributable_under5_deaths[i]),fmt(sums$ihme_malaria_deaths[i])),""),"",
  "![Country estimates](burden/country_vs_ihme.png)","",
  "Both axes use a true log10 scale with identical limits starting at 1,000 deaths. Countries below 1,000 on either axis fall outside the displayed range; [these rows](burden/log10_below_display_minimum.csv) and all underlying estimates are retained. There are 42 estimable countries per year; [excluded country-year rows](burden/log10_plot_exclusions.csv) have unavailable estimates. No pseudocount is added. Titles and subtitles are omitted from the paper's PfPR and country-comparison figures; years and age-band facet labels remain.","",
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
if(settings$regional) md <- c("# Current primary figures and tables", "",
  sprintf("The %s 17-variable PfPR-ACM model uses %s children, %s child-band records and %s deaths in %s surveys and %s countries.",
    if(isTRUE(settings$imputed)) "imputed-covariate (sensitivity)" else "revised",
    fmt(sample$distinct_children),fmt(sample$records),fmt(sample$deaths),fmt(sample$surveys),fmt(sample$countries)), "",
  "Seven separate MAP gamma=2 age-band fits, including survey-region urban percentage, with fixed PfPR/time knots and the declared HIV, UNICEF and available-region substitutions. The raw-source data and HIV imputation were not refitted. All seven models passed the fitted-input and numerical checks.", "",
  "![PfPR curves](pfpr_splines.png)", "", age_table, "",
  "[Table CSV](tables/age_band_results.csv) · [LaTeX source as plain text](tables/age_band_results.latex.txt)", "",
  "[40% to 20% contrast figure](pfpr_40_to_20.png) · [Full contrasts and comparison with the previous adjustment](REPORT.md)", "",
  "![Country estimates](burden/country_vs_ihme.png)", "",
  "Country/IHME axes use identical base-10 scales starting at 1,000 deaths. All country estimates, signed age contributions and missing/support flags remain in the CSVs. Age intervals are conditional on smoothing parameters, exposure and filled covariates; marginal intervals are not summed into total intervals.", "",
  "[Country totals](burden/country_totals.csv) · [Country-age estimates](burden/country_age_estimates.csv) · [Annual comparison](annual_comparison/README.md) · [Nigeria states](nigeria_states/README.md)", "",
  "[Age contributions](burden/deaths_by_age.png) · [DRC synthetic-cohort survival](burden/drc_survival.png) · [Aggregate outcome check](outcome_residuals_by_pfpr.png)", "",
  "[Paper figure manifest/captions](paper_figures/CAPTIONS.md) · [Study flow](study_flow/CAPTION.md) · [Fit diagnostics](fit_diagnostics.csv)", "",
  "These are primary results. Sahel, geographic, period and other sensitivity analyses are separate; they are not refitted by reporting. No TeX files are written. Reproduce with `Rscript R_cbh/primary/run_regional.R --report-only`.")
writeLines(md,file.path(out,if(settings$regional) "RESULTS.md" else "REPORT.md"))
files <- c("pfpr_curves.csv","pfpr_40_to_20_contrasts.csv","pfpr_20_to_zero_contrasts.csv","tables/age_band_results.csv","fit_diagnostics.csv","primary_sample.csv","pfpr_edf.csv",
  "burden/year_summary.csv","burden/country_totals.csv","burden/country_age_estimates.csv","burden/drc_life_table.csv",
  "burden/comparison_with_previous_primary.csv","comparison_with_previous_curves.csv","grouped_outcome_checks.csv","fitted_outcome_checks.csv")
paths <- c(file.path(out,files),"R_cbh/primary/03_report.R","R_cbh/primary/settings.R","R_cbh/reporting/age_band_table.R","R_cbh/reporting/labels.R","R_cbh/reporting/log_axes.R")
cbh_atomic_csv(data.frame(file=paths,md5=vapply(paths,cbh_file_hash,"")),file.path(out,"report_provenance.csv"))
message("Primary MAP figures and report complete")

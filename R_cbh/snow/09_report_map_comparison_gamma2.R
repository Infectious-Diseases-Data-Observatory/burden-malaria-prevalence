#!/usr/bin/env Rscript
source("R_cbh/load_pipeline.R")
library(ggplot2)
out <- "results/cbh/map_snow_gamma2_v1"
writeLines(trimws(readLines(file.path(out,"model_formula.txt")),which="right"),
  file.path(out,"model_formula.txt"))
d <- cbh_read_csv(file.path(out,"pfpr_curves.csv"))
diag <- cbh_read_csv(file.path(out,"fit_diagnostics.csv"))
sm <- cbh_read_csv(file.path(out,"smooth_summaries.csv"))
manifest <- cbh_read_csv(file.path(out,"fit_manifest.csv"))
stopifnot(nrow(diag)==28,all(diag$converged),all(diag$finite_coefficients),all(diag$finite_covariance),
  all(diag$input_verified),all(diag$rank==diag$coefficients))
ages <- cbh_config()$age_bands$age_band
series <- c("map_matched","snow_matched","map_full","snow_full")
labels <- c(map_matched="MAP: matched through 2015",snow_matched="Snow: matched through 2015",
  map_full="MAP: full fitting period",snow_full="Snow: full sample through 2015")
for(a in ages) {
  z <- diag[diag$age_band==a & diag$series %in% series[1:2],]
  stopifnot(nrow(z)==2,length(unique(z$nuisance_hash))==1,length(unique(z$rows))==1,length(unique(z$deaths))==1)
}
cbh_unique(d,c("series","age_band","pfpr_pct"),"MAP/Snow curve grid")
d$label <- factor(labels[d$series],levels=labels)
d$age_label <- factor(d$age_band,levels=ages,labels=ifelse(ages=="<1","<1 month",paste0(ages," months")))
colors <- setNames(c("#215E91","#C26B28","#6C9BC2","#DCA878"),labels)
types <- setNames(c("solid","solid","dashed","dashed"),labels)
plot_curve <- function(z,title,subtitle) ggplot(z[z$within_central_support,],aes(pfpr_pct,log_hazard_ratio,colour=label,fill=label,linetype=label))+
  geom_hline(yintercept=0,colour="grey75",linewidth=.3)+
  geom_vline(xintercept=20,colour="grey85",linewidth=.3)+
  geom_ribbon(aes(ymin=lower_95,ymax=upper_95),alpha=.12,colour=NA,linetype=0)+
  geom_line(linewidth=.9)+facet_wrap(~age_label,ncol=4)+
  scale_colour_manual(values=colors)+scale_fill_manual(values=colors)+scale_linetype_manual(values=types)+
  scale_x_continuous(limits=c(0,100),breaks=seq(0,100,20))+
  labs(title=title,subtitle=subtitle,x="PfPR[2-10] (%)",y="Log mortality hazard ratio relative to PfPR = 20%",
    colour=NULL,fill=NULL,linetype=NULL,caption=paste(
      "Seven separate age-band models; gamma = 2; PfPR cr, k = 5; calendar year cr, k = 6.",
      "Curves cover each source's central 95% exposure range. Shading: conditional pointwise 95% intervals.",
      "Fixed median child HIV imputation; exposure, survey-design and smoothing-parameter uncertainty omitted.",sep="\n"))+
  theme_minimal(base_size=11)+theme(legend.position="bottom",panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold",size=16),strip.text=element_text(face="bold"),plot.caption=element_text(hjust=0,size=9))
plots <- list(
  matched_pfpr_splines=plot_curve(d[d$series %in% series[1:2],],"MAP versus Snow: gamma = 2","Same pre-2016 records, outcomes, confounders and offsets; source-specific exposure knots"),
  full_sample_pfpr_splines=plot_curve(d[d$series %in% series[3:4],],"MAP versus Snow: available fitting samples","Both gamma = 2; this comparison also changes the fitting period and sample"),
  all_pfpr_splines=plot_curve(d,"MAP versus Snow: gamma = 2","Matched comparison and full available samples"))
for(nm in names(plots)) ggsave(file.path(out,paste0(nm,".png")),plots[[nm]],width=14,height=8.5,dpi=180,device=ragg::agg_png,bg="white")
edf <- sm[sm$term=="s(pfpr_pct)",c("series","age_band","edf")]
cbh_atomic_csv(edf,file.path(out,"pfpr_edf.csv"))
hr <- d[d$pfpr_pct==40,c("series","age_band","log_hazard_ratio","standard_error","pfpr_p025","pfpr_p975")]
hr$hazard_ratio_40_to_20 <- exp(-hr$log_hazard_ratio)
hr$lower_95 <- exp(-hr$log_hazard_ratio-1.96*hr$standard_error)
hr$upper_95 <- exp(-hr$log_hazard_ratio+1.96*hr$standard_error)
hr$within_central_support <- hr$pfpr_p025<=20 & hr$pfpr_p975>=40
cbh_atomic_csv(hr,file.path(out,"pfpr_40_to_20_contrasts.csv"))
zero <- d[d$pfpr_pct==0,c("series","age_band","log_hazard_ratio","standard_error","pfpr_min")]
zero$hazard_ratio_20_to_zero <- exp(zero$log_hazard_ratio)
zero$attributable_fraction_at_20pct <- 1-zero$hazard_ratio_20_to_zero
zero$lower_95 <- 1-exp(zero$log_hazard_ratio+1.96*zero$standard_error)
zero$upper_95 <- 1-exp(zero$log_hazard_ratio-1.96*zero$standard_error)
zero$zero_below_observed_support <- zero$pfpr_min>0
cbh_atomic_csv(zero,file.path(out,"pfpr_20_to_zero_contrasts.csv"))
counts <- do.call(rbind,lapply(series,function(s) {
  z <- diag[diag$series==s,]
  data.frame(series=s,child_band_records=sum(z$rows),deaths=sum(z$deaths))
}))
cbh_atomic_csv(counts,file.path(out,"sample_totals.csv"))
fmt <- function(x) format(round(x),big.mark=",",trim=TRUE)
cell <- function(s,a) {z <- hr[hr$series==s & hr$age_band==a,];sprintf("%.3f (%.3f–%.3f)",z$hazard_ratio_40_to_20,z$lower_95,z$upper_95)}
ecell <- function(s,a) sprintf("%.2f",edf$edf[edf$series==s & edf$age_band==a])
zcell <- function(s,a) sprintf("%.1f%%",100*zero$attributable_fraction_at_20pct[zero$series==s & zero$age_band==a])
flags <- diag$fit_id[!is.finite(diag$min_smoothing_hessian_eigenvalue) | diag$min_smoothing_hessian_eigenvalue<=0]
md <- c("# MAP versus Snow, gamma=2","",
  "The primary exposure-source comparison uses identical child-band records through 2015. MAP and Snow fits retain the same outcomes, confounders, covariate scaling, offsets and median child HIV-incidence imputation. PfPR differs by source; each source retains its original exposure knot locations. All fits use cr (PfPR k=5, calendar-year k=6), gamma=2, unweighted binomial cloglog likelihood and separately estimated age-band coefficients and random effects.","",
  "The full-period MAP models used for the previous burden analysis are also refitted at gamma=2. The existing full-sample Snow gamma=2 fits are reused. Comparing those two available samples combines exposure-source and sample/period differences.","",
  "| Series | Child-band records | Deaths |","|---|---:|---:|",
  vapply(seq_len(nrow(counts)),function(i)sprintf("| %s | %s | %s |",labels[counts$series[i]],fmt(counts$child_band_records[i]),fmt(counts$deaths[i])),""),"",
  "![Matched curves](matched_pfpr_splines.png)","",
  "In the matched comparison, the Snow curves have lower effective degrees of freedom in every age band. MAP generally rises more sharply at low prevalence and more often flattens or turns down at high prevalence, especially at older ages. At 20% PfPR, the estimated contrast to zero gives a lower attributable fraction with Snow in every band. Differences over 40% to 20% are less uniform, illustrating why the entire curve matters for the national mortality calculation.","",
  "## Effective degrees of freedom","",
  "| Completed months | MAP matched | Snow matched | MAP full | Snow full |","|---|---:|---:|---:|---:|",
  vapply(ages,function(a)paste0("| ",a," | ",paste(vapply(series,ecell,"",a=a),collapse=" | ")," |"),""),"",
  "## Hazard ratio for PfPR 40% to 20%","",
  "| Completed months | MAP matched | Snow matched | MAP full | Snow full |","|---|---:|---:|---:|---:|",
  vapply(ages,function(a)paste0("| ",a," | ",paste(vapply(series,cell,"",a=a),collapse=" | ")," |"),""),"",
  "Intervals are conditional pointwise 95% intervals. Between-fit differences are descriptive; sampling covariance between fits is not estimated. A 40-to-20 comparison need not rank curves the same way as a current-to-zero burden comparison.","",
  "## Attributable fraction at 20% PfPR relative to zero","",
  "This contrast is closer to the national burden estimand. It depends on the low-exposure part of the spline, and may extrapolate to zero outside observed support; flags and conditional intervals are in the CSV.","",
  "| Completed months | MAP matched | Snow matched | MAP full | Snow full |","|---|---:|---:|---:|---:|",
  vapply(ages,function(a)paste0("| ",a," | ",paste(vapply(series,zcell,"",a=a),collapse=" | ")," |"),""),"",
  "## Checks and limitations","",
  "All selected fits converged, have full rank and finite coefficients/covariance, and reproduce the stored input records. Paired MAP/Snow nuisance inputs are identical, checked over every record. Compact PfPR contrasts and variances match full prediction-matrix calculations at zero and nonzero exposures. Original knot locations are verified. Model hashes and numerical diagnostics are saved.",
  if(length(flags)) paste("Smoothing-Hessian flags:",paste(flags,collapse=", ")) else "All selected smoothing Hessians are positive.",
  "Gamma changes smoothing selection throughout each model. A tighter-tolerance restart is selected only if it improves numerical diagnostics within the same specification. No fREML comparison is used to select between exposure sources.","",
  "Snow is available only through 2015. Any later mortality scenario using these curves requires transport over time and an explicitly specified external exposure. Zero PfPR can lie outside observed exposure support. Fixed Snow/MAP exposures and a single HIV imputation omit their uncertainty; survey design and residual within-child dependence are also omitted.","",
  "- [Full-sample curve comparison](full_sample_pfpr_splines.png) · [All curves](all_pfpr_splines.png)",
  "- [Burden comparison for 2005, 2015 and 2024](burden/REPORT.md)",
  "- [EDF](pfpr_edf.csv) · [40-to-20 contrasts](pfpr_40_to_20_contrasts.csv) · [20-to-zero contrasts](pfpr_20_to_zero_contrasts.csv) · [Diagnostics](fit_diagnostics.csv) · [Fit manifest](fit_manifest.csv)",
  "- Reproduce: Rscript R_cbh/snow/08_fit_map_comparison_gamma2.R; Rscript R_cbh/snow/09_report_map_comparison_gamma2.R.")
writeLines(md,file.path(out,"REPORT.md"))
inputs <- c("R_cbh/snow/09_report_map_comparison_gamma2.R",file.path(out,c("pfpr_curves.csv","fit_diagnostics.csv","smooth_summaries.csv","fit_manifest.csv")))
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"report_provenance.csv"))
message("MAP/Snow spline report complete")

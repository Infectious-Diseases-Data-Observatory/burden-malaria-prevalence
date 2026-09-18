#!/usr/bin/env Rscript
# Plot aggregate Snow-only penalty comparisons without refitting or microdata.
source("R_cbh/load_pipeline.R")
library(ggplot2)
out <- "results/cbh/age_band_snow_2000_2015_v1/penalty_sensitivity"
args <- commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% "--gamma2"))
include_gamma2 <- "--gamma2" %in% args
source_out <- out
d <- cbh_read_csv(file.path(out,"pfpr_curves.csv"))
diag <- cbh_read_csv(file.path(out,"fit_diagnostics.csv"))
sm <- cbh_read_csv(file.path(out,"smooth_summaries.csv"))
manifest <- cbh_read_csv(file.path(out,"fit_manifest.csv"))
if(include_gamma2) {
  extra <- function(name,original) {
    x <- cbh_read_csv(file.path(source_out,"gamma2",name))
    stopifnot(setequal(unique(x$series),c("reference","gamma2")))
    rbind(original,x[x$series=="gamma2",])
  }
  # Verify the additional run reused the exact original reference objects.
  ref2 <- cbh_read_csv(file.path(source_out,"gamma2","fit_manifest.csv"))
  ref1 <- manifest[manifest$series=="reference",]
  ref2 <- ref2[match(ref1$fit_id,ref2$fit_id),]
  stopifnot(identical(ref1$md5,ref2$md5),identical(ref1$input_signature,ref2$input_signature))
  d <- extra("pfpr_curves.csv",d)
  diag <- extra("fit_diagnostics.csv",diag)
  sm <- extra("smooth_summaries.csv",sm)
  manifest <- extra("fit_manifest.csv",manifest)
  out <- file.path(source_out,"gamma2","comparison")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  for(nm in c("pfpr_curves","fit_diagnostics","smooth_summaries","fit_manifest"))
    cbh_atomic_csv(get(switch(nm,pfpr_curves="d",fit_diagnostics="diag",smooth_summaries="sm",fit_manifest="manifest")),file.path(out,paste0(nm,".csv")))
}
n_series <- if(include_gamma2) 4L else 3L
stopifnot(nrow(diag)==7*n_series,nrow(manifest)==7*n_series,all(diag$converged),all(diag$finite_covariance),
  all(diag$input_verified),all(diag$rank==diag$coefficients))
ages <- cbh_config()$age_bands$age_band
for(a in ages) {
  z <- diag[diag$age_band==a,]
  stopifnot(nrow(z)==n_series,length(unique(z$rows))==1,length(unique(z$deaths))==1)
}
cbh_unique(d,c("series","age_band","pfpr_pct"),"Penalty curve grid")
stopifnot(all(abs(d$log_hazard_ratio[d$pfpr_pct==20])<1e-10),
  all(d$standard_error[d$pfpr_pct==20]<1e-10))
labels <- c(reference="Current: cr, gamma = 1",gamma14="cr, gamma = 1.4",cs="PfPR cs, gamma = 1")
if(include_gamma2) labels <- append(labels,c(gamma2="cr, gamma = 2"),after=2)
d$label <- factor(labels[d$series],levels=labels)
d$age_label <- factor(d$age_band,levels=ages,labels=ifelse(ages=="<1","<1 month",paste0(ages," months")))
central <- d[d$within_central_support,]
limits <- range(central$lower_95,central$upper_95)
colors <- setNames(if(include_gamma2) c("#59636B","#CA681E","#7452A5","#207E76") else c("#59636B","#CA681E","#207E76"),labels)
types <- setNames(if(include_gamma2) c("solid","dashed","longdash","dotdash") else c("solid","dashed","dotdash"),labels)
plot_curves <- function(x,title,subtitle) ggplot(x,aes(pfpr_pct,log_hazard_ratio,colour=label,fill=label,linetype=label))+
  geom_hline(yintercept=0,colour="grey75",linewidth=.3)+
  geom_vline(xintercept=20,colour="grey85",linewidth=.3)+
  geom_ribbon(aes(ymin=lower_95,ymax=upper_95),alpha=.08,colour=NA,linetype=0)+
  geom_line(linewidth=.85)+facet_wrap(~age_label,ncol=4)+
  scale_colour_manual(values=colors)+scale_fill_manual(values=colors)+scale_linetype_manual(values=types)+
  scale_x_continuous(limits=c(0,100),breaks=seq(0,100,20))+coord_cartesian(ylim=limits)+
  labs(title=title,subtitle=subtitle,x="Snow PfPR[2-10] (%)",y="Log mortality hazard ratio\nrelative to PfPR = 20%",
    colour=NULL,fill=NULL,linetype=NULL,caption=paste(
      "Identical Snow records through 2015; seven separate age-band models; PfPR k = 5 and calendar-year k = 6.",
      "Shading: pointwise 95% conditional intervals. Curves shown over the same central 95% exposure range within each age band.",
      "Snow prevalence and one median child HIV imputation are fixed. Survey-design and smoothing-parameter uncertainty are omitted.",sep="\n"))+
  theme_minimal(base_size=11)+theme(legend.position="bottom",panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold",size=16),strip.text=element_text(face="bold"),plot.caption=element_text(hjust=0,size=9))
plots <- list(
  all_penalty_splines=plot_curves(central,"Snow prevalence: spline penalisation sensitivities","Separate changes: higher gamma, or a PfPR shrinkage cubic basis"),
  gamma14_splines=plot_curves(central[central$series %in% c("reference","gamma14"),],"Snow prevalence: gamma = 1.4","Both use cr for PfPR and time; gamma changes smoothing selection for the model as a whole"),
  cs_splines=plot_curves(central[central$series %in% c("reference","cs"),],"Snow prevalence: PfPR shrinkage basis","Only the PfPR basis changes from cr to cs; gamma stays at 1 and calendar time stays cr"))
if(include_gamma2) plots$gamma2_splines <- plot_curves(central[central$series %in% c("reference","gamma14","gamma2"),],
  "Snow prevalence: gamma = 1, 1.4 and 2","All use cr for PfPR and time; identical records, covariates and knots")
for(nm in names(plots)) ggsave(file.path(out,paste0(nm,".png")),plots[[nm]],width=14,height=8.5,dpi=180,device=ragg::agg_png,bg="white")
hr <- d[d$pfpr_pct==40,c("series","age_band","log_hazard_ratio","standard_error","pfpr_p025","pfpr_p975")]
hr$hazard_ratio_40_to_20 <- exp(-hr$log_hazard_ratio)
hr$lower_95 <- exp(-hr$log_hazard_ratio-1.96*hr$standard_error)
hr$upper_95 <- exp(-hr$log_hazard_ratio+1.96*hr$standard_error)
hr$within_central_support <- hr$pfpr_p025<=20 & hr$pfpr_p975>=40
cbh_atomic_csv(hr,file.path(out,"pfpr_40_to_20_contrasts.csv"))
edf <- sm[sm$term=="s(pfpr_pct)",c("series","age_band","edf")]
cbh_atomic_csv(edf,file.path(out,"pfpr_edf.csv"))
reference <- d[d$series=="reference",c("age_band","pfpr_pct","log_hazard_ratio")]
names(reference)[3] <- "reference_log_hr"
delta <- merge(d[d$series!="reference" & d$within_central_support,],reference,by=c("age_band","pfpr_pct"))
delta$difference <- delta$log_hazard_ratio-delta$reference_log_hr
differences <- do.call(rbind,lapply(split(delta,paste(delta$series,delta$age_band)),function(z)data.frame(
  series=z$series[1],age_band=z$age_band[1],max_absolute_log_hr_change=max(abs(z$difference)),
  rms_log_hr_change=sqrt(mean(z$difference^2)))))
cbh_atomic_csv(differences,file.path(out,"curve_differences.csv"))
fmt <- function(x) format(round(x),big.mark=",",trim=TRUE)
cell <- function(series,a) {
  z <- hr[hr$series==series & hr$age_band==a,]
  sprintf("%.3f (%.3f–%.3f)",z$hazard_ratio_40_to_20,z$lower_95,z$upper_95)
}
edf_cell <- function(series,a) sprintf("%.2f",edf$edf[edf$series==series & edf$age_band==a])
baseline <- diag[diag$series=="reference",]
flags <- diag$fit_id[is.finite(diag$min_smoothing_hessian_eigenvalue) & diag$min_smoothing_hessian_eigenvalue<=0]
md <- c("# Snow-only spline penalisation sensitivities","",
  sprintf("All analyses use the same %s child-band records and %s deaths from 69 surveys in 31 countries as the full Snow analysis. Surveys and individual interviews after 2015 remain excluded. No MAP estimates or MAP-matched sample restriction enter these fits.",fmt(sum(baseline$rows)),fmt(sum(baseline$deaths))),"",
  "| Series | PfPR basis | Gamma | Calendar-year basis |",
  "|---|---|---:|---|","| Current reference | cr, k=5 | 1 | cr, k=6 |",
  "| Higher gamma | cr, k=5 | 1.4 | cr, k=6 |",if(include_gamma2) "| Higher gamma | cr, k=5 | 2 | cr, k=6 |","| PfPR shrinkage | cs, k=5 | 1 | cr, k=6 |","",
  "These sensitivities change gamma or the PfPR basis separately. Each age band retains its own time function, confounder coefficients and survey/country/region random intercepts. The unweighted binomial cloglog likelihood, band-width offset, covariates and their scaling, and fixed median child HIV-incidence imputation are unchanged. The exact PfPR and calendar-year knots from each reference model are reused.","",
  "Increasing gamma changes smoothing selection across the model, including time and random effects. The cs analysis changes only the PfPR basis: it penalises the otherwise unpenalised linear component as well as curvature, allowing the entire effect to shrink toward zero. Neither change enforces monotonicity. select=FALSE is retained in all fits.","",
  "![All comparisons](all_penalty_splines.png)","",
  "The overall age-specific curve shapes are similar across specifications. Higher gamma generally makes modest changes; it does not remove the curvature in older age bands. The cs fit attenuates the neonatal curve, while the older-age patterns persist. The tables below quantify these descriptive differences.","",
  if(include_gamma2) "With gamma=2, curvature persists at ages 12–47 months. The estimated 40% to 20% hazard ratios move from 0.904 to 0.870 at 12–23 months and from 0.856 to 0.822 at 24–35 months. Thus the additional smoothing does not necessarily attenuate the estimated association: calendar-year and random-effect smoothing are also re-estimated.","",
  "## PfPR effective degrees of freedom","",
  if(include_gamma2) "| Completed months | Reference | Gamma 1.4 | Gamma 2 | cs |" else "| Completed months | Reference | Gamma 1.4 | cs |",
  if(include_gamma2) "|---|---:|---:|---:|---:|" else "|---|---:|---:|---:|",
  vapply(ages,function(a)paste0("| ",a," | ",edf_cell("reference",a)," | ",edf_cell("gamma14",a)," | ",if(include_gamma2) paste0(edf_cell("gamma2",a)," | "),edf_cell("cs",a)," |"),""),"",
  "For cr, EDF near one is approximately linear; cs may shrink below one. More penalisation need not lower every individual term's EDF when all parameters are re-estimated.","",
  "## Hazard ratios for 40% to 20% PfPR","",
  if(include_gamma2) "| Completed months | Reference | Gamma 1.4 | Gamma 2 | cs |" else "| Completed months | Reference | Gamma 1.4 | cs |",
  if(include_gamma2) "|---|---:|---:|---:|---:|" else "|---|---:|---:|---:|",
  vapply(ages,function(a)paste0("| ",a," | ",cell("reference",a)," | ",cell("gamma14",a)," | ",if(include_gamma2) paste0(cell("gamma2",a)," | "),cell("cs",a)," |"),""),"",
  "Intervals are pointwise 95% conditional model intervals. Exposure support flags are in the contrast CSV. Differences between fits are descriptive: the same children appear in each analysis, and covariance between model estimates is not estimated.","",
  "## Validation","",
  paste0(if(include_gamma2) "Twenty-one sensitivity fits" else "Fourteen sensitivity fits", " are compared with seven saved Snow reference fits. All selected fits converged with full rank and finite coefficients/covariances. Every stored outcome, predictor and offset was checked against the reference model frame. Basis classes and original knot locations were verified. Prediction-matrix contrasts match direct link predictions, with zero contrast and uncertainty at 20% PfPR."),
  if(length(flags)) paste("Nonpositive smoothing-Hessian flags remain:",paste(flags,collapse=", ")) else "All selected smoothing Hessians are positive.",
  "Numerically flagged new fits receive a tighter-tolerance restart with the same gamma and basis; a restart is selected only on convergence, Hessian, gradient and fREML criteria within the same specification. Original attempts are retained. The selected-file manifest records the exact fitted objects. fREML criteria are not used to rank the different gamma/basis specifications.","",
  "These intervals condition on Snow point estimates, the HIV imputation and smoothing parameters. They omit source prevalence uncertainty, survey design and residual within-child dependence; source-model limitations from the original Snow report continue to apply.","",
  "## Files","",
  if(include_gamma2) "- [Gamma 1, 1.4 and 2 comparison](gamma2_splines.png)",
  "- [Gamma comparison](gamma14_splines.png) · [cs comparison](cs_splines.png)",
  "- [Curves](pfpr_curves.csv) · [EDF](pfpr_edf.csv) · [40-to-20 contrasts](pfpr_40_to_20_contrasts.csv) · [Curve differences](curve_differences.csv)",
  "- [Diagnostics](fit_diagnostics.csv) · [Selected-fit manifest](fit_manifest.csv) · [Smooth summaries](smooth_summaries.csv)",
  if(include_gamma2) "- Add gamma=2 with R_cbh/snow/06_fit_penalties.R --gamma2 and R_cbh/snow/07_report_penalties.R --gamma2. Existing gamma=1.4 and cs outputs are retained in the parent directory.",
  "- Reproduce with R_cbh/snow/06_fit_penalties.R followed by R_cbh/snow/07_report_penalties.R. Private fits remain under data/derived_cbh/models/age_band_snow_2000_2015_v1/penalty_sensitivity/.")
writeLines(md,file.path(out,"REPORT.md"))
inputs <- c("R_cbh/snow/07_report_penalties.R",file.path(out,c("pfpr_curves.csv","fit_diagnostics.csv","smooth_summaries.csv","fit_manifest.csv")))
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"report_provenance.csv"))
message("Snow penalty report complete")

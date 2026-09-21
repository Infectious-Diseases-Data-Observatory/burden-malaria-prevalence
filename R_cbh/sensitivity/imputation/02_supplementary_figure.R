#!/usr/bin/env Rscript
# Supplementary figure: PfPR splines from the complete-case primary analysis and
# the imputed-covariate sensitivity version, drawn from the saved comparison
# curves only (no refitting). Writes a PNG, a caption and provenance under the
# imputed version's results directory; never edits manuscript or TeX files.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
library(ggplot2)
settings <- cbh_primary_settings("regional_imputed");out <- settings$out
mi_dir <- file.path(out,"multiple_imputation")
inputs <- c(file.path(out,"comparison_pfpr_curves.csv"),file.path(out,"comparison_sample.csv"),
  file.path(out,"comparison_contrasts.csv"),file.path(mi_dir,"pooled_contrasts.csv"),
  file.path(out,"imputation_record_counts.csv"))
stopifnot(file.exists(inputs))
curves <- cbh_read_csv(inputs[1]);sample <- cbh_read_csv(inputs[2])
contrasts <- cbh_read_csv(inputs[3]);pooled <- cbh_read_csv(inputs[4])
counts <- cbh_read_csv(inputs[5])
stopifnot(setequal(unique(curves$iteration),settings$comparison_labels))

# Paper-facing labels replace the internal version identifiers.
paper_labels <- c("Complete-case covariates (primary analysis)","All covariate gaps imputed (sensitivity)")
curves$series_label <- factor(paper_labels[match(curves$iteration,settings$comparison_labels)],levels=paper_labels)
ages <- cbh_config()$age_bands$age_band
curves$age_label <- factor(curves$age_label,levels=ifelse(ages=="<1","<1 month",paste(ages,"months")))
colours <- setNames(c("#B36B39","#215E91"),paper_labels)
paper <- theme_minimal(base_size=18)+theme(panel.grid.minor=element_blank(),
  axis.text=element_text(size=15),strip.text=element_text(size=18,face="bold"),
  legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=16),
  panel.spacing=grid::unit(20,"pt"),plot.margin=margin(12,22,12,12))
central <- curves[curves$within_central_support,]
p <- ggplot(central,aes(pfpr_pct,log_hazard_ratio,colour=series_label,fill=series_label))+
  geom_hline(yintercept=0,colour="grey65",linewidth=.4)+
  geom_vline(xintercept=20,colour="grey80",linewidth=.4)+
  geom_ribbon(aes(ymin=lower_95,ymax=upper_95),alpha=.12,colour=NA)+
  geom_line(linewidth=1)+facet_wrap(~age_label,ncol=4)+
  scale_colour_manual(values=colours)+scale_fill_manual(values=colours)+
  scale_x_continuous(limits=c(0,80),breaks=c(0,20,40,60,80))+
  labs(x="PfPR[2–10] (%)",y="Log hazard ratio\nrelative to PfPR = 20%")+paper
figure <- file.path(out,"sfig_pfpr_splines_imputed_covariates.png")
ggsave(figure,p,width=14,height=8,dpi=220,device=ragg::agg_png,bg="white")

# Caption from saved numbers.
row <- function(m) sample[sample$measure==m,]
fmt <- function(x) format(round(x),big.mark=",")
hr <- function(iter,contrast) {
  d <- contrasts[contrasts$iteration==iter & contrasts$contrast==contrast,]
  d <- d[match(ages,d$age_band),];sprintf("%.2f",d$hazard_ratio)
}
max_share <- max(pooled$between_imputation_share)
stopifnot(counts$records==sample$revised[sample$measure=="Child-band records"])
caption <- c("# Supplementary figure caption: imputed-covariate sensitivity","",
  paste0("Sensitivity of the fitted PfPR–mortality relationship to the treatment of missing covariates. ",
  "Log mortality hazard ratios relative to PfPR[2–10] = 20% by completed-month age band, with pointwise 95% conditional intervals shaded, from the PfPR-ACM model fitted to the complete-case primary sample (",
  fmt(row("Child-band records")$previous)," child–age-band records, ",fmt(row("Deaths")$previous)," deaths, ",row("Surveys")$previous," surveys in ",row("Countries")$previous,
  " countries) and to every MAP-eligible record after imputing all remaining covariate gaps (",
  fmt(row("Child-band records")$revised)," records, ",fmt(row("Deaths")$revised)," deaths, ",row("Surveys")$revised," surveys in ",row("Countries")$revised,
  " countries, including all five Malaria Indicator Surveys). ",
  "Whole-survey gaps in the published regional indicators (wasting, stunting, facility delivery, electricity and the wealth score) were imputed by chained equations at the survey-region level (predictive mean matching, 10 imputations, point value their mean); the missing 2001 World Governance Indicators round by linear interpolation; health expenditure for Zimbabwe 2000–2009 and all countries in 2024 from a generalised additive model on the observed panel; and child HIV incidence for Liberia and São Tomé and Príncipe from the incidence model with a latent adolescent series. ",
  "Observed values were never replaced. The specification, reference knots, basis dimensions and gamma = 2 are those of the primary analysis; covariates were re-standardised on the enlarged sample. ",
  "Each curve is drawn over the central 95% of its own exposure distribution. ",
  "Hazard ratios for PfPR 40% to 20% are ",paste(hr(settings$comparison_labels[1],"40% to 20%"),collapse=", ")," (complete case) and ",
  paste(hr(settings$comparison_labels[2],"40% to 20%"),collapse=", ")," (imputed) for the seven bands from youngest to oldest. ",
  "Refitting the imputed version in each of the ten imputed datasets and pooling with Rubin's rules reproduced these hazard ratios to two decimals; the between-imputation share of interval variance was at most ",
  sprintf("%.1f%%",100*max_share),". The two samples are nested, so the comparison is descriptive."))
writeLines(caption,file.path(out,"CAPTION_sfig_imputed_covariates.md"))
prov <- data.frame(file=c(inputs,figure),md5=vapply(c(inputs,figure),cbh_file_hash,""),row.names=NULL)
cbh_atomic_csv(prov,file.path(out,"sfig_provenance.csv"))
cat("Supplementary imputed-covariate figure saved:",figure,"\n")

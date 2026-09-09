#!/usr/bin/env Rscript
# Report aggregate curves; never load research microdata or refit models.
source("R_cbh/load_pipeline.R")
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% "--retain-later-surveys"))
later <- "--retain-later-surveys" %in% args
id <- if(later) "age_band_snow_2000_2015_retrospective_v1" else "age_band_snow_2000_2015_v1"
out <- file.path("results/cbh",id)
d <- cbh_read_csv(file.path(out,"pfpr_curves_selected.csv"))
diag <- cbh_read_csv(file.path(out,"fit_diagnostics_selected.csv"))
jobs <- cbh_read_csv(file.path(out,"selected_fit_manifest.csv"))
selection <- cbh_read_csv(file.path(out,"selection.csv"))
sample <- cbh_read_csv(file.path(out,"sample_by_age_band.csv"))
stopifnot(nrow(diag)==nrow(jobs),all(diag$converged),all(diag$finite_coefficients),all(diag$finite_covariance),
  all(diag$rank==diag$coefficients))
cbh_unique(d,c("fit_id","age_band","pfpr_pct"),"Reported curve grid")
for(a in unique(d$age_band)) {
  z <- diag[diag$age_band==a,]
  pair <- z[z$series %in% c("map_matched",if("snow_matched" %in% z$series) "snow_matched" else "snow"),]
  stopifnot(nrow(pair)==2L,length(unique(pair$rows))==1L,length(unique(pair$deaths))==1L)
}
full <- cbh_read_csv("results/cbh/age_band_hiv_incidence_shared_time_v3/sensitivity_single_imputation/pfpr_curves.csv")
full <- full[grepl("^age_",full$fit_id),]; full$series <- "map_full"
full$fit_id <- paste0("map_full_",full$fit_id)
curves <- rbind(d,full[names(d)])
ages <- cbh_config()$age_bands$age_band
labels <- c(snow="Snow annual, through 2015",map_matched="MAP, same records",
  snow_matched="Snow, matched records",map_full="MAP, original full period")
colors <- setNames(c("#C35423","#176B87","#D09012","#777777"),labels)
curves$label <- factor(labels[curves$series],levels=labels)
curves$age_label <- factor(curves$age_band,levels=ages,
  labels=ifelse(ages=="<1","<1 month",paste0(ages," months")))
central <- curves[curves$within_central_support,]
range_y <- range(c(central$lower_95,central$upper_95))
style <- function(x,title,subtitle) ggplot(x,aes(pfpr_pct,log_hazard_ratio,colour=label,fill=label))+
  geom_hline(yintercept=0,colour="grey70",linetype="dashed",linewidth=.35)+
  geom_vline(xintercept=20,colour="grey85",linewidth=.3)+
  geom_ribbon(aes(ymin=lower_95,ymax=upper_95),alpha=.10,colour=NA)+
  geom_line(linewidth=.85)+facet_wrap(~age_label,ncol=4)+
  scale_colour_manual(values=colors,drop=TRUE)+scale_fill_manual(values=colors,drop=TRUE)+
  scale_x_continuous(limits=c(0,100),breaks=seq(0,100,20))+coord_cartesian(ylim=range_y)+
  labs(title=title,subtitle=subtitle,x="PfPR[2-10] (%)",y="Log mortality hazard ratio\nrelative to PfPR = 20%",colour=NULL,fill=NULL,
    caption=paste("Separate age-band fits; identical adjustment set and one fixed median child HIV-incidence imputation.",
      "Curves shown over their own central 95% exposure ranges. Shading: pointwise 95% conditional model intervals.",
      "Intervals omit prevalence-estimation, HIV-imputation, survey-design and smoothing-parameter uncertainty.",sep="\n"))+
  theme_minimal(base_size=11)+theme(legend.position="bottom",panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold",size=16),strip.text=element_text(face="bold"),
    plot.caption=element_text(hjust=0,size=9))
matched_snow <- if("snow_matched" %in% d$series) "snow_matched" else "snow"
plots <- list(
  snow_vs_map_splines=style(central[central$series %in% c(matched_snow,"map_matched"),],
    "Alternative prevalence estimates: Snow versus MAP","Same child-band records through 2015; smoothing re-estimated in each fit"),
  snow_map_full_period_splines=style(central,
    "Snow sensitivity in the context of the primary MAP analysis","Matched-period comparison and original full-period MAP curves"))
for(nm in names(plots)) ggsave(file.path(out,paste0(nm,".png")),plots[[nm]],width=14,height=8.5,dpi=180,
  device=ragg::agg_png,bg="white")
hr <- curves[curves$pfpr_pct==40,c("series","age_band","log_hazard_ratio","standard_error",
  "pfpr_min","pfpr_p025","pfpr_p975","pfpr_max")]
hr$hazard_ratio_40_to_20 <- exp(-hr$log_hazard_ratio)
hr$lower_95 <- exp(-hr$log_hazard_ratio-1.96*hr$standard_error)
hr$upper_95 <- exp(-hr$log_hazard_ratio+1.96*hr$standard_error)
hr$contrast_within_central_support <- hr$pfpr_p025<=20 & hr$pfpr_p975>=40
cbh_atomic_csv(hr,file.path(out,"pfpr_40_to_20_contrasts.csv"))
cbh_atomic_csv(unique(curves[c("series","age_band","pfpr_min","pfpr_p025","pfpr_p975","pfpr_max",
  "reference_within_observed_support")]),file.path(out,"pfpr_support.csv"))
delta <- merge(d[d$series==matched_snow,],d[d$series=="map_matched",],by=c("age_band","pfpr_pct"),suffixes=c("_snow","_map"))
delta <- delta[delta$within_central_support_snow & delta$within_central_support_map & delta$pfpr_pct%%.5==0,]
delta$difference <- delta$log_hazard_ratio_snow-delta$log_hazard_ratio_map
dist <- do.call(rbind,lapply(split(delta,delta$age_band),function(z)data.frame(age_band=z$age_band[1],
  overlap_min=min(z$pfpr_pct),overlap_max=max(z$pfpr_pct),grid_points=nrow(z),
  rms_log_hr_difference=sqrt(mean(z$difference^2)),max_absolute_log_hr_difference=max(abs(z$difference)))))
cbh_atomic_csv(dist,file.path(out,"curve_differences.csv"))
fmt <- function(x) format(round(x),big.mark=",",trim=TRUE)
line <- function(series,a) {
  z <- hr[hr$series==series & hr$age_band==a,]
  sprintf("%.2f (%.2f–%.2f)%s",z$hazard_ratio_40_to_20,z$lower_95,z$upper_95,
    ifelse(z$contrast_within_central_support,"","*"))
}
flagged <- diag[is.finite(diag$min_smoothing_hessian_eigenvalue) & diag$min_smoothing_hessian_eigenvalue<=0,]
stability <- cbh_read_csv(file.path(out,"smoothing_stability.csv"),required=FALSE)
stability_text <- if(!nrow(flagged)) "All selected fits have positive smoothing Hessians." else
  paste0("Nonpositive smoothing Hessians were recorded for: ",paste(flagged$fit_id,collapse=", "),
    ". These flags remain visible despite the convergence indicators.")
if(!is.null(stability)) stability_text <- paste(stability_text,
  paste(vapply(seq_len(nrow(stability)),function(i)sprintf(
    "For %s, a tighter-tolerance restart changed the central-support log-hazard-ratio curve by at most %.6f; the 40-to-20 HR was %.6f versus %.6f, and the restarted minimum Hessian eigenvalue was %.6g.",
    stability$fit_id[i],stability$max_central_log_hr_change[i],stability$original_hr40to20[i],
    stability$strict_hr40to20[i],stability$strict_min_hessian[i]),""),collapse=" "),
  paste0("Original fits are preserved. Tighter restarts selected on convergence, positive Hessian, reduced gradient and improved fREML: ",
    paste(stability$fit_id[stability$strict_selected],collapse=", "),
    ". These selected fits supply the reported comparison curves; see smoothing_stability.csv and selected_fit_manifest.csv."))
md <- c("# Snow annual prevalence sensitivity, 2000–2015","",
  sprintf("The alternative analysis includes %s child-band records and %s deaths from %d surveys in %d countries.",
    fmt(sum(sample$rows)),fmt(sum(sample$deaths)),sum(selection$snow_rows>0),length(unique(selection$country[selection$snow_rows>0]))),"",
  "## Exposure and eligibility","",
  "Use the posterior **mean** PfPR2–10 from the annual re-fit of the Snow survey database, not the published five-year Snow estimates. These are microscopy-equivalent percentages at a reference diagnostic/location/seasonal profile. Script 51 intersects each survey's DHS region polygons with Snow polygons in the same country, and weights intersection pieces by GPW 2020 population density × cell area × exact polygon overlap. The population weights remain fixed over time.","",
  if(later) "Retain retrospective bands from later surveys only when the full potential band finished by the end of 2015." else
    "Exclude surveys with metadata years after 2015 and individual interviews after December 2015. Require each full potential band to finish by the end of 2015.",
  "Exposure is joined by survey boundary version, validated region key and **band-entry calendar year**, exactly within 2000–2015. The `year` column in script 51 is the survey year; `period` is the annual exposure year. No 2015 carry-forward or temporal interpolation is used. Existing five-year retrospective eligibility and complete-potential-band rules are retained.","",
  "Exclude regions with less than 50% population coverage, following script 51's explicit threshold. Missing exposures are never replaced by zero. `selection.csv` separates time exclusions, missing Snow exposure, inadequate coverage and incomplete confounders. Existing unavailable dataset shards are listed separately.","",
  "## Model and comparison","",
  "Seven separate binomial complementary-log-log `mgcv::bam` models use the same PfPR cubic spline (k=5), one calendar-time spline per fit (k=6), full-band-years offset, confounders and survey/country/region random intercepts as the primary separate-age analysis. Smoothing parameters and coefficients are re-estimated. Sex, multiple birth, birth order, maternal age/education, wealth, urban residence, log child HIV incidence, log GDP, log health spending and political stability are included; vaccines are excluded. The saved full-sample covariate scaling and posterior-median child HIV imputation are held fixed. The likelihood remains unweighted.","",
  sprintf("Snow models use all eligible Snow records. %s of those records have MAP for the matched comparison; a separate matched Snow fit is added only if this differs from the Snow sample.",fmt(sum(selection$matched_map_rows))),
  "MAP is refitted on exactly the same records as its Snow comparator. Original full-period MAP curves provide context but also differ in period and sample composition. Curve differences are descriptive; no independence is assumed between fits of overlapping records.","",
  "![Snow versus matched MAP](snow_vs_map_splines.png)","",
  "## Hazard ratio for reducing PfPR from 40% to 20%","",
  "| Age (completed months) | Snow, matched | MAP, matched | MAP, full period |",
  "|---|---:|---:|---:|",
  vapply(ages,function(a)paste0("| ",a," | ",line(matched_snow,a)," | ",line("map_matched",a)," | ",line("map_full",a)," |"),""),"",
  "Intervals are pointwise 95% conditional model intervals. An asterisk flags contrasts outside that fit's central 95% exposure range. See the curves and support table; the single 40-to-20 contrast does not summarize all shape differences.","",
  "## Validation and limitations","",
  "All fits passed convergence, full-rank, finite coefficient/covariance and unchanged row-count checks. Prediction-matrix contrasts match direct link predictions, and both contrast and uncertainty are zero at the 20% reference. `fit_diagnostics.csv` records warnings, smoothing Hessians, gradients and fitting times; convergence alone does not establish adequate confounding control or causal identification.",
  stability_text,"",
  "The extraction audit checks the source CSV against its saved posterior draws and recomputes every region mean from weighted draws. It also checks country membership, unique annual keys, fractional-cell population weighting and absence of post-2015 carry-forward.","",
  "The Snow input README reports sparse polygon-year observations, no transmission-limits mask and weak information in 2014–2015. Source MCMC diagnostics include maximum R-hat 1.0595, eight parameters above 1.05 and 697 maximum-tree-depth events. The exposure fit is supplied by the user and is not refitted here. Its uncertainty is not propagated into mortality estimates; the CSV alone cannot yield regional posterior quantiles, which script 51 now leaves missing. The independent-polygon SD is a benchmark, not a guaranteed lower bound.","",
  "Conditional mortality intervals also omit HIV-imputation, survey-design and smoothing-parameter uncertainty. Partial spatial coverage, fixed 2020 population weights and retrospective assignment to residence at interview remain limitations.","",
  "## Outputs","",
  "- [Curve comparison with full-period MAP](snow_map_full_period_splines.png)",
  "- [Selected fitted curves](pfpr_curves_selected.csv), [40-to-20 contrasts](pfpr_40_to_20_contrasts.csv), [support](pfpr_support.csv)",
  "- [Selection](selection.csv), [country sample](sample_by_country.csv), [age-band sample](sample_by_age_band.csv)",
  "- [Region-year join audit](region_year_join_audit.csv), [selected-fit diagnostics](fit_diagnostics_selected.csv), [curve differences](curve_differences.csv)",
  "- [Saved-fit input checks](saved_fit_input_checks.csv); extraction checks are recorded in `extraction_validation.txt` and `aggregation_checks.csv`.",
  "- Exact model data and fits remain under the ignored `data/derived_cbh/models/` directory.")
writeLines(md,file.path(out,"REPORT.md"))
inputs <- c("R_cbh/snow/03_report.R",file.path(out,c("pfpr_curves_selected.csv","fit_diagnostics_selected.csv","selected_fit_manifest.csv","selection.csv")),
  if(!is.null(stability)) file.path(out,"smoothing_stability.csv"),
  "results/cbh/age_band_hiv_incidence_shared_time_v3/sensitivity_single_imputation/pfpr_curves.csv")
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"report_provenance.csv"))
message("Snow sensitivity report complete")

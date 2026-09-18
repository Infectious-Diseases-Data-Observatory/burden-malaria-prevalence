#!/usr/bin/env Rscript
# Figure 6: Nigerian state estimates versus IHME malaria deaths, 2024.
# Plotting only: reads saved state totals; never refits or recomputes exposure.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
source("R_cbh/reporting/log_axes.R")
library(ggplot2)
root <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional"))$out
out <- file.path(root,"nigeria_states")
model_label <- cbh_paper_model_label()
x <- cbh_read_csv(file.path(out,"state_totals_2024.csv"))
recon <- cbh_read_csv(file.path(out,"national_reconciliation_2024.csv"))
stopifnot(nrow(x)==37,all(x$year==2024),all(x$series=="map_full"),
  all(is.finite(x$attributable_under5_deaths)),all(x$attributable_under5_deaths>0),
  all(is.finite(x$ihme_malaria_deaths)),all(x$ihme_malaria_deaths>0))

# Nigeria's six geopolitical zones, for the legend only; the estimates are by state.
zones <- list(
  `North Central`=c("Benue","FCT (Abuja)","Kogi","Kwara","Nasarawa","Niger","Plateau"),
  `North East`=c("Adamawa","Bauchi","Borno","Gombe","Taraba","Yobe"),
  `North West`=c("Jigawa","Kaduna","Kano","Katsina","Kebbi","Sokoto","Zamfara"),
  `South East`=c("Abia","Anambra","Ebonyi","Enugu","Imo"),
  `South South`=c("Akwa Ibom","Bayelsa","Cross River","Delta","Edo","Rivers"),
  `South West`=c("Ekiti","Lagos","Ogun","Ondo","Osun","Oyo"))
zone_of <- setNames(rep(names(zones),lengths(zones)),unlist(zones))
stopifnot(length(zone_of)==37,setequal(names(zone_of),x$state))
x$zone <- factor(zone_of[x$state],levels=names(zones))
x$label <- sub(" [(]Abuja[)]$","",x$state)

# True log10 scales with identical limits on both axes, as in Figure 4.
log_limits <- 10^c(floor(log10(min(x$ihme_malaria_deaths,x$attributable_under5_deaths))),
  ceiling(log10(max(x$ihme_malaria_deaths,x$attributable_under5_deaths))))
log_breaks <- 10^seq(log10(log_limits[1]),log10(log_limits[2]))
log_labels <- function(v) parse(text=paste0("10^",round(log10(v))))
paper <- theme_minimal(base_size=14)+theme(panel.grid.minor=element_blank(),
  text=element_text(size=20),axis.title=element_text(size=20),axis.text=element_text(size=16),
  legend.title=element_blank(),legend.text=element_text(size=16),legend.position="bottom",
  plot.margin=margin(10,24,10,10))
p <- ggplot(x,aes(ihme_malaria_deaths,attributable_under5_deaths,colour=zone))+
  geom_abline(slope=1,intercept=0,colour="grey65",linetype=2)+
  geom_point(size=3.2)+
  ggrepel::geom_text_repel(aes(label=label),size=4.5,seed=20260917,max.overlaps=Inf,
    min.segment.length=0,segment.colour="grey70",show.legend=FALSE)+
  coord_equal(xlim=log_limits,ylim=log_limits,expand=FALSE)+
  scale_x_log10(breaks=log_breaks,labels=log_labels,minor_breaks=cbh_log10_minor_breaks)+
  scale_y_log10(breaks=log_breaks,labels=log_labels,minor_breaks=cbh_log10_minor_breaks)+
  scale_colour_manual(values=c(`North Central`="#E69F00",`North East`="#D55E00",`North West`="#CC79A7",
    `South East`="#009E73",`South South`="#0072B2",`South West`="#56B4E9"))+
  guides(colour=guide_legend(nrow=1,override.aes=list(size=4)))+
  labs(x="IHME malaria deaths before age 5, 2024",y=paste0(model_label," deaths\nbefore age 5, 2024"))+paper+cbh_log10_grid_theme()
stopifnot(is.null(p$labels$title),is.null(p$labels$subtitle),is.null(p$labels$caption),
  identical(p$scales$get_scales("x")$trans$name,"log-10"),identical(p$scales$get_scales("y")$trans$name,"log-10"))
ggsave(file.path(out,"fig6_nigeria_states_vs_ihme.png"),p,width=11,height=11,dpi=240,device=ragg::agg_png,bg="white")

## ---- caption, table and audit ---------------------------------------------------
fmt <- function(z) format(round(z),big.mark=",",trim=TRUE,scientific=FALSE)
ratio_nat <- recon$model_over_ihme_malaria_national
caption <- paste(
  paste0(model_label," malaria-attributable deaths before age 5 versus IHME cause-specific malaria deaths, by Nigerian state (36 states and the Federal Capital Territory) in 2024."),
  "Both axes use a base-10 logarithmic scale with identical limits; the dashed line denotes equality. Colours group states into Nigeria's six geopolitical zones for orientation; all estimates are made state by state.",
  "Model-attributable deaths equal IHME state all-cause deaths in each of seven age bands multiplied by 1 − exp[f_g(0) − f_g(P_state)], summed over bands, where P_state is population-weighted MAP PfPR[2–10] for 2024 and f_g are the seven separate primary MAP gamma=2 age-band effects.",
  "IHME early and late neonatal deaths share the <1-month effect; the IHME 2–4 year group shares one rate and divides deaths and person-time equally across the 24–35, 36–47 and 48–59 month bands. IHME malaria deaths are the state malaria death rate applied to the under-5 person-years implied by the all-cause count and rate.",
  sprintf("Summed over states the model gives %s deaths against %s IHME malaria deaths (ratio %.2f); the state sum differs from the national estimate evaluated at national mean prevalence by %+.1f%%, because the age-band curves are nonlinear.",
    fmt(recon$state_sum_attributable_deaths),fmt(recon$ihme_malaria_deaths_state_sum),ratio_nat,100*recon$relative_difference),
  "State totals are point estimates; conditional age-band intervals are retained in the source table and are not summed. Model-attributable all-cause reductions and IHME cause-specific malaria deaths are different estimands, and the model shares IHME all-cause inputs, so agreement or disagreement is not independent validation. Zero prevalence lies below observed exposure support.")
writeLines(c("# Figure 6 caption","",caption),file.path(out,"CAPTION.md"))
o <- x[order(-x$model_over_ihme_malaria),]
rows <- sprintf("| %s | %s | %.1f | %s | %s | %.1f%% | %s | %.1f%% | %.2f |",o$state,o$zone,o$pfpr_pct,
  fmt(o$ihme_under5_deaths),fmt(o$attributable_under5_deaths),100*o$attributable_fraction,
  fmt(o$ihme_malaria_deaths),100*o$ihme_malaria_share_of_allcause,o$model_over_ihme_malaria)
north <- grepl("^North",o$zone)
writeLines(c("# Figure 6: Nigerian states versus IHME, 2024","",
  "Secondary subnational extrapolation of the primary MAP gamma=2 model to the 37 Nigerian states, using the same convention as the national burden: IHME all-cause deaths by age band multiplied by the age-specific attributable fraction at the state's MAP prevalence. Calculation code is [R_cbh/burden/05_nigeria_state_burden.R](../../../../R_cbh/burden/05_nigeria_state_burden.R); this figure and table are produced by [R_cbh/reporting/08_nigeria_state_comparison.R](../../../../R_cbh/reporting/08_nigeria_state_comparison.R).","",
  "## Input audit","",
  "- Primary model: seven verified saved MAP gamma=2 fits; no refitting.",
  "- IHME all-cause by state and GBD age group, 2024: the seven constructed bands reproduce every state's under-5 total, and their sum over states reproduces the national 2024 age-band inputs used by `primary/02_effects.R` to within 1e-6, confirming the same export vintage.",
  "- IHME malaria by state, 2024: under-5 death rate per 100,000, applied to all-cause implied person-years.",
  "- MAP PfPR[2–10] by state, 2024: `data/pfpr_admin1_ng_cd_2024.csv`, an inherited extraction (MAP admin-1 boundaries, GPW 2020 density weights without cell area). All 37 states lie inside the central 95% of fitted prevalence; none requires extrapolation beyond observed support.",
  "- Source metadata and hashes: [ihme_source.json](ihme_source.json), [input_provenance.csv](input_provenance.csv).","",
  "## Reconciliation with the national estimate","",
  sprintf("State sum %s attributable deaths versus the national estimate %s evaluated at the national mean PfPR of %.2f%% (state person-year-weighted mean %.2f%%): difference %+.2f%%. The national and state exposure extractions use different boundary sources and population weights (see the note in [national_reconciliation_2024.csv](national_reconciliation_2024.csv)). IHME state malaria deaths sum to %s, so the model attributes %.2f times as many under-5 deaths to malaria as IHME's cause-specific estimate.",
    fmt(recon$state_sum_attributable_deaths),fmt(recon$national_attributable_deaths),recon$national_pfpr_pct,
    recon$state_person_year_weighted_pfpr_pct,100*recon$relative_difference,fmt(recon$ihme_malaria_deaths_state_sum),ratio_nat),"",
  "## Pattern of the discrepancy","",
  sprintf("The model exceeds IHME in %d of 37 states. The ratio of model to IHME deaths has median %.2f in the three northern zones and %.2f in the three southern zones; its Spearman correlation with state PfPR is %.2f. IHME's malaria share of all-cause under-5 deaths ranges from %.1f%% to %.1f%% and is highest in the southeast, whereas the model's attributable fraction ranges from %.1f%% to %.1f%% and rises with prevalence. The two methods therefore disagree on geography more than on the national total.",
    sum(o$model_over_ihme_malaria>1),median(o$model_over_ihme_malaria[north]),median(o$model_over_ihme_malaria[!north]),
    cor(o$pfpr_pct,o$model_over_ihme_malaria,method="spearman"),100*min(o$ihme_malaria_share_of_allcause),100*max(o$ihme_malaria_share_of_allcause),
    100*min(o$attributable_fraction),100*max(o$attributable_fraction)),"",
  "## State table","","Sorted by the ratio of model-attributable to IHME malaria deaths. Deaths are 2024 annual counts before age 5.","",
  paste0("| State | Zone | PfPR (%) | IHME all-cause | ",model_label," | Attributable fraction | IHME malaria | IHME malaria share | Model / IHME |"),
  "|---|---|---:|---:|---:|---:|---:|---:|---:|",rows,"",
  "[State totals](state_totals_2024.csv) · [State-by-age estimates](state_age_estimates_2024.csv) · [Reconciliation](national_reconciliation_2024.csv) · [Model specification](model_specification.txt)","",
  "![Figure 6](fig6_nigeria_states_vs_ihme.png)","",caption),file.path(out,"README.md"))
invisible(file.copy("data/external/nigeria_states/source.json",file.path(out,"ihme_source.json"),overwrite=TRUE))
paths <- c(file.path(out,c("state_totals_2024.csv","national_reconciliation_2024.csv")),
  "R_cbh/reporting/08_nigeria_state_comparison.R","R_cbh/reporting/labels.R","R_cbh/reporting/log_axes.R")
cbh_atomic_csv(data.frame(file=paths,md5=vapply(paths,cbh_file_hash,"")),file.path(out,"figure_provenance.csv"))
message("Figure 6 and state comparison table saved")

#!/usr/bin/env Rscript
# Aggregate-only comparison of revised regional adjustment with the saved benchmark.
# Writes CSV, PNG and Markdown only; never edits manuscript or TeX files.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
source("R_cbh/reporting/log_axes.R")
library(ggplot2)
settings <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional_mics"));out <- settings$out
stopifnot(Sys.getenv("CBH_PRIMARY_VERSION","regional_mics") %in% c("regional","regional_imputed","regional_mics","regional_mics_v5"))
formula_path <- file.path(out,"model_formula.txt")
writeLines(trimws(readLines(formula_path),which="right"),formula_path)
ages <- cbh_config()$age_bands$age_band
inputs <- character()
read <- function(name,previous=FALSE) {
  path <- file.path(if(previous) settings$reference else out,name)
  inputs <<- c(inputs,path)
  cbh_read_csv(path)
}
save_plot <- function(p,name,width=14,height=8) ggsave(file.path(out,name),p,
  width=width,height=height,dpi=220,device=ragg::agg_png,bg="white")
write <- function(d,name) cbh_atomic_csv(d,file.path(out,name))
labels <- settings$comparison_labels
colours <- setNames(c("#B36B39","#215E91"),labels)
paper <- theme_minimal(base_size=18)+theme(panel.grid.minor=element_blank(),
  axis.text=element_text(size=15),strip.text=element_text(size=18,face="bold"),
  legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=16),
  panel.spacing=grid::unit(20,"pt"),plot.margin=margin(12,22,12,12))
age_label <- function(x) factor(x,levels=ages,labels=ifelse(ages=="<1","<1 month",paste(ages,"months")))
diag <- read("fit_diagnostics.csv");checks <- read("fitted_outcome_checks.csv")
stopifnot(nrow(diag)==7L,all(diag$converged),all(diag$input_verified),all(diag$gamma==2),
  all(diag$min_smoothing_hessian_eigenvalue>0),all(checks$model_md5_verified))
curves <- do.call(rbind,lapply(c(TRUE,FALSE),function(old) {
  d <- read("pfpr_curves.csv",old);d <- d[d$series=="map_full",]
  d$iteration <- if(old) labels[1] else labels[2];d
}))
curves$iteration <- factor(curves$iteration,levels=labels)
curves$age_label <- age_label(curves$age_band)
write(curves,"comparison_pfpr_curves.csv")
# Both complete curves are exported. Display each central 95% exposure range,
# with zero included separately in the contrast table as an extrapolation.
central <- curves[curves$within_central_support,]
p <- ggplot(central,aes(pfpr_pct,log_hazard_ratio,colour=iteration,fill=iteration))+
  geom_hline(yintercept=0,colour="grey65",linewidth=.4)+
  geom_vline(xintercept=20,colour="grey80",linewidth=.4)+
  geom_ribbon(aes(ymin=lower_95,ymax=upper_95),alpha=.12,colour=NA)+
  geom_line(linewidth=1)+facet_wrap(~age_label,ncol=4)+
  scale_colour_manual(values=colours)+scale_fill_manual(values=colours)+
  scale_x_continuous(limits=c(0,80),breaks=c(0,20,40,60,80))+
  labs(x="PfPR[2–10] (%)",y="Log hazard ratio\nrelative to PfPR = 20%")+paper
save_plot(p,"comparison_pfpr_splines.png")

contrasts <- do.call(rbind,lapply(c(TRUE,FALSE),function(old) {
  a <- read("pfpr_40_to_20_contrasts.csv",old)
  b <- read("pfpr_20_to_zero_contrasts.csv",old)
  a <- a[match(ages,a$age_band),];b <- b[match(ages,b$age_band),]
  rbind(data.frame(age_band=ages,contrast="40% to 20%",hazard_ratio=a$hazard_ratio_40_to_20,
    lower_95=a$lower_95,upper_95=a$upper_95,zero_extrapolation=FALSE),
    data.frame(age_band=ages,contrast="20% to 0%",hazard_ratio=b$hazard_ratio_20_to_zero,
      lower_95=exp(b$log_hazard_ratio-1.96*b$standard_error),
      upper_95=exp(b$log_hazard_ratio+1.96*b$standard_error),zero_extrapolation=b$zero_below_observed_support)) |>
    transform(iteration=if(old) labels[1] else labels[2])
}))
contrasts$reduction_pct <- 100*(1-contrasts$hazard_ratio)
contrasts$reduction_lower_95 <- 100*(1-contrasts$upper_95)
contrasts$reduction_upper_95 <- 100*(1-contrasts$lower_95)
write(contrasts,"comparison_contrasts.csv")
contrasts$age_label <- factor(contrasts$age_band,levels=rev(ages))
contrasts$iteration <- factor(contrasts$iteration,levels=labels)
p <- ggplot(contrasts,aes(age_label,hazard_ratio,colour=iteration))+
  geom_hline(yintercept=1,colour="grey60",linetype=2)+
  geom_errorbar(aes(ymin=lower_95,ymax=upper_95),position=position_dodge(width=.5),width=.15)+
  geom_point(position=position_dodge(width=.5),size=2.5)+coord_flip()+
  facet_wrap(~contrast,nrow=1)+scale_y_log10(minor_breaks=cbh_log10_minor_breaks)+scale_colour_manual(values=colours)+
  labs(x="Age (months)",y="Mortality hazard ratio (95% interval)")+paper+cbh_log10_grid_theme()
save_plot(p,"comparison_contrasts.png",12,7)

sample <- read("prepared_sample.csv");old_sample <- read("primary_sample.csv",TRUE)
old_coverage <- read("survey_coverage.csv",TRUE)
old_counts <- c(old_sample$records,old_sample$distinct_children,old_sample$deaths,
  sum(old_coverage$regions),old_sample$surveys,old_sample$countries)
sample_comparison <- data.frame(measure=c("Child-band records","Children","Deaths","Survey-regions","Surveys","Countries"),
  previous=old_counts,revised=as.numeric(sample[1,c("records","children","deaths","regions","surveys","countries")]))
sample_comparison$change <- sample_comparison$revised-sample_comparison$previous
write(sample_comparison,"comparison_sample.csv")
edf <- read("pfpr_edf.csv");old_edf <- read("pfpr_edf.csv",TRUE)
edf_comparison <- data.frame(age_band=ages,previous_edf=old_edf$edf[match(ages,old_edf$age_band)],
  revised_edf=edf$edf[match(ages,edf$age_band)])
write(edf_comparison,"comparison_edf.csv")
age_results <- diag[match(ages,diag$age_band),c("age_band","rows","deaths")]
age_results$share_of_observed_under5_deaths_pct <- 100*age_results$deaths/sum(age_results$deaths)
age_results$pfpr_edf <- edf$edf[match(ages,edf$age_band)]
for(contrast in c("40% to 20%","20% to 0%")) {
  z <- contrasts[contrasts$iteration==labels[2] & contrasts$contrast==contrast,]
  z <- z[match(ages,z$age_band),]
  prefix <- if(contrast=="40% to 20%") "hr_40_to_20" else "hr_20_to_zero"
  for(v in c("hazard_ratio","lower_95","upper_95"))
    age_results[[paste0(prefix,if(v=="hazard_ratio") "" else paste0("_",v))]] <- z[[v]]
}
stopifnot(abs(sum(age_results$share_of_observed_under5_deaths_pct)-100)<1e-10,
  sum(age_results$deaths)==sample$deaths,sum(age_results$rows)==sample$records)
write(age_results,"primary_age_band_results.csv")

country <- read("burden/country_totals.csv");previous <- read("burden/country_totals.csv",TRUE)
j <- match(cbh_key(country,c("year","iso3")),cbh_key(previous,c("year","iso3")))
stopifnot(!anyNA(j),identical(country$pfpr_pct,previous$pfpr_pct[j]),
  identical(country$ihme_under5_deaths,previous$ihme_under5_deaths[j]),
  identical(country$ihme_malaria_deaths,previous$ihme_malaria_deaths[j]))
country$previous_deaths <- previous$attributable_under5_deaths[j]
country$change_deaths <- country$attributable_under5_deaths-country$previous_deaths
country$change_pct <- 100*country$change_deaths/country$previous_deaths
write(country,"burden/comparison_country_totals.csv")
year <- read("burden/year_summary.csv");old_year <- read("burden/year_summary.csv",TRUE)
j <- match(year$year,old_year$year)
year$previous_deaths <- old_year$attributable_under5_deaths[j]
year$change_deaths <- year$attributable_under5_deaths-year$previous_deaths
year$change_pct <- 100*year$change_deaths/year$previous_deaths
write(year,"burden/comparison_year_totals.csv")
bands <- read("burden/country_age_estimates.csv");old_bands <- read("burden/country_age_estimates.csv",TRUE)
by_age <- do.call(rbind,lapply(settings$years,function(y)do.call(rbind,lapply(ages,function(a) {
  z <- bands[bands$year==y & bands$age_band==a & is.finite(bands$attributable_deaths),]
  b <- old_bands[old_bands$year==y & old_bands$age_band==a & old_bands$iso3 %in% z$iso3,]
  stopifnot(nrow(z)==42,nrow(b)==42,setequal(z$iso3,b$iso3))
  data.frame(year=y,age_band=a,previous_deaths=sum(b$attributable_deaths),
    revised_deaths=sum(z$attributable_deaths),change_deaths=sum(z$attributable_deaths)-sum(b$attributable_deaths))
}))))
write(by_age,"burden/comparison_deaths_by_age.csv")
plot_rows <- country[is.finite(country$previous_deaths) & country$previous_deaths>0 &
  is.finite(country$attributable_under5_deaths) & country$attributable_under5_deaths>0,]
write(country[!cbh_key(country,c("year","iso3")) %in% cbh_key(plot_rows,c("year","iso3")),],
  "burden/comparison_log10_exclusions.csv")
limits <- 10^c(floor(log10(min(plot_rows$previous_deaths,plot_rows$attributable_under5_deaths))),
  ceiling(log10(max(plot_rows$previous_deaths,plot_rows$attributable_under5_deaths))))
log_labels <- function(x)parse(text=paste0("10^",round(log10(x))))
p <- ggplot(plot_rows,aes(previous_deaths,attributable_under5_deaths))+
  geom_abline(slope=1,intercept=0,linetype=2,colour="grey60")+
  geom_point(colour=colours[2],size=2.5,alpha=.8)+
  ggrepel::geom_text_repel(data=plot_rows[plot_rows$iso3 %in% c("COD","NGA","AGO","TZA","UGA"),],
    aes(label=iso3),size=4.5,seed=20260917,max.overlaps=Inf)+
  facet_wrap(~year,nrow=1)+scale_x_log10(labels=log_labels,minor_breaks=cbh_log10_minor_breaks)+scale_y_log10(labels=log_labels,minor_breaks=cbh_log10_minor_breaks)+
  coord_equal(xlim=limits,ylim=limits)+labs(x=paste0(cbh_paper_model_label(),": previous adjustment\nDeaths before age 5"),
    y=paste0(cbh_paper_model_label(),": revised adjustment\nDeaths before age 5"))+paper+cbh_log10_grid_theme()
save_plot(p,"burden/comparison_country_deaths.png",15,6)

fmt <- function(x)format(round(x),big.mark=",",scientific=FALSE,trim=TRUE)
hrfmt <- function(x) sprintf("%.3f (%.3f–%.3f)",x$hazard_ratio,x$lower_95,x$upper_95)
contrast_rows <- unlist(lapply(c("40% to 20%","20% to 0%"),function(c) {
  a <- contrasts[contrasts$contrast==c & contrasts$iteration==labels[1],]
  b <- contrasts[contrasts$contrast==c & contrasts$iteration==labels[2],]
  c(paste0("### PfPR ",c),"","| Age (months) | Previous HR (95% interval) | Revised HR (95% interval) |",
    "|---|---:|---:|",paste0("| ",ages," | ",hrfmt(a)," | ",hrfmt(b)," |"),"")
}))
old40 <- contrasts[contrasts$contrast=="40% to 20%" & contrasts$iteration==labels[1],]
new40 <- contrasts[contrasts$contrast=="40% to 20%" & contrasts$iteration==labels[2],]
largest <- which.max(abs(new40$hazard_ratio-old40$hazard_ratio))
old0 <- contrasts[contrasts$contrast=="20% to 0%" & contrasts$iteration==labels[1],]
new0 <- contrasts[contrasts$contrast=="20% to 0%" & contrasts$iteration==labels[2],]
older <- match(ages[4:6],old0$age_band)
interpretation <- c(sprintf("For PfPR 40%% to 20%%, the largest absolute change in the point-estimate HR is %.3f, at %s months (%.3f previously; %.3f revised).",
  abs(new40$hazard_ratio[largest]-old40$hazard_ratio[largest]),new40$age_band[largest],
  old40$hazard_ratio[largest],new40$hazard_ratio[largest]),"",
  paste0("For PfPR 20% to zero, the predicted mortality reductions at ",
    paste(ages[4:6],collapse=", ")," months change from ",
    paste(sprintf("%.1f%%",old0$reduction_pct[older]),collapse=", ")," to ",
    paste(sprintf("%.1f%%",new0$reduction_pct[older]),collapse=", "),", respectively."),"")
if(settings$imputed) {
  ic <- read("imputation_record_counts.csv")
  title <- "# PfPR-ACM model: imputed-covariate sensitivity"
  what_changed <- c(
    sprintf("Seven separate MAP age-band models, gamma=2, with the same 17-variable specification as the current primary (`%s`), fitted after imputing every remaining covariate gap so that all %s MAP-eligible child-band records are retained. All seven converged, have finite covariance and full rank, and pass the positive smoothing-Hessian and exact fitted-input checks.",basename(settings$reference),fmt(sample$records)),"",
    "## What changed","",
    sprintf("The complete-case primary excludes records whose survey-region or country-year covariates are unavailable after the UNICEF vaccination and within-survey region-mean fallbacks. Here those gaps are imputed instead (plan section 2.4): whole-survey gaps in wasting, stunting, facility delivery, electricity and the wealth score by chained-equation multiple imputation at the survey-region level (predictive mean matching, 10 imputations, point value = their mean); the missing 2001 WGI political-stability round by interpolation; health expenditure for Zimbabwe 2000–2009 and every country's 2024 from a GAM on the observed panel; and child HIV incidence for Liberia and São Tomé and Príncipe from the extended incidence model with a latent adolescent series. In this sample %s records carry at least one model-imputed regional covariate, %s an interpolated political-stability value, %s a modelled health-expenditure value and %s an HIV incidence imputed without an adolescent series.",
      fmt(ic$regional_any_model_imputed),fmt(ic$political_stability_interpolated),fmt(ic$health_expenditure_imputed),fmt(ic$hiv_no_adolescent_series)),"",
    "MAP band-entry exposure, age bands, full-band offset, unweighted binomial/cloglog likelihood, reference spline knots, cr basis dimensions and gamma=2 are held fixed. Confounders are re-scaled on the enlarged sample. Each age has its own time spline, covariate coefficients and random effects. Imputation uncertainty is not reflected in these intervals; the separate multiple-imputation check propagates it.","",
    "| Sample | Complete case | Imputed | Change |","|---|---:|---:|---:|",
    paste0("| ",sample_comparison$measure," | ",fmt(sample_comparison$previous)," | ",fmt(sample_comparison$revised)," | ",fmt(sample_comparison$change)," |"),"",
    "The imputed sample contains the complete-case sample entirely, plus the previously excluded records, 25 of them whole surveys including all five MIS surveys. Differences combine the added records and the imputed covariate values; they are not a test of either alone.","")
} else if(isTRUE(settings$liberia)) {
  lbr_check <- cbh_read_csv("results/cbh/hiv_incidence/liberia_aidsinfo/definition_check_country_years.csv")
  title <- "# Primary PfPR-ACM model: DHS and MICS surveys, with Liberia"
  what_changed <- c(
    sprintf("Seven separate MAP age-band models, gamma=2, with the same 17-variable specification and sample rules as `%s`, now including Liberia. All seven converged, have finite covariance and full rank, and pass the positive smoothing-Hessian and exact fitted-input checks.",basename(settings$reference)),"",
    "## What changed","",
    sprintf("Liberia had no child HIV incidence series in the UNAIDS/UNICEF workbook, so its surveys failed complete-case selection. Its child rate is now derived from the UNAIDS new HIV infections among children aged 0-14 published on AIDSinfo (epidemiological estimates 2026), divided by IHME-implied under-5 person-years less children aged 0-4 living with HIV; this reproduces the workbook's published child rates (median derived/published %.2f across %d country-years; [check](../hiv_incidence/liberia_aidsinfo/REPORT.md)). Liberia's 2007, 2013 and 2019-20 DHS enter; its 2009 MIS still lacks anthropometry. Every other country's HIV value is unchanged.",median(lbr_check$ratio),nrow(lbr_check)),"",
    "MAP band-entry exposure, age bands, full-band offset, unweighted binomial/cloglog likelihood, reference spline knots, cr basis dimensions and gamma=2 are held fixed. Confounders are re-scaled on the enlarged sample.","",
    "| Sample | DHS and MICS (v5) | With Liberia (v7) | Change |","|---|---:|---:|---:|",
    paste0("| ",sample_comparison$measure," | ",fmt(sample_comparison$previous)," | ",fmt(sample_comparison$revised)," | ",fmt(sample_comparison$change)," |"),"",
    "The v7 sample contains the v5 sample entirely. Differences reflect the added Liberia records and the rescaling of confounders on the larger sample.","")
} else if(settings$mics) {
  title <- "# Primary PfPR-ACM model: DHS and MICS surveys"
  what_changed <- c(
    sprintf("Seven separate MAP age-band models, gamma=2, with the same 17-variable specification as `%s`, fitted to the DHS complete-case sample plus every MICS survey with a complete birth history that passes the same complete-case selection. All seven converged, have finite covariance and full rank, and pass the positive smoothing-Hessian and exact fitted-input checks.",basename(settings$reference)),"",
    "## What changed","",
    "MICS birth histories were converted to the DHS Births Recode layout (R_mics/07_make_recodes.R) and passed through the same child age-band builder, with annual regional MAP PfPR extracted by the DHS method on DHS boundary polygons (new admin-1 polygons for Guinea-Bissau and the Central African Republic). All 13 regional covariates were computed from the MICS microdata (R_mics/09_regional_covariates.R); Guinea 2016, Comoros 2022 and Chad 2019 use the national WUENIC DTP3 and measles estimates because their recall doses are unusable. Complete-case selection then keeps only MICS surveys with anthropometry and complete national series, as for DHS. The DHS part reproduces the previous sample exactly.","",
    "MAP band-entry exposure, fixed median child HIV incidence imputation, age bands, full-band offset, unweighted binomial/cloglog likelihood, reference spline knots, cr basis dimensions and gamma=2 are held fixed. Confounders are re-scaled on the enlarged sample.","",
    "| Sample | DHS only | DHS and MICS | Change |","|---|---:|---:|---:|",
    paste0("| ",sample_comparison$measure," | ",fmt(sample_comparison$previous)," | ",fmt(sample_comparison$revised)," | ",fmt(sample_comparison$change)," |"),"",
    "The combined sample contains the DHS sample entirely. Differences reflect the added MICS records and the rescaling of confounders on the larger sample; MICS covariate definitions differ from DHS in documented ways (education level converted to years, facility delivery for the last birth in two years).","")
} else {
  title <- "# Primary PfPR-ACM model: revised regional adjustment"
  what_changed <- c(
  "Seven separate MAP age-band models, gamma=2, fitted with the revised 17-variable specification. All seven converged, have finite covariance and full rank, and pass the positive smoothing-Hessian and exact fitted-input checks. The saved PfPR components reproduce full model prediction contrasts and variances.","",
  "## What changed", "",
  "The previous fitted benchmark used 18 regional/annual covariates. The new model uses 13 survey-region summaries and four national annual covariates: sex, multiple births and birth order are removed; maternal age at each birth is replaced with age at first birth; wasting and stunting are added. Urban percentage, DTP3/measles coverage, facility delivery, short birth interval, water, sanitation and electricity are included. Hib3, PCV, rotavirus and exclusive breastfeeding are excluded. Regional gaps use the declared UNICEF and within-survey available-region fallback; whole-survey gaps remain excluded.","",
  "MAP band-entry exposure, fixed median child HIV incidence imputation, age bands, full-band offset, unweighted binomial/cloglog likelihood, reference spline knots, cr basis dimensions and gamma=2 are held fixed. Confounders are scaled on the new selected sample. Each age has its own time spline, covariate coefficients and random effects.","",
  "| Sample | Previous | Revised | Change |","|---|---:|---:|---:|",
  paste0("| ",sample_comparison$measure," | ",fmt(sample_comparison$previous)," | ",fmt(sample_comparison$revised)," | ",fmt(sample_comparison$change)," |"),"",
  "The availability audit confirms that the samples share 5,465,305 child-band records. The revision loses 214,812 previous records and recovers none. Differences therefore combine changes in adjustment definitions, added covariates, imputation and sample selection. They do not isolate any one of these changes; that would require additional fits on a common sample.","")
}
report <- c(title, "", what_changed,
  "## PfPR effects","","![PfPR spline comparison](comparison_pfpr_splines.png)","",
  "Curves are log hazard ratios relative to PfPR=20%; each is displayed over its own central 95% exposure range. Shading is a conditional 95% interval. The full 0–100% curves and support flags remain in comparison_pfpr_curves.csv. Intervals condition on smoothing parameters, exposure, one HIV imputation and other filled covariates. Overlapping-sample estimates are dependent; no formal test or interval for the difference between iterations is implied.","",
  interpretation,contrast_rows,
  "Zero PfPR is just below observed support in every band; the 20% to 0% contrasts involve extrapolation. HR below one indicates lower mortality at the lower prevalence. Conditional intervals are calculated from the joint covariance of both evaluation points.","",
  "## National attributable mortality","",
  "The same population-weighted national MAP values and IHME all-cause baselines are used in both iterations, for the same 42 estimable countries. These are predicted all-cause reductions under zero PfPR. Country-age contributions remain signed, with missing countries retained in source tables; marginal age-specific intervals are not summed into total intervals.","",
  "| Year | Previous deaths | Revised deaths | Change |","|---|---:|---:|---:|",
  paste0("| ",year$year," | ",fmt(year$previous_deaths)," | ",fmt(year$attributable_under5_deaths)," | ",sprintf("%+.1f%%",year$change_pct)," |"),"",
  "![Country estimate comparison](burden/comparison_country_deaths.png)","",
  "Country and age-specific changes are in burden/comparison_country_totals.csv and burden/comparison_deaths_by_age.csv. Identical national PfPR and IHME baselines were verified. The existing equal-person-time assumption for the IHME 2–4-year group is retained.","",
  "## Reproduction and output scope","",
  "Run `Rscript R_cbh/primary/run_regional.R` for fresh fits, `--resume` to reuse only verified fit caches, or `--report-only` to recalculate effects, diagnostics and comparisons from the saved fits. No raw-source extraction, HIV refitting, supplementary fitting, manuscript editing or TeX generation is performed.","",
  if(settings$imputed) sprintf("The comparator is the current complete-case primary in %s, which remains the primary analysis. This imputed-covariate version is a sensitivity analysis; it does not feed the paper figures, tables or burden estimates.",settings$reference)
  else if(settings$mics) sprintf("The comparator is the DHS-only primary in %s, preserved unchanged.",settings$reference)
  else "The previous 18-variable fits are preserved in primary_map_regional18_gamma2_v2. Current primary reporting, including annual 2000–2024 comparisons and Nigerian state estimates, is indexed in RESULTS.md and paper_figures/CAPTIONS.md. Subgroup and exposure sensitivity fits remain historical until separately refitted; they are outside this reporting refresh.","",
  "Numerical checks and provenance: fit_diagnostics.csv, fit_manifest.csv, fit_input_provenance.csv, fitted_outcome_checks.csv, comparison_edf.csv, and comparison_provenance.csv. In-sample outcome checks are not external validation or survey influence analyses.")
writeLines(report,file.path(out,"REPORT.md"))
inputs <- unique(c(inputs,formula_path,"R_cbh/primary/05_compare_regional.R","R_cbh/primary/settings.R",
  "R_cbh/reporting/labels.R","R_cbh/reporting/log_axes.R"))
write(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),"comparison_provenance.csv")
message("Revised primary comparison complete: ",out,"/REPORT.md")

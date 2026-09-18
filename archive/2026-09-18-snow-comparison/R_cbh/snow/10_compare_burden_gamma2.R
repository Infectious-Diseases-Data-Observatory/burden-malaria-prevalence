#!/usr/bin/env Rscript
# Hold national MAP exposure and IHME inputs fixed to compare fitted curves.
source("R_cbh/load_pipeline.R")
library(mgcv)
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(identical(args,"--exposure=map"))
root <- "results/cbh/map_snow_gamma2_v1"
out <- file.path(root,"burden")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
component_path <- "data/derived_cbh/models/map_snow_gamma2_v1/pfpr_components.rds"
components <- readRDS(component_path)
manifest <- cbh_read_csv(file.path(root,"fit_manifest.csv"))
diag <- cbh_read_csv(file.path(root,"fit_diagnostics.csv"))
stopifnot(length(components)==28,nrow(manifest)==28,all(diag$converged),all(diag$input_verified))
for(piece in components) stopifnot(identical(piece$model_md5,manifest$md5[match(piece$fit_id,manifest$fit_id)]))
ages <- cbh_config()$age_bands$age_band
years <- c(2005L,2015L,2024L)
series <- c("map_matched","snow_matched","map_full","snow_full")
labels <- c(map_matched="MAP matched, gamma 2",snow_matched="Snow matched, gamma 2",
  map_full="MAP full, gamma 2",snow_full="Snow full, gamma 2",map_previous="Previous MAP full, gamma 1")
inputs <- c(component_path,file.path(root,c("fit_manifest.csv","fit_diagnostics.csv")))
all_age <- all_totals <- list()
for(year in years) {
  olddir <- file.path("results/cbh/age_band_separate_v1",paste0("country_burden_",year))
  age_path <- file.path(olddir,paste0("country_age_attributable_",year,".csv"))
  compare_path <- file.path(olddir,paste0("model_vs_ihme_malaria_",year,".csv"))
  nat_path <- file.path(olddir,paste0("national_pfpr_",year,".csv"))
  inputs <- c(inputs,age_path,compare_path,nat_path)
  old <- cbh_read_csv(age_path);previous <- cbh_read_csv(compare_path);nat <- cbh_read_csv(nat_path)
  cbh_unique(old,c("iso3","age_band"),"Original age inputs")
  stopifnot(nrow(old)==45*7,nrow(previous)==45,all(old$exposure_year==year),
    identical(old$pfpr_pct,nat$pfpr_pct[match(old$iso3,nat$iso3)]))
  keep <- c("iso3","age_band","Location","pfpr_pct","ihme_deaths","implied_person_years","ihme_rate_per100000",
    "age_mapping","status","map_population_coverage_within_raster","map_coverage_below_95pct","population_weight_year")
  bands <- old[keep]
  for(s in series) {
    res <- do.call(rbind,lapply(ages,function(a) {
      z <- bands[bands$age_band==a,]
      piece <- components[[paste0(s,"_age_",match(a,ages))]]
      use <- is.finite(z$pfpr_pct)
      nd <- data.frame(pfpr_pct=z$pfpr_pct[use]);zero <- nd;zero$pfpr_pct <- 0
      L <- PredictMat(piece$smooth,zero)-PredictMat(piece$smooth,nd)
      variance <- rowSums((L%*%piece$covariance)*L)
      stopifnot(all(is.finite(variance)),min(variance)>=-1e-9)
      z$log_hr_zero_vs_current <- z$log_hr_se <- NA_real_
      z$log_hr_zero_vs_current[use] <- drop(L%*%piece$coef)
      z$log_hr_se[use] <- sqrt(pmax(variance,0))
      z$hr_zero_vs_current <- exp(z$log_hr_zero_vs_current)
      z$hr_lower_95 <- exp(z$log_hr_zero_vs_current-1.96*z$log_hr_se)
      z$hr_upper_95 <- exp(z$log_hr_zero_vs_current+1.96*z$log_hr_se)
      z$attributable_fraction <- 1-z$hr_zero_vs_current
      z$af_lower_95 <- 1-z$hr_upper_95;z$af_upper_95 <- 1-z$hr_lower_95
      for(unit in c("deaths","rate_per100000")) {
        base <- z[[paste0("ihme_",unit)]]
        z[[paste0("counterfactual_",unit)]] <- base*z$hr_zero_vs_current
        z[[paste0("attributable_",unit)]] <- base*z$attributable_fraction
        z[[paste0("attributable_",unit,"_lower_95")]] <- base*z$af_lower_95
        z[[paste0("attributable_",unit,"_upper_95")]] <- base*z$af_upper_95
        stopifnot(max(abs(base-z[[paste0("counterfactual_",unit)]]-z[[paste0("attributable_",unit)]]),na.rm=TRUE)<1e-7)
      }
      z$zero_below_observed_support <- piece$support[1]>0
      z$current_pfpr_outside_central95 <- z$pfpr_pct<piece$support[2] | z$pfpr_pct>piece$support[3]
      z$current_pfpr_outside_observed_support <- z$pfpr_pct<piece$support[1] | z$pfpr_pct>piece$support[4]
      z$negative_attributable_estimate <- z$attributable_fraction<0
      z$series <- s;z$year <- year;z$exposure_source <- "MAP"
      z$temporal_transport_from_pre2016_fit <- s!="map_full" && year>2015
      z
    }))
    all_age[[paste(s,year)]] <- res
    total <- do.call(rbind,lapply(split(res,res$iso3),function(z) {
      j <- match(z$iso3[1],previous$iso3)
      stopifnot(nrow(z)==7,abs(sum(z$ihme_deaths)-previous$ihme_under5_deaths[j])<1e-7)
      data.frame(iso3=z$iso3[1],country=z$Location[1],year=year,series=s,pfpr_pct=z$pfpr_pct[1],
        exposure_source="MAP",status=z$status[1],ihme_under5_deaths=sum(z$ihme_deaths),
        attributable_under5_deaths=sum(z$attributable_deaths),counterfactual_under5_deaths=sum(z$counterfactual_deaths),
        attributable_fraction=sum(z$attributable_deaths)/sum(z$ihme_deaths),
        ihme_malaria_deaths=previous$ihme_malaria_deaths[j],
        previous_map_gamma1_deaths=previous$attributable_under5_deaths[j],
        map_population_coverage=z$map_population_coverage_within_raster[1],
        any_current_pfpr_outside_observed_support=any(z$current_pfpr_outside_observed_support),
        any_zero_below_observed_support=any(z$zero_below_observed_support),
        temporal_transport_from_pre2016_fit=any(z$temporal_transport_from_pre2016_fit))
    }))
    all_totals[[paste(s,year)]] <- total
  }
}
r <- do.call(rbind,all_age);totals <- do.call(rbind,all_totals)
cbh_unique(r,c("series","year","iso3","age_band"),"New country age estimates")
cbh_unique(totals,c("series","year","iso3"),"New country totals")
stopifnot(nrow(r)==4*3*45*7,nrow(totals)==4*3*45,
  all(table(totals$series[is.finite(totals$attributable_under5_deaths)],totals$year[is.finite(totals$attributable_under5_deaths)])==42))
cbh_atomic_csv(r,file.path(out,"country_age_estimates.csv"))
cbh_atomic_csv(totals,file.path(out,"country_totals.csv"))
wide <- reshape(totals[c("iso3","country","year","series","attributable_under5_deaths")],
  idvar=c("iso3","country","year"),timevar="series",direction="wide")
names(wide) <- sub("attributable_under5_deaths[.]","",names(wide))
meta <- totals[totals$series=="map_full",]
j <- match(cbh_key(wide,c("iso3","year")),cbh_key(meta,c("iso3","year")))
wide$pfpr_pct <- meta$pfpr_pct[j]
wide$previous_map_gamma1 <- meta$previous_map_gamma1_deaths[j]
wide$ihme_malaria <- meta$ihme_malaria_deaths[j]
wide$matched_snow_minus_map <- wide$snow_matched-wide$map_matched
wide$matched_snow_vs_map_percent <- 100*(wide$snow_matched/wide$map_matched-1)
wide$full_snow_minus_map <- wide$snow_full-wide$map_full
wide$full_snow_vs_map_percent <- 100*(wide$snow_full/wide$map_full-1)
wide$map_gamma2_minus_previous <- wide$map_full-wide$previous_map_gamma1
wide$map_gamma2_vs_previous_percent <- 100*(wide$map_full/wide$previous_map_gamma1-1)
wide <- wide[order(wide$year,wide$iso3),]
cbh_atomic_csv(wide,file.path(out,"country_comparisons.csv"))
matched <- wide[is.finite(wide$map_full) & is.finite(wide$snow_matched),]
sums <- do.call(rbind,lapply(years,function(y) {
  z <- matched[matched$year==y,]
  data.frame(year=y,countries=nrow(z),previous_map_gamma1=sum(z$previous_map_gamma1),
    map_full_gamma2=sum(z$map_full),snow_full_gamma2=sum(z$snow_full),
    map_matched_gamma2=sum(z$map_matched),snow_matched_gamma2=sum(z$snow_matched),ihme_malaria=sum(z$ihme_malaria),
    matched_snow_vs_map_percent=100*(sum(z$snow_matched)/sum(z$map_matched)-1),
    full_snow_vs_map_percent=100*(sum(z$snow_full)/sum(z$map_full)-1),
    map_gamma2_vs_previous_percent=100*(sum(z$map_full)/sum(z$previous_map_gamma1)-1))
}))
cbh_atomic_csv(sums,file.path(out,"year_summary.csv"))
# Retain negative values; no truncation and no unsupported total confidence interval.
scatter <- function(z,x,y,title,subtitle) {
  special <- z[z$iso3 %in% c("COD","NGA","AGO","UGA","TZA"),]
  ggplot(z,aes(.data[[x]],.data[[y]]))+
    geom_abline(slope=1,intercept=0,colour="grey65",linetype=2)+geom_point(colour="#236B80",alpha=.8,size=2)+
    geom_text(data=special,aes(label=iso3),hjust=-.1,vjust=-.5,size=3,check_overlap=TRUE)+
    facet_wrap(~year,nrow=1)+coord_equal()+
    scale_x_continuous(labels=scales::label_number(scale_cut=scales::cut_short_scale()),expand=expansion(mult=c(.04,.15)))+
    scale_y_continuous(labels=scales::label_number(scale_cut=scales::cut_short_scale()),expand=expansion(mult=c(.04,.15)))+
    labs(title=title,subtitle=subtitle,x="MAP-fitted malaria-attributable deaths",y="Snow-fitted malaria-attributable deaths",
      caption="Same annual national MAP prevalence and IHME all-cause inputs for both curves. Dashed line: equality. Country totals are point estimates.")+
    theme_minimal(base_size=11)+theme(plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0))
}
ps <- list(
  matched_country_comparison=scatter(matched,"map_matched","snow_matched","MAP versus Snow mortality estimates: gamma = 2","Both fitted on identical records through 2015; both evaluated at annual national MAP prevalence"),
  full_sample_country_comparison=scatter(matched,"map_full","snow_full","MAP versus Snow mortality estimates: available fitting samples","MAP uses the full fitting period; Snow is fitted through 2015; national exposure is MAP in both"))
for(nm in names(ps)) ggsave(file.path(out,paste0(nm,".png")),ps[[nm]],width=13,height=5,dpi=180,device=ragg::agg_png,bg="white")
long <- rbind(data.frame(matched[c("country","iso3","year")],deaths=matched$map_full,source="MAP full, gamma 2"),
  data.frame(matched[c("country","iso3","year")],deaths=matched$snow_full,source="Snow full, gamma 2"),
  data.frame(matched[c("country","iso3","year")],deaths=matched$ihme_malaria,source="IHME malaria"))
long$country <- factor(long$country,levels=rev(sort(unique(long$country))))
p <- ggplot(long,aes(deaths,country,colour=source,shape=source))+
  geom_point(position=position_dodge(width=.5),size=1.9)+facet_wrap(~year,nrow=1)+
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=100),breaks=c(0,1000,10000,100000),
    labels=scales::label_number(scale_cut=scales::cut_short_scale()))+
  scale_colour_manual(values=c("MAP full, gamma 2"="#215E91","Snow full, gamma 2"="#C26B28","IHME malaria"="#59636B"))+
  labs(title="Country malaria-attributable deaths: MAP versus Snow, gamma = 2",subtitle="Available fitting samples; annual national MAP prevalence used for both models",
    x="Deaths (pseudo-log scale)",y=NULL,colour=NULL,shape=NULL,
    caption="42 countries per year with MAP coverage. IHME malaria is a cause-specific comparator. Model estimates are signed all-cause reductions under zero PfPR.")+
  theme_minimal(base_size=10)+theme(legend.position="bottom",panel.grid.minor=element_blank(),plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0))
ggsave(file.path(out,"country_deaths_all_years.png"),p,width=14,height=12,dpi=180,device=ragg::agg_png,bg="white")
fmt <- function(x) format(round(x),big.mark=",",trim=TRUE)
missing <- unique(wide$country[!is.finite(wide$map_full)])
drc <- wide[wide$iso3=="COD",]
md <- c("# MAP versus Snow: mortality comparison at gamma=2","",
  "Both sets of fitted curves are evaluated at the SAME national population-weighted MAP PfPR in 2005, 2015 and 2024. This holds the scenario exposure fixed and compares fitted relationships. Snow has no 2024 prevalence estimate: the Snow-fitted 2024 result is a transport scenario using MAP exposure, not a source-specific Snow estimate. Applying a Snow-fitted relationship to MAP exposure also assumes the two prevalence scales can be used interchangeably for this scenario.","",
  "The observed IHME all-cause rates and death counts remain fixed. For country c and age g, HR = exp[f_g(0)-f_g(P_c)], counterfactual mortality = IHME mortality × HR, and malaria-attributable mortality = IHME mortality × (1−HR). This is a calibrated attributable-mortality calculation, not an independent prediction of baseline all-cause mortality.","",
  "## Full available fitting samples","",
  "MAP full uses the previous primary model's full fitting sample at gamma=2; Snow full uses its sample through 2015. Their difference includes both exposure source and fitting sample/period. The previous MAP gamma=1 totals are included to show the smoothing update.","",
  "| Year | Countries | Previous MAP gamma 1 | MAP full gamma 2 | Snow full gamma 2 | Snow vs MAP | IHME malaria |",
  "|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(sums)),function(i)sprintf("| %d | %d | %s | %s | %s | %+.1f%% | %s |",sums$year[i],sums$countries[i],fmt(sums$previous_map_gamma1[i]),fmt(sums$map_full_gamma2[i]),fmt(sums$snow_full_gamma2[i]),sums$full_snow_vs_map_percent[i],fmt(sums$ihme_malaria[i])),""),"",
  "## Same-record exposure-source comparison","",
  "MAP and Snow are fitted on identical child-band records through 2015, with matching outcomes, confounders and offsets. This is the cleaner comparison of fitted exposure relationships. Source-specific original PfPR knots are retained; all other basis choices and gamma=2 are the same.","",
  "| Year | MAP matched gamma 2 | Snow matched gamma 2 | Snow vs MAP |","|---|---:|---:|---:|",
  vapply(seq_len(nrow(sums)),function(i)sprintf("| %d | %s | %s | %+.1f%% |",sums$year[i],fmt(sums$map_matched_gamma2[i]),fmt(sums$snow_matched_gamma2[i]),sums$matched_snow_vs_map_percent[i]),""),"",
  "Snow gives lower totals even after matching the fitting records, so the difference is not solely explained by the shorter fitting period. In the matched comparison, Snow estimates are lower in 40 of 42 countries in 2005 and all 42 in 2015 and 2024. With full available samples, Snow is lower in all 42 countries in each year. Increasing full-period MAP gamma from 1 to 2 raises the aggregate estimates by approximately 9–10%.","",
  "![Matched country comparison](matched_country_comparison.png)","",
  "## DRC example","",
  "| Year | MAP full | Snow full | MAP matched | Snow matched |","|---|---:|---:|---:|---:|",
  vapply(seq_len(nrow(drc)),function(i)sprintf("| %d | %s | %s | %s | %s |",drc$year[i],fmt(drc$map_full[i]),fmt(drc$snow_full[i]),fmt(drc$map_matched[i]),fmt(drc$snow_matched[i])),""),"",
  "## Interpretation and checks","",
  "The national PfPR values and IHME inputs are reused from the reviewed earlier calculations and checked against their saved tables. Population weighting remains GPW 2020 density × cell area with fractional country overlap. The IHME 2–4 rate is shared across the three model bands, and deaths/person-time are split equally. Every new age-band calculation satisfies observed = counterfactual + attributable deaths/rates. Matched model frames and compact spline predictions are verified in the fitting stage.","",
  paste("All 45 countries remain in the CSVs; estimates are missing for",paste(missing,collapse=", "),"because national MAP exposure is unavailable. Summaries use the same 42 countries in all series. Missing values are not replaced by zero; negative attributable estimates are retained."),"",
  "Age-band intervals use within-model covariance, conditional on fitted smoothing parameters, fixed HIV imputation, prevalence and IHME baselines. Cross-age and cross-model sampling covariance is unavailable, so country-total intervals and significance tests of between-model differences are not constructed. Source, survey-design, residual within-child, temporal-transport and IHME allocation uncertainty are omitted. Zero-PfPR and current-exposure support flags accompany every age-band estimate; incomplete geographic coverage remains relevant. IHME malaria deaths and the modeled reduction in all-cause deaths are different estimands, and export release alignment remains unverified.","",
  "- [Every country's comparison and changes](country_comparisons.csv) · [Country totals and counterfactuals](country_totals.csv)",
  "- [Age-band attributable rates, counts, intervals and support flags](country_age_estimates.csv)",
  "- [All-country figure](country_deaths_all_years.png) · [Full-sample scatter](full_sample_country_comparison.png) · [Year totals](year_summary.csv)",
  "- [Spline comparison and model specification](../REPORT.md)",
  "- Reproduce after fitting: Rscript R_cbh/snow/10_compare_burden_gamma2.R --exposure=map.")
writeLines(md,file.path(out,"REPORT.md"))
inputs <- c(inputs,"R_cbh/snow/10_compare_burden_gamma2.R")
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"provenance.csv"))
message("MAP/Snow country burden comparison complete")

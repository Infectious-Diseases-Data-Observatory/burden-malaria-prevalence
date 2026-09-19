#!/usr/bin/env Rscript
# Figure 4 (combined): A = 2024 country comparison with IHME; B = Nigerian
# states versus IHME, 2024; C = annual under-five malaria mortality, 2004–2024.
# The main figure expresses every panel as deaths per 1,000 under-five
# child-years; the death-count version is a supplementary figure.
# Plotting only: reads saved aggregates written by the burden/reporting stages.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
source("R_cbh/reporting/log_axes.R")
library(ggplot2)
library(patchwork)
root <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional"))$out
out <- file.path(root,"burden_comparison")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
model_label <- cbh_paper_model_label()
read <- function(p) cbh_read_csv(file.path(root,p))

# Display limits for the supplementary count version: both axes cut at these
# values. The lower country limit of 100 keeps 36 of 42 countries in view.
country_limits <- c(100,2.5e5)
state_limits <- c(200,3e4)
comma <- scales::label_comma()
paper <- theme_minimal(base_size=14)+theme(panel.grid.minor=element_blank(),
  text=element_text(size=16),axis.title=element_text(size=16),axis.text=element_text(size=13),
  legend.title=element_blank(),legend.text=element_text(size=13),legend.position="bottom",
  plot.tag=element_text(size=22,face="bold"),plot.margin=margin(8,16,8,8))

## ---- A: countries, 2024 -------------------------------------------------------
totals <- read("burden/country_totals.csv")
cty <- totals[totals$year==2024 & totals$series=="map_full",]
cty <- cty[is.finite(cty$ihme_malaria_deaths) & cty$ihme_malaria_deaths>0 &
  is.finite(cty$attributable_under5_deaths) & cty$attributable_under5_deaths>0,]
stopifnot(nrow(cty)>0,!anyDuplicated(cty$iso3))
outside <- function(d,x,y,lim) d[d[[x]]<lim[1] | d[[x]]>lim[2] | d[[y]]<lim[1] | d[[y]]>lim[2],]
clipped_a <- outside(cty,"ihme_malaria_deaths","attributable_under5_deaths",country_limits)$iso3
clipped <- data.frame(panel=rep("A",length(clipped_a)),unit=clipped_a)
pA <- ggplot(cty,aes(ihme_malaria_deaths,attributable_under5_deaths))+
  geom_abline(slope=1,intercept=0,colour="grey65",linetype=2)+
  geom_point(colour="#215E91",size=2.6,alpha=.85)+
  ggrepel::geom_text_repel(data=cty[cty$iso3 %in% c("COD","NGA","AGO","UGA","TZA") |
      abs(log(cty$attributable_under5_deaths/cty$ihme_malaria_deaths))>=.5 |
      abs(cty$attributable_under5_deaths-cty$ihme_malaria_deaths)>=5000,],
    aes(label=iso3),size=3.8,seed=20260915,max.overlaps=Inf,min.segment.length=0,segment.colour="grey60")+
  coord_equal(xlim=country_limits,ylim=country_limits,expand=FALSE)+
  scale_x_log10(breaks=c(1e2,1e3,1e4,1e5),labels=comma,minor_breaks=cbh_log10_minor_breaks)+
  scale_y_log10(breaks=c(1e2,1e3,1e4,1e5),labels=comma,minor_breaks=cbh_log10_minor_breaks)+
  labs(x="IHME malaria deaths before age 5, 2024",y=paste0(model_label," deaths\nbefore age 5, 2024"))+
  paper+cbh_log10_grid_theme()

## ---- B: Nigerian states, 2024 --------------------------------------------------
x <- read("nigeria_states/state_totals_2024.csv")
stopifnot(nrow(x)==37,all(x$year==2024),all(x$series=="map_full"),
  all(x$attributable_under5_deaths>0),all(x$ihme_malaria_deaths>0))
zones <- list(
  `North Central`=c("Benue","FCT (Abuja)","Kogi","Kwara","Nasarawa","Niger","Plateau"),
  `North East`=c("Adamawa","Bauchi","Borno","Gombe","Taraba","Yobe"),
  `North West`=c("Jigawa","Kaduna","Kano","Katsina","Kebbi","Sokoto","Zamfara"),
  `South East`=c("Abia","Anambra","Ebonyi","Enugu","Imo"),
  `South South`=c("Akwa Ibom","Bayelsa","Cross River","Delta","Edo","Rivers"),
  `South West`=c("Ekiti","Lagos","Ogun","Ondo","Osun","Oyo"))
zone_of <- setNames(rep(names(zones),lengths(zones)),unlist(zones))
stopifnot(setequal(names(zone_of),x$state))
x$zone <- factor(zone_of[x$state],levels=names(zones))
x$label <- sub(" [(]Abuja[)]$","",x$state)
clipped_b <- outside(x,"ihme_malaria_deaths","attributable_under5_deaths",state_limits)$state
clipped <- rbind(clipped,data.frame(panel=rep("B",length(clipped_b)),unit=clipped_b))
pB <- ggplot(x,aes(ihme_malaria_deaths,attributable_under5_deaths,colour=zone))+
  geom_abline(slope=1,intercept=0,colour="grey65",linetype=2)+
  geom_point(size=2.6)+
  ggrepel::geom_text_repel(aes(label=label),size=3.3,seed=20260917,max.overlaps=Inf,
    min.segment.length=0,segment.colour="grey70",show.legend=FALSE)+
  coord_equal(xlim=state_limits,ylim=state_limits,expand=FALSE)+
  scale_x_log10(breaks=c(1e3,1e4),labels=comma,minor_breaks=cbh_log10_minor_breaks)+
  scale_y_log10(breaks=c(1e3,1e4),labels=comma,minor_breaks=cbh_log10_minor_breaks)+
  scale_colour_manual(values=c(`North Central`="#E69F00",`North East`="#D55E00",`North West`="#CC79A7",
    `South East`="#009E73",`South South`="#0072B2",`South West`="#56B4E9"))+
  guides(colour=guide_legend(nrow=2,override.aes=list(size=3.5)))+
  labs(x="IHME malaria deaths before age 5, 2024",y=paste0(model_label," deaths\nbefore age 5, 2024"))+
  paper+cbh_log10_grid_theme()

## ---- C: annual mortality, 2004–2024 --------------------------------------------
long <- read("annual_comparison/figure5_data.csv")
sources <- c(model_label,"IHME","UN IGME")
stopifnot(setequal(unique(long$source),sources),all(long$countries==42),
  identical(sort(unique(long$year)),2004:2024))
long$source <- factor(long$source,levels=sources)
pC <- ggplot(long,aes(year,rate_per100000,colour=source,linetype=source))+
  geom_line(linewidth=1.15)+
  scale_colour_manual(values=setNames(c("#16747C","#253746","#CC6A30"),sources))+
  scale_linetype_manual(values=setNames(c("solid","longdash","dotdash"),sources))+
  scale_x_continuous(breaks=c(2004,2009,2014,2019,2024),limits=c(2004,2024),
    expand=expansion(mult=c(.015,.025)))+
  scale_y_continuous(limits=c(0,NA),labels=comma,expand=expansion(mult=c(0,.05)))+
  labs(x="Year",y="Under-five malaria deaths\nper 100,000 child-years")+
  paper+theme(legend.key.width=grid::unit(1.2,"cm"))

## ---- assemble -------------------------------------------------------------------
for(p in list(pA,pB,pC)) stopifnot(is.null(p$labels$title),is.null(p$labels$subtitle))
for(p in list(pA,pB)) stopifnot(identical(p$scales$get_scales("x")$trans$name,"log-10"),
  identical(p$scales$get_scales("y")$trans$name,"log-10"))
fig <- (pA|pB)/pC+plot_layout(heights=c(1,.7))+plot_annotation(tag_levels="A")
ggsave(file.path(out,"sfig_burden_comparison_counts.png"),fig,width=14,height=12.5,dpi=240,
  device=ragg::agg_png,bg="white")
cbh_atomic_csv(clipped,file.path(out,"points_outside_display_limits.csv"))

## ---- main Figure 4: every panel in deaths per 1,000 under-five child-years -------
# Same three comparisons expressed as mortality rates, so that panels A and B are not
# dominated by population size. Linear axes with identical limits and an equality line.
ann <- read("annual_comparison/country_estimates_2004_2024.csv")
cr <- ann[ann$year==2024 & ann$complete & is.finite(ann$model_rate_per100000) &
  is.finite(ann$ihme_malaria_rate_per100000),]
stopifnot(nrow(cr)==42,setequal(cr$iso3,cty$iso3),!anyDuplicated(cr$iso3))
cr$model_rate <- cr$model_rate_per100000/100; cr$ihme_rate <- cr$ihme_malaria_rate_per100000/100
x$model_rate <- x$attributable_rate_per100000/100; x$ihme_rate <- x$ihme_malaria_rate_per100000/100
long$rate_per1000 <- long$rate_per100000/100
rate_limits <- function(...) c(0,ceiling(max(...)))
lim_a <- rate_limits(cr$model_rate,cr$ihme_rate); lim_b <- rate_limits(x$model_rate,x$ihme_rate)
rA <- ggplot(cr,aes(ihme_rate,model_rate))+
  geom_abline(slope=1,intercept=0,colour="grey65",linetype=2)+
  geom_point(colour="#215E91",size=2.6,alpha=.85)+
  # Label the large-burden countries and every country whose two rates differ by
  # at least a factor of 1.65 or by at least 1 death per 1,000 child-years.
  ggrepel::geom_text_repel(data=cr[cr$iso3 %in% c("COD","NGA","AGO","UGA","TZA") |
      abs(log(cr$model_rate/cr$ihme_rate))>=.5 | abs(cr$model_rate-cr$ihme_rate)>=1,],
    aes(label=iso3),size=3.8,seed=20260915,max.overlaps=Inf,min.segment.length=0,segment.colour="grey60",
    box.padding=.35)+
  coord_equal(xlim=lim_a,ylim=lim_a,expand=FALSE)+
  labs(x="IHME malaria deaths per 1,000 child-years\nbefore age 5, 2024",
    y=paste0(model_label," deaths per 1,000\nchild-years before age 5, 2024"))+paper
rB <- ggplot(x,aes(ihme_rate,model_rate,colour=zone))+
  geom_abline(slope=1,intercept=0,colour="grey65",linetype=2)+
  geom_point(size=2.6)+
  ggrepel::geom_text_repel(aes(label=label),size=3.3,seed=20260917,max.overlaps=Inf,
    min.segment.length=0,segment.colour="grey70",show.legend=FALSE)+
  coord_equal(xlim=lim_b,ylim=lim_b,expand=FALSE)+
  scale_colour_manual(values=c(`North Central`="#E69F00",`North East`="#D55E00",`North West`="#CC79A7",
    `South East`="#009E73",`South South`="#0072B2",`South West`="#56B4E9"))+
  guides(colour=guide_legend(nrow=2,override.aes=list(size=3.5)))+
  labs(x="IHME malaria deaths per 1,000 child-years\nbefore age 5, 2024",
    y=paste0(model_label," deaths per 1,000\nchild-years before age 5, 2024"))+paper
rC <- ggplot(long,aes(year,rate_per1000,colour=source,linetype=source))+
  geom_line(linewidth=1.15)+
  scale_colour_manual(values=setNames(c("#16747C","#253746","#CC6A30"),sources))+
  scale_linetype_manual(values=setNames(c("solid","longdash","dotdash"),sources))+
  scale_x_continuous(breaks=c(2004,2009,2014,2019,2024),limits=c(2004,2024),
    expand=expansion(mult=c(.015,.025)))+
  scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.05)))+
  labs(x="Year",y="Under-five malaria deaths\nper 1,000 child-years")+
  paper+theme(legend.key.width=grid::unit(1.2,"cm"))
for(p in list(rA,rB,rC)) stopifnot(is.null(p$labels$title),is.null(p$labels$subtitle))
rates <- (rA|rB)/rC+plot_layout(heights=c(1,.7))+plot_annotation(tag_levels="A")
ggsave(file.path(out,"fig4_burden_comparison.png"),rates,width=14,height=12.5,dpi=240,
  device=ragg::agg_png,bg="white")
if(file.exists(file.path(out,"fig4_burden_comparison_rates.png"))) file.remove(file.path(out,"fig4_burden_comparison_rates.png"))
cbh_atomic_csv(cr[order(cr$iso3),c("iso3","country","under5_person_years","model_deaths","ihme_malaria_deaths","model_rate","ihme_rate")],
  file.path(out,"country_rates_2024.csv"))

## ---- caption and provenance ------------------------------------------------------
recon <- read("nigeria_states/national_reconciliation_2024.csv")
fmt <- function(z) format(round(z),big.mark=",",trim=TRUE,scientific=FALSE)
definition <- paste0("In all panels, ",model_label," deaths equal IHME all-cause deaths in each of seven age bands multiplied by 1 − exp[f_g(0) − f_g(P)], where P is population-weighted MAP PfPR[2–10] for the country or state and year and f_g are the seven fitted age-band effects; IHME early and late neonatal deaths share the <1-month effect and the IHME 2–4 year group is divided equally across the three oldest bands.")
caveats <- "All values are point estimates; joint uncertainty across countries and age bands is not available and marginal intervals are not summed. Model-attributable all-cause reductions and cause-specific malaria deaths are different estimands, and the model shares IHME all-cause inputs, so agreement is not independent validation."
caption <- paste(
  paste0("Malaria mortality before age 5 from the ",model_label," compared with IHME and UN IGME cause-specific malaria estimates, expressed as deaths per 1,000 under-five child-years."),
  sprintf("(A) National rates for 2024 in the %d countries with both estimates available; each point is a country and the dashed line denotes equality. Each country's model and IHME deaths are divided by the same under-five person-years implied by the IHME all-cause death count and rate.",nrow(cr)),
  sprintf("(B) The same comparison for the 36 Nigerian states and the Federal Capital Territory in 2024, using the state person-years implied by the IHME all-cause inputs; colours group states into Nigeria's six geopolitical zones for orientation only. Summed over states the model gives %s deaths against %s IHME malaria deaths (ratio %.2f).",
    fmt(recon$state_sum_attributable_deaths),fmt(recon$ihme_malaria_deaths_state_sum),recon$model_over_ihme_malaria_national),
  "(C) Annual under-five malaria mortality, 2004–2024, pooled across the same 42 countries covered by the national burden analysis. For each source the rate is the sum of national under-five malaria deaths divided by the sum of under-five person-years implied by the IHME all-cause death counts and rates, so the three series share one denominator. The UN IGME series is the CA-CODE 2026 release.",
  "Axes in A and B are linear with identical limits. Within each country or state the model and comparator share a denominator, so their ratio is the same as for death counts; the supplementary count version of this figure shows the same comparisons as numbers of deaths.",
  definition,caveats)
writeLines(c("# Figure 4 caption","",caption),file.path(out,"CAPTION.md"))
count_caption <- paste(
  paste0("Count version of Figure 4: malaria-attributable deaths before age 5 from the ",model_label," compared with IHME and UN IGME cause-specific malaria death estimates, as numbers of deaths."),
  sprintf("(A) National estimates for 2024 in the %d countries with both estimates available; each point is a country and the dashed line denotes equality. Both axes use a base-10 logarithmic scale from %s to %s deaths; %d countries with fewer than %s deaths on either axis lie outside the displayed range and are retained in the source table.",
    nrow(cty),fmt(country_limits[1]),fmt(country_limits[2]),sum(clipped$panel=="A"),fmt(country_limits[1])),
  sprintf("(B) The same comparison for the 36 Nigerian states and the Federal Capital Territory in 2024, with axes from %s to %s deaths; colours group states into Nigeria's six geopolitical zones for orientation only.",
    fmt(state_limits[1]),fmt(state_limits[2])),
  "(C) Annual under-five malaria mortality, 2004–2024, per 100,000 child-years, pooled across the same 42 countries on the common IHME-implied person-year denominator.",
  definition,caveats)
writeLines(c("# Supplementary figure caption: count version of Figure 4","",count_caption),file.path(out,"CAPTION_counts.md"))
if(file.exists(file.path(out,"CAPTION_rates.md"))) file.remove(file.path(out,"CAPTION_rates.md"))
paths <- c(file.path(root,c("burden/country_totals.csv","nigeria_states/state_totals_2024.csv",
  "nigeria_states/national_reconciliation_2024.csv","annual_comparison/figure5_data.csv",
  "annual_comparison/country_estimates_2004_2024.csv")),
  "R_cbh/reporting/11_burden_comparison_figure.R","R_cbh/reporting/labels.R","R_cbh/reporting/log_axes.R")
cbh_atomic_csv(data.frame(file=paths,md5=vapply(paths,cbh_file_hash,"")),file.path(out,"figure_provenance.csv"))
message("Combined Figure 4 saved: ",file.path(out,"fig4_burden_comparison.png"))

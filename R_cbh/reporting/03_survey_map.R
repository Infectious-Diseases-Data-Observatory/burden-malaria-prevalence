#!/usr/bin/env Rscript
# Figure 1: current primary survey coverage, adapted from archived 28_survey_map.R.
# Uses aggregate selected-survey counts and existing local outlines; no downloads.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
library(sf)
library(ggplot2)
library(patchwork)
st <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional_mics")); root <- st$out; mics <- isTRUE(st$mics)
out <- file.path(root,"survey_map")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
paths <- c(coverage=file.path(root,"survey_coverage.csv"),
  registry=cbh_config()$registry,selection=file.path(root,"study_flow/survey_selection.csv"))
x <- cbh_read_csv(paths[["coverage"]]);r <- cbh_read_csv(paths[["registry"]]);sel <- cbh_read_csv(paths[["selection"]])
if(mics) {
  # MICS surveys with a complete birth history join the registry; the 46 without one are an
  # up-front ineligibility shown in the study flow, not on this timeline.
  paths <- c(paths,mics_registry="data/derived_mics/survey_registry_mics.csv")
  cols <- c("svkey","iso3","year","SurveyType","boundary_file")
  r <- rbind(r[cols],cbh_read_csv(paths[["mics_registry"]])[cols])
}
cbh_unique(x,"survey","Primary survey coverage")
cbh_unique(r,"svkey","Survey registry")
sample <- cbh_read_csv(file.path(root,"primary_sample.csv"))
stopifnot(nrow(x)==sample$surveys,length(unique(x$country))==sample$countries,
  setequal(x$survey,sel$survey[sel$complete_case_rows>0]))
j <- match(x$survey,r$svkey)
stopifnot(!anyNA(j),all(x$country==r$iso3[j]))
x$year <- r$year[j];x$type <- r$SurveyType[j]
stopifnot(!anyNA(x),all(x$type %in% c("DHS","MIS","MICS")))
summary <- do.call(rbind,lapply(split(x,x$country),function(d)
  data.frame(country=d$country[1],surveys=nrow(d),first_year=min(d$year),last_year=max(d$year),survey_regions=sum(d$regions))))

# Every registry survey, with its fate in the primary selection. Surveys absent from
# the selection ledger were never processed (no MAP geography); processed surveys with
# no complete-case rows lost every record to covariate availability (or, if it arose,
# to missing MAP), as recorded per survey in the study-flow ledger.
all <- data.frame(survey=r$svkey,country=r$iso3,year=r$year,type=r$SurveyType)
k <- match(all$survey,sel$survey)
all$eligible_rows <- sel$eligible_rows[k];all$pfpr_available_rows <- sel$pfpr_available_rows[k]
all$complete_case_rows <- sel$complete_case_rows[k]
all$regions <- x$regions[match(all$survey,x$survey)]
all$status <- ifelse(all$survey %in% x$survey,"Included in primary analysis",
  ifelse(is.na(k) | all$pfpr_available_rows==0,"Excluded: malaria free, no MAP prevalence (Lesotho)",
    "Excluded: required covariates unavailable"))
stopifnot(all(!is.na(all$regions)==(all$status=="Included in primary analysis")),
  sum(all$status=="Included in primary analysis")==nrow(x),!anyNA(all$year))
status_levels <- c("Included in primary analysis","Excluded: malaria free, no MAP prevalence (Lesotho)","Excluded: required covariates unavailable")
all$status <- factor(all$status,levels=status_levels)

# Most recent available boundary in each country with any registry survey.
shape_paths <- character()
outline <- lapply(unique(all$country),function(iso) {
  rr <- r[r$iso3==iso,];rr <- rr[order(-rr$year),]
  for(i in seq_len(nrow(rr))) {
    path <- rr$boundary_file[i]
    if(is.na(path) || !nzchar(path) || !file.exists(path)) next
    shape <- sf::st_make_valid(readRDS(path))
    if(!nrow(shape)) next
    geom <- sf::st_union(sf::st_geometry(sf::st_transform(shape,4326)))
    shape_paths <<- c(shape_paths,path)
    return(sf::st_sf(country=iso,geometry=geom))
  }
  if(iso %in% summary$country) stop("Missing local country outline for ",iso)
  message("No local outline for excluded country ",iso,"; it appears on the timeline only")
  NULL
})
outline <- do.call(rbind,outline[!vapply(outline,is.null,logical(1))])
outline <- merge(outline,summary,by="country",all.x=TRUE)
outline$surveys[is.na(outline$surveys)] <- 0L
centres <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(outline))))
outline$lon <- centres[,1];outline$lat <- centres[,2]
outline$fill <- ifelse(outline$surveys>0,outline$surveys,NA_integer_)
map <- ggplot(outline)+geom_sf(aes(fill=fill),colour="white",linewidth=.3)+
  ggrepel::geom_text_repel(data=sf::st_drop_geometry(outline),aes(lon,lat,label=country),
    size=4.5,seed=20260915,max.overlaps=Inf,min.segment.length=0,box.padding=.15,
    segment.colour="grey60",segment.size=.25)+
  scale_fill_gradient(low="#DCEAF3",high="#12557A",name="Included surveys",breaks=scales::breaks_pretty(4),
    na.value="grey90")+
  guides(fill=guide_colourbar(barwidth=grid::unit(1.9,"in"),barheight=grid::unit(.2,"in")))+
  theme_minimal(base_size=20)+theme(axis.title=element_blank(),axis.text=element_blank(),
    panel.grid=element_blank(),legend.position="bottom",
    legend.title=element_text(size=18),legend.text=element_text(size=16))
# Countries without an outline (none expected) sort to the top of the timeline.
lat_order <- outline$country[order(outline$lat)]
all$country_label <- factor(all$country,levels=c(lat_order,setdiff(unique(all$country),lat_order)))
disp_levels <- if(mics) c("Included: DHS","Included: MICS",status_levels[2:3]) else status_levels
all$display <- if(mics) ifelse(all$status==status_levels[1],ifelse(all$type=="MICS","Included: MICS","Included: DHS"),as.character(all$status)) else as.character(all$status)
all$display <- factor(all$display,levels=disp_levels)
shape_values <- if(mics) c(16,17,5,4) else c(16,5,4)
colour_values <- if(mics) c("#12557A","#C07A12","#B3261E","#B3261E") else c("#12557A","#B3261E","#B3261E")
included <- all[all$status==status_levels[1],];excluded <- all[all$status!=status_levels[1],]
timeline <- ggplot(all,aes(year,country_label))+
  geom_line(aes(group=country_label),colour="grey85",linewidth=.4)+
  geom_point(data=included,aes(size=regions,shape=display,colour=display),alpha=.9)+
  geom_point(data=excluded,aes(shape=display,colour=display),size=2.6,stroke=1.1)+
  scale_size_continuous(range=c(1.2,4),name="Regions")+
  scale_shape_manual(values=setNames(shape_values,disp_levels),name=NULL,drop=FALSE)+
  scale_colour_manual(values=setNames(colour_values,disp_levels),name=NULL,drop=FALSE)+
  scale_x_continuous(breaks=seq(2000,2025,5))+
  labs(x="Survey year",y=NULL)+theme_minimal(base_size=20)+
  guides(size=guide_legend(order=1,nrow=1),
    shape=guide_legend(order=2,ncol=1,override.aes=list(size=3.2)),colour=guide_legend(order=2,ncol=1))+
  theme(panel.grid.minor=element_blank(),legend.position="bottom",axis.text=element_text(size=16),
    axis.title=element_text(size=20),legend.title=element_text(size=18),legend.text=element_text(size=15),
    legend.box="vertical",legend.spacing.y=grid::unit(0,"pt"),plot.margin=margin(8,26,8,8))
p <- map+timeline+plot_layout(widths=c(1,1.15))
ggsave(file.path(out,"survey_map_and_timing.png"),p,width=13,height=10.5,dpi=300,device=ragg::agg_png,bg="white")
x$country_label <- NULL;all$country_label <- NULL;all$display <- NULL
cbh_atomic_csv(x,file.path(out,"survey_coverage.csv"))
cbh_atomic_csv(all,file.path(out,"survey_timeline_all.csv"))
cbh_atomic_csv(summary,file.path(out,"country_summary.csv"))
inputs <- unique(c(unname(paths),file.path(root,"primary_sample.csv"),shape_paths,"R_cbh/reporting/03_survey_map.R","R_cbh/primary/settings.R"))
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"provenance.csv"))
n_map <- sum(all$status==status_levels[2]);n_cov <- sum(all$status==status_levels[3])
mis_dropped <- sum(all$type=="MIS" & all$status!=status_levels[1]);mis_total <- sum(all$type=="MIS")
if(mics) {
  nd <- sum(all$type!="MICS"); nm <- sum(all$type=="MICS"); inc_d <- sum(included$type=="DHS"); inc_m <- sum(included$type=="MICS")
  writeLines(c("# Figure 1 caption","",
    sprintf("Geographic coverage and timing of the %s DHS and MIS surveys and %s MICS surveys with complete birth histories (%s countries). Filled circles are the %s DHS surveys and filled triangles the %s MICS surveys that contribute to the primary analysis (%s countries in total), with point area indicating the number of survey regions contributing analysis records; country shading indicates the number of included surveys, and grey countries have none. Open diamonds mark %s surveys in Lesotho, which is malaria free, so MAP prevalence cannot be assigned; crosses mark %s surveys excluded because a required covariate was unavailable for every region after the declared substitutions (all %s MIS surveys, and MICS surveys without anthropometry, without a child HIV incidence series, or whose band-entry years precede a national series). The 46 MICS surveys without a complete birth history are not shown.",
      nd,nm,length(unique(all$country)),inc_d,inc_m,nrow(summary),n_map,n_cov,mis_total),"",
    sprintf("Region counts are the union across the seven fitted age-band samples; the figure contains %s included survey–region pairs, which are not counts of geographically distinct regions across survey years. Country outlines are dissolved from the most recent available boundary file in each country (DHS files, or the analysis-region polygons built for MICS). Survey year describes fieldwork, not the calendar year assigned to each child's band entry. Per-survey exclusion reasons are in survey_timeline_all.csv.",sum(x$regions))),
    file.path(out,"CAPTION.md"))
} else writeLines(c("# Figure 1 caption","",
  sprintf("Geographic coverage and timing of the %s DHS and MIS surveys with complete birth histories in the survey registry (%s countries). Filled circles are the %s surveys in %s countries that contribute to the primary analysis, with point area indicating the number of survey regions contributing analysis records; country shading indicates the number of included surveys, and grey countries have none. Open diamonds mark %s surveys excluded because MAP prevalence could not be assigned to their regions; crosses mark %s surveys excluded because a required covariate was unavailable for every region after the declared substitutions (%s of the %s MIS surveys fall in this group, so the analysed sample is DHS only).",
    nrow(all),length(unique(all$country)),nrow(x),nrow(summary),n_map,n_cov,mis_dropped,mis_total),"",
  sprintf("Region counts are the union across the seven fitted age-band samples; the figure contains %s included survey–region pairs, which are not counts of geographically distinct regions across survey years. Country outlines are dissolved from the most recent available DHS boundary file in each country. Survey year describes fieldwork, not the calendar year assigned to each child's band entry. Per-survey exclusion reasons are in survey_timeline_all.csv.",sum(x$regions))),
  file.path(out,"CAPTION.md"))
message("Generated survey map: ",nrow(x)," included of ",nrow(all)," registry surveys; ",n_map," excluded for MAP, ",n_cov," for covariates")

#!/usr/bin/env Rscript
# Figure 1: current primary survey coverage, adapted from archived 28_survey_map.R.
# Uses aggregate selected-survey counts and existing local outlines; no downloads.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
library(sf)
library(ggplot2)
library(patchwork)
root <- cbh_primary_settings()$out
out <- file.path(root,"survey_map")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
paths <- c(coverage=file.path(root,"survey_coverage.csv"),
  registry=cbh_config()$registry,selection=file.path(root,"study_flow/survey_selection.csv"))
x <- cbh_read_csv(paths[["coverage"]]);r <- cbh_read_csv(paths[["registry"]]);sel <- cbh_read_csv(paths[["selection"]])
cbh_unique(x,"survey","Primary survey coverage")
cbh_unique(r,"svkey","Survey registry")
stopifnot(nrow(x)==105,length(unique(x$country))==34,setequal(x$survey,sel$survey[sel$complete_case_rows>0]))
j <- match(x$survey,r$svkey)
stopifnot(!anyNA(j),all(x$country==r$iso3[j]))
x$year <- r$year[j];x$type <- r$SurveyType[j]
stopifnot(!anyNA(x),all(x$type %in% c("DHS","MIS")))
summary <- do.call(rbind,lapply(split(x,x$country),function(d)
  data.frame(country=d$country[1],surveys=nrow(d),first_year=min(d$year),last_year=max(d$year),survey_regions=sum(d$regions))))
# Most recent available boundary among contributing surveys in each country.
shape_paths <- character()
outline <- lapply(summary$country,function(iso) {
  rr <- r[r$iso3==iso & r$svkey %in% x$survey,];rr <- rr[order(-rr$year),]
  for(i in seq_len(nrow(rr))) {
    path <- rr$boundary_file[i]
    if(!file.exists(path)) next
    shape <- sf::st_make_valid(readRDS(path))
    if(!nrow(shape)) next
    geom <- sf::st_union(sf::st_geometry(sf::st_transform(shape,4326)))
    shape_paths <<- c(shape_paths,path)
    return(sf::st_sf(country=iso,geometry=geom))
  }
  stop("Missing local country outline for ",iso)
})
outline <- do.call(rbind,outline)
outline <- merge(outline,summary,by="country")
centres <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(outline))))
outline$lon <- centres[,1];outline$lat <- centres[,2]
map <- ggplot(outline)+geom_sf(aes(fill=surveys),colour="white",linewidth=.3)+
  ggrepel::geom_text_repel(data=sf::st_drop_geometry(outline),aes(lon,lat,label=country),
    size=4.5,seed=20260915,max.overlaps=Inf,min.segment.length=0,box.padding=.15,
    segment.colour="grey60",segment.size=.25)+
  scale_fill_gradient(low="#DCEAF3",high="#12557A",name="Surveys",breaks=scales::breaks_pretty(4))+
  guides(fill=guide_colourbar(barwidth=grid::unit(1.9,"in"),barheight=grid::unit(.2,"in")))+
  theme_minimal(base_size=20)+theme(axis.title=element_blank(),axis.text=element_blank(),
    panel.grid=element_blank(),legend.position="bottom",
    legend.title=element_text(size=18),legend.text=element_text(size=16))
x$country_label <- factor(x$country,levels=outline$country[order(outline$lat)])
timeline <- ggplot(x,aes(year,country_label))+
  geom_line(aes(group=country_label),colour="grey85",linewidth=.4)+
  geom_point(aes(size=regions,shape=type),colour="#12557A",alpha=.85)+
  scale_size_continuous(range=c(1.2,4),name="Regions")+
  scale_shape_manual(values=c(DHS=16,MIS=17),name="Type")+
  scale_x_continuous(breaks=seq(2000,2025,5))+
  labs(x="Survey year",y=NULL)+theme_minimal(base_size=20)+
  guides(size=guide_legend(order=1,nrow=1),shape=guide_legend(order=2,nrow=1))+
  theme(panel.grid.minor=element_blank(),legend.position="bottom",axis.text=element_text(size=16),
    axis.title=element_text(size=20),legend.title=element_text(size=18),legend.text=element_text(size=16),
    legend.box="vertical",legend.spacing.y=grid::unit(0,"pt"),plot.margin=margin(8,26,8,8))
p <- map+timeline+plot_layout(widths=c(1,1.15))
ggsave(file.path(out,"survey_map_and_timing.png"),p,width=13,height=10,dpi=300,device=ragg::agg_png,bg="white")
x$country_label <- NULL
cbh_atomic_csv(x,file.path(out,"survey_coverage.csv"))
cbh_atomic_csv(summary,file.path(out,"country_summary.csv"))
inputs <- unique(c(unname(paths),shape_paths,"R_cbh/reporting/03_survey_map.R","R_cbh/primary/settings.R"))
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"provenance.csv"))
writeLines(c("# Figure 1 caption","",
  sprintf("Geographic coverage and timing of the %s surveys in %s countries contributing to the primary MAP analysis. Country shading indicates the number of included surveys. Timeline point area indicates the number of regions contributing analysis records in each survey; symbols distinguish DHS and MIS. Region counts are the union across the seven fitted age-band samples. The figure contains %s survey–region pairs; these are not counts of geographically distinct regions across survey years.",nrow(x),nrow(summary),sum(x$regions)),"",
  "Country outlines are dissolved from the most recent available DHS boundary files among the included surveys. Survey year describes fieldwork, not the calendar year assigned to each child's band entry. The historical figure showed 120 surveys/36 countries; this version shows the 105-survey/34-country primary complete-case sample."),file.path(out,"CAPTION.md"))
message("Generated current primary survey map: ",nrow(x)," surveys, ",nrow(summary)," countries")

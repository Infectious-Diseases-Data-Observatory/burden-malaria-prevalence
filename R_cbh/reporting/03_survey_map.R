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
  # MICS surveys with a complete birth history join the registry; those without one (counted
  # from the MICS inventory) are an up-front ineligibility shown in the study flow, not on this timeline.
  paths <- c(paths,mics_registry="data/derived_mics/survey_registry_mics.csv",mics_inventory="results/mics_inventory/survey_inventory.csv")
  cols <- c("svkey","iso3","year","SurveyType","boundary_file")
  r <- rbind(r[cols],cbh_read_csv(paths[["mics_registry"]])[cols])
  mics_inv <- cbh_read_csv(paths[["mics_inventory"]])
  stopifnot(is.logical(mics_inv$has_bh),!anyNA(mics_inv$has_bh),sum(mics_inv$has_bh)==sum(r$SurveyType=="MICS"))
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
# Mothers (complete birth histories) contributing at least one child-band record to the primary
# analysis, per included survey: children in the prepared sample matched to their shard's mother id.
analysed <- unique(as.character(readRDS(st$data)$data$child_id))
shard_dirs <- c(cbh_config()$output_dir,if(mics) st$mics_output_dir)
mothers <- do.call(rbind,lapply(shard_dirs,function(dir) {
  m <- readRDS(file.path(dir,"manifest.rds"))$manifest; m <- m[m$survey %in% x$survey,]
  do.call(rbind,lapply(seq_len(nrow(m)),function(i) { o <- readRDS(file.path(dir,m$file[i])); stopifnot(identical(o$signature,m$signature[i]))
    d <- o$data[as.character(o$data$child_id) %in% analysed,,drop=FALSE]
    data.frame(survey=m$survey[i],children=length(unique(d$child_id)),mothers=length(unique(d$mother_id))) })) }))
cbh_unique(mothers,"survey","Mothers by survey")
stopifnot(setequal(mothers$survey,x$survey),sum(mothers$children)==sample$distinct_children,all(mothers$mothers>0))
x$mothers <- mothers$mothers[match(x$survey,mothers$survey)]
summary <- do.call(rbind,lapply(split(x,x$country),function(d)
  data.frame(country=d$country[1],surveys=nrow(d),first_year=min(d$year),last_year=max(d$year),survey_regions=sum(d$regions),mothers=sum(d$mothers))))

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
  scale_fill_gradient(low="#DCEAF3",high="#12557A",name="Included surveys",breaks=scales::breaks_pretty(4),
    na.value="grey90")+
  guides(fill=guide_colourbar(barwidth=grid::unit(1.9,"in"),barheight=grid::unit(.2,"in")))+
  theme_minimal(base_size=20)+theme(axis.title=element_blank(),axis.text=element_blank(),
    panel.grid=element_blank(),legend.position="bottom",
    legend.title=element_text(size=18),legend.text=element_text(size=16))
# Timeline rows grouped by UN M49 sub-region (Middle Africa labelled Central Africa), countries in
# alphabetical order of ISO3 code within each region, with the number of mothers on the right.
m49 <- list(`West Africa`=c("BEN","BFA","CIV","CPV","GHA","GIN","GMB","GNB","LBR","MLI","MRT","NER","NGA","SEN","SLE","TGO"),
  `Central Africa`=c("AGO","CAF","CMR","COD","COG","GAB","GNQ","STP","TCD"),
  `East Africa`=c("BDI","COM","DJI","ERI","ETH","KEN","MDG","MOZ","MUS","MWI","RWA","SOM","SSD","TZA","UGA","ZMB","ZWE"),
  `Southern Africa`=c("BWA","LSO","NAM","SWZ","ZAF"))
region_of <- setNames(rep(names(m49),lengths(m49)),unlist(m49))
stopifnot(all(all$country %in% names(region_of)))
all$region_group <- factor(region_of[all$country],levels=names(m49))
# Full country names (as in Figure 4), with DRC, CAF and RC abbreviated (decided 24 September 2026);
# alphabetical by name within each region.
country_names <- c(AGO="Angola",BDI="Burundi",BEN="Benin",BFA="Burkina Faso",BWA="Botswana",CAF="CAF",CIV="Côte d'Ivoire",
  CMR="Cameroon",COD="DRC",COG="RC",COM="Comoros",CPV="Cabo Verde",DJI="Djibouti",ERI="Eritrea",ETH="Ethiopia",GAB="Gabon",
  GHA="Ghana",GIN="Guinea",GMB="The Gambia",GNB="Guinea-Bissau",GNQ="Equatorial Guinea",KEN="Kenya",LBR="Liberia",LSO="Lesotho",
  MDG="Madagascar",MLI="Mali",MOZ="Mozambique",MRT="Mauritania",MUS="Mauritius",MWI="Malawi",NAM="Namibia",NER="Niger",
  NGA="Nigeria",RWA="Rwanda",SEN="Senegal",SLE="Sierra Leone",SOM="Somalia",SSD="South Sudan",STP="São Tomé and Príncipe",
  SWZ="Eswatini",TCD="Chad",TGO="Togo",TZA="Tanzania",UGA="Uganda",ZAF="South Africa",ZMB="Zambia",ZWE="Zimbabwe")
stopifnot(all(all$country %in% names(country_names)))
all$country_name <- unname(country_names[all$country])
row_order <- unique(all[order(all$region_group,all$country_name),"country_name"])
all$country_label <- factor(all$country_name,levels=rev(row_order))
levels(all$region_group) <- sub(" ","\n",levels(all$region_group))
n_rows <- unique(all[c("country","country_label","region_group")])
n_rows$mothers <- summary$mothers[match(n_rows$country,summary$country)]; n_rows$mothers[is.na(n_rows$mothers)] <- 0L
n_rows$label <- sprintf("(n=%s)",format(n_rows$mothers,big.mark=",",trim=TRUE))
# Shape gives the programme (DHS and MIS circles, MICS triangles); fill gives inclusion (filled
# included, open excluded, whatever the exclusion reason; reasons stay in survey_timeline_all.csv).
all$programme <- factor(ifelse(all$type=="MICS","MICS","DHS/MIS"),levels=c("DHS/MIS","MICS"))
all$inclusion <- factor(ifelse(all$status==status_levels[1],"Included","Excluded"),levels=c("Included","Excluded"))
ink <- "#12557A"
# Where a DHS or MIS survey shares the country-year, the MICS symbol is lifted slightly so
# that neither symbol hides the other.
all$nudge <- mics & all$type=="MICS" & ave(all$type!="MICS",all$country,all$year,FUN=any)
included <- all[all$status==status_levels[1],];excluded <- all[all$status!=status_levels[1],]
survey_points <- function(d,mapping,...) list(geom_point(data=d[!d$nudge,],mapping=mapping,...),
  geom_point(data=d[d$nudge,],mapping=mapping,position=position_nudge(y=.3),...))
timeline <- ggplot(all,aes(year,country_label))+
  geom_line(aes(group=country_label),colour="grey85",linewidth=.4)+
  geom_text(data=n_rows,aes(x=2025.9,y=country_label,label=label),hjust=0,size=4.3,colour="grey25")+
  survey_points(all,aes(shape=programme,fill=inclusion),colour=ink,size=3.4,stroke=1)+
  scale_shape_manual(values=c(`DHS/MIS`=21,MICS=24),name=NULL,drop=FALSE)+
  scale_fill_manual(values=c(Included=ink,Excluded="white"),name=NULL,drop=FALSE)+
  scale_x_continuous(breaks=seq(2000,2025,5))+
  coord_cartesian(xlim=c(1999.3,2025.2),clip="off")+
  facet_grid(region_group~.,scales="free_y",space="free_y",switch="y")+
  labs(x="Survey year",y=NULL)+theme_minimal(base_size=20)+
  guides(shape=guide_legend(order=1,nrow=1,override.aes=list(fill=ink,size=4)),
    fill=guide_legend(order=2,nrow=1,override.aes=list(shape=21,size=4)))+
  theme(panel.grid.minor=element_blank(),legend.position="bottom",axis.text.x=element_text(size=17),axis.text.y.left=element_text(size=15),
    axis.title=element_text(size=20),legend.text=element_text(size=18),
    legend.box="horizontal",legend.spacing.x=grid::unit(14,"pt"),strip.placement="outside",
    strip.text.y.left=element_text(size=15,face="bold",angle=0,hjust=1),panel.spacing.y=grid::unit(10,"pt"),plot.margin=margin(8,120,8,8))
p <- map+timeline+plot_layout(widths=c(1,1.45))
ggsave(file.path(out,"survey_map_and_timing.png"),p,width=16,height=12,dpi=300,device=ragg::agg_png,bg="white")
nudged <- unique(paste(all$country[all$nudge],all$year[all$nudge]))
x$country_label <- NULL;all$country_label <- NULL;all$country_name <- NULL;all$programme <- NULL;all$inclusion <- NULL;all$nudge <- NULL;all$region_group <- NULL
cbh_atomic_csv(x,file.path(out,"survey_coverage.csv"))
cbh_atomic_csv(all,file.path(out,"survey_timeline_all.csv"))
cbh_atomic_csv(summary,file.path(out,"country_summary.csv"))
cbh_atomic_csv(mothers,file.path(out,"mothers_by_survey.csv"))
inputs <- unique(c(unname(paths),file.path(root,"primary_sample.csv"),st$data,file.path(shard_dirs,"manifest.rds"),shape_paths,
  "R_cbh/reporting/03_survey_map.R","R_cbh/primary/settings.R"))
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"provenance.csv"))
n_map <- sum(all$status==status_levels[2]);n_cov <- sum(all$status==status_levels[3])
mis_dropped <- sum(all$type=="MIS" & all$status!=status_levels[1]);mis_total <- sum(all$type=="MIS")
if(mics) {
  nd <- sum(all$type!="MICS"); nm <- sum(all$type=="MICS"); inc_d <- sum(included$type=="DHS"); inc_m <- sum(included$type=="MICS")
  cov_type <- all$type[all$status==status_levels[3]]
  cov_mis <- sum(cov_type=="MIS"); cov_mics <- sum(cov_type=="MICS")
  cov_text <- paste0(sprintf("%s DHS, ",sum(cov_type=="DHS")),
    if(cov_mis==mis_total) sprintf("all %s MIS",mis_total) else sprintf("%s of the %s MIS",cov_mis,mis_total),
    sprintf(" and %s MICS surveys",cov_mics),
    if(cov_mics) "; the MICS surveys lack anthropometry or a child HIV incidence series, or their band-entry years precede a national series" else "")
  nudge_text <- if(length(nudged)) sprintf(" Where a MICS survey shares a country and year with a DHS or MIS survey (%s), the MICS symbol is drawn slightly above the line.",paste(nudged,collapse=", ")) else ""
  writeLines(c("# Figure 1 caption","",
    sprintf("Geographic coverage and timing of the %s DHS and MIS surveys and %s MICS surveys with complete birth histories (%s countries). Circles are DHS and MIS surveys and triangles MICS surveys; filled symbols are the %s DHS and %s MICS surveys that contribute to the primary analysis (%s countries), and open symbols the %s excluded surveys: %s in Lesotho, which is malaria free, so MAP prevalence cannot be assigned, and %s for which a required covariate was unavailable for every region after the declared substitutions (%s). Country shading indicates the number of included surveys; grey countries have none. Timeline rows are grouped by UN sub-region (Central Africa: UN Middle Africa), and n is the number of mothers (complete birth histories) contributing at least one child to the primary analysis, summed over the country's included surveys (%s in total).%s The %s MICS surveys without a complete birth history are not shown.",
      nd,nm,length(unique(all$country)),inc_d,inc_m,nrow(summary),n_map+n_cov,n_map,n_cov,cov_text,format(sum(mothers$mothers),big.mark=","),nudge_text,sum(!mics_inv$has_bh)),"",
    "Country outlines are dissolved from the most recent available boundary file in each country (DHS files, or the analysis-region polygons built for MICS). Survey year describes fieldwork, not the calendar year assigned to each child's band entry. Per-survey exclusion reasons are in survey_timeline_all.csv."),
    file.path(out,"CAPTION.md"))
} else writeLines(c("# Figure 1 caption","",
  sprintf("Geographic coverage and timing of the %s DHS and MIS surveys with complete birth histories in the survey registry (%s countries). Filled circles are the %s surveys in %s countries that contribute to the primary analysis and open circles the %s excluded surveys: %s because MAP prevalence could not be assigned to their regions and %s because a required covariate was unavailable for every region after the declared substitutions (%s of the %s MIS surveys fall in this group, so the analysed sample is DHS only). Country shading indicates the number of included surveys; grey countries have none.",
    nrow(all),length(unique(all$country)),nrow(x),nrow(summary),n_map+n_cov,n_map,n_cov,mis_dropped,mis_total),"",
  "Country outlines are dissolved from the most recent available DHS boundary file in each country. Survey year describes fieldwork, not the calendar year assigned to each child's band entry. Per-survey exclusion reasons are in survey_timeline_all.csv."),
  file.path(out,"CAPTION.md"))
message("Generated survey map: ",nrow(x)," included of ",nrow(all)," registry surveys; ",n_map," excluded for MAP, ",n_cov," for covariates")

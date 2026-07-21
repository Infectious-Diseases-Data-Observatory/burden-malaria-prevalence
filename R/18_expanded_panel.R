# =============================================================================
# 18_expanded_panel.R — EXPANDED DHS/MIS panel, 2000-2024, MAP-for-all prevalence.
# Downloads every SSA births-recode (mortality) survey from 2000+, computes
# region-level all-cause child mortality (DHS.rates::chmort by v024), assigns
# each region a MAP PfPR2-10 (pop-weighted at region + survey year), adds
# GDP/DTP3/urban, and refits the primary Method-2 model (linear nb-GAM, country
# random slope, prevalence x year interaction). Extends the DHS panel back to
# ~2000, so the temporal slope is estimated from data rather than extrapolated.
#
# Resumable: per-survey rows cached in data/panel_cache/, boundaries in
# data/dhs_boundaries/ (with retry), MAP rasters in data/map_annual/.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(rdhs); library(DHS.rates); library(malariaAtlas); library(terra); library(sf); library(mgcv)})
options(rappdir_permission = TRUE)
BDIR <- file.path(DATA,"dhs_boundaries"); MDIR <- file.path(DATA,"map_annual"); PDIR <- file.path(DATA,"panel_cache")
for (dd in c(BDIR,MDIR,PDIR)) dir.create(dd, showWarnings = FALSE)
ext <- matrix(c(-18,-35,52,38), nrow=2)
den0 <- terra::rast(GPW_TIF)
gdp <- read.csv(file.path(DATA,"wb_gdp_pc.csv")); dtp <- read.csv(file.path(DATA,"wb_dtp3.csv"))
nrst <- function(panel,iso,yr,col){ s<-panel[panel$iso3==iso & is.finite(panel[[col]]),]; if(!nrow(s)) NA_real_ else s[[col]][which.min(abs(s$year-yr))] }

## ---- 1. enumerate BR recodes, SSA, 2000+ ------------------------------------
ds <- dhs_datasets(fileFormat="FL")
br <- ds[ds$FileType=="Births Recode",]
br$iso3 <- countrycode::countrycode(br$CountryName,"country.name","iso3c",warn=FALSE)
br$reg  <- countrycode::countrycode(br$iso3,"iso3c","region",warn=FALSE)
br$yr   <- as.integer(br$SurveyYear)
sv <- br[!is.na(br$reg) & br$reg=="Sub-Saharan Africa" & br$yr>=2000 & br$SurveyType %in% c("DHS","MIS","AIS"),]
sv$svkey <- paste0(substr(sv$FileName,1,2), substr(sv$FileName,5,8))
cat(sprintf("[1] BR surveys SSA 2000+: %d (%d countries)\n", nrow(sv), length(unique(sv$iso3))))

## ---- 2. download BR recodes (cached in rdhs) --------------------------------
cat("[2] downloading BR recodes (cached) ...\n")
paths <- tryCatch(get_datasets(dataset_filenames=sv$FileName, download_option="rds", reformat=TRUE, clear_cache=FALSE),
                  error=function(e){cat("  get_datasets error:",conditionMessage(e),"\n"); list()})
cat(sprintf("    downloaded/cached: %d of %d\n", sum(!vapply(paths,is.null,TRUE)), nrow(sv)))

## ---- helpers: MAP raster (cache) + boundary (cache+retry) -------------------
get_map <- function(yr){ f<-file.path(MDIR,sprintf("pfpr2_10_%d.tif",yr))
  if(!file.exists(f)){ dl<-file.path(tempdir(),paste0("m",yr)); dir.create(dl,showWarnings=FALSE)
    invisible(tryCatch(getRaster(dataset_id="Malaria__202508_Global_Pf_Parasite_Rate",year=yr,extent=ext,file_path=dl),error=function(e)NULL))
    tif<-list.files(dl,pattern="\\.tiff?$",full.names=TRUE); if(!length(tif))return(NULL)
    r<-terra::rast(tif); r<-r[[grep("_1$",names(r))[1]]]; terra::writeRaster(r,f,overwrite=TRUE) }
  terra::rast(f) }
get_bnd <- function(sid){ f<-file.path(BDIR,paste0(sid,".rds"))
  if(file.exists(f)) return(readRDS(f))
  for(a in 1:3){ b<-tryCatch(download_boundaries(surveyId=sid,method="sf")[[1]],error=function(e)NULL)
    if(!is.null(b)){ saveRDS(b,f); return(b) } }
  NULL }

## ---- 3. per-survey region rows (cached) -------------------------------------
process <- function(i){
  s<-sv[i,]; cf<-file.path(PDIR,paste0(s$svkey,".rds")); if(file.exists(cf)) return(readRDS(cf))
  p<-paths[[sub("\\..*$","",s$FileName)]]; if(is.null(p)) return(NULL)   # get_datasets keys by name w/o extension
  br<-tryCatch(readRDS(p),error=function(e)NULL); if(is.null(br)||!"v024"%in%names(br)) return(NULL)
  mo<-tryCatch(mort_by_region(br,"v024"),error=function(e)NULL); if(is.null(mo)||!nrow(mo)) return(NULL)
  # % urban by region (v025 weighted by v005)
  w<-as.numeric(br$v005)/1e6; reg<-as.character(br$v024); urb<-tolower(as.character(br$v025))=="urban"
  ok<-!is.na(reg); pu<-100*tapply(w[ok]*urb[ok],reg[ok],sum)/tapply(w[ok],reg[ok],sum)
  murb<-data.frame(regkey=rkey(names(pu)), pct_urban=as.numeric(pu))
  bd<-get_bnd(s$SurveyId); if(is.null(bd)||!"DHSREGEN"%in%names(bd)) return(NULL)
  pf<-get_map(s$yr); if(is.null(pf)) return(NULL)
  v<-terra::makeValid(terra::vect(sf::st_make_valid(bd)))
  den<-terra::resample(terra::crop(den0,pf),pf,method="bilinear"); wt<-terra::mask(den,pf)
  an<-terra::extract(pf*wt,v,fun=sum,na.rm=TRUE,ID=FALSE)[[1]]; aw<-terra::extract(wt,v,fun=sum,na.rm=TRUE,ID=FALSE)[[1]]
  bmap<-data.frame(regkey=rkey(bd$DHSREGEN), pfpr2_10=100*an/aw)
  r<-merge(mo,bmap,by="regkey"); r<-merge(r,murb,by="regkey",all.x=TRUE)
  if(!nrow(r)) return(NULL)
  r$svkey<-s$svkey; r$iso3<-s$iso3; r$year<-s$yr; r$pfpr10<-r$pfpr2_10/10
  r$log_gdp<-log(nrst(gdp,s$iso3,s$yr,"gdp_pc")); r$dtp3<-nrst(dtp,s$iso3,s$yr,"dtp3")
  saveRDS(r,cf); r
}
cat("[3] building per-survey region rows ...\n")
rows<-list()
for(i in seq_len(nrow(sv))){ r<-tryCatch(process(i),error=function(e){message(sv$svkey[i],": ",conditionMessage(e));NULL})
  rows[[sv$svkey[i]]]<-r; cat(sprintf("  %-8s %-11s %s\n", sv$svkey[i], sv$SurveyId[i], if(is.null(r))"-" else paste(nrow(r),"reg"))) }
d<-do.call(rbind,rows)
d$year_c<-d$year-round(mean(d$year,na.rm=TRUE))
d<-d[is.finite(d$u5mr)&d$u5mr>5&is.finite(d$m1mo5y)&d$m1mo5y>0&is.finite(d$pfpr2_10),]
write.csv(d,file.path(RESULTS,"component2_region_data_expanded.csv"),row.names=FALSE)
cat(sprintf("[3] EXPANDED PANEL: %d region-years, %d surveys, %d countries, years %d-%d\n",
    nrow(d),length(unique(d$svkey)),length(unique(d$iso3)),min(d$year),max(d$year)))

## ---- 4. refit primary model on the expanded panel ---------------------------
fitd<-d[complete.cases(d[,c("m1mo5y","pfpr10","dtp3","log_gdp","pct_urban","year_c","iso3","svkey")]) &
          is.finite(d$exposure)&d$exposure>0&d$pfpr2_10>=1,]
fitd$deaths<-round(fitd$m1mo5y/1000*fitd$exposure); fitd$country<-factor(fitd$iso3); fitd$svkey<-factor(fitd$svkey)
cat(sprintf("[4] model sample: %d regions, %d countries, years %d-%d\n",nrow(fitd),nlevels(fitd$country),min(fitd$year),max(fitd$year)))
m_lin <- gam(deaths ~ pfpr10 + dtp3 + log_gdp + pct_urban + s(year_c) + s(country,bs="re") + s(country,pfpr10,bs="re") + offset(log(exposure)),
             family=nb(),method="REML",data=fitd)
m_int <- gam(deaths ~ pfpr10 + pfpr10:year_c + dtp3 + log_gdp + pct_urban + s(year_c) + s(country,bs="re") + s(country,pfpr10,bs="re") + offset(log(exposure)),
             family=nb(),method="REML",data=fitd)
b<-summary(m_lin)$p.table["pfpr10",]; bi<-summary(m_int)$p.table
cat(sprintf("[4] MAIN effect: %+.1f%% per +10 PfPR pts (95%% CI %+.1f to %+.1f)\n",(exp(b[1])-1)*100,(exp(b[1]-1.96*b[2])-1)*100,(exp(b[1]+1.96*b[2])-1)*100))
cat(sprintf("[4] prevalence x year interaction: %.5f/yr (p=%.3f)\n",bi["pfpr10:year_c","Estimate"],bi["pfpr10:year_c","Pr(>|z|)"]))
saveRDS(list(m_lin=m_lin,m_int=m_int,center=round(mean(d$year,na.rm=TRUE))),file.path(RESULTS,"expanded_models.rds"))
cat("saved: results/component2_region_data_expanded.csv + expanded_models.rds\n")
cat("[DONE] expanded panel + refit complete\n")

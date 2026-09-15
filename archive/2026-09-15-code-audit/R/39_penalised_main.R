# =============================================================================
# 39_penalised_main.R — regenerate the MAIN-analysis headline figures from the
# RIDGE-PENALISED model (R/33 pn/nn), now the primary specification. The
# confounder block enters as a single ridge-penalised matrix term G; predictions
# set G = 0 (all covariates at their standardised mean). Regenerates:
#   paper_fig1.png (scatter + spline + AF), malaria_spline_shape.png,
#   paper_fig4_temporal.png, country_af_curves_postneonatal.png, paper_fig3_rct.png
# and prints the penalised-model numbers for the text. Copies fig1/3/4 to Overleaf.
# =============================================================================
source("R/00_utils.R"); suppressMessages({library(mgcv); library(ggplot2); library(patchwork); library(malariaAtlas)})
P <- readRDS(file.path(RESULTS, "penalized_models.rds")); m <- P$pn; mn <- P$nn
Gcn <- colnames(m$model$G); addG <- function(nd) { nd$G <- matrix(0, nrow(nd), length(Gcn), dimnames = list(NULL, Gcn)); nd }
EX <- c("s(country)", "s(country,pfpr10)")
nd0 <- function(p10) addG(data.frame(pfpr10 = p10, year_c = 0, exposure = 1, country = levels(m$model$country)[1]))
d <- read.csv(file.path(RESULTS, "component2_region_data_full.csv"), stringsAsFactors = FALSE)
fitd <- d[is.finite(d$exposure) & d$exposure > 0 & d$pfpr2_10 >= 1 & is.finite(d$m1mo5y) & d$m1mo5y > 0, ]
fitd$country <- factor(fitd$iso3)
sterm <- function(mod, p10) { nd <- nd0(p10); as.numeric(predict(mod, nd, type = "terms", terms = "s(pfpr10)")) }
af1 <- function(mod, p) { e1 <- sterm(mod, 0.1); 100 * pmax(1 - exp(-(sterm(mod, p/10) - e1)), 0) }
cat(sprintf("[penalised main] AF vs 1%%: 10%%=%.0f%%  30%%=%.0f%%  50%%=%.0f%%\n", af1(m,10), af1(m,30), af1(m,50)))

## ===== FIG 1 =====
grid <- data.frame(pfpr2_10 = exp(seq(log(1), log(max(fitd$pfpr2_10)), length.out = 220))); grid$pfpr10 <- grid$pfpr2_10/10
ndg <- nd0(grid$pfpr10)
pr <- predict(m, ndg, type = "link", se.fit = TRUE, exclude = EX)
grid$rate <- 1000*exp(pr$fit); grid$lo <- 1000*exp(pr$fit-1.96*pr$se.fit); grid$hi <- 1000*exp(pr$fit+1.96*pr$se.fit)
Xg <- predict(m, ndg, type = "lpmatrix"); X1 <- predict(m, nd0(0.1), type = "lpmatrix")
si <- grep("s\\(pfpr10\\)", colnames(Xg)); dX <- sweep(Xg[, si, drop=FALSE], 2, X1[1, si])
eta <- as.numeric(dX %*% coef(m)[si]); se <- sqrt(rowSums((dX %*% vcov(m)[si, si]) * dX))
grid$af <- 100*pmax(1-exp(-eta),0); grid$aflo <- 100*pmax(1-exp(-(eta-1.96*se)),0); grid$afhi <- 100*pmax(1-exp(-(eta+1.96*se)),0)
GMIN<-element_line(colour="grey92",linewidth=0.3); GMAJ<-element_line(colour="grey85",linewidth=0.4)
XSC<-scale_x_log10(breaks=c(1,2,5,10,20,50,80),expand=expansion(mult=c(0.02,0.03))); XLAB<-expression(italic(Pf)*"PR"[2-10]*" (%), MAP model, log scale")
pA <- ggplot() + geom_point(data=fitd,aes(pfpr2_10,m1mo5y,size=exposure),colour="grey45",alpha=0.35) +
  geom_ribbon(data=grid,aes(pfpr2_10,ymin=lo,ymax=hi),fill="#08519c",alpha=0.18)+geom_line(data=grid,aes(pfpr2_10,rate),colour="#08519c",linewidth=1.1)+
  scale_size_area(max_size=6,name="Births in\nsurvey-region",breaks=c(500,2000,5000),labels=scales::comma)+
  XSC+scale_y_log10(breaks=c(10,20,50,100,200))+coord_cartesian(ylim=c(10,200))+
  labs(x=XLAB,y="All-cause post-neonatal mortality\n(1 month-5 years, per 1,000 live births, log scale)")+
  theme_bw(base_size=12)+theme(panel.grid.minor=GMIN,panel.grid.major=GMAJ,legend.position=c(0.99,0.02),legend.justification=c(1,0),
    legend.background=element_rect(fill=scales::alpha("white",0.7),colour=NA),legend.key.size=unit(4,"mm"),legend.title=element_text(size=8),legend.text=element_text(size=7))
pB <- ggplot(grid,aes(pfpr2_10,af))+geom_ribbon(aes(ymin=aflo,ymax=afhi),fill="#d73027",alpha=0.16)+geom_line(colour="#d73027",linewidth=1.1)+
  geom_hline(yintercept=0,linetype="dotted",colour="grey55")+XSC+
  labs(x=XLAB,y="Share of post-neonatal deaths attributable\nto malaria (%, vs 1% PfPR counterfactual)")+
  theme_bw(base_size=12)+theme(panel.grid.minor=GMIN,panel.grid.major=GMAJ)
ggsave(file.path(RESULTS,"paper_fig1.png"),(pA+pB)+plot_annotation(tag_levels="A")&theme(plot.tag=element_text(face="bold",size=14)),width=13,height=5.8,dpi=320)

## ===== SPLINE (post-neonatal pn vs neonatal nn) =====
sp <- function(mod,lab){ pr<-predict(mod,ndg,type="terms",terms="s(pfpr10)",se.fit=TRUE)
  data.frame(pfpr2_10=grid$pfpr2_10,outcome=lab,fit=as.numeric(pr$fit),lo=as.numeric(pr$fit)-1.96*as.numeric(pr$se.fit),hi=as.numeric(pr$fit)+1.96*as.numeric(pr$se.fit)) }
L <- rbind(sp(m,"Post-neonatal (1mo-5y)"), sp(mn,"Neonatal (<1mo)")); cols<-c("Post-neonatal (1mo-5y)"="#08519c","Neonatal (<1mo)"="#d73027")
ps <- ggplot(L,aes(pfpr2_10,fit,colour=outcome,fill=outcome))+geom_hline(yintercept=0,linetype="dotted",colour="grey55")+
  geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.16,colour=NA)+geom_line(linewidth=1.1)+
  geom_rug(data=data.frame(x=fitd$pfpr2_10),aes(x=x),inherit.aes=FALSE,alpha=0.2,sides="b")+
  scale_colour_manual(values=cols,name=NULL)+scale_fill_manual(values=cols,name=NULL)+
  labs(x=expression(italic(Pf)*"PR"[2-10]~"(%)"),y="s(PfPR) partial effect on mortality (log scale)")+
  theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position="top")
ggsave(file.path(RESULTS,"malaria_spline_shape.png"),ps,width=8,height=5.8,dpi=300)

## ===== FIG 4 temporal (burden under penalised pn) =====
sm <- m$smooth[[which(vapply(m$smooth,function(s)s$label,"")=="s(country,pfpr10)")]]
bsl <- setNames(as.numeric(coef(m)[sm$first.para:sm$last.para]), levels(m$model$country)); f1 <- sterm(m,0.1)
af_ct <- function(iso,pf){ bc<-unname(bsl[iso]);bc[is.na(bc)]<-0; pmax(ifelse(pf>=1,1-exp(-((sterm(m,pf/10)-f1)+bc*(pf-1)/10)),0),0) }
ts<-read.csv(file.path(DATA,"wb_mortality_timeseries.csv")); pcy<-read.csv(file.path(DATA,"pfpr_by_country_year.csv"))
dz<-merge(ts[,c("iso3","year","allcause_1mo5y")],pcy,by=c("iso3","year")); dz$mal<-mapply(af_ct,dz$iso3,dz$pfpr_pct)*dz$allcause_1mo5y
agg<-aggregate(mal~year,dz,sum); o<-function(y)agg$mal[agg$year==y]
U5<-0.75; who<-data.frame(year=2000:2024,point=c(804,815,786,758,751,712,723,704,666,673,650,618,578,556,550,550,546,548,552,545,598,577,573,567,579)*1000)
ih<-tryCatch(read.csv(file.path(RESULTS,"ihme_u5_deaths_ssa_timeseries.csv")),error=function(e)NULL)
OPN<-"Prevalence/all-cause mortality (penalised model)";WH<-"WHO — African-region under-5";IH<-"IHME/GBD — SSA under-5"
pl<-rbind(data.frame(year=agg$year,series=OPN,deaths=agg$mal),data.frame(year=who$year,series=WH,deaths=who$point*U5))
if(!is.null(ih))pl<-rbind(pl,data.frame(year=ih$year,series=IH,deaths=ih$point))
pl<-pl[pl$year<=2024,];lev<-c(OPN,WH,IH);pl$series<-factor(pl$series,levels=lev)
p4<-ggplot(pl,aes(year,deaths/1000,colour=series))+geom_line(aes(linetype=series),linewidth=1)+geom_point(size=1.1)+
  scale_colour_manual(values=setNames(c("#d73027","grey35","#238b45"),lev),name=NULL)+scale_linetype_manual(values=setNames(c("solid","22","44"),lev),name=NULL)+
  scale_x_continuous(breaks=c(2000,2005,2010,2015,2020,2024))+expand_limits(y=0)+labs(x=NULL,y="Malaria-attributable child deaths (thousands/yr)")+
  theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="top",axis.title=element_text(size=13),axis.text=element_text(size=11),legend.text=element_text(size=9))
ggsave(file.path(RESULTS,"paper_fig4_temporal.png"),p4,width=11.5,height=6.4,dpi=300)
cat(sprintf("[penalised main] temporal: %.0fk (2000) -> %.0fk (2024)  %+.0f%%\n",o(2000)/1e3,o(2024)/1e3,100*(o(2024)/o(2000)-1)))

## ===== copy to Overleaf =====
FIG<-path.expand("~/Dropbox/Apps/Overleaf/malaria-burden-reassessment/figures")
for(f in c("fig1_methods.png","fig1_survey_gam.png")) file.copy(file.path(RESULTS,"paper_fig1.png"),file.path(FIG,f),overwrite=TRUE)
file.copy(file.path(RESULTS,"paper_fig4_temporal.png"),file.path(FIG,"fig4_temporal.png"),overwrite=TRUE)
cat("saved: paper_fig1, malaria_spline_shape, paper_fig4_temporal ; copied fig1/fig4 to Overleaf\n")

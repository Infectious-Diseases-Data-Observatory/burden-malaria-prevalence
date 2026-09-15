# =============================================================================
# 29_paper_figures.R — rebuild the Overleaf paper figures on the NEW analysis:
# expanded 2000-2024 panel (936 region-years, MAP-for-all PfPR2-10), region-level
# health-access covariates, and the combined best-fitting model m_comb
# (nonlinear s(pfpr10) + region health + country RE/slope; R/25).
#
#   paper_fig1        -> figures/fig1_methods.png     (scatter + spline; AF vs 1%)
#   paper_fig2_negctl -> figures/fig2_negcontrol.png  (U5 / post-neonatal / neonatal effect)
#   paper_fig3_rct    -> figures/fig3_rct.png         (RCT triangulation; Panel B from m_comb)
#   paper_fig4_temporal -> figures/fig4_temporal.png  (burden under m_comb vs WHO/IHME)
# fig_ihme_share_spline (Method-1 GBD) is unaffected by the DHS panel — re-run R/14.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2); library(patchwork); library(malariaAtlas)})
CM <- readRDS(file.path(RESULTS, "combined_models.rds")); m <- CM$m_comb
EX <- c("s(country)", "s(country,pfpr10)")

d <- read.csv(file.path(RESULTS, "component2_region_data_full.csv"), stringsAsFactors = FALSE)  # imputed panel (main analysis)
d$nnmr <- d$u5mr - d$m1mo5y
need <- c("m1mo5y","nnmr","u5mr","pfpr10","dtp3","dtp3_reg","facility","educ_yrs","wealth_q",
          "log_gdp","pct_urban","year_c","iso3","svkey")
fitd <- d[complete.cases(d[, need]) & is.finite(d$exposure) & d$exposure > 0 &
            d$pfpr2_10 >= 1 & d$m1mo5y > 0 & d$nnmr > 0, ]
fitd$country <- factor(fitd$iso3)
cmeans <- lapply(fitd[, c("dtp3_reg","facility","educ_yrs","wealth_q","log_gdp","pct_urban")], mean)
nd0 <- function(p10) data.frame(c(list(pfpr10 = p10, year_c = 0), cmeans,
                                  exposure = 1, country = levels(fitd$country)[1]))

## ===== FIG 1: scatter + m_comb spline (A); attributable fraction vs 1% (B) ====
grid <- data.frame(pfpr2_10 = exp(seq(log(1), log(max(fitd$pfpr2_10)), length.out = 220)))
grid$pfpr10 <- grid$pfpr2_10 / 10
ndg <- nd0(grid$pfpr10)
pr  <- predict(m, ndg, type = "link", se.fit = TRUE, exclude = EX)
grid$rate <- 1000 * exp(pr$fit); grid$lo <- 1000*exp(pr$fit-1.96*pr$se.fit); grid$hi <- 1000*exp(pr$fit+1.96*pr$se.fit)
# AF vs 1%: smooth-only difference eta(p)-eta(1) with delta-method CI via lpmatrix
Xg <- predict(m, ndg, type = "lpmatrix"); X1 <- predict(m, nd0(0.1), type = "lpmatrix")
si <- grep("s\\(pfpr10\\)", colnames(Xg)); dX <- sweep(Xg[, si, drop=FALSE], 2, X1[1, si])
eta <- as.numeric(dX %*% coef(m)[si]); se <- sqrt(rowSums((dX %*% vcov(m)[si, si]) * dX))
grid$af <- 100*pmax(1-exp(-eta),0); grid$af_lo <- 100*pmax(1-exp(-(eta-1.96*se)),0); grid$af_hi <- 100*pmax(1-exp(-(eta+1.96*se)),0)
GMIN <- element_line(colour="grey92", linewidth=0.3); GMAJ <- element_line(colour="grey85", linewidth=0.4)
XSC <- scale_x_log10(breaks=c(1,2,5,10,20,50,80), expand=expansion(mult=c(0.02,0.03)))
XLAB <- expression(italic(Pf)*"PR"[2-10]*" (%), MAP model, log scale"); BLU <- "#08519c"
pA <- ggplot() +
  geom_point(data=fitd, aes(pfpr2_10, m1mo5y, size=exposure), colour="grey45", alpha=0.35) +
  geom_ribbon(data=grid, aes(pfpr2_10, ymin=lo, ymax=hi), fill=BLU, alpha=0.18) +
  geom_line(data=grid, aes(pfpr2_10, rate), colour=BLU, linewidth=1.1) +
  scale_size_area(max_size=6, name="Births in\nsurvey-region", breaks=c(500,2000,5000), labels=scales::comma) +
  XSC + scale_y_log10(breaks=c(10,20,50,100,200)) + coord_cartesian(ylim=c(10,200)) +
  labs(x=XLAB, y="All-cause post-neonatal mortality\n(1 month-5 years, per 1,000 live births, log scale)") +
  theme_bw(base_size=12) + theme(panel.grid.minor=GMIN, panel.grid.major=GMAJ,
    legend.position=c(0.99,0.02), legend.justification=c(1,0),
    legend.background=element_rect(fill=scales::alpha("white",0.7), colour=NA),
    legend.key.size=unit(4,"mm"), legend.title=element_text(size=8), legend.text=element_text(size=7))
pB <- ggplot(grid, aes(pfpr2_10, af)) +
  geom_ribbon(aes(ymin=af_lo, ymax=af_hi), fill="#d73027", alpha=0.16) +
  geom_line(colour="#d73027", linewidth=1.1) + geom_hline(yintercept=0, linetype="dotted", colour="grey55") +
  XSC + labs(x=XLAB, y="Share of post-neonatal deaths attributable\nto malaria (%, vs 1% PfPR counterfactual)") +
  theme_bw(base_size=12) + theme(panel.grid.minor=GMIN, panel.grid.major=GMAJ)
ggsave(file.path(RESULTS,"paper_fig1.png"), (pA+pB)+plot_annotation(tag_levels="A") &
  theme(plot.tag=element_text(face="bold", size=14)), width=13, height=5.8, dpi=320)
cat("fig1 done: n =", nrow(fitd), "region-years; AF@30% =", round(grid$af[which.min(abs(grid$pfpr2_10-30))],1), "%\n")

## ===== FIG 2: negative control — effect across U5 / post-neonatal / neonatal ==
RE  <- "s(country,bs=\"re\") + s(country,pfpr10,bs=\"re\") + offset(log(exposure))"
RH  <- "pfpr10 + dtp3_reg + facility + educ_yrs + wealth_q + log_gdp + pct_urban + s(year_c)"
efit <- function(rate) { fitd$deaths <- round(fitd[[rate]]/1000*fitd$exposure)
  mm <- gam(as.formula(paste("deaths ~", RH, "+", RE)), family=nb(), method="REML", data=fitd)
  b <- unname(summary(mm)$p.table["pfpr10", ])
  data.frame(pct=(exp(b[1])-1)*100, lo=(exp(b[1]-1.96*b[2])-1)*100, hi=(exp(b[1]+1.96*b[2])-1)*100, p=b[4]) }
nc <- do.call(rbind, Map(function(r,l){x<-efit(r);x$outcome<-l;x},
  c("u5mr","m1mo5y","nnmr"), c("All under-5","1mo-5y (post-neonatal)","Neonatal (control)")))
nc$outcome <- factor(nc$outcome, levels=c("Neonatal (control)","1mo-5y (post-neonatal)","All under-5"))
write.csv(nc, file.path(RESULTS,"paper_fig2_negcontrol.csv"), row.names=FALSE)
pnc <- ggplot(nc, aes(pct, outcome)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey50") +
  geom_linerange(aes(xmin=lo, xmax=hi)) + geom_point(size=3, colour="#08519c") +
  labs(x="% change in mortality per +10 PfPR2-10 points (95% CI)", y=NULL) +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(),
    axis.title=element_text(size=14), axis.text=element_text(size=13))
ggsave(file.path(RESULTS,"paper_fig2_negcontrol.png"), pnc, width=9, height=3.8, dpi=300)
cat(sprintf("fig2 done: U5 %+.1f%%, post-neo %+.1f%%, neonatal %+.1f%% (p=%.3f)\n",
    nc$pct[nc$outcome=="All under-5"], nc$pct[grepl("post",nc$outcome)],
    nc$pct[grepl("Neonatal",nc$outcome)], nc$p[grepl("Neonatal",nc$outcome)]))

## ===== FIG 3: RCT triangulation — Panel A (data), Panel B (predicted from m_comb) =
rct <- data.frame(
  study=c(paste0("D'Alessandro z",1:5),"Habluetzel 1997","Phillips-Howard 2003","Nevill 1996"),
  comparison=c(rep("untreated net",5),"no curtain","no net","no net"),
  prev_treated=c(28.2,22.8,15.9,43.2,71.0,83,54,11.9), prev_control=c(36.7,33.0,25.8,53.3,44.8,91,66,25.1),
  prev_amin=c(rep(1,5),0.5,0.25,0.08), prev_amax=c(rep(5,5),5,3,0.92),
  d_t=c(4,13,18,35,66,NA,NA,NA), cy_t=c(1056,1472,1908,1939,3079,NA,NA,NA),
  d_c=c(10,25,33,32,83,NA,NA,NA), cy_c=c(1008,1374,1717,1676,4632,NA,NA,NA),
  mort_rr_direct=c(rep(NA,5),0.85,0.846,0.70), mort_logse_direct=c(rep(NA,5),0.099,0.044,0.142),
  stringsAsFactors=FALSE)
rct$mort_rr <- ifelse(!is.na(rct$d_t),(rct$d_t/rct$cy_t)/(rct$d_c/rct$cy_c),rct$mort_rr_direct)
rct$mort_logse <- ifelse(!is.na(rct$d_t),sqrt(1/rct$d_t+1/rct$d_c),rct$mort_logse_direct)
std <- function(p,a,b) 100*as.numeric(malariaAtlas::convertPrevalence(p/100,a,b,2,10))
rct$pfpr_t <- mapply(std,rct$prev_treated,rct$prev_amin,rct$prev_amax)
rct$pfpr_c <- mapply(std,rct$prev_control,rct$prev_amin,rct$prev_amax)
rct$dPfPR <- rct$pfpr_c-rct$pfpr_t; rct$obs_red <- (1-rct$mort_rr)*100
rct$obs_lo <- (1-exp(log(rct$mort_rr)+1.96*rct$mort_logse))*100
rct$obs_hi <- (1-exp(log(rct$mort_rr)-1.96*rct$mort_logse))*100
# predicted reduction from m_comb nonlinear response: RR = exp(f(p_t) - f(p_c))
fsm <- function(p) as.numeric(predict(m, nd0(p/10), type="terms", terms="s(pfpr10)"))
rct$pred <- (1-exp(fsm(rct$pfpr_t)-fsm(rct$pfpr_c)))*100
lab <- sub("Phillips-Howard","P-Howard",sub("D'Alessandro ","DA ",sub(" 199.| 2003","",rct$study)))
rct$lab <- lab                                   # hand-tuned label offsets (row order DA z1-5, Hab, P-H, Nevill)
rct$nax_A <- c(0,2.0,-2.0,0,4.5,-1.8,1.8,-2.4); rct$nay_A <- c(8,8,-8,-8,8,8,-8,7)
rct$nax_B <- c(-4.5,4.5,4.5,4.5,4.5,-4.5,4.5,-5); rct$nay_B <- c(0,5,-5,0,0,0,2,0)
pA <- ggplot(rct, aes(dPfPR, obs_red)) +
  geom_hline(yintercept=0, colour="grey80") + geom_vline(xintercept=0, colour="grey80") +
  geom_smooth(method="lm", mapping=aes(weight=1/mort_logse^2), se=TRUE, colour="grey30", fill="grey85", linewidth=0.7) +
  geom_errorbar(aes(ymin=obs_lo, ymax=obs_hi), width=0.7, alpha=0.35) +
  geom_point(aes(colour=comparison), size=3) +
  geom_text(aes(x=dPfPR+nax_A, y=obs_red+nay_A, label=lab), size=3.0) +
  scale_colour_manual(values=c("no net"="#08519c","no curtain"="#41ab5d","untreated net"="#d73027"), name=NULL) +
  labs(x="Prevalence reduction (PfPR2-10 points, control - intervention)", y="Observed U5 mortality reduction (%)") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(), legend.position="top",
    axis.title=element_text(size=13), axis.text=element_text(size=11), legend.text=element_text(size=11))
xr <- range(rct$pred,0)+c(-4,6); yr <- range(rct$obs_lo,rct$obs_hi,0)+c(-6,10)
pB <- ggplot(rct, aes(pred, obs_red)) +
  geom_abline(slope=1, intercept=0, linetype="dotted", colour="grey50") +
  geom_hline(yintercept=0, colour="grey85") + geom_vline(xintercept=0, colour="grey85") +
  geom_errorbar(aes(ymin=obs_lo, ymax=obs_hi), width=0.9, alpha=0.35) +
  geom_point(aes(colour=comparison), size=3) + geom_text(aes(x=pred+nax_B, y=obs_red+nay_B, label=lab), size=3.0, colour="grey25") +
  scale_colour_manual(values=c("no net"="#08519c","no curtain"="#41ab5d","untreated net"="#d73027"), name=NULL) +
  coord_cartesian(xlim=xr, ylim=yr) +
  labs(x="Predicted U5 mortality reduction (%, m_comb)", y="Observed U5 mortality reduction (%)") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(), legend.position="top",
    axis.title=element_text(size=13), axis.text=element_text(size=11), legend.text=element_text(size=11))
ggsave(file.path(RESULTS,"paper_fig3_rct.png"), (pA+pB)+patchwork::plot_layout(widths=c(1.35,1)) +
  patchwork::plot_annotation(tag_levels="A") & theme(plot.tag=element_text(face="bold", size=16)),
  width=14, height=6.4, dpi=300)
cat("fig3 done\n")

## ===== FIG 4: temporal burden under m_comb vs WHO / IHME =====================
sm <- m$smooth[[which(vapply(m$smooth,function(s)s$label,"")=="s(country,pfpr10)")]]
bsl <- setNames(as.numeric(coef(m)[sm$first.para:sm$last.para]), levels(fitd$country))
f1 <- fsm(1)
af_ct <- function(iso, pf) { bc <- unname(bsl[iso]); bc[is.na(bc)] <- 0
  pmax(ifelse(pf>=1, 1-exp(-((fsm(pf)-f1) + bc*(pf-1)/10)), 0), 0) }
ts  <- read.csv(file.path(DATA,"wb_mortality_timeseries.csv"))
pcy <- read.csv(file.path(DATA,"pfpr_by_country_year.csv"))
dz <- merge(ts[,c("iso3","year","allcause_1mo5y")], pcy, by=c("iso3","year"))
dz$mal <- mapply(af_ct, dz$iso3, dz$pfpr_pct) * dz$allcause_1mo5y
agg <- aggregate(mal ~ year, dz, sum); o <- function(y) agg$mal[agg$year==y]
U5 <- 0.75
who <- data.frame(year=2000:2024, point=c(804,815,786,758,751,712,723,704,666,673,650,618,578,556,550,550,546,548,552,545,598,577,573,567,579)*1000)
ih  <- tryCatch(read.csv(file.path(RESULTS,"ihme_u5_deaths_ssa_timeseries.csv")), error=function(e) NULL)
OPN <- "Prevalence/all-cause mortality (combined model)"; WH <- "WHO — African-region under-5"; IH <- "IHME/GBD — SSA under-5"
pl <- rbind(data.frame(year=agg$year, series=OPN, deaths=agg$mal),
            data.frame(year=who$year, series=WH, deaths=who$point*U5))
if(!is.null(ih)) pl <- rbind(pl, data.frame(year=ih$year, series=IH, deaths=ih$point))
pl <- pl[pl$year<=2024,]; lev <- c(OPN,WH,IH); pl$series <- factor(pl$series, levels=lev)
p4 <- ggplot(pl, aes(year, deaths/1000, colour=series)) +
  geom_line(aes(linetype=series), linewidth=1) + geom_point(size=1.1) +
  scale_colour_manual(values=setNames(c("#d73027","grey35","#238b45"),lev), name=NULL) +
  scale_linetype_manual(values=setNames(c("solid","22","44"),lev), name=NULL) +
  scale_x_continuous(breaks=c(2000,2005,2010,2015,2020,2024)) + expand_limits(y=0) +
  labs(x=NULL, y="Malaria-attributable child deaths (thousands/yr)") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(), legend.position="top",
    axis.title=element_text(size=13), axis.text=element_text(size=11), legend.text=element_text(size=9))
ggsave(file.path(RESULTS,"paper_fig4_temporal.png"), p4, width=11.5, height=6.4, dpi=300)
cat(sprintf("fig4 done: %.0fk (2000) -> %.0fk (2024)  %+.0f%%\n", o(2000)/1e3, o(2024)/1e3, 100*(o(2024)/o(2000)-1)))
cat("ALL FIGURES REBUILT\n")

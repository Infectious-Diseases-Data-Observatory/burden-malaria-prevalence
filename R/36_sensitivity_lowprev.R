# =============================================================================
# 36_sensitivity_lowprev.R — SENSITIVITY: include region-years with PfPR2-10 < 1%.
# The primary analysis restricts to PfPR2-10 >= 1% (C2_MIN_PFPR). Here we refit the
# combined model (m_comb spec: nonlinear s(pfpr10) + region-health covars + country
# RE/slope) on ALL region-years and compare the malaria effect and spline shape to
# the >=1% primary, for post-neonatal and neonatal. AF still referenced to 1%.
# =============================================================================
source("R/00_utils.R"); suppressMessages({library(mgcv); library(ggplot2)})
d <- read.csv(file.path(RESULTS, "component2_region_data_expanded_health.csv"), stringsAsFactors = FALSE)
d$nnmr <- d$u5mr - d$m1mo5y
need <- c("pfpr10","dtp3_reg","facility","educ_yrs","wealth_q","log_gdp","pct_urban","year_c","iso3","svkey")
base <- complete.cases(d[, need]) & is.finite(d$exposure) & d$exposure > 0
cat(sprintf("region-years: total with covariates=%d ; PfPR>=1%%=%d ; PfPR<1%%=%d\n",
            sum(base), sum(base & d$pfpr2_10 >= 1), sum(base & d$pfpr2_10 < 1)))
cat(sprintf("PfPR<1%% added set: prevalence range %.2f-%.2f%%, mean post-neonatal %.1f/1000 (vs %.1f overall)\n",
            min(d$pfpr2_10[base & d$pfpr2_10<1]), max(d$pfpr2_10[base & d$pfpr2_10<1]),
            mean(d$m1mo5y[base & d$pfpr2_10<1], na.rm=TRUE), mean(d$m1mo5y[base], na.rm=TRUE)))

RE  <- "s(country,bs=\"re\") + s(country,pfpr10,bs=\"re\") + offset(log(exposure))"
COV <- "dtp3_reg + facility + educ_yrs + wealth_q + log_gdp + pct_urban"
fit <- function(rate, lowincl) { dd <- d[base & is.finite(d[[rate]]) & d[[rate]] > 0 &
                                          (lowincl | d$pfpr2_10 >= 1), ]
  dd$deaths <- round(dd[[rate]]/1000*dd$exposure); dd$country <- factor(dd$iso3)
  gam(as.formula(paste("deaths ~ s(pfpr10) +", COV, "+ s(year_c) +", RE)),
      family = nb(), method = "REML", data = dd) }
fsm <- function(m,p){ nd <- m$model[1,]; nd$pfpr10 <- p/10; as.numeric(predict(m, nd, type="terms", terms="s(pfpr10)")) }
af30 <- function(m) 100*(1-exp(-(fsm(m,30)-fsm(m,1))))
sm   <- function(m) summary(m)$s.table["s(pfpr10)", c("edf","p-value")]

res <- data.frame()
for (oc in c("m1mo5y","nnmr")) { lab <- ifelse(oc=="m1mo5y","post-neonatal","neonatal")
  for (li in c(FALSE, TRUE)) { m <- fit(oc, li)
    res <- rbind(res, data.frame(outcome=lab, sample=ifelse(li,"all (incl <1%)",">=1% (primary)"),
      n=nrow(m$model), edf=round(sm(m)[1],2), p=signif(sm(m)[2],2), af30=round(af30(m),1)))
    if (oc=="m1mo5y") assign(ifelse(li,"pn_all","pn_pri"), m) } }
cat("\n=== combined model: malaria effect with vs without sub-1% region-years ===\n"); print(res, row.names=FALSE)

## ---- spline overlay: post-neonatal, >=1% vs all ----------------------------
sp <- function(mod, lab){ mf <- mod$model
  g <- data.frame(pfpr10 = seq(min(mf$pfpr10), max(mf$pfpr10), length.out=220))
  nd <- data.frame(pfpr10=g$pfpr10, dtp3_reg=mean(mf$dtp3_reg), facility=mean(mf$facility),
    educ_yrs=mean(mf$educ_yrs), wealth_q=mean(mf$wealth_q), log_gdp=mean(mf$log_gdp),
    pct_urban=mean(mf$pct_urban), year_c=0, exposure=1, country=levels(mf$country)[1])
  pr <- predict(mod, nd, type="terms", terms="s(pfpr10)", se.fit=TRUE)
  data.frame(pfpr2_10=g$pfpr10*10, sample=lab, fit=as.numeric(pr$fit),
             lo=as.numeric(pr$fit)-1.96*as.numeric(pr$se.fit), hi=as.numeric(pr$fit)+1.96*as.numeric(pr$se.fit)) }
L <- rbind(sp(pn_pri, ">=1% (primary)"), sp(pn_all, "all (incl <1%)"))
cols <- c(">=1% (primary)"="#08519c", "all (incl <1%)"="#d73027")
p <- ggplot(L, aes(pfpr2_10, fit, colour=sample, fill=sample)) +
  geom_hline(yintercept=0, linetype="dotted", colour="grey55") +
  geom_vline(xintercept=1, linetype="dashed", colour="grey70") +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, colour=NA) + geom_line(linewidth=1.1) +
  scale_colour_manual(values=cols, name="Sample") + scale_fill_manual(values=cols, name="Sample") +
  scale_x_log10(breaks=c(0.2,0.5,1,2,5,10,20,50,80)) +
  labs(x=expression(italic(Pf)*"PR"[2-10]~"(%), log scale"),
       y="s(PfPR) partial effect on post-neonatal mortality (log)") +
  theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(), legend.position="top")
ggsave(file.path(RESULTS,"sensitivity_lowprev_spline.png"), p, width=8.5, height=5.8, dpi=300)
write.csv(res, file.path(RESULTS,"sensitivity_lowprev.csv"), row.names=FALSE)
cat("\nsaved: results/sensitivity_lowprev_spline.png + sensitivity_lowprev.csv\n")

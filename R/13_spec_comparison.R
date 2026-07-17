# =============================================================================
# 13_spec_comparison.R — Method-2 specification comparison (backs SText 1)
# Four fixed-effect functional forms for prevalence, all with the same nb() GAM
# scaffold (dtp3 + log_gdp + pct_urban + s(year_c) + country random intercept &
# slope + log-exposure offset). Fitted on two samples so within-sample AIC is
# valid: (i) full; (ii) the analysis sample PfPR2-10 >= C2_MIN_PFPR.
#   M1 linear in PfPR2-10               <- PRIMARY (00_utils::fit_c2_primary)
#   M2 penalised spline s(PfPR2-10)
#   M3 linear in log PfPR2-10 (power law)
#   M4 penalised spline s(log PfPR2-10)
# On the analysis sample M1 minimises AIC and M2 collapses onto it (edf=1); the
# log-scale forms fit marginally worse and the power law clearly worst.
# Outputs: results/component2_spec_comparison.csv + results/sfig1_spec_comparison.png
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})

d <- read.csv(file.path(RESULTS, "component2_region_data.csv"), stringsAsFactors = FALSE)
covs <- c("dtp3", "log_gdp", "pct_urban")
d0 <- d[complete.cases(d[, c("m1mo5y", covs, "pfpr2_10", "year_c", "country", "svkey")]) &
          is.finite(d$exposure) & d$exposure > 0 & d$m1mo5y > 0, ]

base <- deaths ~ dtp3 + log_gdp + pct_urban + s(year_c) + offset(log(exposure))
fit4 <- function(dat) {
  dat$deaths   <- round(dat$m1mo5y / 1000 * dat$exposure)
  dat$country  <- factor(dat$country); dat$svkey <- factor(dat$svkey)
  dat$log_pfpr <- log(dat$pfpr2_10)
  list(
    M1 = gam(update(base, . ~ . + pfpr10          + s(country, bs="re") + s(country, pfpr10, bs="re")),   family=nb(), method="REML", data=dat),
    M2 = gam(update(base, . ~ . + s(pfpr10, k=4)  + s(country, bs="re") + s(country, pfpr10, bs="re")),   family=nb(), method="REML", data=dat),
    M3 = gam(update(base, . ~ . + log_pfpr          + s(country, bs="re") + s(country, log_pfpr, bs="re")), family=nb(), method="REML", data=dat),
    M4 = gam(update(base, . ~ . + s(log_pfpr, k=4)  + s(country, bs="re") + s(country, log_pfpr, bs="re")), family=nb(), method="REML", data=dat))
}
row1 <- function(m, nm, form) {
  s  <- summary(m)
  sm <- setdiff(grep("pfpr", rownames(s$s.table), value=TRUE),
                grep("country", rownames(s$s.table), value=TRUE))   # fixed prevalence smooth only
  edf <- if (length(sm)) s$s.table[sm[1], "edf"] else NA_real_
  eff <- if (nm == "M1") { b <- s$p.table["pfpr10", ];   sprintf("%+.1f%% /+10pt", (exp(b[1])-1)*100) }
         else if (nm == "M3") { b <- s$p.table["log_pfpr", ]; sprintf("%+.1f%% /doubling", (2^b[1]-1)*100) }
         else sprintf("edf=%.2f", edf)
  data.frame(model=nm, form=form, AIC=round(AIC(m),1), dev=round(100*s$dev.expl,1),
             edf=round(edf,2), effect=eff, stringsAsFactors=FALSE)
}
report <- function(dat, lab) {
  ms <- fit4(dat)
  forms <- c(M1="linear in PfPR", M2="spline s(PfPR)", M3="linear in log PfPR (power law)", M4="spline s(log PfPR)")
  tab <- do.call(rbind, Map(function(m, nm) row1(m, nm, forms[nm]), ms, names(ms)))
  tab$dAIC <- round(tab$AIC - min(tab$AIC), 1); tab$sample <- lab; tab$n <- nrow(dat)
  cat(sprintf("\n== %s (n=%d, %d countries) ==\n", lab, nrow(dat), nlevels(factor(dat$country))))
  print(tab[order(tab$AIC), c("model","form","AIC","dAIC","dev","edf","effect")], row.names=FALSE)
  tab
}
full <- report(d0, "Full sample")
tr1  <- report(d0[d0$pfpr2_10 >= C2_MIN_PFPR, ], sprintf("Analysis sample PfPR2-10 >= %g%%", C2_MIN_PFPR))
out <- rbind(full, tr1)[, c("sample","n","model","form","AIC","dAIC","dev","edf","effect")]
write.csv(out, file.path(RESULTS, "component2_spec_comparison.csv"), row.names=FALSE)

## figure: the four candidate fits on the analysis sample (log-log axes)
dd <- d0[d0$pfpr2_10 >= C2_MIN_PFPR, ]; dd$country <- factor(dd$country); dd$svkey <- factor(dd$svkey)
dd$deaths <- round(dd$m1mo5y/1000 * dd$exposure); dd$log_pfpr <- log(dd$pfpr2_10)
ms <- fit4(dd)
re_of <- function(m) grep("country", vapply(m$smooth, function(s) s$label, ""), value=TRUE)
gp <- data.frame(pfpr2_10 = exp(seq(log(C2_MIN_PFPR), log(max(dd$pfpr2_10)), length.out=200)))
predr <- function(m) { g <- data.frame(pfpr10=gp$pfpr2_10/10, log_pfpr=log(gp$pfpr2_10), dtp3=mean(dd$dtp3),
   log_gdp=mean(dd$log_gdp), pct_urban=mean(dd$pct_urban), year_c=0, exposure=1, country=dd$country[1])
   1000*exp(as.numeric(predict(m, g, type="link", exclude=re_of(m)))) }
cur <- data.frame(pfpr2_10 = rep(gp$pfpr2_10, 4),
  rate = c(predr(ms$M1), predr(ms$M2), predr(ms$M3), predr(ms$M4)),
  spec = rep(c("M1 linear PfPR [primary]","M2 spline s(PfPR)","M3 power law","M4 spline s(log PfPR)"), each=200))
p <- ggplot() +
  geom_point(data=dd, aes(pfpr2_10, m1mo5y, size=exposure), colour="grey60", alpha=0.3) +
  geom_line(data=cur, aes(pfpr2_10, rate, colour=spec, linetype=spec), linewidth=1) +
  scale_size_area(max_size=4.5, guide="none") +
  scale_x_log10(breaks=c(1,2,5,10,20,50,80)) + scale_y_log10(breaks=c(10,20,50,100,200)) +
  coord_cartesian(ylim=c(10,200)) +
  scale_colour_manual(values=c("M1 linear PfPR [primary]"="#d73027","M2 spline s(PfPR)"="#08519c",
     "M3 power law"="#984ea3","M4 spline s(log PfPR)"="#66c2a5"), name=NULL) +
  scale_linetype_manual(values=c("M1 linear PfPR [primary]"="solid","M2 spline s(PfPR)"="dashed",
     "M3 power law"="longdash","M4 spline s(log PfPR)"="dotted"), name=NULL) +
  labs(x=expression("Age-standardised "*italic(Pf)*"PR"[2-10]*" (%), log scale"),
       y="Post-neonatal mortality (per 1,000 lb, log scale)") +
  theme_bw(base_size=11) + theme(legend.position=c(0.99,0.02), legend.justification=c(1,0),
     legend.background=element_rect(fill=scales::alpha("white",0.75), colour=NA), legend.key.width=unit(9,"mm"))
ggsave(file.path(RESULTS, "sfig1_spec_comparison.png"), p, width=8, height=6, dpi=300)
cat("\nsaved: results/component2_spec_comparison.csv + results/sfig1_spec_comparison.png\n")

# =============================================================================
# 33_penalized_model.R — RIDGE-PENALISED version of the full combined model.
# The confounder block (standardised, logit for proportions) is placed under a
# single ridge penalty (mgcv paraPen, identity S; REML picks the shrinkage), so
# collinear predictors are shrunk jointly instead of sign-flipping. The malaria
# term s(pfpr10) and country RE/slope are left free. MEASLES dropped (DPT3 already
# captures immunisation); health-exp %GDP, exclusive BF, national log-GDP dropped
# as before. Fit for post-neonatal and neonatal; compare to the unpenalised full
# model (R/32).
# =============================================================================
source("R/00_utils.R"); suppressMessages({library(mgcv); library(ggplot2)})
d <- read.csv(file.path(RESULTS, "component2_region_data_full.csv"), stringsAsFactors = FALSE)
d$nnmr <- d$u5mr - d$m1mo5y

prop_vars <- c("dtp3_reg","facility","stunting","underweight","wasting","birth_int",
               "imp_water","imp_sanit","elec_dhs","pct_urban")                 # logit  (measles removed)
cont_vars <- c("educ_yrs","wealth_q","mage1","log_hexp_pc")
labs <- c(dtp3_reg="DPT3 coverage", facility="Facility delivery", stunting="Stunting",
  underweight="Underweight", wasting="Wasting", birth_int="Birth interval <24mo",
  imp_water="Improved water", imp_sanit="Improved sanitation", elec_dhs="Electricity (DHS)",
  pct_urban="% urban", educ_yrs="Maternal education", wealth_q="Wealth quintile",
  mage1="Age at first birth", log_hexp_pc="Health exp/cap (log)")
lgt <- function(x){ p <- pmin(pmax(x/100, 0.005), 0.995); log(p/(1-p)) }
for (v in prop_vars) d[[paste0(v,"_z")]] <- as.numeric(scale(lgt(d[[v]])))
for (v in cont_vars) d[[paste0(v,"_z")]] <- as.numeric(scale(d[[v]]))
Zc <- paste0(c(prop_vars, cont_vars), "_z")

fitd <- d[is.finite(d$exposure) & d$exposure > 0 & d$pfpr2_10 >= 1 &
            is.finite(d$m1mo5y) & d$m1mo5y > 0 & is.finite(d$nnmr) & d$nnmr > 0, ]
fitd$country <- factor(fitd$iso3)
G <- as.matrix(fitd[, Zc]); colnames(G) <- sub("_z$","",Zc); fitd$G <- G       # confounder block (ridge-penalised)
cat(sprintf("fit sample: %d region-years, %d countries; ridge block = %d covariates (measles removed)\n",
            nrow(fitd), nlevels(fitd$country), ncol(G)))

RE <- "s(country,bs=\"re\") + s(country,pfpr10,bs=\"re\") + offset(log(exposure))"
fitp <- function(rate){ fitd$deaths <- round(fitd[[rate]]/1000*fitd$exposure)
  gam(as.formula(paste("deaths ~ G + s(pfpr10) + s(year_c) +", RE)),
      family=nb(), method="REML", paraPen=list(G=list(diag(ncol(G)))), data=fitd) }
pn <- fitp("m1mo5y"); nn <- fitp("nnmr")

## ---- PfPR survival + shrinkage ---------------------------------------------
fsm <- function(m,p){ nd <- fitd[1,]; nd$pfpr10 <- p/10; as.numeric(predict(m, nd, type="terms", terms="s(pfpr10)")) }
af30 <- function(m) 100*(1-exp(-(fsm(m,30)-fsm(m,1))))
edf_G <- function(m) sum(m$edf[grep("^G", names(m$edf))])
sp <- function(m) summary(m)$s.table["s(pfpr10)", c("edf","p-value")]
cat(sprintf("\n[post-neonatal] s(pfpr10) edf=%.2f p=%.2g ; AF@30%%=%.1f%% ; ridge block edf=%.1f/%d (shrinkage)\n",
            sp(pn)[1], sp(pn)[2], af30(pn), edf_G(pn), ncol(G)))
cat(sprintf("[neonatal]      s(pfpr10) edf=%.2f p=%.2g ; AF@30%%=%.1f%% ; ridge block edf=%.1f/%d\n",
            sp(nn)[1], sp(nn)[2], af30(nn), edf_G(nn), ncol(G)))

## ---- penalised per-SD coefficients + forest --------------------------------
gcoef <- function(m, oc){ cf <- coef(m); V <- vcov(m); gi <- grep("^G", names(cf))
  nm <- sub("^G","",names(cf)[gi]); est <- unname(cf[gi]); se <- sqrt(diag(V)[gi])
  data.frame(outcome=oc, var=labs[nm], pct=(exp(est)-1)*100,
             lo=(exp(est-1.96*se)-1)*100, hi=(exp(est+1.96*se)-1)*100, row.names=NULL) }
CT <- rbind(gcoef(pn,"Post-neonatal (1mo-5y)"), gcoef(nn,"Neonatal (<1mo)"))
ord <- gcoef(pn,"x"); CT$var <- factor(CT$var, levels=ord$var[order(ord$pct)])
write.csv(CT, file.path(RESULTS,"penalized_model_coefficients.csv"), row.names=FALSE)
cat("\n=== penalised post-neonatal — % change per +1 SD (ridge-shrunk) ===\n")
pc <- CT[CT$outcome=="Post-neonatal (1mo-5y)",]; pc <- pc[order(-abs(pc$pct)),]
print(data.frame(predictor=as.character(pc$var), pct_per_SD=round(pc$pct,1),
                 CI=sprintf("%.1f,%.1f",pc$lo,pc$hi)), row.names=FALSE)
p <- ggplot(CT, aes(pct, var, colour=outcome)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey55") +
  geom_linerange(aes(xmin=lo, xmax=hi), position=position_dodge(width=0.55), linewidth=0.6) +
  geom_point(position=position_dodge(width=0.55), size=2) +
  scale_colour_manual(values=c("Post-neonatal (1mo-5y)"="#d73027","Neonatal (<1mo)"="#4575b4"), name=NULL) +
  labs(x="Ridge-penalised % change in mortality per +1 SD (95% CI)", y=NULL) +
  theme_bw(base_size=11) + theme(legend.position="top", panel.grid.minor=element_blank())
ggsave(file.path(RESULTS,"penalized_model_forest.png"), p, width=9, height=6, dpi=200)

## ---- compare to unpenalised full model (R/32) ------------------------------
uf <- tryCatch(readRDS(file.path(RESULTS,"full_models.rds"))$pn_full, error=function(e) NULL)
if (!is.null(uf)) cat(sprintf("\nAIC: unpenalised full %.1f  vs  ridge-penalised %.1f\n", AIC(uf), AIC(pn)))
saveRDS(list(pn=pn, nn=nn), file.path(RESULTS,"penalized_models.rds"))
cat("saved: results/penalized_model_forest.png + penalized_model_coefficients.csv + penalized_models.rds\n")

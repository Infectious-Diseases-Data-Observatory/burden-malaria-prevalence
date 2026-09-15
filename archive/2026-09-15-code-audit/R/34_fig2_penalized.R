# =============================================================================
# 34_fig2_penalized.R — rebuild PAPER FIGURE 2 (negative control) from the
# RIDGE-PENALISED model. Same confounder block as R/33 (14 covariates, measles
# removed, ridge paraPen), but PfPR entered LINEARLY so each outcome gives one
# comparable "% change per +10 PfPR pts" for the forest. Fit for all-U5,
# post-neonatal and neonatal; neonatal is the negative control.
# Writes results/paper_fig2_negcontrol.{png,csv} and copies to Overleaf.
# =============================================================================
source("R/00_utils.R"); suppressMessages({library(mgcv); library(ggplot2)})
d <- read.csv(file.path(RESULTS, "component2_region_data_full.csv"), stringsAsFactors = FALSE)
d$nnmr <- d$u5mr - d$m1mo5y

prop_vars <- c("dtp3_reg","facility","stunting","underweight","wasting","birth_int",
               "imp_water","imp_sanit","elec_dhs","pct_urban")            # measles removed
cont_vars <- c("educ_yrs","wealth_q","mage1","log_hexp_pc")
lgt <- function(x){ p <- pmin(pmax(x/100, 0.005), 0.995); log(p/(1-p)) }
for (v in prop_vars) d[[paste0(v,"_z")]] <- as.numeric(scale(lgt(d[[v]])))
for (v in cont_vars) d[[paste0(v,"_z")]] <- as.numeric(scale(d[[v]]))
Zc <- paste0(c(prop_vars, cont_vars), "_z")

fitd <- d[is.finite(d$exposure) & d$exposure > 0 & d$pfpr2_10 >= 1 &
            is.finite(d$m1mo5y) & d$m1mo5y > 0 & is.finite(d$nnmr) & d$nnmr > 0, ]
fitd$country <- factor(fitd$iso3)
G <- as.matrix(fitd[, Zc]); colnames(G) <- sub("_z$","",Zc); fitd$G <- G

# ridge-penalised confounder block; PfPR LINEAR (free), so its coef = per +10 pts
fit <- function(rate){ fitd$deaths <- round(fitd[[rate]]/1000*fitd$exposure)
  gam(deaths ~ pfpr10 + G + s(year_c) + s(country,bs="re") + s(country,pfpr10,bs="re") + offset(log(exposure)),
      family=nb(), method="REML", paraPen=list(G=list(diag(ncol(G)))), data=fitd) }
eff <- function(rate, lab){ m <- fit(rate); b <- unname(summary(m)$p.table["pfpr10", ])
  data.frame(outcome=lab, pct=(exp(b[1])-1)*100, lo=(exp(b[1]-1.96*b[2])-1)*100,
             hi=(exp(b[1]+1.96*b[2])-1)*100, p=b[4]) }
nc <- rbind(eff("u5mr","All under-5"), eff("m1mo5y","1mo-5y (post-neonatal)"), eff("nnmr","Neonatal (control)"))
nc$outcome <- factor(nc$outcome, levels=c("Neonatal (control)","1mo-5y (post-neonatal)","All under-5"))
write.csv(nc, file.path(RESULTS, "paper_fig2_negcontrol.csv"), row.names=FALSE)
cat("=== negative control (ridge-penalised, linear PfPR) — % change per +10 PfPR2-10 pts ===\n")
print(data.frame(outcome=as.character(nc$outcome), pct=round(nc$pct,1),
                 CI=sprintf("%.1f,%.1f",nc$lo,nc$hi), p=signif(nc$p,2)), row.names=FALSE)

pnc <- ggplot(nc, aes(pct, outcome)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey50") +
  geom_linerange(aes(xmin=lo, xmax=hi)) + geom_point(size=3, colour="#08519c") +
  labs(x="% change in mortality per +10 PfPR2-10 points (95% CI, fully adjusted + penalised)", y=NULL) +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(),
    axis.title=element_text(size=13), axis.text=element_text(size=13))
ggsave(file.path(RESULTS,"paper_fig2_negcontrol.png"), pnc, width=9, height=3.8, dpi=300)
FIG <- path.expand("~/Dropbox/Apps/Overleaf/malaria-burden-reassessment/figures/fig2_negcontrol.png")
file.copy(file.path(RESULTS,"paper_fig2_negcontrol.png"), FIG, overwrite=TRUE)
cat("\nsaved: results/paper_fig2_negcontrol.png  + copied to Overleaf\n")

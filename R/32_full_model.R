# =============================================================================
# 32_full_model.R — fold the FULL predictor set into the combined mortality model
# (nonlinear s(pfpr10) + country RE/slope), for post-neonatal AND neonatal.
# Per user: use imputed data; LOGIT-transform proportions; DROP health exp %GDP
# and exclusive BF; KEEP health exp per capita (log). Also drop national log-GDP
# (collinear with / superseded by log health-exp-pc). All covariates z-scored so
# coefficients read as % change in mortality per +1 SD (comparable forest plot).
# =============================================================================
source("R/00_utils.R"); suppressMessages({library(mgcv); library(ggplot2)})
d <- read.csv(file.path(RESULTS, "component2_region_data_full.csv"), stringsAsFactors = FALSE)
d$nnmr <- d$u5mr - d$m1mo5y

## ---- transform + standardise predictors ------------------------------------
prop_vars <- c("dtp3_reg","measles","facility","stunting","underweight","wasting",
               "birth_int","imp_water","imp_sanit","elec_dhs","pct_urban")   # logit
cont_vars <- c("educ_yrs","wealth_q","mage1","log_hexp_pc")                   # as-is (log_hexp_pc already log)
lgt <- function(x){ p <- pmin(pmax(x/100, 0.005), 0.995); log(p/(1-p)) }
labs <- c(dtp3_reg="DPT3 coverage", measles="Measles coverage", facility="Facility delivery",
  stunting="Stunting", underweight="Underweight", wasting="Wasting", birth_int="Birth interval <24mo",
  imp_water="Improved water", imp_sanit="Improved sanitation", elec_dhs="Electricity (DHS)",
  pct_urban="% urban", educ_yrs="Maternal education", wealth_q="Wealth quintile",
  mage1="Age at first birth", log_hexp_pc="Health exp/cap (log)")
for (v in prop_vars) d[[paste0(v,"_z")]] <- as.numeric(scale(lgt(d[[v]])))
for (v in cont_vars) d[[paste0(v,"_z")]] <- as.numeric(scale(d[[v]]))
Z <- paste0(c(prop_vars, cont_vars), "_z")

## ---- shared fit sample (imputed -> covariates complete) --------------------
fitd <- d[is.finite(d$exposure) & d$exposure > 0 & d$pfpr2_10 >= 1 &
            is.finite(d$m1mo5y) & d$m1mo5y > 0 & is.finite(d$nnmr) & d$nnmr > 0, ]
fitd$country <- factor(fitd$iso3)
cat(sprintf("fit sample: %d region-years, %d countries, %d surveys\n",
            nrow(fitd), nlevels(fitd$country), length(unique(fitd$svkey))))

RE  <- "s(country,bs=\"re\") + s(country,pfpr10,bs=\"re\") + offset(log(exposure))"
fit <- function(rate, rhs){ fitd$deaths <- round(fitd[[rate]]/1000*fitd$exposure)
  gam(as.formula(paste("deaths ~", rhs, "+ s(year_c) +", RE)), family=nb(), method="REML", data=fitd) }
CORE <- "dtp3_reg_z + facility_z + educ_yrs_z + wealth_q_z + pct_urban_z"     # reduced (combined m_comb covars)
FULL <- paste(Z, collapse=" + ")
pn_red  <- fit("m1mo5y", paste("s(pfpr10) +", CORE));  pn_full <- fit("m1mo5y", paste("s(pfpr10) +", FULL))
nn_full <- fit("nnmr",   paste("s(pfpr10) +", FULL))

## ---- PfPR survival + fit improvement ---------------------------------------
sp <- function(m) summary(m)$s.table["s(pfpr10)", c("edf","p-value")]
fsm <- function(m,p){ nd <- fitd[1,]; nd$pfpr10 <- p/10
  as.numeric(predict(m, nd, type="terms", terms="s(pfpr10)")) }
af30 <- function(m) 100*(1-exp(-(fsm(m,30)-fsm(m,1))))
cat(sprintf("\n[post-neonatal] reduced AIC %.1f  ->  FULL AIC %.1f  (dAIC %+.1f)\n",
            AIC(pn_red), AIC(pn_full), AIC(pn_full)-AIC(pn_red)))
cat(sprintf("[post-neonatal] s(pfpr10): edf=%.2f p=%.2g ; AF@30%% = %.1f%%  (was ~40%% pre-full-adjust)\n",
            sp(pn_full)[1], sp(pn_full)[2], af30(pn_full)))
cat(sprintf("[neonatal]      s(pfpr10): edf=%.2f p=%.2g ; AF@30%% = %.1f%%\n",
            sp(nn_full)[1], sp(nn_full)[2], af30(nn_full)))

## ---- coefficient table (per +1 SD % change) + forest -----------------------
ctab <- function(m, oc){ pt <- summary(m)$p.table; rn <- intersect(Z, rownames(pt))
  data.frame(outcome=oc, var=labs[sub("_z$","",rn)], est=pt[rn,"Estimate"], se=pt[rn,"Std. Error"],
             p=pt[rn,"Pr(>|z|)"], row.names=NULL) }
CT <- rbind(ctab(pn_full,"Post-neonatal (1mo-5y)"), ctab(nn_full,"Neonatal (<1mo)"))
CT$pct <- (exp(CT$est)-1)*100; CT$lo <- (exp(CT$est-1.96*CT$se)-1)*100; CT$hi <- (exp(CT$est+1.96*CT$se)-1)*100
ord <- ctab(pn_full,"x"); ord <- ord$var[order(ord$est)]
CT$var <- factor(CT$var, levels=ord)
write.csv(CT[,c("outcome","var","pct","lo","hi","p")], file.path(RESULTS,"full_model_coefficients.csv"), row.names=FALSE)
cat("\n=== FULL post-neonatal model — adjusted % change in mortality per +1 SD ===\n")
pc <- CT[CT$outcome=="Post-neonatal (1mo-5y)",]; pc <- pc[order(-abs(pc$pct)),]
print(data.frame(predictor=as.character(pc$var), pct_per_SD=round(pc$pct,1),
                 CI=sprintf("%.1f,%.1f",pc$lo,pc$hi), p=signif(pc$p,2)), row.names=FALSE)
p <- ggplot(CT, aes(pct, var, colour=outcome)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey55") +
  geom_linerange(aes(xmin=lo, xmax=hi), position=position_dodge(width=0.55), linewidth=0.6) +
  geom_point(position=position_dodge(width=0.55), size=2) +
  scale_colour_manual(values=c("Post-neonatal (1mo-5y)"="#d73027","Neonatal (<1mo)"="#4575b4"), name=NULL) +
  labs(x="Adjusted % change in mortality per +1 SD (95% CI)", y=NULL) +
  theme_bw(base_size=11) + theme(legend.position="top", panel.grid.minor=element_blank())
ggsave(file.path(RESULTS,"full_model_forest.png"), p, width=9, height=6.5, dpi=200)
saveRDS(list(pn_full=pn_full, nn_full=nn_full, pn_red=pn_red), file.path(RESULTS,"full_models.rds"))
cat("\nsaved: results/full_model_forest.png + full_model_coefficients.csv + full_models.rds\n")

# =============================================================================
# 08_country_malaria_deaths.R — malaria-attributable child deaths for EVERY
# country with national PfPR2-10 > 10%, by both methods, at national level, for
# TWO denominators: all under-5, and under-5 EXCLUDING neonates (1mo-5y).
#
#   Component 1 (share): deaths = share_hat(PfPR, GDP, DTP3)/100 x D
#   Component 2 (AF):    deaths = (1 - exp(-beta*PfPR/10)) x D
#   D(U5)    = U5MR/1000     x births         (all-cause under-5 deaths)
#   D(1mo5y) = (U5MR-NNMR)/1000 x births      (all-cause post-neonatal deaths)
# with the outcome-matched coefficients (share_u5 / share_1mo5y ; u5mr / m1mo5y).
#
# (The NG/DRC check in script 07 showed subnational aggregation changes the
# total by <1%, so national PfPR is used as the representative scale here.)
#
# Uncertainty: ONLY the malaria coefficient (C1 pfpr_pct slope; C2 pfpr10 beta),
# everything else fixed. The coefficient is shared across countries, so the CI on
# a POOLED total is the total re-evaluated at the coefficient's CI bounds.
# =============================================================================
source("R/00_utils.R")
suppressMessages(library(lme4))

## ---- national inputs; keep PfPR2-10 > 10% -----------------------------------
c1 <- read.csv(file.path(RESULTS, "component1_country_data.csv"), stringsAsFactors = FALSE)
c1$log_gdp <- log(c1$gdp_pc)
hi <- c1[is.finite(c1$pfpr_pct) & c1$pfpr_pct > 10 &
         is.finite(c1$u5mr) & is.finite(c1$m_1mo_5y) & is.finite(c1$births) &
         is.finite(c1$log_gdp) & is.finite(c1$dtp3), ]
hi$allcause_u5 <- hi$u5mr      / 1000 * hi$births      # all under-5
hi$allcause_pn <- hi$m_1mo_5y  / 1000 * hi$births      # post-neonatal (1mo-5y)

## ---- Component 1 share models (per outcome) ---------------------------------
fit_c1 <- function(y) { m <- lm(reformulate(c("pfpr_pct", "log_gdp", "dtp3"), y), data = c1)
  list(b = coef(m), se = summary(m)$coefficients["pfpr_pct", "Std. Error"]) }
c1u <- fit_c1("share_u5"); c1p <- fit_c1("share_1mo5y")
share_hat <- function(b, p, lg, dt, bpf) b["(Intercept)"] + bpf * p + b["log_gdp"] * lg + b["dtp3"] * dt

## ---- Component 2 mixed models (per outcome) ---------------------------------
d2 <- read.csv(file.path(RESULTS, "component2_region_data.csv"), stringsAsFactors = FALSE)
covs <- c("pfpr10", "dtp3", "log_gdp", "pct_urban", "year_c")
fit_c2 <- function(y) { dd <- d2[complete.cases(d2[, c(y, covs, "country", "svkey")]), ]
  m <- lmer(reformulate(c(covs, "(1 + pfpr10 || country)"), response = paste0("log(", y, ")")),
            data = dd, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  list(beta = fixef(m)["pfpr10"], se = sqrt(vcov(m)["pfpr10", "pfpr10"])) }
c2u <- fit_c2("u5mr"); c2p <- fit_c2("m1mo5y")
AF <- function(p, b) 1 - exp(-b * p / 10)

## ---- per-country estimates (point + coefficient-only 95% CI) ----------------
est_c1 <- function(mod, D) { s <- function(bpf) share_hat(mod$b, hi$pfpr_pct, hi$log_gdp, hi$dtp3, bpf) / 100 * D
  data.frame(x = s(mod$b["pfpr_pct"]), lo = s(mod$b["pfpr_pct"] - 1.96*mod$se), hi = s(mod$b["pfpr_pct"] + 1.96*mod$se)) }
est_c2 <- function(mod, D) data.frame(x = AF(hi$pfpr_pct, mod$beta) * D,
  lo = AF(hi$pfpr_pct, mod$beta - 1.96*mod$se) * D, hi = AF(hi$pfpr_pct, mod$beta + 1.96*mod$se) * D)

C1U <- est_c1(c1u, hi$allcause_u5); C2U <- est_c2(c2u, hi$allcause_u5)   # all under-5
C1P <- est_c1(c1p, hi$allcause_pn); C2P <- est_c2(c2p, hi$allcause_pn)   # 1mo-5y
res <- data.frame(country = hi$country, pfpr_pct = hi$pfpr_pct,
  allcause_u5 = hi$allcause_u5, c1_u5 = C1U$x, c1_u5_lo = C1U$lo, c1_u5_hi = C1U$hi,
  c2_u5 = C2U$x, c2_u5_lo = C2U$lo, c2_u5_hi = C2U$hi,
  allcause_pn = hi$allcause_pn, c1_pn = C1P$x, c1_pn_lo = C1P$lo, c1_pn_hi = C1P$hi,
  c2_pn = C2P$x, c2_pn_lo = C2P$lo, c2_pn_hi = C2P$hi, ihme = hi$deaths)
res <- res[order(-res$c2_u5), ]

## ---- pooled total (coefficient shared -> sum of bounds) ---------------------
num <- setdiff(names(res), c("country", "pfpr_pct"))
tot <- as.data.frame(as.list(colSums(res[, num]))); tot$country <- sprintf("TOTAL (%d countries)", nrow(res)); tot$pfpr_pct <- NA
res <- rbind(res, tot[, names(res)])
write.csv(res, file.path(RESULTS, "country_malaria_deaths_gt10.csv"), row.names = FALSE)

## ---- print ------------------------------------------------------------------
cat(sprintf("C1 slopes: U5 %.3f, 1mo-5y %.3f pp/pt | C2 beta: U5 %.4f (+%.1f%%/10pt), 1mo-5y %.4f (+%.1f%%/10pt)\n\n",
  c1u$b["pfpr_pct"], c1p$b["pfpr_pct"], c2u$beta, (exp(c2u$beta)-1)*100, c2p$beta, (exp(c2p$beta)-1)*100))
k <- function(x) format(round(x), big.mark = ",")
out <- data.frame(Country = res$country, PfPR = ifelse(is.na(res$pfpr_pct), "", sprintf("%.1f%%", res$pfpr_pct)),
  U5_all = k(res$allcause_u5), U5_C1 = k(res$c1_u5), U5_C2 = k(res$c2_u5),
  PN_all = k(res$allcause_pn), PN_C1 = k(res$c1_pn), PN_C2 = k(res$c2_pn),
  IHME = k(res$ihme), check.names = FALSE)
cat("Malaria-attributable deaths (point). U5 = all under-5; PN = under-5 excl. neonates (1mo-5y):\n\n")
print(out, row.names = FALSE)
tl <- res[nrow(res), ]
cat(sprintf("\nPOOLED 95%% CI (coefficient only):\n  U5    -> C1 %s [%s-%s] ; C2 %s [%s-%s]\n  1mo-5y-> C1 %s [%s-%s] ; C2 %s [%s-%s]\n",
  k(tl$c1_u5), k(tl$c1_u5_lo), k(tl$c1_u5_hi), k(tl$c2_u5), k(tl$c2_u5_lo), k(tl$c2_u5_hi),
  k(tl$c1_pn), k(tl$c1_pn_lo), k(tl$c1_pn_hi), k(tl$c2_pn), k(tl$c2_pn_lo), k(tl$c2_pn_hi)))
cat("\nsaved: results/country_malaria_deaths_gt10.csv\n")

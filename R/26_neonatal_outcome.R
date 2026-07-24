# =============================================================================
# 26_neonatal_outcome.R — refit the R/25 model set to NEONATAL mortality (<1mo)
# as a specificity check. Malaria's direct child toll is post-neonatal; the
# neonatal link is only indirect (placental malaria -> low birthweight/preterm),
# so a much weaker PfPR effect on neonatal than post-neonatal supports a causal
# reading rather than shared confounding.
#
# Neonatal rate recovered from the panel: nnmr = u5mr - m1mo5y (the identity used
# in mort_by_region). Same 2x2 as R/25 (linear/nonlinear PfPR x national/region-
# health covars), and a head-to-head linear effect vs post-neonatal on ONE shared
# sample. nb GAM, country RE + random slope, offset(log exposure), REML.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})

d <- read.csv(file.path(RESULTS, "component2_region_data_full.csv"), stringsAsFactors = FALSE)  # imputed panel (main analysis)
d$nnmr <- d$u5mr - d$m1mo5y                                        # neonatal = U5MR - post-neonatal
need <- c("m1mo5y","nnmr","pfpr10","dtp3","dtp3_reg","facility","educ_yrs","wealth_q",
          "log_gdp","pct_urban","year_c","iso3","svkey")
fitd <- d[complete.cases(d[, need]) & is.finite(d$exposure) & d$exposure > 0 &
            d$pfpr2_10 >= 1 & d$m1mo5y > 0 & d$nnmr > 0, ]
fitd$deaths_pn <- round(fitd$m1mo5y / 1000 * fitd$exposure)        # post-neonatal death count
fitd$deaths_nn <- round(fitd$nnmr  / 1000 * fitd$exposure)        # neonatal death count
fitd$country <- factor(fitd$iso3)
cat(sprintf("shared sample: %d region-years, %d countries, %d surveys\n",
            nrow(fitd), nlevels(fitd$country), length(unique(fitd$svkey))))
cat(sprintf("mean rates: post-neonatal %.1f/1000, neonatal %.1f/1000\n",
            mean(fitd$m1mo5y), mean(fitd$nnmr)))

RE  <- "s(country,bs=\"re\") + s(country,pfpr10,bs=\"re\") + offset(log(exposure))"
NAT <- "dtp3 + log_gdp + pct_urban"; RH <- "dtp3_reg + facility + educ_yrs + wealth_q + log_gdp + pct_urban"
fit <- function(outcome, core) gam(as.formula(paste(outcome, "~", core, "+ s(year_c) +", RE)),
                                   family = nb(), method = "REML", data = fitd)

## ---- 2x2 AIC ladder for the NEONATAL outcome --------------------------------
n_base <- fit("deaths_nn", paste("pfpr10 +", NAT))
n_full <- fit("deaths_nn", paste("pfpr10 +", RH))
n_addN <- fit("deaths_nn", paste("s(pfpr10) +", NAT))
n_comb <- fit("deaths_nn", paste("s(pfpr10) +", RH))
tab <- data.frame(
  model = c("n_base (linear, national)", "n_full (linear, region-health)",
            "n_addN (nonlinear, national)", "n_comb (nonlinear, region-health)"),
  edf = round(c(sum(n_base$edf), sum(n_full$edf), sum(n_addN$edf), sum(n_comb$edf)), 1),
  AIC = round(c(AIC(n_base), AIC(n_full), AIC(n_addN), AIC(n_comb)), 1))
tab$dAIC <- round(tab$AIC - min(tab$AIC), 1)
cat("\n=== NEONATAL 2x2 fit comparison (lower AIC better) ===\n"); print(tab, row.names = FALSE)
st <- summary(n_comb)$s.table
cat(sprintf("[n_comb] s(pfpr10): edf = %.2f, p = %.3g\n", st["s(pfpr10)","edf"], st["s(pfpr10)","p-value"]))

## ---- head-to-head: linear region-health PfPR effect, neonatal vs post-neo ---
p_full <- fit("deaths_pn", paste("pfpr10 +", RH))
eff <- function(m) { b <- unname(summary(m)$p.table["pfpr10", ])
  c(pct = (exp(b[1])-1)*100, lo = (exp(b[1]-1.96*b[2])-1)*100, hi = (exp(b[1]+1.96*b[2])-1)*100, p = b[4]) }
ep <- eff(p_full); en <- eff(n_full)
cat("\n=== PfPR effect per +10 PfPR2-10 pts (linear, region-health adjusted, SAME sample) ===\n")
cat(sprintf("  post-neonatal (1mo-5y): %+.1f%% (%+.1f to %+.1f)  p=%.2g\n", ep["pct"], ep["lo"], ep["hi"], ep["p"]))
cat(sprintf("  neonatal (<1mo)       : %+.1f%% (%+.1f to %+.1f)  p=%.2g\n", en["pct"], en["lo"], en["hi"], en["p"]))

## ---- AF vs 1% implied by the neonatal nonlinear term ------------------------
EX <- c("s(country)", "s(country,pfpr10)")
nd <- function(p10) data.frame(pfpr10 = p10, year_c = 0, dtp3_reg = mean(fitd$dtp3_reg),
  facility = mean(fitd$facility), educ_yrs = mean(fitd$educ_yrs), wealth_q = mean(fitd$wealth_q),
  log_gdp = mean(fitd$log_gdp), pct_urban = mean(fitd$pct_urban), exposure = 1, country = levels(fitd$country)[1])
e1 <- predict(n_comb, nd(0.1), type = "link", exclude = EX)
af <- function(p) as.numeric(1 - exp(-(predict(n_comb, nd(p/10), type = "link", exclude = EX) - e1)))
cat("\nNeonatal AF vs 1% (n_comb):\n")
for (p in c(5,10,20,30,50)) cat(sprintf("  PfPR2-10=%2d%%: AF=%.1f%%\n", p, 100*af(p)))

## ---- figure: neonatal vs post-neonatal fitted response ----------------------
g <- data.frame(pfpr2_10 = seq(1, 70, 0.5)); g$pfpr10 <- g$pfpr2_10/10
grid <- cbind(g, nd(g$pfpr10)[, setdiff(names(nd(0.1)), "pfpr10")]); grid$exposure <- 1000
p_comb <- fit("deaths_pn", paste("s(pfpr10) +", RH))
pl <- rbind(data.frame(pfpr2_10 = grid$pfpr2_10, outcome = "Post-neonatal (1mo-5y)",
                       rate = predict(p_comb, grid, type = "response", exclude = EX)),
            data.frame(pfpr2_10 = grid$pfpr2_10, outcome = "Neonatal (<1mo)",
                       rate = predict(n_comb, grid, type = "response", exclude = EX)))
p <- ggplot(pl, aes(pfpr2_10, rate, colour = outcome)) + geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("Post-neonatal (1mo-5y)" = "#d73027", "Neonatal (<1mo)" = "#4575b4"), name = NULL) +
  labs(x = expression(PfPR[2-10]~"(%)"), y = "Predicted deaths / 1000 (combined model)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave(file.path(RESULTS, "neonatal_vs_postneonatal_response.png"), p, width = 7.5, height = 6, dpi = 300)
saveRDS(list(n_base=n_base, n_full=n_full, n_addN=n_addN, n_comb=n_comb, tab=tab),
        file.path(RESULTS, "neonatal_models.rds"))
cat("\nsaved: results/neonatal_vs_postneonatal_response.png + neonatal_models.rds\n")

# =============================================================================
# 25_combined_model.R — the definitive best-fitting Method-2 model: NONLINEAR
# prevalence response s(pfpr10) AND region-level health-access covariates, i.e.
# the winners of both axes (R/20 form + R/23 covariates) combined.
#
# Fit a 2x2 on ONE sample (complete region-health covariates, PfPR2-10 >= 1%) so
# AIC is comparable: {linear, nonlinear PfPR} x {national, region-health covars}.
#   m_base : pfpr10     + national dtp3/gdp/urban          (R/23 baseline)
#   m_full : pfpr10     + region health + gdp/urban        (R/23 winner, covar axis)
#   m_addN : s(pfpr10)  + national dtp3/gdp/urban          (R/20 winner, form axis)
#   m_comb : s(pfpr10)  + region health + gdp/urban        (COMBINED)
# All: nb GAM, country random intercept + random slope on pfpr10, offset(log exposure), REML.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})

d <- read.csv(file.path(RESULTS, "component2_region_data_full.csv"), stringsAsFactors = FALSE)  # imputed panel: main analysis on all 921 region-years
need <- c("m1mo5y","pfpr10","dtp3","dtp3_reg","facility","educ_yrs","wealth_q",
          "log_gdp","pct_urban","year_c","iso3","svkey")
fitd <- d[complete.cases(d[, need]) & is.finite(d$exposure) & d$exposure > 0 & d$pfpr2_10 >= 1, ]
fitd$deaths  <- round(fitd$m1mo5y / 1000 * fitd$exposure)
fitd$country <- factor(fitd$iso3)
cat(sprintf("fit sample: %d region-years, %d countries, %d surveys, years %d-%d\n",
            nrow(fitd), nlevels(fitd$country), length(unique(fitd$svkey)), min(fitd$year), max(fitd$year)))

RE  <- "s(country,bs=\"re\") + s(country,pfpr10,bs=\"re\") + offset(log(exposure))"
fit <- function(core) gam(as.formula(paste("deaths ~", core, "+ s(year_c) +", RE)),
                          family = nb(), method = "REML", data = fitd)
NAT <- "dtp3 + log_gdp + pct_urban"; RH <- "dtp3_reg + facility + educ_yrs + wealth_q + log_gdp + pct_urban"
m_base <- fit(paste("pfpr10 +", NAT))
m_full <- fit(paste("pfpr10 +", RH))
m_addN <- fit(paste("s(pfpr10) +", NAT))
m_comb <- fit(paste("s(pfpr10) +", RH))

## ---- 2x2 AIC ladder ---------------------------------------------------------
tab <- data.frame(
  model  = c("m_base  (linear PfPR, national)", "m_full  (linear PfPR, region-health)",
             "m_addN  (nonlinear PfPR, national)", "m_comb  (nonlinear PfPR, region-health)"),
  pfpr   = c("linear","linear","s(pfpr10)","s(pfpr10)"),
  covars = c("national","region-health","national","region-health"),
  edf    = round(c(sum(m_base$edf), sum(m_full$edf), sum(m_addN$edf), sum(m_comb$edf)), 1),
  AIC    = round(c(AIC(m_base), AIC(m_full), AIC(m_addN), AIC(m_comb)), 1))
tab$dAIC <- round(tab$AIC - min(tab$AIC), 1)
cat("\n=== 2x2 fit comparison (same sample; lower AIC better) ===\n"); print(tab, row.names = FALSE)

## ---- combined model: prevalence term + covariate coefficients ---------------
st <- summary(m_comb)$s.table
cat(sprintf("\n[m_comb] s(pfpr10): edf = %.2f, p = %.2g\n",
            st["s(pfpr10)","edf"], st["s(pfpr10)","p-value"]))
cat("[m_comb] covariate coefficients (log-mortality scale):\n")
print(round(summary(m_comb)$p.table[c("dtp3_reg","facility","educ_yrs","wealth_q","log_gdp","pct_urban"),
                                    c("Estimate","Std. Error","Pr(>|z|)")], 4))

## ---- attributable fraction vs 1% implied by the nonlinear term --------------
EX <- c("s(country)", "s(country,pfpr10)")
nd <- function(p10) data.frame(pfpr10 = p10, year_c = 0, dtp3_reg = mean(fitd$dtp3_reg),
  facility = mean(fitd$facility), educ_yrs = mean(fitd$educ_yrs), wealth_q = mean(fitd$wealth_q),
  log_gdp = mean(fitd$log_gdp), pct_urban = mean(fitd$pct_urban), exposure = 1, country = levels(fitd$country)[1])
e1  <- predict(m_comb, nd(0.1), type = "link", exclude = EX)
af  <- function(p) as.numeric(1 - exp(-(predict(m_comb, nd(p/10), type = "link", exclude = EX) - e1)))
# same, national-covariate nonlinear model, to see if adjustment changes the response shape
e1n <- predict(m_addN, within(nd(0.1), {dtp3 <- mean(fitd$dtp3)}), type = "link", exclude = EX)
afN <- function(p) as.numeric(1 - exp(-(predict(m_addN, within(nd(p/10), {dtp3 <- mean(fitd$dtp3)}),
                                                 type = "link", exclude = EX) - e1n)))
cat("\nAF vs 1% (post-neonatal), combined vs national-covariate:\n")
for (p in c(5,10,20,30,50)) cat(sprintf("  PfPR2-10=%2d%%:  m_comb=%.1f%%   m_addN=%.1f%%\n", p, 100*af(p), 100*afN(p)))

## ---- figure: adjusted prevalence-response (m_comb) --------------------------
g <- data.frame(pfpr2_10 = seq(1, 70, 0.5)); g$pfpr10 <- g$pfpr2_10/10
gg <- cbind(g, nd(g$pfpr10)[, setdiff(names(nd(0.1)), c("pfpr10"))]); gg$exposure <- 1000
gg$rate <- predict(m_comb, gg, type = "response", exclude = EX)
p <- ggplot(gg, aes(pfpr2_10, rate)) +
  geom_point(data = fitd, aes(pfpr2_10, m1mo5y), inherit.aes = FALSE, alpha = 0.15, size = 0.8, colour = "grey45") +
  geom_line(linewidth = 1.1, colour = "#d73027") +
  labs(x = expression(PfPR[2-10]~"(%)"),
       y = "Predicted post-neonatal deaths / 1000 (combined model)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(RESULTS, "combined_model_response.png"), p, width = 7.5, height = 6, dpi = 300)
saveRDS(list(m_base=m_base, m_full=m_full, m_addN=m_addN, m_comb=m_comb, tab=tab,
             center = round(mean(d$year, na.rm = TRUE))), file.path(RESULTS, "combined_models.rds"))
cat("\nsaved: results/combined_model_response.png + combined_models.rds\n")

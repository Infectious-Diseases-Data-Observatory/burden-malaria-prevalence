# =============================================================================
# 27_random_slope_test.R — does the country random SLOPE on prevalence
# s(country,pfpr10,bs="re") improve fit, for post-neonatal AND neonatal? And for
# neonatal, is a nonlinear s(pfpr10) any better than a linear pfpr10 term?
#
# The random slope is a RANDOM-effect term, so AIC with vs without it (identical
# fixed effects) is a valid comparison; also report the between-country slope SD
# via gam.vcomp. All on the same 847-region-year sample; nb GAM, REML, country
# random intercept + region-health covariates + s(year_c) + offset(log exposure).
# =============================================================================
source("R/00_utils.R")
suppressMessages(library(mgcv))

d <- read.csv(file.path(RESULTS, "component2_region_data_expanded_health.csv"), stringsAsFactors = FALSE)
d$nnmr <- d$u5mr - d$m1mo5y
need <- c("m1mo5y","nnmr","pfpr10","dtp3_reg","facility","educ_yrs","wealth_q",
          "log_gdp","pct_urban","year_c","iso3","svkey")
fitd <- d[complete.cases(d[, need]) & is.finite(d$exposure) & d$exposure > 0 &
            d$pfpr2_10 >= 1 & d$m1mo5y > 0 & d$nnmr > 0, ]
fitd$deaths_pn <- round(fitd$m1mo5y / 1000 * fitd$exposure)
fitd$deaths_nn <- round(fitd$nnmr  / 1000 * fitd$exposure)
fitd$country   <- factor(fitd$iso3)
cat(sprintf("sample: %d region-years, %d countries\n\n", nrow(fitd), nlevels(fitd$country)))

RH  <- "dtp3_reg + facility + educ_yrs + wealth_q + log_gdp + pct_urban"
mk <- function(outcome, pf, slope) {
  terms <- c(pf, RH, "s(year_c)", "s(country,bs=\"re\")",
             if (slope) "s(country,pfpr10,bs=\"re\")", "offset(log(exposure))")
  gam(as.formula(paste(outcome, "~", paste(terms, collapse = " + "))),
      family = nb(), method = "REML", data = fitd)
}
slope_sd <- function(m) {                                   # between-country SD of the prevalence slope
  vc <- tryCatch(gam.vcomp(m, rescale = FALSE), error = function(e) NULL)
  if (is.null(vc)) return(NA_real_)
  r <- rownames(as.data.frame(vc)); i <- grep("country,pfpr10", r)
  if (!length(i)) NA_real_ else as.data.frame(vc)[i, 1]
}

## ---- Q1: random slope WITH vs WITHOUT, both outcomes (nonlinear s(pfpr10)) ---
pn_rs <- mk("deaths_pn","s(pfpr10)",TRUE);  pn_no <- mk("deaths_pn","s(pfpr10)",FALSE)
nn_rs <- mk("deaths_nn","s(pfpr10)",TRUE);  nn_no <- mk("deaths_nn","s(pfpr10)",FALSE)
q1 <- data.frame(
  outcome = c("post-neonatal","neonatal"),
  AIC_with_slope    = round(c(AIC(pn_rs), AIC(nn_rs)), 1),
  AIC_without_slope = round(c(AIC(pn_no), AIC(nn_no)), 1))
q1$dAIC_drop_slope <- round(q1$AIC_without_slope - q1$AIC_with_slope, 1)   # >0 => slope helps
q1$slope_SD_per10pts <- round(c(slope_sd(pn_rs), slope_sd(nn_rs)), 4)
cat("=== Q1: country random slope s(country,pfpr10) — with vs without ===\n")
cat("(dAIC_drop_slope > 0 means removing the slope WORSENS fit, i.e. the slope helps)\n")
print(q1, row.names = FALSE)
cat(sprintf("\n[post-neo] s(country,pfpr10) approx p = %.3g ; [neonatal] p = %.3g\n",
            summary(pn_rs)$s.table["s(country,pfpr10)","p-value"],
            summary(nn_rs)$s.table["s(country,pfpr10)","p-value"]))

## ---- Q2: neonatal — nonlinear s(pfpr10) vs linear pfpr10 --------------------
nl_rs <- mk("deaths_nn","pfpr10",TRUE);  nl_no <- mk("deaths_nn","pfpr10",FALSE)
q2 <- data.frame(
  model = c("s(pfpr10) + random slope", "linear pfpr10 + random slope",
            "s(pfpr10), no random slope", "linear pfpr10, no random slope"),
  edf = round(c(sum(nn_rs$edf), sum(nl_rs$edf), sum(nn_no$edf), sum(nl_no$edf)), 1),
  AIC = round(c(AIC(nn_rs), AIC(nl_rs), AIC(nn_no), AIC(nl_no)), 1))
q2$dAIC <- round(q2$AIC - min(q2$AIC), 1)
cat("\n=== Q2: NEONATAL — nonlinear s(pfpr10) vs linear pfpr10 ===\n")
print(q2, row.names = FALSE)
b <- unname(summary(nl_rs)$p.table["pfpr10", ])
cat(sprintf("\nneonatal linear PfPR effect: %+.1f%% per +10 pts (%.1f to %.1f), p=%.3g ; s(pfpr10) edf=%.2f\n",
            (exp(b[1])-1)*100, (exp(b[1]-1.96*b[2])-1)*100, (exp(b[1]+1.96*b[2])-1)*100, b[4],
            summary(nn_rs)$s.table["s(pfpr10)","edf"]))

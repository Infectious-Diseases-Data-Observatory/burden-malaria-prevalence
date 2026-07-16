# =============================================================================
# 06_attributable_fraction.R — malaria-attributable fraction of child mortality
# implied by the Component 2 DHS models.
#
# AF(p) = 1 - exp(-(eta(p) - eta(0))): the proportion of child deaths at
# prevalence p that would be averted if transmission fell to zero, holding
# covariates fixed (a counterfactual attributable fraction, assuming the
# adjusted log-rate association is causal). Malaria deaths per 1,000 = N x AF(p).
#
# Two functional forms x two outcomes:
#   linear  — closed form 1 - exp(-beta * p/10) from the log-rate LMM coefficient
#   spline  — 1 - exp(-(s_hat(p) - s_hat(0))) from the nb-GAM smooth s(pfpr10)
#   outcomes: U5MR (5q0) and 1mo-5y (neonatal-excluded)
# No-stunting specification, full sample (600 survey-regions, 23 countries).
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(lme4); library(mgcv); library(ggplot2); library(patchwork)})

d <- read.csv(file.path(RESULTS, "component2_region_data.csv"), stringsAsFactors = FALSE)
covs   <- c("pfpr10", "dtp3", "log_gdp", "pct_urban", "year_c")   # stunting dropped
grid_p <- seq(0, 80, by = 0.5)                                    # PfPR2-10 (%), within observed 0-81%
LAB    <- c(u5mr = "All under-5 (U5MR)", m1mo5y = "1mo-5y (neonatal excluded)")

## ---- linear AF from the LMM log-rate coefficient ----------------------------
af_linear <- function(outcome) {
  dd <- d[complete.cases(d[, c(outcome, covs, "country", "svkey")]), ]
  m  <- lmer(reformulate(c(covs, "(1 + pfpr10 || country)"), response = paste0("log(", outcome, ")")),
             data = dd, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  b  <- fixef(m)["pfpr10"]; se <- sqrt(vcov(m)["pfpr10", "pfpr10"])
  data.frame(pfpr = grid_p,
             af = 1 - exp(-b * grid_p / 10),
             lo = 1 - exp(-(b - 1.96 * se) * grid_p / 10),
             hi = 1 - exp(-(b + 1.96 * se) * grid_p / 10),
             method = "linear (LMM)", outcome = LAB[[outcome]])
}

## ---- nonlinear AF from the GAM smooth, contrast vs p=0 (delta-method CI) -----
af_spline <- function(outcome) {
  dg <- d[complete.cases(d[, c(outcome, covs, "country", "svkey")]) & is.finite(d$exposure) & d$exposure > 0, ]
  dg$deaths <- round(dg[[outcome]] / 1000 * dg$exposure); dg$country <- factor(dg$country)
  m <- mgcv::gam(deaths ~ s(pfpr10, k = 4) + dtp3 + log_gdp + pct_urban + s(year_c) +
                   s(country, bs = "re") + s(country, pfpr10, bs = "re") + offset(log(exposure)),
                 family = mgcv::nb(), method = "REML", data = dg)   # k=4: low-df smooth (max ~3 edf)
  nd <- data.frame(pfpr10 = grid_p / 10, dtp3 = mean(dg$dtp3), log_gdp = mean(dg$log_gdp),
                   pct_urban = mean(dg$pct_urban), year_c = 0, exposure = 1, country = dg$country[1])
  X   <- predict(m, nd, type = "lpmatrix"); V <- vcov(m); bb <- coef(m)
  idx <- grep("^s\\(pfpr10\\)", colnames(X))                      # population prevalence smooth only
  X0  <- X[1, ]                                                   # reference row = p=0
  est <- se <- numeric(nrow(X))
  for (i in seq_len(nrow(X))) { dv <- X[i, ] - X0; dv[-idx] <- 0
    est[i] <- sum(dv * bb); se[i] <- sqrt(drop(t(dv) %*% V %*% dv)) }
  edf <- summary(m)$s.table["s(pfpr10)", "edf"]
  attr_df <- data.frame(pfpr = grid_p,
             af = 1 - exp(-est),
             lo = 1 - exp(-(est - 1.96 * se)),
             hi = 1 - exp(-(est + 1.96 * se)),
             method = sprintf("spline (GAM, edf=%.1f)", edf), outcome = LAB[[outcome]])
  attr_df
}

af <- do.call(rbind, lapply(c("u5mr", "m1mo5y"), function(o) rbind(af_linear(o), af_spline(o))))
af$method_grp <- ifelse(grepl("^linear", af$method), "linear (LMM)", "spline (GAM)")
af$outcome <- factor(af$outcome, levels = LAB)

write.csv(af, file.path(RESULTS, "attributable_fraction.csv"), row.names = FALSE)

## ---- figure -----------------------------------------------------------------
cols <- c("linear (LMM)" = "#2166ac", "spline (GAM)" = "#b2182b")
p <- ggplot(af, aes(pfpr, 100 * af, colour = method_grp, fill = method_grp)) +
  geom_ribbon(aes(ymin = 100 * lo, ymax = 100 * hi), alpha = 0.13, colour = NA) +
  geom_line(aes(linetype = method_grp), linewidth = 1) +
  facet_wrap(~outcome) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_linetype_manual(values = c("linear (LMM)" = "22", "spline (GAM)" = "solid"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 80, 20)) +
  labs(x = expression("Age-standardised "*PfPR[2-10]*" (%)"),
       y = "Malaria-attributable fraction of child deaths (%)",
       title = "Model-implied malaria-attributable fraction of child mortality (Component 2 DHS models)",
       subtitle = "AF(p) = 1 - exp(-(η(p) - η(0))): proportion of child deaths averted if transmission fell to zero, covariates fixed.",
       caption = "Linear = LMM log-rate coefficient; spline = nb-GAM smooth s(PfPR). Shaded = 95% CI (prevalence-effect sampling error only). No-stunting spec, 600 survey-regions.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave(file.path(RESULTS, "attributable_fraction.png"), p, width = 11, height = 5.5, dpi = 140)

## ---- console summary at round prevalences -----------------------------------
show <- af[af$pfpr %in% c(10, 20, 30, 40, 50), c("outcome", "method_grp", "pfpr", "af", "lo", "hi")]
show[, c("af", "lo", "hi")] <- round(100 * show[, c("af", "lo", "hi")], 1)
cat("Malaria-attributable fraction of child deaths (%), by outcome/method:\n"); print(show, row.names = FALSE)
cat("\nsaved: results/attributable_fraction.png + attributable_fraction.csv\n")

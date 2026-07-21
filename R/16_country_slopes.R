# =============================================================================
# 16_country_slopes.R — caterpillar plot of country-specific prevalence slopes
# from the primary Method-2 model (M1: linear nb-GAM, PfPR2-10 >= 1%,
# post-neonatal outcome). Each country's slope = fixed pfpr10 effect + its
# random-slope deviation; 95% shrinkage intervals from the Bayesian posterior
# covariance (Vp). Ordered by effect size.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv); library(ggplot2)})

d <- read.csv(file.path(RESULTS, "component2_region_data.csv"), stringsAsFactors = FALSE)
f <- fit_c2_primary(d, "m1mo5y")                 # M1 primary model
m <- f$m; dd <- f$dd; b <- f$beta
cf <- coef(m); V <- vcov(m)                       # Vp: posterior covariance (gives shrinkage SEs)

## country slope_c = pfpr10 (fixed) + random-slope deviation_c; SE from a contrast
ki  <- which(vapply(m$smooth, function(s) s$label, "") == "s(country,pfpr10)")
sm  <- m$smooth[[ki]]; idx <- sm$first.para:sm$last.para; lv <- levels(dd$country)
bp  <- which(names(cf) == "pfpr10")
df  <- do.call(rbind, lapply(seq_along(idx), function(i) {
  v <- numeric(length(cf)); v[bp] <- 1; v[idx[i]] <- 1
  est <- sum(v * cf); se <- sqrt(drop(t(v) %*% V %*% v))
  data.frame(country = lv[i], est = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
}))
pc <- function(x) (exp(x) - 1) * 100
df$pct <- pc(df$est); df$pct_lo <- pc(df$lo); df$pct_hi <- pc(df$hi)
df <- df[order(df$pct), ]; df$country <- factor(df$country, levels = df$country)
pop <- pc(b)
write.csv(df, file.path(RESULTS, "component2_country_slopes_m1.csv"), row.names = FALSE)

p <- ggplot(df, aes(pct, country)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = pop, linetype = "dashed", colour = "#d73027", linewidth = 0.6) +
  geom_linerange(aes(xmin = pct_lo, xmax = pct_hi), orientation = "y", colour = "grey55", linewidth = 0.5) +
  geom_point(colour = "#08519c", size = 2.3) +
  labs(x = "% change in post-neonatal mortality per +10 PfPR2-10 points (95% CI)", y = NULL,
       caption = sprintf("Dashed red = population estimate (%+.1f%%); grey line = no effect. Blue points = country slopes with 95%% shrinkage intervals.", pop)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "grey93"),
        axis.title = element_text(size = 13), axis.text.y = element_text(size = 10))
ggsave(file.path(RESULTS, "fig_country_slopes.png"), p, width = 8, height = 6.8, dpi = 320)
cat(sprintf("saved: results/fig_country_slopes.png; %d countries, range %+.1f%% to %+.1f%%, all positive=%s\n",
            nrow(df), min(df$pct), max(df$pct), all(df$pct > 0)))
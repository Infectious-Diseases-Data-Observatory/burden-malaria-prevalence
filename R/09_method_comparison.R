# =============================================================================
# 09_method_comparison.R — overlay the two methods on ONE set of axes.
# Both estimate the same quantity: malaria as a % of child deaths vs PfPR2-10.
#   Component 1 (country level): IHME/IGME malaria SHARE per country (points) +
#     fitted share model (share ~ PfPR + log GDP + DTP3) at mean covariates.
#   Component 2 (DHS):           attributable fraction AF = 1 - exp(-(s(PfPR)))
#     from the log-rate model (the spline curve produced by script 06).
# Faceted by outcome: all under-5, and 1mo-5y (neonatal-excluded).
# =============================================================================
source("R/00_utils.R")
suppressMessages(library(ggplot2))

OUT   <- c(share_u5 = "All under-5", share_1mo5y = "1mo-5y (neonatal excluded)")
AFMAP <- c("All under-5 (U5MR)" = "All under-5", "1mo-5y (neonatal excluded)" = "1mo-5y (neonatal excluded)")
LVL   <- c("All under-5", "1mo-5y (neonatal excluded)")
C1LAB <- "Component 1 — country share (IHME/IGME)"
C2LAB <- "Component 2 — DHS attributable fraction"

## ---- Component 1: country points + adjusted fitted line (mean covariates) ---
c1 <- read.csv(file.path(RESULTS, "component1_country_data.csv"), stringsAsFactors = FALSE)
c1$log_gdp <- log(c1$gdp_pc)
pts <- do.call(rbind, lapply(names(OUT), function(o)
  data.frame(pfpr = c1$pfpr_pct, share = c1[[o]], outcome = OUT[[o]])))
c1line <- do.call(rbind, lapply(names(OUT), function(o) {
  m  <- fit_c1_share(c1, o)                                    # shared spec (00_utils.R)
  g  <- data.frame(pfpr_pct = seq(min(c1$pfpr_pct), max(c1$pfpr_pct), length = 120),
                   log_gdp = mean(c1$log_gdp), dtp3 = mean(c1$dtp3))
  pr <- predict(m, g, se.fit = TRUE)
  data.frame(pfpr = g$pfpr_pct, y = pr$fit, lo = pr$fit - 1.96*pr$se.fit, hi = pr$fit + 1.96*pr$se.fit,
             method = C1LAB, outcome = OUT[[o]]) }))

## ---- Component 2: attributable-fraction spline (from script 06) -------------
af <- read.csv(file.path(RESULTS, "attributable_fraction.csv"), stringsAsFactors = FALSE)
af <- af[grepl("spline", af$method), ]
c2line <- data.frame(pfpr = af$pfpr, y = 100*af$af, lo = 100*af$lo, hi = 100*af$hi,
                     method = C2LAB, outcome = AFMAP[af$outcome])

lines <- rbind(c1line, c2line)
lines$outcome <- factor(lines$outcome, levels = LVL); pts$outcome <- factor(pts$outcome, levels = LVL)
lines$method <- factor(lines$method, levels = c(C1LAB, C2LAB))

## ---- figure -----------------------------------------------------------------
cols <- setNames(c("#2c7fb8", "#d73027"), c(C1LAB, C2LAB))
p <- ggplot() +
  geom_point(data = pts, aes(pfpr, share), colour = "grey65", size = 1.4, alpha = 0.6) +
  geom_ribbon(data = lines, aes(pfpr, ymin = lo, ymax = hi, fill = method), alpha = 0.15) +
  geom_line(data = lines, aes(pfpr, y, colour = method), linewidth = 1) +
  facet_wrap(~outcome) +
  scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, NA)) +
  labs(x = expression(PfPR[2-10]*" (%)"), y = "Malaria as % of child deaths",
       title = "Two methods for malaria's share of child deaths vs prevalence",
       subtitle = "Grey points = Component 1 country shares. Lines at mean covariates, shaded = 95% CI. C1 drawn over its national-prevalence data range.",
       caption = "C1: share ~ PfPR + log GDP + DTP3 (country level, IHME/IGME). C2: AF = 1 - exp(-(s(PfPR))) from the DHS log-rate model (script 06 spline).") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave(file.path(RESULTS, "method_comparison.png"), p, width = 12, height = 5.8, dpi = 140)
cat("saved: results/method_comparison.png\n")

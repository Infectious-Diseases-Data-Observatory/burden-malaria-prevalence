# =============================================================================
# 31_bivariate_plots.R — bivariate relationships between each candidate predictor
# and the two outcomes: NEONATAL (<1mo) and POST-NEONATAL (1mo-5y) all-cause
# mortality, region-level. Uses OBSERVED values only (imputed cells dropped via
# the *_imp flags from R/30). Mortality on a log y-axis (the model is
# multiplicative). Also prints/saves a Spearman-correlation summary.
# =============================================================================
source("R/00_utils.R"); suppressMessages(library(ggplot2))
d <- read.csv(file.path(RESULTS, "component2_region_data_full.csv"), stringsAsFactors = FALSE)
d$nnmr <- d$u5mr - d$m1mo5y

# variable -> (display label, x-axis transform). Proportions on natural %; skewed on log.
V <- list(
  dtp3_reg    = list("DPT3 coverage (%)",          identity),
  measles     = list("Measles coverage (%)",       identity),
  facility    = list("Facility delivery (%)",       identity),
  excl_bf     = list("Exclusive BF <6mo (%)",       identity),
  stunting    = list("Stunting HAZ<-2 (%)",         identity),
  underweight = list("Underweight WAZ<-2 (%)",      identity),
  wasting     = list("Wasting WHZ<-2 (%)",          identity),
  educ_yrs    = list("Maternal educ (years)",       identity),
  wealth_q    = list("Wealth quintile (1-5)",       identity),
  birth_int   = list("Birth interval <24mo (%)",    identity),
  mage1       = list("Age at first birth (yrs)",    identity),
  imp_water   = list("Improved water (%)",          identity),
  imp_sanit   = list("Improved sanitation (%)",     identity),
  elec_dhs    = list("Electricity, DHS region (%)", identity),
  hexp_gdp    = list("Health exp (% GDP)",          identity),
  log_hexp_pc = list("Health exp per cap (log USD)", identity))

## ---- long frame: observed covariate x each outcome -------------------------
mk <- function(v) { lab <- V[[v]][[1]]; tr <- V[[v]][[2]]
  obs <- !d[[paste0(v, "_imp")]] & is.finite(d[[v]])
  x <- tr(d[[v]][obs])
  rbind(data.frame(var = lab, x = x, outcome = "Post-neonatal (1mo-5y)", mort = d$m1mo5y[obs]),
        data.frame(var = lab, x = x, outcome = "Neonatal (<1mo)",        mort = d$nnmr[obs])) }
L <- do.call(rbind, lapply(names(V), mk))
L <- L[is.finite(L$mort) & L$mort > 0, ]
L$var <- factor(L$var, levels = vapply(V, `[[`, "", 1))

p <- ggplot(L, aes(x, mort, colour = outcome)) +
  geom_point(alpha = 0.18, size = 0.6) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  facet_wrap(~ var, scales = "free_x", ncol = 4) +
  scale_y_log10(breaks = c(5, 10, 20, 50, 100, 200)) +
  scale_colour_manual(values = c("Post-neonatal (1mo-5y)" = "#d73027", "Neonatal (<1mo)" = "#4575b4"), name = NULL) +
  labs(x = NULL, y = "All-cause mortality (per 1,000 live births, log scale)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9))
ggsave(file.path(RESULTS, "bivariate_predictors.png"), p, width = 13, height = 12, dpi = 200)

## ---- Spearman correlation summary (observed) -------------------------------
sp <- function(v) { obs <- !d[[paste0(v,"_imp")]] & is.finite(d[[v]])
  c(n = sum(obs & is.finite(d$m1mo5y) & is.finite(d$nnmr)),
    post_neo = suppressWarnings(cor(d[[v]][obs], d$m1mo5y[obs], method = "spearman", use = "complete.obs")),
    neonatal = suppressWarnings(cor(d[[v]][obs], d$nnmr[obs],   method = "spearman", use = "complete.obs"))) }
tab <- data.frame(variable = vapply(V, `[[`, "", 1), t(sapply(names(V), sp)))
tab[, c("post_neo","neonatal")] <- round(tab[, c("post_neo","neonatal")], 2)
tab <- tab[order(-abs(tab$post_neo)), ]
write.csv(tab, file.path(RESULTS, "bivariate_correlations.csv"), row.names = FALSE)
cat("=== Spearman correlation with mortality (region-level, observed) ===\n")
print(tab, row.names = FALSE)
cat("\nsaved: results/bivariate_predictors.png + bivariate_correlations.csv\n")

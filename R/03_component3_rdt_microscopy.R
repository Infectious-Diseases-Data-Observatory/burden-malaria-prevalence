# =============================================================================
# 03_component3_rdt_microscopy.R — RDT vs MICROSCOPY (minor component)
# Across all DHS/MIS survey-regions with both tests, relate microscopy to RDT
# prevalence and derive a consensus conversion factor. Also builds the shared
# per-region prevalence + covariate table that Component 2 consumes.
# Runs BEFORE Component 2.
# =============================================================================
source("R/00_utils.R")
suppressMessages(library(ggplot2))
DHS <- file.path(DATA, "dhs")

## ---- DHS 2-letter code -> country / iso3 ------------------------------------
univ <- read.csv(file.path(DATA, "ssa_usable_dhs_surveys.csv"), stringsAsFactors = FALSE)
cc_iso <- unique(data.frame(cc = univ$DHS_CountryCode, country = univ$CountryName,
                            iso3 = countrycode::countrycode(univ$CountryName, "country.name", "iso3c", warn = FALSE),
                            stringsAsFactors = FALSE))

## ---- per-region prevalence + covariates across all PR recodes (cached) ------
PREV <- file.path(DATA, "dhs_prevalence_by_region.csv")
if (!file.exists(PREV)) {
  message("Extracting per-region prevalence + covariates from PR recodes ...")
  bases <- toupper(sub("[.]rds$", "", list.files(DHS, pattern = "rds$")))
  rows <- list()
  for (b in bases[substr(bases, 3, 4) == "PR"]) {
    cc  <- substr(b, 1, 2)
    rvp <- if (cc == "NG") "shstate" else "hv024"
    pr  <- readRDS(file.path(DHS, paste0(b, ".rds")))
    tab <- tryCatch(prev_cov_by_region(pr, rvp), error = function(e) NULL)
    if (is.null(tab)) next
    tab$survey  <- b; tab$cc <- cc
    tab$country <- cc_iso$country[match(cc, cc_iso$cc)]
    tab$iso3    <- cc_iso$iso3[match(cc, cc_iso$cc)]
    tab$year    <- tryCatch(1900L + floor((median(as.numeric(pr$hv008), na.rm = TRUE) - 1) / 12), error = function(e) NA)
    rows[[b]]   <- tab
    cat(sprintf("  %s (%s, %s): %d regions, RDT %.1f-%.1f%%\n", b, tab$country[1], tab$year[1],
                nrow(tab), min(tab$rdt, na.rm = TRUE), max(tab$rdt, na.rm = TRUE)))
  }
  write.csv(do.call(rbind, rows), PREV, row.names = FALSE)
}
prev <- read.csv(PREV, stringsAsFactors = FALSE)

## ---- RDT vs microscopy conversion -------------------------------------------
b2 <- prev[is.finite(prev$rdt) & is.finite(prev$mic) & prev$rdt > 0, ]
fit0 <- lm(mic ~ 0 + rdt, data = b2)    # through-origin: microscopy = k x RDT
fit1 <- lm(mic ~ rdt,     data = b2)    # with intercept
k <- unname(coef(fit0)[1])
conv <- data.frame(
  model     = c("through_origin", "linear"),
  slope     = c(coef(fit0)[1], coef(fit1)[2]),
  intercept = c(0, coef(fit1)[1]),
  pearson_r = cor(b2$rdt, b2$mic),
  n_regions = nrow(b2))
write.csv(conv, file.path(RESULTS, "rdt_microscopy_conversion.csv"), row.names = FALSE)
cat(sprintf("\nConversion (n=%d regions, r=%.2f): microscopy ~= %.3f x RDT  (RDT over-detects; k<1)\n",
            nrow(b2), cor(b2$rdt, b2$mic), k))

## ---- figure -----------------------------------------------------------------
lim <- max(b2$rdt, b2$mic, na.rm = TRUE)
p <- ggplot(b2, aes(rdt, mic)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_abline(slope = k, intercept = 0, colour = "#d73027", linewidth = 0.8) +
  geom_point(aes(colour = country), size = 1.6, alpha = 0.85) +
  annotate("text", x = lim * 0.05, y = lim * 0.95, hjust = 0, fontface = "bold", size = 3.6,
           label = sprintf("microscopy = %.2f x RDT\n(r = %.2f, n = %d regions)", k, cor(b2$rdt, b2$mic), nrow(b2))) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2))) +
  coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
  labs(x = "RDT-positive, children 6-59mo (%)", y = "Microscopy-positive, children 6-59mo (%)",
       title = "Component 3 — RDT vs microscopy parasite prevalence (DHS/MIS survey-regions)",
       subtitle = sprintf("Red = through-origin conversion (microscopy = %.2f x RDT); dashed = 1:1.", k),
       caption = "All DHS/MIS survey-regions with both tests. Conversion used in Component 2 to impute microscopy where only RDT exists.") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
    legend.key.size = unit(3.3, "mm"), legend.text = element_text(size = 6.5), legend.title = element_text(size = 8))
ggsave(file.path(RESULTS, "component3_rdt_vs_microscopy.png"), p, width = 10, height = 7.5, dpi = 140)
cat("saved: results/component3_rdt_vs_microscopy.png + rdt_microscopy_conversion.csv\n")

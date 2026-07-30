# =============================================================================
# 23_deprivation_proxy.R — is PfPR2-10 a stronger proxy for deprivation in West
# Africa than elsewhere, and did that entanglement change over time?
#
# MOTIVATION. Script 20 found the prevalence-mortality slope steepening in West
# Africa and flattening in East & Southern Africa. One explanation is confounding
# rather than causation: if higher-transmission regions are also systematically
# poorer, and more so in West Africa, then the West African PfPR coefficient
# carries more residual socioeconomic confounding. If the entanglement also GREW
# over time it could generate the apparent steepening.
#
# METHOD. Correlations are computed WITHIN SURVEY: both PfPR and each marker are
# residualised on survey fixed effects, so the question is whether, inside one
# country at one moment, the higher-transmission regions are the worse-off ones.
# Between-country and between-era differences cannot contribute.
#
# FINDING. The entanglement is about twice as strong in West Africa (composite
# r = 0.49) as in East & Southern Africa (r = 0.27), a difference that is unlikely
# to be chance (Fisher z = 3.46, p = 0.0005), and it holds across five of the six
# markers. But it did NOT grow in West Africa (0.49, 0.57, 0.43 across three eras)
# while it did rise in East & Southern Africa (0.19, 0.36, 0.31). So deprivation
# confounding plausibly inflates the West African slope at ALL times - a level
# effect relevant to the regional attributable fraction - but does not explain the
# steepening, which needs another mechanism.
#
# Outputs (results/dhs_rebuild/):
#   deprivation_proxy_by_marker.csv   within-survey r for each marker and region
#   deprivation_proxy_composite.csv   composite index, region comparison and eras
# =============================================================================
source("R_dhs/00_config.R")

WEST <- c("BEN","BFA","CIV","GHA","GIN","GMB","LBR","MLI","MRT","NER","NGA","SEN","SLE","TGO")
CENTRAL <- c("AGO","CMR","COD","COG","GAB","TCD")
EAST <- c("BDI","COM","ETH","KEN","MDG","MOZ","MWI","NAM","RWA","SWZ","UGA","ZMB","ZWE")
MARKERS <- c(pct_urban = "Urban residence", educ_yrs = "Maternal education",
             facility = "Facility delivery", imp_water = "Improved water",
             imp_sanit = "Improved sanitation", elec_dhs = "Household electricity")
GROUPS <- c("West Africa", "Central Africa", "East & Southern")

analysis <- read_analysis_data()
analysis <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
analysis$region_group <- ifelse(analysis$iso3 %in% WEST, "West Africa",
  ifelse(analysis$iso3 %in% CENTRAL, "Central Africa",
    ifelse(analysis$iso3 %in% EAST, "East & Southern", NA)))
analysis <- analysis[!is.na(analysis$region_group), , drop = FALSE]

# residualise on survey fixed effects, so only within-survey variation remains
residualise <- function(v, survey) {
  out <- rep(NA_real_, length(v)); ok <- is.finite(v)
  if (sum(ok) < 10) return(out)
  out[ok] <- residuals(lm(v[ok] ~ factor(survey[ok])))
  out
}
within_cor <- function(z, column) {
  if (!column %in% names(z)) return(NA_real_)
  a <- residualise(z$pfpr2_10, z$svkey)
  b <- residualise(as.numeric(z[[column]]), z$svkey)
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 20) return(NA_real_)
  cor(a[ok], b[ok])
}

marker_res <- do.call(rbind, lapply(GROUPS, function(g) {
  z <- analysis[analysis$region_group == g, , drop = FALSE]
  data.frame(region_group = g, n = nrow(z), surveys = length(unique(z$svkey)),
             marker = unname(MARKERS), variable = names(MARKERS),
             correlation = vapply(names(MARKERS), function(v) within_cor(z, v), 0),
             row.names = NULL)
}))
write.csv(marker_res, file.path(RESULTS_DIR, "deprivation_proxy_by_marker.csv"), row.names = FALSE)
cat("=== Within-survey correlation of PfPR2-10 with each deprivation marker ===\n")
cat("    negative means higher transmission goes with worse conditions\n\n")
wide <- reshape(marker_res[, c("region_group", "marker", "correlation")],
                idvar = "region_group", timevar = "marker", direction = "wide")
names(wide) <- sub("^correlation\\.", "", names(wide))
print(within(wide, for (nm in unname(MARKERS)) assign(nm, round(get(nm), 2))), row.names = FALSE)

## ---- composite deprivation index -------------------------------------------
mat <- analysis[, names(MARKERS)]
ok <- complete.cases(mat)
pc <- prcomp(scale(mat[ok, ]))
analysis$deprivation <- NA_real_
analysis$deprivation[ok] <- pc$x[, 1]
# orient so that higher = more deprived
if (cor(analysis$deprivation[ok], analysis$pct_urban[ok]) > 0) {
  analysis$deprivation <- -analysis$deprivation
}
cat(sprintf("\nComposite index: PC1 explains %.0f%% of the variance in the six markers.\n",
            100 * summary(pc)$importance[2, 1]))

fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
composite <- do.call(rbind, lapply(GROUPS, function(g) {
  z <- analysis[analysis$region_group == g & is.finite(analysis$deprivation), , drop = FALSE]
  r <- within_cor(z, "deprivation")
  data.frame(region_group = g, n = nrow(z), correlation = r, z = fisher_z(r))
}))
cat("\n=== Composite deprivation index versus PfPR2-10, within survey ===\n")
print(within(composite, { correlation <- round(correlation, 3); z <- round(z, 3) }), row.names = FALSE)
w <- composite[composite$region_group == "West Africa", ]
e <- composite[composite$region_group == "East & Southern", ]
zdiff <- (w$z - e$z) / sqrt(1 / (w$n - 3) + 1 / (e$n - 3))
cat(sprintf("West Africa versus East & Southern: Fisher z = %.2f, p = %.4g\n",
            zdiff, 2 * pnorm(-abs(zdiff))))

## ---- did the entanglement change over time? --------------------------------
analysis$era <- cut(analysis$year, c(1999.5, 2008.5, 2016.5, 2024.5),
                    labels = c("2000-2008", "2009-2016", "2017-2024"))
era_res <- do.call(rbind, lapply(c("West Africa", "East & Southern"), function(g) {
  do.call(rbind, lapply(levels(analysis$era), function(e) {
    z <- analysis[analysis$region_group == g & analysis$era == e &
                    is.finite(analysis$deprivation), , drop = FALSE]
    if (nrow(z) < 40) return(NULL)
    data.frame(region_group = g, era = e, n = nrow(z),
               surveys = length(unique(z$svkey)), correlation = within_cor(z, "deprivation"))
  }))
}))
write.csv(rbind(cbind(composite, era = "all"),
                cbind(era_res[, c("region_group", "n", "correlation")],
                      z = fisher_z(era_res$correlation), era = era_res$era)),
          file.path(RESULTS_DIR, "deprivation_proxy_composite.csv"), row.names = FALSE)
cat("\n=== Has the entanglement changed over time? ===\n")
print(within(era_res, correlation <- round(correlation, 3)), row.names = FALSE)
cat("\nThe entanglement is persistently high in West Africa but not rising, so it\n",
    "plausibly inflates the West African slope at all times without explaining the\n",
    "steepening reported in script 20.\n", sep = "")

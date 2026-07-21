# =============================================================================
# 17_region_map_check.R — CHECKPOINT before expanding the DHS panel.
# For the existing 44 survey-regions, extract MAP PfPR2-10 pop-weighted at each
# survey's admin-1 region + survey year, and compare to the survey-MEASURED
# prevalence (pfpr2_10). Strong region-level agreement => MAP is a valid
# prevalence substitute, so we can extend the panel to 2000 with MAP.
#
# Boundaries: rdhs::download_boundaries (DHS SDR); matched to component2 regions
# by normalised name. MAP: annual PfPR2-10 rasters (getRaster), pop-weighted with
# GPW. Both cached under data/ (git-ignored).
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(rdhs); library(malariaAtlas); library(terra); library(sf); library(ggplot2)})
options(rappdir_permission = TRUE)

BDIR <- file.path(DATA, "dhs_boundaries"); dir.create(BDIR, showWarnings = FALSE)
MDIR <- file.path(DATA, "map_annual");     dir.create(MDIR, showWarnings = FALSE)

d  <- read.csv(file.path(RESULTS, "component2_region_data.csv"), stringsAsFactors = FALSE)
sv <- unique(d[, c("svkey","survey","iso3","year")])
ds <- dhs_datasets(fileFormat = "FL"); ds$fn <- sub("\\..*$", "", ds$FileName)
sv$SurveyId <- ds$SurveyId[match(sv$survey, ds$fn)]

## ---- 1. MAP annual PfPR2-10 rasters for the survey years (cached) -----------
ext <- matrix(c(-18, -35, 52, 38), nrow = 2)                       # Africa bbox
get_map_year <- function(yr) {
  f <- file.path(MDIR, sprintf("pfpr2_10_%d.tif", yr))
  if (!file.exists(f)) {
    dl <- file.path(tempdir(), paste0("map", yr)); dir.create(dl, showWarnings = FALSE)
    invisible(tryCatch(getRaster(dataset_id = "Malaria__202508_Global_Pf_Parasite_Rate",
                                 year = yr, extent = ext, file_path = dl), error = function(e) NULL))
    tif <- list.files(dl, pattern = "\\.tiff?$", full.names = TRUE)
    if (!length(tif)) return(NULL)
    r <- terra::rast(tif); r <- r[[grep("_1$", names(r))[1]]]      # band 1 = mean estimate
    terra::writeRaster(r, f, overwrite = TRUE)
  }
  terra::rast(f)
}
years <- sort(unique(sv$year))
cat("MAP rasters for years:", paste(range(years), collapse = "-"), "\n")
den0 <- terra::rast(GPW_TIF)                                        # GPW density (pop weight)

## ---- 2. per survey: boundaries + pop-weighted MAP per region ----------------
region_map <- function(i) {
  s <- sv[i, ]
  bf <- file.path(BDIR, paste0(s$SurveyId, ".rds"))
  bd <- if (file.exists(bf)) readRDS(bf) else tryCatch({
    b <- download_boundaries(surveyId = s$SurveyId, method = "sf")[[1]]; saveRDS(b, bf); b
  }, error = function(e) NULL)
  if (is.null(bd) || !"DHSREGEN" %in% names(bd)) return(NULL)
  pf <- get_map_year(s$year); if (is.null(pf)) return(NULL)
  v  <- terra::makeValid(terra::vect(sf::st_make_valid(bd)))
  den <- terra::resample(terra::crop(den0, pf), pf, method = "bilinear")
  w   <- terra::mask(den, pf)
  an  <- terra::extract(pf * w, v, fun = sum, na.rm = TRUE, ID = FALSE)[[1]]
  aw  <- terra::extract(w,      v, fun = sum, na.rm = TRUE, ID = FALSE)[[1]]
  data.frame(svkey = s$svkey, regkey = rkey(bd$DHSREGEN), map_pfpr = 100 * an / aw,
             stringsAsFactors = FALSE)
}
rows <- list()
for (i in seq_len(nrow(sv))) {
  r <- tryCatch(region_map(i), error = function(e) { message(sv$svkey[i], ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) rows[[sv$svkey[i]]] <- r
  cat(sprintf("  %-8s %s: %s\n", sv$svkey[i], sv$SurveyId[i], if (is.null(r)) "FAILED" else paste(nrow(r), "regions")))
}
mp <- do.call(rbind, rows)

## ---- 3. merge with survey-measured prevalence & compare ---------------------
cmp <- merge(d[, c("svkey","regkey","region","iso3","year","pfpr2_10","exposure")], mp, by = c("svkey","regkey"))
cmp <- cmp[is.finite(cmp$pfpr2_10) & is.finite(cmp$map_pfpr), ]
write.csv(cmp, file.path(RESULTS, "region_map_vs_survey.csv"), row.names = FALSE)
cat(sprintf("\nMatched region-year pairs: %d (from %d surveys)\n", nrow(cmp), length(unique(cmp$svkey))))
cat(sprintf("Pearson r = %.3f, Spearman = %.3f; mean survey=%.1f%% MAP=%.1f%% diff=%+.1f pts\n",
            cor(cmp$pfpr2_10, cmp$map_pfpr), cor(cmp$pfpr2_10, cmp$map_pfpr, method = "spearman"),
            mean(cmp$pfpr2_10), mean(cmp$map_pfpr), mean(cmp$pfpr2_10 - cmp$map_pfpr)))
f <- lm(pfpr2_10 ~ map_pfpr, cmp); cat(sprintf("survey = %.2f + %.2f*MAP (R2=%.2f)\n", coef(f)[1], coef(f)[2], summary(f)$r.squared))

p <- ggplot(cmp, aes(map_pfpr, pfpr2_10)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(size = exposure), alpha = 0.4, colour = "#08519c") +
  scale_size_area(max_size = 4, guide = "none") +
  labs(x = "MAP region PfPR2-10 (%)", y = "Survey-measured region PfPR2-10 (%)",
       title = "Region-level agreement: MAP vs survey-measured prevalence") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(RESULTS, "region_map_vs_survey.png"), p, width = 7, height = 6.5, dpi = 300)
cat("saved: results/region_map_vs_survey.{png,csv}\n")

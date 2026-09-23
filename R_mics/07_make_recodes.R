#!/usr/bin/env Rscript
# Convert each MICS birth history into a DHS-shaped births file for the child age-band
# builder (R_cbh/R/child_bands.R). Only the fields the builder uses for eligibility, deaths,
# keys and region are essential: v008, b3, b5, b6/b7, v005, v001-v003 and v024.
# Codings are converted to DHS conventions:
#   b5  MICS 1 = alive, 2 = dead  ->  DHS 1 = alive, 0 = dead (otherwise every death is lost)
#   b0  MICS 1 = single, 2 = multiple -> DHS 0 = single, 1 = multiple
#   b6  100 * unit + number (unit 1 days, 2 months, 3 years), as in DHS
#   v005 women's weight x 1e6
#   v024 labelled with the analysis-region names used in the boundary polygons
suppressPackageStartupMessages({ library(haven); library(data.table) })
setup <- fread("data/derived_mics/survey_setup.csv")
rmap <- fread("data/derived_mics/region_map.csv")
root <- "data/MICS_extracted"; out <- "data/derived_mics/recodes"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
num <- function(x) as.numeric(zap_missing(zap_labels(x)))
NR <- 0L
nv <- function(x) if (is.null(x)) rep(NA_real_, NR) else num(x)   # absent field -> missing, not an error
getc <- function(d, ...) { low <- tolower(names(d)); for (a in tolower(c(...))) { j <- match(a, low); if (!is.na(j)) return(d[[j]]) }; NULL }
cmc_from <- function(mo, yr) { if (is.null(mo) || is.null(yr)) return(NULL); mo <- num(mo); yr <- num(yr)
  yr <- ifelse(yr < 100, yr + 1900, yr); ifelse(mo %in% 1:12 & yr > 1900 & yr < 2030, (yr - 1900) * 12 + mo, NA) }
log <- list()
for (i in seq_len(nrow(setup))) {
  s <- setup[i]; bh <- read_sav(file.path(root, s$folder, "bh.sav")); NR <- nrow(bh)
  wm <- read_sav(file.path(root, s$folder, "wm.sav"))
  # Woman keys: cluster, household, line number (MICS names vary by round).
  bk <- list(getc(bh, "HH1", "WM1"), getc(bh, "HH2", "WM2"), getc(bh, "LN", "WM4", "WM3"))
  wk <- list(getc(wm, "HH1", "WM1"), getc(wm, "HH2", "WM2"), getc(wm, "LN", "WM4", "WM3"))
  if (any(vapply(bk, is.null, TRUE))) {           # Mozambique 2008 uses a combined reference number
    ref <- getc(bh, "XX1"); stopifnot(!is.null(ref)); bk <- list(num(ref), 0, num(getc(bh, "MEMID")))
    wk <- list(num(getc(wm, "XX1")), 0, num(getc(wm, "MEMID", "LN")))
  }
  bkey <- paste(num(bk[[1]]), num(bk[[2]]), num(bk[[3]]))
  wkey <- if (!any(vapply(wk, is.null, TRUE))) paste(num(wk[[1]]), num(wk[[2]]), num(wk[[3]])) else NULL
  from_wm <- function(...) { v <- getc(wm, ...); if (is.null(v) || is.null(wkey)) return(NULL); v[match(bkey, wkey)] }
  pick <- function(...) { v <- getc(bh, ...); if (!is.null(v)) v else from_wm(...) }
  # Dates.
  dob <- getc(bh, "BH4C", "CCDOB"); dob <- if (!is.null(dob)) num(dob) else cmc_from(getc(bh, "BH4MC", "BH4M", "HN5_MES"), getc(bh, "BH4YC", "BH4Y", "HN5_ANO"))
  doi <- pick("WDOI", "WM6C", "CMCDOIW"); doi <- if (!is.null(doi)) num(doi) else cmc_from(pick("WM6M"), pick("WM6Y"))
  wdob <- pick("WDOB", "WM8C"); wdob <- if (!is.null(wdob)) num(wdob) else rep(NA_real_, nrow(bh))
  # Survival and age at death.
  alive <- num(getc(bh, "BH5", "HN6")); b5 <- ifelse(alive == 1, 1, ifelse(alive == 2, 0, NA))
  u <- nv(pick("BH9UC", "BH9U", "BH9A", "HN10_U")); n <- nv(pick("BH9NC", "BH9N", "BH9B", "HN10_I"))
  b6 <- ifelse(b5 == 0 & u %in% 1:3 & is.finite(n) & n < 97, 100 * u + n, ifelse(b5 == 0, 999, NA))
  b7 <- getc(bh, "BH9C"); b7 <- if (!is.null(b7)) num(b7) else
    ifelse(b5 == 0 & u == 1 & n < 97, floor(n / 30.4375), ifelse(b5 == 0 & u == 2 & n < 97, n, ifelse(b5 == 0 & u == 3 & n < 97, 12 * n, NA)))
  b7[b5 != 0 | is.na(b5)] <- NA
  # Weight, urban, wealth, strata.
  wt <- nv(pick("wmweight", "WMWEIGHT")); urban <- nv(pick("HH6")); wealth <- nv(pick("windex5", "WLTHIND5"))
  if (all(is.na(urban)) && grepl("Dakar", s$folder)) urban <- rep(1, NR)   # Dakar city survey: urban by design
  strat <- pick("stratum", "HH7_stratum", "WMSTRAT", "HHSTRAT")
  # Region: map MICS labels to analysis regions.
  m <- rmap[svkey == s$svkey]
  if (startsWith(s$region_var, "CONST:")) lab <- rep(sub("CONST:", "", s$region_var), nrow(bh)) else {
    rv <- if (startsWith(s$region_var, "WM:")) from_wm(sub("^WM:", "", s$region_var)) else pick(s$region_var)
    lab <- as.character(as_factor(rv)) }
  ar <- m$analysis_region[match(lab, m$mics_label)]
  levels_ar <- sort(unique(m$analysis_region))
  v024 <- labelled(match(ar, levels_ar), setNames(seq_along(levels_ar), levels_ar))
  # Birth order and preceding interval from the birth dates within each mother.
  ord <- order(bkey, dob, seq_len(nrow(bh)))
  bord <- integer(nrow(bh)); b11 <- rep(NA_real_, nrow(bh))
  grp <- bkey[ord]; d <- dob[ord]
  bord[ord] <- ave(seq_along(grp), grp, FUN = seq_along)
  prev <- ave(d, grp, FUN = function(z) c(NA, head(z, -1))); b11[ord] <- d - prev
  twin <- nv(getc(bh, "BH2", "HN3")); b0 <- ifelse(twin == 2, 1, 0)
  sex <- nv(getc(bh, "BH3", "HN4"))
  bidx <- ave(seq_len(nrow(bh)), bkey, FUN = seq_along)
  st <- if (!is.null(strat)) num(strat) else as.numeric(factor(paste(match(ar, levels_ar), urban)))
  rec <- data.frame(caseid = bkey, v001 = num(bk[[1]]), v002 = num(bk[[2]]), v003 = num(bk[[3]]),
    v005 = round(wt * 1e6), v008 = doi, v011 = wdob, v021 = num(bk[[1]]), v022 = st, v023 = st,
    v025 = urban, v133 = NA_real_, v190 = wealth, bidx = bidx, bord = bord, b0 = b0, b3 = dob, b4 = sex,
    b5 = b5, b6 = b6, b7 = b7, b10 = NA_real_, b11 = b11, b13 = NA_real_)
  rec$v024 <- v024
  saveRDS(rec, file.path(out, paste0(s$svkey, ".rds")))
  log[[s$svkey]] <- data.table(svkey = s$svkey, births = nrow(rec), deaths = sum(rec$b5 == 0, na.rm = TRUE),
    missing_region = sum(is.na(ar)), missing_dob = sum(is.na(rec$b3)), missing_doi = sum(is.na(rec$v008)),
    missing_weight = sum(!is.finite(rec$v005) | rec$v005 <= 0), regions = length(levels_ar),
    interview_year = as.integer(1900 + floor((median(rec$v008, na.rm = TRUE) - 1) / 12)))
}
L <- rbindlist(log); fwrite(L, "data/derived_mics/recode_log.csv")
print(L[, .(svkey, births, deaths, missing_region, missing_dob, missing_doi, missing_weight, regions)], nrows = 60)

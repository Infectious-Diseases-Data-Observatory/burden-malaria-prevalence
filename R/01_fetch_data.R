# =============================================================================
# 01_fetch_data.R — fetch/prepare all raw inputs into data/ (idempotent: each
# block is skipped if its output already exists). Run from the repo root.
#
# Access-controlled inputs a collaborator must supply themselves (see README):
#   * data/ihme_malaria_u5_deaths_by_country.csv  — IHME/GBD export (free login)
#   * data/dhs/*.rds                              — DHS microdata (free DHS login)
# =============================================================================
source("R/00_utils.R")
dir.create(file.path(DATA, "pop"), showWarnings = FALSE, recursive = TRUE)

## ---- 1. MAP PfPR2-10 raster (Africa, 2024) ---------------------------------
if (!file.exists(PFPR_TIF)) {
  message("Downloading MAP PfPR2-10 raster (2024)...")
  suppressMessages(library(malariaAtlas))
  dl <- file.path(tempdir(), "map_dl"); dir.create(dl, showWarnings = FALSE)
  ext <- matrix(c(-18, -35, 52, 38), nrow = 2)                       # Africa bbox
  invisible(tryCatch(getRaster(dataset_id = "Malaria__202508_Global_Pf_Parasite_Rate",
                               year = MAP_YEAR, extent = ext, file_path = dl),
                     error = function(e) message("note: ", conditionMessage(e))))
  tif <- list.files(dl, pattern = "\\.tiff?$", full.names = TRUE)
  tif <- tif[which.max(file.info(tif)$mtime)]
  terra::writeRaster(terra::rast(tif)[[1]], PFPR_TIF, overwrite = TRUE)   # band 1 = mean estimate
}

## ---- 2. GPW population-density raster (2020, 2.5 arc-min) -------------------
if (!file.exists(GPW_TIF)) {
  message("Downloading GPW population raster...")
  suppressMessages(library(geodata))
  geodata::population(year = 2020, res = 2.5, path = file.path(DATA, "pop"))
  hit <- list.files(file.path(DATA, "pop"), pattern = "gpw.*2020.*2\\.5m\\.tif$",
                    recursive = TRUE, full.names = TRUE)[1]
  if (!is.na(hit) && normalizePath(hit) != normalizePath(GPW_TIF)) file.copy(hit, GPW_TIF, overwrite = TRUE)
}

## ---- 3. Country boundaries (admin-0, African countries) --------------------
ADM0 <- file.path(DATA, "africa_admin0.rds")
if (!file.exists(ADM0)) {
  message("Building admin-0 boundaries...")
  suppressMessages(library(geodata))
  w <- geodata::world(path = file.path(DATA, "pop"))
  w$iso <- w$GID_0
  saveRDS(sf::st_as_sf(w[, "iso"]), ADM0)
}

## ---- 4. IGME child mortality (U5MR, neonatal, infant) via World Bank -------
IGME <- file.path(DATA, "igme_mortality_by_country.csv")
if (!file.exists(IGME)) {
  message("Fetching IGME mortality (World Bank)...")
  u5 <- wb_latest("SH.DYN.MORT");    names(u5)[2:3]  <- c("u5mr_year", "u5mr")     # per 1,000 live births
  nn <- wb_latest("SH.DYN.NMRT")[, c("iso3","value")]; names(nn)[2] <- "nmr"       # neonatal
  im <- wb_latest("SP.DYN.IMRT.IN")[, c("iso3","value")]; names(im)[2] <- "imr"    # infant
  m  <- merge(merge(u5, nn, by = "iso3"), im, by = "iso3")
  m$m_1mo_5y <- m$u5mr - m$nmr                                                      # deaths 1mo-5y / 1,000 lb
  write.csv(m, IGME, row.names = FALSE)
}

## ---- 5. Live births (crude birth rate x population) ------------------------
BIRTHS <- file.path(DATA, "wb_livebirths_by_country.csv")
if (!file.exists(BIRTHS)) {
  message("Fetching live births (World Bank)...")
  cbrt <- wb_latest("SP.DYN.CBRT.IN"); pop <- wb_latest("SP.POP.TOTL")
  b <- merge(cbrt[, c("iso3","value")], pop[, c("iso3","value")], by = "iso3")
  b$births <- b$value.x / 1000 * b$value.y
  write.csv(b[, c("iso3","births")], BIRTHS, row.names = FALSE)
}

## ---- 6. GDP per capita (country x year) — NEW covariate --------------------
GDP <- file.path(DATA, "wb_gdp_pc.csv")
if (!file.exists(GDP)) {
  message("Fetching GDP per capita (World Bank NY.GDP.PCAP.CD)...")
  g <- wb_fetch("NY.GDP.PCAP.CD"); names(g)[3] <- "gdp_pc"
  write.csv(g, GDP, row.names = FALSE)
}

## ---- 7. DTP3 vaccine coverage (country x year) — NEW covariate -------------
DTP3 <- file.path(DATA, "wb_dtp3.csv")
if (!file.exists(DTP3)) {
  message("Fetching DTP3 coverage (World Bank SH.IMM.IDPT)...")
  v <- wb_fetch("SH.IMM.IDPT"); names(v)[3] <- "dtp3"
  write.csv(v, DTP3, row.names = FALSE)
}

## ---- 8. DHS usable-survey universe (malaria biomarker + births recode) -----
UNIV <- file.path(DATA, "ssa_usable_dhs_surveys.csv")
if (!file.exists(UNIV)) {
  message("Building DHS usable-survey universe...")
  api <- function(p) fromJSON(paste0("https://api.dhsprogram.com/rest/dhs/", p))
  pg  <- function(base) { f <- api(base); o <- list(f$Data)
    if (f$TotalPages > 1) for (p in 2:f$TotalPages) o[[p]] <- api(paste0(base, "&page=", p))$Data
    do.call(rbind, o) }
  mal <- pg("data?indicatorIds=ML_PMAL_C_RDT&perpage=3000&f=json")   # surveys measuring child malaria (RDT)
  mal$iso3   <- countrycode(mal$CountryName, "country.name", "iso3c", warn = FALSE)
  mal$region <- countrycode(mal$iso3, "iso3c", "region", warn = FALSE)
  mal <- mal[!is.na(mal$region) & mal$region == "Sub-Saharan Africa", ]
  ccs <- unique(mal$DHS_CountryCode)                                 # DHS 2-letter country codes
  ds  <- pg(sprintf("datasets?countryIds=%s&perpage=3000&f=json", paste(ccs, collapse = ",")))
  ds$recode <- toupper(substr(sub("[.].*$", "", ds$FileName), 3, 4))
  msurv  <- unique(mal[, c("DHS_CountryCode","CountryName","SurveyId","SurveyYear","SurveyType")])
  usable <- msurv[msurv$SurveyId %in% unique(ds$SurveyId[ds$recode == "BR"]), ]   # + has births recode
  write.csv(usable[order(usable$CountryName, usable$SurveyYear), ], UNIV, row.names = FALSE)
}

## ---- 9. DHS microdata (PR + BR flat recodes) via rdhs -----------------------
# Requires a configured rdhs login cached in ~/.rdhs.json (see README). Downloads
# only the files not already present. NOTE: do NOT call set_rdhs_config() here
# (it triggers a password prompt); library(rdhs) auto-loads the cached config.
fetch_dhs_microdata <- function() {
  suppressMessages(library(rdhs)); options(rappdir_permission = TRUE)
  usable <- read.csv(UNIV, stringsAsFactors = FALSE)
  api <- function(p) fromJSON(paste0("https://api.dhsprogram.com/rest/dhs/", p))
  ds  <- api(sprintf("datasets?countryIds=%s&perpage=3000&f=json",
                     paste(unique(usable$DHS_CountryCode), collapse = ",")))$Data
  fn  <- toupper(sub("[.].*$", "", ds$FileName))
  want <- fn[ds$SurveyId %in% usable$SurveyId & substr(fn,3,4) %in% c("PR","BR") & substr(fn,7,8) == "FL"]
  have <- toupper(sub("[.]rds$", "", list.files(file.path(DATA,"dhs"), pattern = "rds$")))
  todo <- setdiff(unique(want), have)
  if (!length(todo)) { message("DHS microdata: all present."); return(invisible()) }
  message("Downloading ", length(todo), " DHS recode files...")
  paths <- get_datasets(dataset_filenames = todo, download_option = "rds", reformat = TRUE, clear_cache = FALSE)
  for (nm in names(paths)) if (is.character(paths[[nm]]) && file.exists(paths[[nm]]))
    file.copy(paths[[nm]], file.path(DATA, "dhs", paste0(nm, ".rds")), overwrite = TRUE)
}
if (!dir.exists(file.path(DATA, "dhs")) ||
    length(list.files(file.path(DATA, "dhs"), pattern = "BR.*rds$")) < 40) {
  tryCatch(fetch_dhs_microdata(),
           error = function(e) message("DHS download skipped (", conditionMessage(e),
                                        "). Configure rdhs and re-run — see README."))
}

## ---- check the two user-supplied, access-controlled inputs ------------------
for (f in c(file.path(DATA, "ihme_malaria_u5_deaths_by_country.csv"))) {
  if (!file.exists(f)) message("MISSING (user must supply, see README): ", f)
}
message("01_fetch_data.R complete.")

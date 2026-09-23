#!/usr/bin/env Rscript
# Map each MICS survey's region labels to boundary polygons. Writes
# data/derived_mics/region_map_auto.csv (one row per survey x MICS label) for review.
suppressPackageStartupMessages({ library(haven); library(sf); library(data.table) })
sf_use_s2(FALSE)
source("R_cbh/R/geography.R")   # cbh_rkey()
snow <- "/Users/jameswatson/Documents/Claude Projects/Prevalence Model/Snow Prevalence Data/shape files/Africa_New_Admin.shp"
# survey folder, svkey, region variable (in bh.sav unless noted), boundary source
S <- fread(text = "
folder|svkey|region_var|boundary
BEN_2014_MICS5_v01_M|MC_BEN2014|HH7|BJ2017DHS
BEN_2021_MICS6_v01_M|MC_BEN2021|HH7|BJ2017DHS
CAF_2018_MICS6_v01_M|MC_CAF2018|HH7|CAF_analysis
CIV_2016_MICS5_v01_M|MC_CIV2016|HH7|CI2012DHS
CMR_2014_MICS5_v01_M|MC_CMR2014|HH7|CM2018DHS
COD_2017_MICS6_v01_M|MC_COD2017|HH7|CD2023DHS
COG_2014_MICS5_v01_M|MC_COG2014|HH7|CG2011DHS
COM_2022_MICS6_v01_M|MC_COM2022|HH7|KM2012DHS
GHA_2011_MICS4_v01_M|MC_GHA2011|HH7|GH2014DHS
GHA_2017_MICS6_v01_M|MC_GHA2017|HH7|GH2014DHS
GIN_2016_MICS5_v01_M|MC_GIN2016|HH7|GN2018DHS
GMB_2018_MICS6_v01_M|MC_GMB2018|HH7|GM2019DHS
GNB_2014_MICS5_v01_M|MC_GNB2014|HH7|GNB_analysis
GNB_2018_MICS6_v01_M|MC_GNB2018|HH7|GNB_analysis
KEN(Bungoma County)_2013_MICS5_v01_M|MC_KEN2013BUN|CONST:Bungoma|KE2022DHS
KEN(Kakamega County)_2013_MICS5_v01_M|MC_KEN2013KAK|CONST:Kakamega|KE2022DHS
KEN(Mombasa Informal Settlements)_2009_MICS4_v01_M|MC_KEN2009MOM|CONST:Mombasa|KE2022DHS
KEN(Nyanza Province)_2011_MICS4_v01_M|MC_KEN2011NYA|HH7|KE2022DHS
KEN(Turkana County)_2013_MICS5_v01_M|MC_KEN2013TUR|CONST:Turkana|KE2022DHS
LSO_2018_MICS6_v01_M|MC_LSO2018|HH7A|LS2014DHS
MDG_2018_MICS6_v01_M|MC_MDG2018|HH7|MD2008DHS
MDG(South)_2012_MICS4_v01_M|MC_MDG2012S|HH7|MD2008DHS
MLI_2015_MICS5_v01_M|MC_MLI2015|HH7|ML2018DHS
MOZ_2008_MICS3_v01_M|MC_MOZ2008|HH7|MZ2011DHS
MRT_2011_MICS4_v01_M|MC_MRT2011|HH7|MR2020DHS
MRT_2015_MICS5_v01_M|MC_MRT2015|HH7|MR2020DHS
MWI_2006_MICS3_v01_M|MC_MWI2006|HHREG|MW2015DHS
MWI_2013_MICS5_v01_M|MC_MWI2013|WM:region|MW2015DHS
MWI_2019_MICS6_v01_M|MC_MWI2019|HH7|MW2015DHS
NGA_2016_MICS5_v01_M|MC_NGA2016|Zone|NG2018DHS
NGA_2021_MICS6_v01_M|MC_NGA2021|zone|NG2018DHS
SEN(Dakar City)_2015_MICS5_v01_M|MC_SEN2015DKR|CONST:Dakar & Thies|SNOW:SEN
SLE_2017_MICS6_v01_M|MC_SLE2017|HH7|SL2013DHS
SOM_2006_MICS3_v01_M|MC_SOM2006|HH7|SNOW:SOM
SOM(Northeast Zone)_2011_MICS4_v01_M|MC_SOM2011NE|HH7|SNOW:SOM
SOM(Somaliland)_2011_MICS4_v01_M|MC_SOM2011SL|HH7|SNOW:SOM
SSD_2010_MICS4_v01_M|MC_SSD2010|HH7|SNOW:SSD
STP_2014_MICS5_v01_M|MC_STP2014|HH7|ST2008DHS
STP_2019_MICS6_v01_M|MC_STP2019|HH7|ST2008DHS
SWZ_2010_MICS4_v01_M|MC_SWZ2010|HH7|SZ2006DHS
SWZ_2014_MICS5_v01_M|MC_SWZ2014|HH7|SZ2006DHS
SWZ_2021_MICS6_v01_M|MC_SWZ2021|HH7|SZ2006DHS
TCD_2019_MICS6_v01_M|MC_TCD2019|HH7|TD2014DHS
TGO_2017_MICS6_v01_M|MC_TGO2017|HH7|TG2013DHS
ZWE_2009_MICS3_v01_M|MC_ZWE2009|hh7|ZW2015DHS
ZWE_2014_MICS5_v01_M|MC_ZWE2014|HH7|ZW2015DHS
ZWE_2019_MICS6_v01_M|MC_ZWE2019|HH7|ZW2015DHS
", sep = "|")
stopifnot(nrow(S) == 47, !anyDuplicated(S$svkey), all(grepl("^[A-Za-z0-9_-]+$", S$svkey)))
fwrite(S, "data/derived_mics/survey_setup.csv")

boundary_names <- function(b) {
  if (b %in% c("CAF_analysis", "GNB_analysis"))
    return(readRDS(sprintf("data/derived_mics/boundaries/%s_analysis_regions.rds", sub("_analysis", "", b)))$analysis_region)
  if (startsWith(b, "SNOW:")) { a <- st_read(snow, quiet = TRUE); return(a$AFR_Admin[a$Country_ID == sub("SNOW:", "", b)]) }
  as.data.frame(readRDS(sprintf("data/dhs_boundaries/%s.rds", b)))$DHSREGEN
}
mics_labels <- function(folder, var) {
  if (startsWith(var, "CONST:")) return(sub("CONST:", "", var))
  f <- if (startsWith(var, "WM:")) "wm.sav" else "bh.sav"; v <- sub("^WM:", "", var)
  d <- read_sav(file.path("data/MICS_extracted", folder, f), col_select = any_of(v))
  if (!ncol(d)) d <- read_sav(file.path("data/MICS_extracted", folder, "wm.sav"), col_select = any_of(v))
  stopifnot(ncol(d) == 1)
  sort(unique(as.character(as_factor(d[[1]]))))
}
rows <- list()
for (i in seq_len(nrow(S))) {
  bn <- boundary_names(S$boundary[i]); lk <- mics_labels(S$folder[i], S$region_var[i])
  bk <- cbh_rkey(bn)
  for (l in lk) {
    j <- match(cbh_rkey(l), bk)
    rows[[length(rows) + 1]] <- data.table(svkey = S$svkey[i], boundary = S$boundary[i], mics_label = l,
      mics_key = cbh_rkey(l), polygon = if (is.na(j)) NA_character_ else bn[j], method = if (is.na(j)) "unmatched" else "exact")
  }
}
m <- rbindlist(rows)
fwrite(m, "data/derived_mics/region_map_auto.csv")
cat("labels:", nrow(m), "| exact:", sum(m$method == "exact"), "| unmatched:", sum(m$method == "unmatched"), "\n")

# Apply the reviewed overrides (translations, spellings and merges) to the unmatched labels and
# write the final map: one row per survey x MICS label, with its analysis region and the boundary
# polygon(s) whose union defines that region.
ov <- fread("R_mics/config/region_overrides.csv", encoding = "UTF-8")
k <- paste(m$svkey, m$mics_label); ko <- paste(ov$svkey, ov$mics_label)
stopifnot(!anyDuplicated(ko), all(ko %in% k))
j <- match(k, ko)
final <- data.table(svkey = m$svkey, boundary = m$boundary, mics_label = m$mics_label,
  analysis_region = ifelse(is.na(j), m$polygon, ov$analysis_region[j]),
  polygons = ifelse(is.na(j), m$polygon, ov$polygons[j]),
  method = ifelse(is.na(j), m$method, "override"), note = ifelse(is.na(j), "", ov$note[j]))
if (any(is.na(final$analysis_region))) stop("Unresolved MICS region labels: ",
  paste(final[is.na(analysis_region), paste(svkey, mics_label)], collapse = "; "))
fwrite(final, "data/derived_mics/region_map.csv")
cat("final map:", nrow(final), "labels |", sum(final$method == "override"), "by override\n")

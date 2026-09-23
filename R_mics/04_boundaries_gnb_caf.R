#!/usr/bin/env Rscript
# Analysis-region polygons for the two MICS countries with no DHS boundary file.
# Source: the local Africa admin-1 shapefile from the Snow prevalence project (path in
# R_mics/config/paths.R; override with MICS_ADMIN1_SHP). It merges
# each capital with a neighbour (Bissau with Biombo; Bangui with Ombella-M'Poko), so the
# corresponding MICS regions are merged into one analysis region:
#   Guinea-Bissau: 9 MICS regions -> 8 analysis regions (SAB + Biombo combined).
#   CAR: 7 MICS regions -> 6 analysis regions (Region 1 + Region 7/Bangui combined).
# CAR's MICS regions are unnamed ("Region 1".."Region 7"). Their prefecture composition follows
# the country's seven health regions; Region 7 is 100% urban in the MICS data, consistent with Bangui.
suppressPackageStartupMessages({ library(sf); library(dplyr) })
sf_use_s2(FALSE)
source("R_mics/config/paths.R")   # mics_admin1_shp
src <- mics_admin1_shp
stopifnot(file.exists(src))
a <- st_make_valid(st_read(src, quiet = TRUE))
out <- "data/derived_mics/boundaries"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

gnb_map <- c("Tombali" = "Tombali", "Quinara" = "Quinara", "Oio" = "Oio", "Bolama/bijagos" = "Bolama/Bijagos",
             "Bafata" = "Bafata", "Gabu" = "Gabu", "Cacheu" = "Cacheu",
             "Sector Autonomo De Bissau & Biombo" = "SAB and Biombo")
g <- a[a$Country_ID == "GNB", ]; stopifnot(nrow(g) == 8, setequal(g$AFR_Admin, names(gnb_map)))
g$analysis_region <- unname(gnb_map[g$AFR_Admin])
g <- g %>% group_by(analysis_region) %>% summarise(n_units = n(), .groups = "drop") %>% mutate(iso3 = "GNB")

caf_map <- c("Ombella-Mpoko & Bangui" = "Regions 1 and 7", "Lobaye" = "Regions 1 and 7",
             "Mambere-kadei" = "Region 2", "Nana Mambere" = "Region 2", "Sangha Mbaere" = "Region 2",
             "Ouham-pende" = "Region 3", "Ouham" = "Region 3",
             "Kemo" = "Region 4", "Nana Grebizi" = "Region 4", "Ouaka" = "Region 4",
             "Hautte-kotto" = "Region 5", "Bamingui-bangora" = "Region 5", "Vakaga" = "Region 5",
             "Basse Kotto" = "Region 6", "Mbomou" = "Region 6", "Haut-mboumou" = "Region 6")
k <- a[a$Country_ID == "CAF", ]; stopifnot(nrow(k) == 16, setequal(k$AFR_Admin, names(caf_map)))
k$analysis_region <- unname(caf_map[k$AFR_Admin])
k <- k %>% group_by(analysis_region) %>% summarise(n_units = n(), .groups = "drop") %>% mutate(iso3 = "CAF")

for (x in list(g, k)) {
  x <- st_make_valid(x); stopifnot(all(st_is_valid(x)))
  saveRDS(x, file.path(out, paste0(x$iso3[1], "_analysis_regions.rds")))
  cat(x$iso3[1], ":", nrow(x), "analysis regions |", paste(x$analysis_region, collapse = " | "), "\n")
}
# Which MICS HH7 labels belong to each analysis region (used when building region keys).
key <- rbind(
  data.frame(iso3 = "GNB", mics_region = c("Tombali","Quinara","Oio","Bolama/Bijagós","Bafatá","Gabú","Cacheu","SAB","Biombo"),
             analysis_region = c("Tombali","Quinara","Oio","Bolama/Bijagos","Bafata","Gabu","Cacheu","SAB and Biombo","SAB and Biombo")),
  data.frame(iso3 = "CAF", mics_region = paste("Région", 1:7),
             analysis_region = c("Regions 1 and 7", paste("Region", 2:6), "Regions 1 and 7")))
write.csv(key, file.path(out, "mics_region_to_analysis_region_gnb_caf.csv"), row.names = FALSE)

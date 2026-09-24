#!/usr/bin/env Rscript
# Serial maps of SMC coverage by admin-1 unit, 2015-2022, for the SMC countries in the
# analysis. Confirmed units (district- or region-level records) are shaded by the share
# of the unit's DHS/MIS clusters in areas with SMC by that year; units whose switch-on is
# known only from national-scope records are shown separately once the national
# programme has started (the codebook: coverage of the area is not recorded).
# Polygons: Nigerian states and Burkina Faso regions from the Snow admin-1 shapefile;
# other countries from the DHS boundary file whose regions best match the units.
suppressPackageStartupMessages({ library(sf); library(data.table); library(ggplot2) })
source("R_cbh/load_pipeline.R")
source("R_cbh/sensitivity/smc/settings.R")
source("R_mics/config/paths.R")
sf_use_s2(FALSE)
st <- cbh_smc_settings(); years <- 2015:2022
units <- fread(file.path(st$private, "unit_smc.csv")); coverage <- fread(file.path(st$private, "unit_smc_coverage.csv"))
countries <- sort(unique(units$iso3))
snow <- st_make_valid(st_read(mics_admin1_shp, quiet = TRUE))
snow_alias <- c(fct = "federalcapitalterritory", fctabuja = "federalcapitalterritory", nassarawa = "nasarawa",
  bugishu = "bugisuelgonfrom2024", busoga = "eastcentralbusoga", central1 = "central1southbuganda", central2 = "central2northbuganda")
key <- function(x) { k <- cbh_rkey(x); ifelse(k %in% names(snow_alias), snow_alias[k], k) }
registry <- fread("data/derived_dhs/survey_registry.csv")
polys <- lapply(countries, function(iso) {
  keys <- units$regkey[units$iso3 == iso]
  cand <- list()
  s <- snow[snow$Country_ID == iso, ]; cand[["snow"]] <- st_sf(regkey = key(s$AFR_Admin), geometry = st_geometry(s))
  for (f in registry$boundary_file[registry$iso3 == iso & file.exists(registry$boundary_file)]) {
    b <- st_make_valid(readRDS(f)); cand[[basename(f)]] <- st_sf(regkey = key(b$DHSREGEN), geometry = st_geometry(b))
  }
  score <- vapply(cand, function(p) sum(unique(p$regkey) %in% keys), 0)
  best <- cand[[which.max(score)]]
  best <- aggregate(best, list(regkey = best$regkey), function(x) x[1])[, c("regkey", "geometry")]
  best$iso3 <- iso; best$source <- names(cand)[which.max(score)]
  message(iso, ": ", best$source[1], ", ", sum(best$regkey %in% keys), " of ", nrow(best), " polygons matched to ", length(keys), " units")
  best
})
polys <- do.call(rbind, polys)
st_crs(polys) <- 4326
# Status by unit and year.
grid <- CJ(i = seq_len(nrow(polys)), year = years)
grid[, `:=`(iso3 = polys$iso3[i], regkey = polys$regkey[i])]
u <- match(paste(grid$iso3, grid$regkey), paste(units$iso3, units$regkey))
cv <- match(paste(grid$iso3, grid$regkey, grid$year), paste(coverage$iso3, coverage$regkey, coverage$year))
grid[, `:=`(basis = units$status_basis[u], national_year = units$first_year_main[u], cov = coverage$coverage_confirmed[cv])]
bins <- c("No SMC", "1–24%", "25–49%", "50–74%", "75–100%")
grid[, status := fcase(
  is.na(u), "No SMC record for this area",
  basis == "national" & !is.na(national_year) & year >= national_year, "National programme, area coverage not recorded",
  basis == "national", "No SMC",
  is.na(cov) | cov <= 0, "No SMC",
  cov < .25, "1–24%", cov < .5, "25–49%", cov < .75, "50–74%", default = "75–100%")]
levels_status <- c(bins, "National programme, area coverage not recorded", "No SMC record for this area")
grid[, status := factor(status, levels = levels_status)]
map <- cbind(polys[grid$i, c("iso3", "regkey")], grid[, .(year, status)])
map <- st_as_sf(map)
background <- snow[!snow$Country_ID %in% countries, ]
bbox <- st_bbox(polys)
outline <- aggregate(polys["iso3"], list(country = polys$iso3), function(x) x[1])
fill <- setNames(c("#F4F1EA", "#CFE8D8", "#8CC7A1", "#3E9A68", "#0B5D36", "#E3A857", "#D9D9D9"), levels_status)
p <- ggplot() +
  geom_sf(data = background, fill = "grey93", colour = NA) +
  geom_sf(data = map, aes(fill = status), colour = "white", linewidth = .08) +
  geom_sf(data = outline, fill = NA, colour = "grey25", linewidth = .25) +
  facet_wrap(~year, ncol = 3) +
  scale_fill_manual(values = fill, drop = FALSE, name = "Share of DHS clusters in areas\nwith SMC (confirmed records)") +
  coord_sf(xlim = c(bbox["xmin"] - .5, bbox["xmax"] + .5), ylim = c(bbox["ymin"] - .5, bbox["ymax"] + .5), expand = FALSE) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_void(base_size = 15) + theme(strip.text = element_text(face = "bold", size = 15, margin = margin(4, 0, 4, 0)),
    legend.position = "inside", legend.position.inside = c(.835, .16), legend.title = element_text(size = 14),
    legend.text = element_text(size = 13), legend.key.size = grid::unit(18, "pt"),
    panel.spacing = grid::unit(8, "pt"), plot.background = element_rect(fill = "white", colour = NA))
figure <- file.path(st$out, "smc_coverage_maps_2015_2022.png")
ggsave(figure, p, width = 16, height = 10.5, dpi = 220, device = ragg::agg_png, bg = "white")
tab <- grid[, .(units = .N), by = .(year, iso3, status)][order(year, iso3, status)]
cbh_atomic_csv(as.data.frame(tab), file.path(st$out, "smc_coverage_map_status.csv"))
unmatched <- units[!paste(iso3, regkey) %in% paste(polys$iso3, polys$regkey), .(iso3, regkey, status_basis, first_year_main)]
cbh_atomic_csv(as.data.frame(unmatched), file.path(st$out, "smc_units_without_map_polygon.csv"))
writeLines(c("# SMC coverage maps, 2015–2022", "",
  paste0("Admin-1 units of the ", length(countries), " SMC countries in the cluster file (", paste(countries, collapse = ", "), "): Nigerian states and Burkina Faso regions (Snow admin-1 shapefile); other countries' DHS survey regions (the DHS boundary file whose regions best match the units). ",
  "Confirmed units (district- or region-level campaign records: Nigeria, Burkina Faso, Uganda, Cameroon North and Far North) are shaded by the share of the unit's DHS/MIS clusters, pooled over survey rounds, located in areas where SMC had started by that year. ",
  "Units whose switch-on is known only from national-scope records (Mali, Niger, Chad, Côte d'Ivoire, the rest of Cameroon) are shown in amber from the national programme's first year: the country ran SMC, but which areas were covered is not recorded. ",
  "Grey areas inside the SMC countries have no cluster in the file. Other countries with SMC programmes (for example Senegal, The Gambia, Guinea, Ghana) are not in the file and are not shaded. Cluster shares approximate, but are not, population coverage; campaign cycles and eligible ages are not shown."),
  "", sprintf("Units without a polygon: %d (listed in smc_units_without_map_polygon.csv).", nrow(unmatched))),
  file.path(st$out, "CAPTION_smc_coverage_maps.md"))
message("SMC coverage maps written: ", figure)

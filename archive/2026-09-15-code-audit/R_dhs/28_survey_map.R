# =============================================================================
# 28_survey_map.R
# Where the surveys are and when they were fielded.
#
# Country outlines are dissolved from the DHS boundary files the analysis
# already uses, so the map needs no external basemap package and shows exactly
# the geography the panel is built on.
#
# Outputs
#   results/dhs_rebuild/figure10_survey_map.png
#   results/dhs_rebuild/survey_map_country_summary.csv
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("sf", "ggplot2", "patchwork"))

analysis <- read_analysis_data()
registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
registry <- registry[registry$svkey %in% analysis$svkey, , drop = FALSE]

surveys <- unique(analysis[, c("svkey", "iso3", "year")])
surveys$type <- registry$SurveyType[match(surveys$svkey, registry$svkey)]
surveys$type[is.na(surveys$type)] <- "DHS"

regions_per_survey <- as.data.frame(table(analysis$svkey),
                                    stringsAsFactors = FALSE)
names(regions_per_survey) <- c("svkey", "regions")
surveys$regions <- regions_per_survey$regions[
  match(surveys$svkey, regions_per_survey$svkey)
]

## ---- country outlines, dissolved from the DHS boundaries --------------------
# The most recent survey's boundary gives the fullest coverage of the country.
country_outline <- function(iso3) {
  rows <- registry[registry$iso3 == iso3, , drop = FALSE]
  rows <- rows[order(-rows$year), , drop = FALSE]
  for (i in seq_len(nrow(rows))) {
    path <- file.path(BOUNDARY_DIR, paste0(rows$SurveyId[i], ".rds"))
    if (!file.exists(path)) next
    shape <- tryCatch(sf::st_make_valid(readRDS(path)), error = function(e) NULL)
    if (is.null(shape) || !nrow(shape)) next
    dissolved <- tryCatch(sf::st_union(sf::st_geometry(shape)),
                          error = function(e) NULL)
    if (is.null(dissolved)) next
    return(sf::st_sf(iso3 = iso3, geometry = sf::st_sfc(dissolved),
                     crs = sf::st_crs(shape)))
  }
  NULL
}

message("Dissolving country outlines from DHS boundaries ...")
outlines <- lapply(sort(unique(surveys$iso3)), country_outline)
outlines <- outlines[!vapply(outlines, is.null, logical(1))]
outlines <- do.call(rbind, lapply(outlines, function(x) {
  sf::st_transform(x, 4326)
}))

summary_by_country <- do.call(rbind, lapply(split(surveys, surveys$iso3),
  function(d) data.frame(
    iso3 = d$iso3[1], surveys = nrow(d),
    first_year = min(d$year), last_year = max(d$year),
    region_years = sum(d$regions), stringsAsFactors = FALSE
  )))
outlines <- merge(outlines, summary_by_country, by = "iso3")

centroids <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(
  sf::st_geometry(outlines)
)))
outlines$lon <- centroids[, 1]
outlines$lat <- centroids[, 2]

write.csv(summary_by_country[order(-summary_by_country$surveys), ],
          file.path(RESULTS_DIR, "survey_map_country_summary.csv"),
          row.names = FALSE)

## ---- panel A: the map -------------------------------------------------------
map_panel <- ggplot2::ggplot(outlines) +
  ggplot2::geom_sf(ggplot2::aes(fill = surveys), colour = "white",
                   linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(x = lon, y = lat, label = iso3),
    size = 2.1, colour = "grey15"
  ) +
  ggplot2::scale_fill_gradient(
    low = "#DCEAF3", high = "#12557A", name = "Surveys",
    breaks = scales::pretty_breaks(4)
  ) +
  ggplot2::labs(
    title = "Countries contributing to the panel",
    subtitle = paste0(nrow(surveys), " surveys, ", nrow(outlines),
                      " countries, ", nrow(analysis), " survey regions")
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank()
  )

## ---- panel B: when each survey was fielded ---------------------------------
# Countries ordered north to south so the timeline echoes the map.
order_by_latitude <- outlines$iso3[order(outlines$lat)]
surveys$country <- factor(surveys$iso3, levels = order_by_latitude)

timeline_panel <- ggplot2::ggplot(
  surveys, ggplot2::aes(x = year, y = country)
) +
  ggplot2::geom_line(ggplot2::aes(group = country), colour = "grey85",
                     linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(size = regions, shape = type),
                      colour = "#12557A", alpha = 0.85) +
  ggplot2::scale_size_continuous(range = c(0.8, 3.2), name = "Regions") +
  ggplot2::scale_shape_manual(values = c(DHS = 16, MIS = 17), name = "Type") +
  ggplot2::scale_x_continuous(breaks = seq(2000, 2025, 5)) +
  ggplot2::labs(title = "When each survey was fielded", x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

combined <- map_panel + timeline_panel +
  patchwork::plot_layout(widths = c(1, 1.15))

output <- file.path(RESULTS_DIR, "figure10_survey_map.png")
ggplot2::ggsave(output, combined, width = 11, height = 7.2, dpi = 200)
message("Wrote ", basename(output), " and survey_map_country_summary.csv")
message(nrow(surveys), " surveys across ", nrow(outlines), " countries, ",
        min(surveys$year), "-", max(surveys$year))

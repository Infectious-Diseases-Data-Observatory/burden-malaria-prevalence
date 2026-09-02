# =============================================================================
# 11_study_flow.R — analysis flow diagram for the DHS/MIS malaria-prevalence
# vs child-mortality study.
#
# Two things are shown together:
#   (1) the inclusion funnel, from the DHS survey registry through the assembled
#       survey-region panel to the shared-outcome analysis sample, with every
#       exclusion counted; and
#   (2) the external, non-DHS data that flow into the panel and covariate
#       block: MAP PfPR2-10 rasters, UNICEF WUENIC vaccine coverage, derived
#       child HIV prevalence (UNAIDS numbers / World Bank child population) and
#       World Bank economic covariates.
#
# Every number on the diagram is read from the pipeline's own outputs (the
# survey coverage audit, the region merge table, the inclusion counts, the
# analysis dataset, the covariate catalogue), so the figure cannot drift from
# the results the way a hardcoded version did.
# =============================================================================
source("R_dhs/00_config.R")
required_packages("ggplot2")

## ---- counts --------------------------------------------------------------------
audit <- read.csv(file.path(RESULTS_DIR, "survey_coverage_audit.csv"),
                  stringsAsFactors = FALSE)
merge_quality <- read.csv(file.path(RESULTS_DIR, "region_merge_quality.csv"),
                          stringsAsFactors = FALSE)
inclusion <- read.csv(file.path(RESULTS_DIR, "analysis_inclusion_counts.csv"),
                      stringsAsFactors = FALSE)
analysis <- read_analysis_data()
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
map_regions <- read.csv(MAP_REGION_CSV, stringsAsFactors = FALSE)

count <- function(stage) {
  v <- inclusion$region_years[inclusion$stage == stage]
  if (!length(v)) stop("Missing inclusion stage: ", stage)
  v
}
n_registry <- nrow(audit)
n_registry_countries <- length(unique(audit$iso3))
n_types <- table(audit$SurveyType)
excluded <- audit[!as.logical(audit$in_analysis), , drop = FALSE]
excluded_causes <- table(excluded$omission_cause)
cause_labels <- c(no_map_series = "no MAP PfPR series (no transmission)",
                  no_boundary = "no admin-1 boundary file",
                  no_recode = "recode not available",
                  no_key_overlap = "region keys do not join")
excluded_text <- paste(sprintf("%d %s", excluded_causes,
                               ifelse(names(excluded_causes) %in% names(cause_labels),
                                      cause_labels[names(excluded_causes)],
                                      names(excluded_causes))),
                       collapse = "\n")
excluded_countries <- paste(sort(unique(excluded$iso3)), collapse = ", ")

n_panel <- count("assembled survey-region panel")
n_panel_surveys <- length(unique(analysis$svkey))
n_panel_countries <- length(unique(analysis$iso3))
n_map_regions <- nrow(map_regions)
n_unmerged <- sum(merge_quality$unmatched_boundary)
n_country_filter <- count("country mean PfPR >1%")
n_region_filter <- count("region PfPR >=1%")
n_main <- count("main shared-outcome sample")
n_complete <- count("complete-case eligible-covariate sensitivity")
main <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
n_main_surveys <- length(unique(main$svkey))
n_main_countries <- length(unique(main$iso3))
n_covariates <- sum(as.logical(catalog$included_in_main))
n_prevalence_excluded <- n_panel - n_region_filter
n_zero_death_excluded <- n_region_filter - n_main
fmt <- function(x) format(x, big.mark = ",")

## ---- drawing helpers -------------------------------------------------------------
# palette: main funnel (blue), external sources (orange), exclusions (red),
# downstream burden (green), sensitivity and model comparison (grey)
MB <- "#08519c"; MF <- "#deebf7"
OB <- "#d95f0e"; OF <- "#fdd0a2"
EB <- "#a50f15"; EF <- "#fee0d2"
GB <- "#238b45"; GF <- "#e5f5e0"
SB <- "#525252"; SF <- "#f0f0f0"

box <- function(xc, yc, label, fill, border, hw, hh, size = 4.1) {
  list(ggplot2::geom_rect(ggplot2::aes(xmin = xc - hw, xmax = xc + hw,
                                       ymin = yc - hh, ymax = yc + hh),
                          fill = fill, colour = border, linewidth = 0.6),
       ggplot2::annotate("text", x = xc, y = yc, label = label, size = size,
                         lineheight = 0.95))
}
arr <- function(x0, y0, x1, y1, colour = "grey25") {
  ggplot2::annotate("segment", x = x0, y = y0, xend = x1, yend = y1,
                    arrow = ggplot2::arrow(length = ggplot2::unit(2.6, "mm"),
                                           type = "closed"),
                    linewidth = 0.55, colour = colour)
}

## ---- the diagram -----------------------------------------------------------------
# Left column (funnel) and right column (exclusions, comparison, burden). Text
# sizes are chosen so the longest line in each box fits at the saved width.
LX <- 6.3; LW <- 5.7      # funnel boxes span 0.6 to 12.0
RX <- 16.3; RW <- 3.3     # side boxes span 13.0 to 19.6
p <- ggplot2::ggplot() +
  # ---- external data sources (top band) ----
  box(3.3, 13.4,
      sprintf("DHS/MIS Births Recodes, 2000-2024 (rdhs)\n%d surveys (%d DHS, %d MIS)\n%d countries",
              n_registry, n_types[["DHS"]], n_types[["MIS"]], n_registry_countries),
      MF, MB, 2.7, 0.95, size = 3.7) +
  box(9.7, 13.4,
      "MAP annual PfPR2-10 rasters (malariaAtlas)\nx GPW population density, weighted\nto each admin-1 survey region",
      OF, OB, 2.9, 0.95, size = 3.7) +
  box(16.3, 13.4,
      "National series (country-year)\nUNICEF WUENIC vaccines; UNAIDS 0-14 HIV\ndivided by World Bank 0-14 population;\nWorld Bank GDP, health spending, stability",
      OF, OB, 3.3, 0.95, size = 3.5) +

  # ---- inclusion funnel ----
  box(LX, 10.0,
      sprintf(paste0(
        "Assembled survey-region panel\n",
        "admin-1 mortality over the %d months before interview (DHS synthetic-cohort life table)\n",
        "regional PfPR2-10 in the survey year; %d covariates\n",
        "%s region-years; %d surveys; %d countries"),
        CHMORT_PERIOD, n_covariates, fmt(n_panel), n_panel_surveys, n_panel_countries),
      MF, MB, LW, 1.15, size = 3.7) +
  box(RX, 10.0,
      sprintf("Excluded: %d surveys (%s)\n%s\n\n%d of %s MAP regions not merged:\nMali 2012 regions not surveyed,\nUganda 2016 Central labelling ambiguous",
              nrow(excluded), excluded_countries, excluded_text,
              n_unmerged, fmt(n_map_regions)),
      EF, EB, RW, 1.15, size = 3.4) +

  box(LX, 6.6,
      sprintf(paste0(
        "Primary sample\n",
        "country-mean PfPR2-10 above 1%% and regional PfPR2-10 at least 1%%: %s region-years\n",
        "at least one post-neonatal and one neonatal death in the window: %s region-years\n",
        "%d surveys; %d countries"),
        fmt(n_region_filter), fmt(n_main), n_main_surveys, n_main_countries),
      MF, MB, LW, 1.05, size = 3.7) +
  box(RX, 6.6,
      sprintf("Excluded: %d region-years with PfPR2-10 below 1%%\n(retained in a sensitivity analysis)\n%d region-years with no post-neonatal or\nno neonatal death within %d months",
              n_prevalence_excluded, n_zero_death_excluded, CHMORT_PERIOD),
      EF, EB, RW, 0.95, size = 3.4) +

  box(LX, 3.3,
      sprintf(paste0(
        "MAIN ANALYSIS: ridge-penalised negative-binomial GAM (mgcv)\n",
        "%d standardised covariates in one ridge block; country random intercept; offset\n",
        "five prevalence-by-time structures compared by AIC\n",
        "post-neonatal mortality primary, neonatal the negative control; attributable fraction vs 1%%"),
        n_covariates),
      MF, MB, LW, 1.05, size = 3.7) +
  box(RX, 3.3,
      "Bayesian refits (brms / Stan)\nt2(PfPR, year) surface, smoothness estimated\nadditive vs linear-in-time vs full surface:\nPSIS-LOO and survey-grouped 10-fold CV\nera x region subgroups; neonatal control;\nneonatal rate as a covariate",
      SF, SB, RW, 1.05, size = 3.4) +

  box(LX, 0.5,
      sprintf(paste0(
        "Sensitivity analyses\n",
        "complete cases (%s); include PfPR2-10 below 1%%; restrict to 5-40%%; covariate blocks;\n",
        "log-Gaussian likelihood; country PfPR slope; mortality window 12-60 months\n",
        "x prevalence lag 0-2 years; timing of fieldwork"),
        fmt(n_complete)),
      SF, SB, LW, 1.0, size = 3.7) +
  box(RX, 0.5,
      "National burden extrapolation\npopulation-average attributable fraction x\nnational MAP PfPR x IGME all-cause deaths\n2000-2024 trend under each time structure;\n2024 against IHME (GBD) and WHO (WMR 2025)",
      GF, GB, RW, 1.0, size = 3.4) +

  # ---- arrows: sources into the panel ----
  arr(3.3, 13.4 - 0.95, 4.4, 10.0 + 1.15) +
  arr(9.7, 13.4 - 0.95, 7.2, 10.0 + 1.15) +
  arr(16.3, 13.4 - 0.95, 10.4, 10.0 + 1.15) +
  # ---- funnel arrows ----
  arr(LX, 10.0 - 1.15, LX, 6.6 + 1.05) +
  arr(LX, 6.6 - 1.05, LX, 3.3 + 1.05) +
  arr(LX, 3.3 - 1.05, LX, 0.5 + 1.0) +
  # ---- exclusions, comparison and burden branches ----
  arr(LX + LW, 10.0, RX - RW, 10.0) +
  arr(LX + LW, 6.6, RX - RW, 6.6) +
  arr(LX + LW, 3.3, RX - RW, 3.3) +
  arr(LX + LW, 3.3 - 0.7, RX - RW, 0.5 + 0.5) +

  ggplot2::coord_cartesian(xlim = c(0, 20), ylim = c(-0.7, 14.6)) +
  ggplot2::theme_void()

output <- file.path(RESULTS_DIR, "study_flow_diagram.png")
ggplot2::ggsave(output, p, width = 13, height = 8.6, dpi = 320, bg = "white")
cat("saved:", output, "\n")
cat(sprintf(paste(
  "registry %d surveys -> panel %s region-years (%d surveys, %d countries) ->",
  "prevalence filters %s -> shared-outcome sample %s (%d surveys, %d countries);",
  "%d covariates; %d MAP regions unmerged\n"),
  n_registry, fmt(n_panel), n_panel_surveys, n_panel_countries,
  fmt(n_region_filter), fmt(n_main), n_main_surveys, n_main_countries,
  n_covariates, n_unmerged))

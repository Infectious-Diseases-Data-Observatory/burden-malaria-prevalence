# =============================================================================
# 11_study_flow.R — analysis flow diagram for the DHS/MIS malaria-prevalence
# vs child-mortality study. (Ported from the legacy R/37_study_flow.R.)
#
# Two things are shown together:
#   (1) the survey-inclusion funnel (enumerate -> assembled panel 936 ->
#       primary sample 921 -> main analysis), and
#   (2) the external, non-DHS data that flow INTO the panel and covariate
#       block: MAP PfPR2-10 rasters, UNICEF WUENIC vaccine coverage, and
#       derived child HIV prevalence (UNAIDS numbers / World Bank child
#       population), plus World Bank economic covariates.
#
# Self-contained: needs only ggplot2. Writes the figure to the DHS-rebuild
# results directory so it sits alongside the other manuscript figures.
# =============================================================================
suppressMessages(library(ggplot2))

OUT_DIR <- file.path("results", "dhs_rebuild")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# palette: main funnel (blue), external sources (orange), exclusions (red),
# downstream burden (green), sensitivity (grey)
MB <- "#08519c"; MF <- "#deebf7"      # main
OB <- "#d95f0e"; OF <- "#fdd0a2"      # outside/external data
EB <- "#a50f15"; EF <- "#fee0d2"      # excluded
GB <- "#238b45"; GF <- "#e5f5e0"      # burden
SB <- "#525252"; SF <- "#f0f0f0"      # sensitivity

box <- function(xc, yc, label, fill, border, hw, hh, size = 3.5)
  list(geom_rect(aes(xmin = xc - hw, xmax = xc + hw, ymin = yc - hh, ymax = yc + hh),
                 fill = fill, colour = border, linewidth = 0.6),
       annotate("text", x = xc, y = yc, label = label, size = size, lineheight = 0.95))
arr <- function(x0, y0, x1, y1, colour = "grey25")
  annotate("segment", x = x0, y = y0, xend = x1, yend = y1,
           arrow = arrow(length = unit(2.6, "mm"), type = "closed"),
           linewidth = 0.55, colour = colour)

p <- ggplot() +
  # ---- external data sources (top band) ----
  box(3.3, 13.2,
      "DHS/MIS Births Recodes\n2000–2024 (rdhs)\n124 surveys · 37 countries",
      MF, MB, 2.7, 1.0) +
  box(10.0, 13.2,
      "MAP annual PfPR₂₋₁₀ rasters\n(malariaAtlas)\n× GPW population density",
      OF, OB, 2.6, 1.0) +
  box(16.7, 13.2,
      "National series (country–year)\nUNICEF WUENIC vaccines ·\nUNAIDS 0–14 HIV ÷ World Bank 0–14 pop ·\nWorld Bank GDP & health expenditure",
      OF, OB, 3.1, 1.0, size = 3.2) +

  # ---- inclusion funnel ----
  box(6.7, 9.6,
      "Assembled survey-region panel\nadmin-1 child mortality (DHS synthetic-cohort life table)\n+ population-weighted regional PfPR₂₋₁₀ + regional & national covariates\n936 region-years · 106 surveys · 35 countries",
      MF, MB, 4.5, 1.15) +
  box(15.3, 9.6,
      "Excluded: 18 surveys\nno usable admin-1 boundary\nor region-key merge",
      EF, EB, 2.4, 0.85, size = 3.2) +

  box(6.7, 6.5,
      "Primary sample: PfPR₂₋₁₀ ≥ 1% (country mean > 1%)\n921 region-years · 105 surveys · 34 countries",
      MF, MB, 4.5, 0.85) +
  box(15.3, 6.5,
      "Excluded: 15 region-years\nPfPR₂₋₁₀ < 1%\n(retained in sensitivity)",
      EF, EB, 2.4, 0.8, size = 3.2) +

  box(6.7, 3.3,
      "MAIN ANALYSIS — ridge-penalised negative-binomial GAM\n17 covariates (regional DHS + national WUENIC, child HIV, World Bank)\ncountry random intercept & random PfPR slope · log-exposure offset",
      MF, MB, 4.5, 1.05) +
  box(15.3, 3.3,
      "National burden extrapolation\nselected model × national MAP PfPR (2024)\n× IGME all-cause mortality × live births\n→ compare with IHME (GBD 2025) & WHO (WMR 2025)",
      GF, GB, 3.3, 1.05, size = 3.1) +

  box(6.7, 0.8,
      "Sensitivity analyses (malaria–mortality association essentially unchanged):\ncomplete-case covariates (859) · include PfPR₂₋₁₀ < 1% · restrict 5–40% ·\nremove HIV/vaccine blocks · log-Gaussian likelihood · alternative response forms",
      SF, SB, 4.5, 0.85, size = 3.1) +

  # ---- arrows: sources into the panel ----
  arr(3.3, 13.2 - 1.0, 4.6, 9.6 + 1.15) +
  arr(10.0, 13.2 - 1.0, 7.4, 9.6 + 1.15) +
  arr(16.7, 13.2 - 1.0, 10.2, 9.6 + 1.15) +
  # ---- funnel arrows ----
  arr(6.7, 9.6 - 1.15, 6.7, 6.5 + 0.85) +
  arr(6.7, 6.5 - 0.85, 6.7, 3.3 + 1.05) +
  arr(6.7, 3.3 - 1.05, 6.7, 0.8 + 0.85) +
  # ---- exclusions + burden branch ----
  arr(6.7 + 4.5, 9.6, 15.3 - 2.4, 9.6) +
  arr(6.7 + 4.5, 6.5, 15.3 - 2.4, 6.5) +
  arr(6.7 + 4.5, 3.3, 15.3 - 3.3, 3.3) +

  coord_cartesian(xlim = c(0, 20), ylim = c(-0.3, 14.5)) +
  theme_void()

ggsave(file.path(OUT_DIR, "study_flow_diagram.png"), p,
       width = 14, height = 9.2, dpi = 300, bg = "white")
cat("saved:", file.path(OUT_DIR, "study_flow_diagram.png"), "\n")

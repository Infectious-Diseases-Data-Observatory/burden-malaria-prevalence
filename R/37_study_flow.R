# =============================================================================
# 37_study_flow.R — PRISMA-style study-inclusion flow diagram for the DHS/MIS
# analysis. Funnel: enumerate -> panel (936) -> [exclude PfPR<1%] -> primary
# sample (921) -> MAIN ANALYSIS. The SENSITIVITY branch tees off the PANEL (936),
# BEFORE the <1% exclusion (it re-includes sub-1% region-years); main + sensitivity
# sit on the same bottom row. Counts from survey_inclusion.csv + the fit samples.
# =============================================================================
source("R/00_utils.R"); suppressMessages(library(ggplot2))

box <- function(xc, yc, label, fill, border, hw, hh, size)
  list(geom_rect(aes(xmin = xc-hw, xmax = xc+hw, ymin = yc-hh, ymax = yc+hh),
                 fill = fill, colour = border, linewidth = 0.6),
       annotate("text", x = xc, y = yc, label = label, size = size, lineheight = 0.98))
arr <- function(x0,y0,x1,y1) annotate("segment", x=x0, y=y0, xend=x1, yend=y1,
       arrow = arrow(length = unit(2.6,"mm"), type="closed"), linewidth = 0.55, colour = "grey30")
seg <- function(x0,y0,x1,y1) annotate("segment", x=x0, y=y0, xend=x1, yend=y1, linewidth = 0.55, colour = "grey30")

MB<-"#08519c"; MF<-"#deebf7"; EB<-"#a50f15"; EF<-"#fee0d2"; SB<-"#238b45"; SF<-"#e5f5e0"; xc<-2.8
p <- ggplot() +
  box(xc, 12.2, "Enumerated: sub-Saharan Africa\nDHS/MIS/AIS Births Recodes 2000–2024\n128 surveys, 37 countries", MF, MB, 2.5, 0.9, 4.0) +
  box(xc,  9.5, "Panel: admin-1 child mortality +\nmalaria prevalence (MAP PfPR₂₋₁₀)\n106 surveys, 35 countries · 936 region-years", MF, MB, 2.5, 0.9, 4.0) +
  box(xc,  6.8, "Primary sample: PfPR₂₋₁₀ ≥ 1%\n921 region-years\n105 surveys, 34 countries", MF, MB, 2.5, 0.82, 4.0) +
  box(xc,  3.0, "MAIN ANALYSIS\n921 region-years · 105 surveys, 34 countries\nridge-penalised covariate adjustment\n(missing covariates imputed)", MF, MB, 2.5, 1.0, 3.8) +
  box(10.5, 3.0, "SENSITIVITY ANALYSES\n(malaria effect essentially unchanged)\n• complete-case covariates: 847 region-years\n• include PfPR₂₋₁₀ < 1%\n• parsimonious region-health model\n• alternative prevalence-response forms", SF, SB, 3.2, 1.05, 3.8) +
  box(8.5, 10.9, "Excluded: 22 surveys —\nno usable admin-1 boundary\nor region-key merge (older\nzone schemes, continuous-DHS\nfiles, malaria-free countries)", EF, EB, 2.15, 1.15, 4.0) +
  box(8.5,  7.9, "Excluded: 15 region-years\nwith PfPR₂₋₁₀ < 1%\n(retained in sensitivity)", EF, EB, 2.15, 0.82, 4.0) +
  arr(xc, 12.2-0.9, xc, 9.5+0.9) + arr(xc, 9.5-0.9, xc, 6.8+0.82) + arr(xc, 6.8-0.82, xc, 3.0+1.0) +
  arr(xc, 10.9, 8.5-2.15, 10.9) + arr(xc, 7.9, 8.5-2.15, 7.9) +
  seg(xc+2.5, 9.5, 11.8, 9.5) + arr(11.8, 9.5, 11.8, 3.0+1.05) +      # sensitivity branch off the PANEL (936)
  coord_cartesian(xlim = c(0, 14.5), ylim = c(1.5, 13.2)) + theme_void()
ggsave(file.path(RESULTS, "study_flow_diagram.png"), p, width = 13.5, height = 8.5, dpi = 300, bg = "white")
cat("saved: results/study_flow_diagram.png\n")

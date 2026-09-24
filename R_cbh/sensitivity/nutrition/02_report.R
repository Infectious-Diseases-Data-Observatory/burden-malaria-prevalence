#!/usr/bin/env Rscript
# Supplementary figure and tables for the no-nutrition sensitivity: overlay the
# primary PfPR curves with those from the refit that drops regional wasting and
# stunting. Reads saved aggregates only; no refitting, no TeX writes.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE); stopifnot(all(args %in% "--dhs-only"))
dhs_only <- "--dhs-only" %in% args
st <- list(out = file.path("results/cbh", if (dhs_only) "nutrition_adjustment_map_gamma2_v1" else "nutrition_adjustment_dhsmics_map_gamma2_v3"))
primary <- cbh_primary_settings(if (dhs_only) "regional" else "regional_mics")
labels <- c(full = "All 17 covariates (primary analysis)",
            drop = "Wasting and stunting removed (sensitivity)")
colours <- setNames(c("#B36B39", "#215E91"), labels)
ages <- cbh_config()$age_bands$age_band
age_label <- function(x) factor(x, levels = ages,
  labels = ifelse(ages == "<1", "<1 month", paste(ages, "months")))

inputs <- c(file.path(primary$out, "pfpr_curves.csv"), file.path(st$out, "pfpr_curves.csv"),
            file.path(primary$out, "fit_diagnostics.csv"), file.path(st$out, "fit_diagnostics.csv"),
            file.path(st$out, "dropped_covariates.csv"))
stopifnot(file.exists(inputs))
pc <- cbh_read_csv(inputs[1]); pc <- pc[pc$series == "map_full", ]
nc <- cbh_read_csv(inputs[2])
pd <- cbh_read_csv(inputs[3]); nd <- cbh_read_csv(inputs[4])
stopifnot(nrow(pd) == 7L, nrow(nd) == 7L, all(nd$converged), all(nd$gamma == 2),
          all(nd$input_verified), all(nd$min_smoothing_hessian_eigenvalue > 0),
          identical(sort(pd$rows), sort(nd$rows)), identical(sort(pd$deaths), sort(nd$deaths)))
pc$series_label <- labels[["full"]]; nc$series_label <- labels[["drop"]]
curves <- rbind(pc[, intersect(names(pc), names(nc))], nc[, intersect(names(pc), names(nc))])
curves$series_label <- factor(c(rep(labels[["full"]], nrow(pc)), rep(labels[["drop"]], nrow(nc))),
                              levels = labels)
curves$age_label <- age_label(curves$age_band)
cbh_atomic_csv(curves, file.path(st$out, "comparison_pfpr_curves.csv"))

paper <- theme_minimal(base_size = 18) + theme(panel.grid.minor = element_blank(),
  axis.text = element_text(size = 15), strip.text = element_text(size = 18, face = "bold"),
  legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16),
  panel.spacing = grid::unit(20, "pt"), plot.margin = margin(12, 22, 12, 12))
central <- curves[curves$within_central_support, ]
p <- ggplot(central, aes(pfpr_pct, log_hazard_ratio, colour = series_label, fill = series_label)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = .4) +
  geom_vline(xintercept = 20, colour = "grey80", linewidth = .4) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = .12, colour = NA) +
  geom_line(linewidth = 1) + facet_wrap(~age_label, ncol = 4) +
  scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
  scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) +
  labs(x = "PfPR[2-10] (%)", y = "Log hazard ratio\nrelative to PfPR = 20%") + paper
figure <- file.path(st$out, "sfig_pfpr_splines_without_nutrition.png")
ggsave(figure, p, width = 14, height = 8, dpi = 220, device = ragg::agg_png, bg = "white")

# Contrasts on the shared grid: 40% to 20% and 20% to 0%.
contrast <- function(d, iteration) do.call(rbind, lapply(ages, function(a) {
  z <- d[d$age_band == a, ]
  g <- function(v) z[which.min(abs(z$pfpr_pct - v)), ]
  hi <- g(40); lo <- g(20); ze <- g(0)
  data.frame(age_band = a, iteration,
    hr_40_to_20 = exp(lo$log_hazard_ratio - hi$log_hazard_ratio),
    hr_20_to_0  = exp(ze$log_hazard_ratio - lo$log_hazard_ratio))
}))
tab <- rbind(contrast(pc, labels[["full"]]), contrast(nc, labels[["drop"]]))
cbh_atomic_csv(tab, file.path(st$out, "comparison_contrasts.csv"))
w <- reshape(tab, idvar = "age_band", timevar = "iteration", direction = "wide")
w <- w[match(ages, w$age_band), ]
delta <- max(abs(w[[2]] - w[[4]]))
fmt <- function(x) sprintf("%.2f", x)
md <- c("# No-nutrition sensitivity: hazard ratios for PfPR 40% to 20%", "",
  "| Age (months) | All 17 covariates | Wasting and stunting removed |", "|---|---:|---:|")
for (i in seq_len(nrow(w)))
  md <- c(md, sprintf("| %s | %s | %s |", w$age_band[i], fmt(w[[2]][i]), fmt(w[[4]][i])))
md <- c(md, "", sprintf("Largest absolute change in the 40%% to 20%% hazard ratio: %.3f.", delta))
writeLines(md, file.path(st$out, "CONTRASTS.md"))

caption <- c("# Supplementary figure caption: adjustment without wasting and stunting", "",
  paste0("Sensitivity of the fitted PfPR-mortality relationship to removing the two nutritional covariates. ",
  "Log mortality hazard ratios relative to PfPR[2-10] = 20% by completed-month age band, with pointwise 95% conditional intervals shaded, from the PfPR-ACM model with all 17 covariates (the primary analysis) and from a refit of the same seven models on the identical sample with regional wasting and stunting prevalence removed, leaving 15 covariates. ",
  "Wasting and stunting are measured in surviving children at the time of the survey and may lie on the causal pathway from malaria to death, so adjusting for them risks removing part of the association of interest. ",
  "Sample, exposure, reference knots, basis dimensions, gamma = 2, offset and random-effect structure are unchanged; only the adjustment set differs. ",
  "Each curve is drawn over the central 95% of its own exposure distribution. ",
  sprintf("The largest absolute change in the hazard ratio for a reduction in PfPR[2-10] from 40%% to 20%% is %.3f, at %s months. ",
          delta, w$age_band[which.max(abs(w[[2]] - w[[4]]))]),
  "The two models are fitted to the same records, so the comparison is descriptive."))
writeLines(caption, file.path(st$out, "CAPTION.md"))
prov <- c(inputs, figure)
cbh_atomic_csv(data.frame(file = prov, md5 = vapply(prov, cbh_file_hash, "")),
               file.path(st$out, "report_provenance.csv"))
cat("No-nutrition sensitivity figure and tables saved:", figure, "\n")
cat(sprintf("Largest 40%%->20%% hazard-ratio change: %.3f\n", delta))

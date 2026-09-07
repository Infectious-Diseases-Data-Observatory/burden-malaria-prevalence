#!/usr/bin/env Rscript
# Plot saved contrasts and fitted PfPR curves; no model refitting.
# The private fitted object is read locally; only aggregate predictions are exported.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
out <- file.path("results/cbh", cbh_trial_spec()$id)
d <- cbh_read_csv(file.path(out, "pfpr_40_to_20_contrasts.csv"))
age <- cbh_config()$age_bands$age_band
stopifnot(nrow(d) == length(age), setequal(d$age_band, age))
d <- d[match(age, d$age_band), ]
stopifnot(all(is.finite(as.matrix(d[c("hazard_ratio", "lower_95", "upper_95")]))), all(d$lower_95 > 0))
path <- file.path(out, "pfpr_40_to_20_by_age.png")
png(path, width = 1500, height = 1100, res = 180, bg = "white")
par(mar = c(5.2, 6.4, 4.7, 1.5), family = "sans", las = 1)
y <- rev(seq_len(nrow(d)))
xlim <- range(c(d$lower_95, d$upper_95, 1))
xlim <- exp(log(xlim) + c(-.12, .12) * max(.5, diff(log(xlim))))
plot(NA, xlim = xlim, ylim = c(.5, nrow(d) + .5), log = "x", yaxt = "n", bty = "n",
     xlab = "Mortality hazard ratio: PfPR 40% to 20%", ylab = "", cex.lab = .95)
abline(h = y, col = "#EEEEEE", lwd = .8)
abline(v = 1, col = "#777777", lty = 2)
segments(d$lower_95, y, d$upper_95, y, col = "#176B87", lwd = 2)
segments(d$lower_95, y - .08, d$lower_95, y + .08, col = "#176B87", lwd = 1.5)
segments(d$upper_95, y - .08, d$upper_95, y + .08, col = "#176B87", lwd = 1.5)
points(d$hazard_ratio, y, pch = 19, col = "#176B87", cex = 1.15)
axis(2, at = y, labels = ifelse(d$age_band == "<1", "<1 month", paste0(d$age_band, " months")),
     tick = FALSE, cex.axis = .9)
title("Age-specific PfPR associations", adj = 0, line = 2.6, cex.main = 1.2)
mtext("Exploratory complete-case fit; adjusted, unweighted cloglog model", side = 3,
      line = 1.15, adj = 0, cex = .78, col = "#444444")
mtext("Bars: pointwise 95% model-based intervals; survey-design uncertainty not included.",
      side = 1, line = 3.8, cex = .67, col = "#444444")
dev.off()
message("Saved ", path)

# Differences in linear predictors isolate the age-specific PfPR smooth.
# Holding all other predictors fixed also cancels intercepts, offsets and REs.
saved <- readRDS(file.path("data/derived_cbh/models", cbh_trial_spec()$id, "fit.rds"))
fit <- saved$fit
stopifnot(isTRUE(fit$converged), all(is.finite(fit$Vp)))
reference <- 20
curves <- lapply(age, function(band) {
  rows <- which(fit$model$age_band == band)
  exposure <- fit$model$pfpr_pct[rows]
  support <- quantile(exposure, c(0, .025, .975, 1), names = FALSE)
  grid <- sort(unique(c(seq(support[1], support[4], length.out = 401),
                        support, reference, 40)))
  stopifnot(reference >= support[1], reference <= support[4])
  template <- fit$model[rows[1], , drop = FALSE]
  template$band_years <- exp(template[["offset(log(band_years))"]])
  new <- template[rep(1L, length(grid)), , drop = FALSE]
  new$pfpr_pct <- grid
  ref <- new
  ref$pfpr_pct <- reference
  # Use exact prediction matrices: the discretized predictor can fail on a
  # single-band template with constant nuisance covariates.
  L <- predict(fit, newdata = new, type = "lpmatrix", discrete = FALSE) -
       predict(fit, newdata = ref, type = "lpmatrix", discrete = FALSE)
  estimate <- drop(L %*% coef(fit))
  variance <- rowSums((L %*% fit$Vp) * L)
  stopifnot(all(is.finite(variance)), min(variance) > -1e-8)
  se <- sqrt(pmax(variance, 0))
  direct <- predict(fit, newdata = new, type = "link", discrete = FALSE) -
            predict(fit, newdata = ref, type = "link", discrete = FALSE)
  stopifnot(max(abs(estimate - direct)) < 1e-7,
            abs(estimate[grid == reference]) < 1e-10,
            se[grid == reference] < 1e-10)
  # This is the reverse of the previously reported 40 -> 20 contrast.
  previous <- d[d$age_band == band, ]
  check <- exp(c(estimate[grid == 40],
                 estimate[grid == 40] - 1.96 * se[grid == 40],
                 estimate[grid == 40] + 1.96 * se[grid == 40]))
  stopifnot(max(abs(check - 1 / unlist(previous[c("hazard_ratio", "upper_95", "lower_95")],
                                              use.names = FALSE))) < 1e-7)
  data.frame(age_band = band, pfpr_pct = grid, reference_pfpr_pct = reference,
             log_hazard_ratio = estimate, standard_error = se,
             lower_95 = estimate - 1.96 * se, upper_95 = estimate + 1.96 * se,
             central_95_support = grid >= support[2] & grid <= support[3],
             pfpr_p025 = support[2], pfpr_p975 = support[3],
             fit_signature = saved$signature)
})
curves <- do.call(rbind, curves)
cbh_atomic_csv(curves, file.path(out, "pfpr_spline_curves.csv"))
curves$age_label <- factor(curves$age_band, levels = age,
                          labels = ifelse(age == "<1", "<1 month", paste0(age, " months")))

library(ggplot2)
plot_curves <- function(full_range = FALSE) {
  pdat <- if (full_range) curves else curves[curves$central_95_support, ]
  p <- ggplot(pdat, aes(pfpr_pct, log_hazard_ratio)) +
    geom_hline(yintercept = 0, colour = "#87939C", linetype = "dashed", linewidth = .4) +
    geom_vline(xintercept = reference, colour = "#CCD2D6", linewidth = .4) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#176B87", alpha = .17) +
    geom_line(colour = "#176B87", linewidth = .8) +
    facet_wrap(~age_label, ncol = 4) +
    scale_x_continuous(breaks = seq(0, 100, 20), expand = expansion(mult = .025)) +
    labs(title = "Estimated PfPR splines by age band",
         subtitle = paste0("Adjusted complete-case model | ",
                           if (full_range) "Full observed exposure range" else "Central 95% of exposure values in each age band"),
         x = expression(paste("PfPR"["2-10"], " at band entry (%)")),
         y = "Log mortality hazard ratio relative to 20% PfPR",
         caption = paste0("Shading: pointwise 95% model-based intervals. All curves are centered at 20% PfPR, where the contrast and its uncertainty are zero.\n",
                          "Common axes across panels. Unweighted exploratory fit; spline-basis and smoothing-optimization checks remain open.",
                          if (full_range) "\nDotted vertical lines mark the 2.5th and 97.5th exposure percentiles." else "")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#EEF0F2"),
          strip.text = element_text(face = "bold", hjust = 0, size = 12),
          plot.title = element_text(face = "bold", size = 19),
          plot.subtitle = element_text(colour = "#52616B", margin = margin(b = 15)),
          plot.caption = element_text(hjust = 0, colour = "#52616B", size = 9, lineheight = 1.3),
          panel.spacing = grid::unit(1.2, "lines"),
          axis.title = element_text(size = 12), plot.margin = margin(16, 20, 12, 16))
  if (full_range) {
    bounds <- unique(curves[c("age_label", "pfpr_p025", "pfpr_p975")])
    p <- p + geom_vline(data = bounds, aes(xintercept = pfpr_p025),
                        linetype = "dotted", colour = "#87939C", linewidth = .4) +
             geom_vline(data = bounds, aes(xintercept = pfpr_p975),
                        linetype = "dotted", colour = "#87939C", linewidth = .4)
  }
  name <- if (full_range) "pfpr_splines_full_range.png" else "pfpr_splines_by_age.png"
  ggsave(file.path(out, name), p, device = ragg::agg_png,
         width = 14, height = 8, units = "in", dpi = 160, bg = "white")
  message("Saved ", file.path(out, name))
}
plot_curves()
plot_curves(full_range = TRUE)
message("Verified seven curves against direct link predictions, the 20% reference and saved 40-to-20 contrasts.")

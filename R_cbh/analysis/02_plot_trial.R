#!/usr/bin/env Rscript
# Plot saved aggregate contrasts; no model fitting or child-level data reads.
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

# =============================================================================
# 45_person_time_data_and_curves.R — the six-band dose-response with the data
# it is fitted to shown above it.
#
# Top row: one point per survey region, MAP PfPR2-10 (person-time-weighted over
# the region's windows) against the band's death rate per 1,000 child-years over
# the five 12-month windows, with each panel on its own scale because the
# neonatal rate is an order of magnitude above the rest; point size is the
# child-years behind the rate, and each panel states the band's total deaths and
# share of all under-5 deaths. Bottom row: the smooth prevalence effect per
# band from script 42 (hazard ratio against 0% prevalence). Both rows read
# script 42's outputs; nothing is refitted here.
#
# Outputs
#   results/dhs_rebuild/figure27_person_time_data_and_curves.png     bands split at 3 months
#   results/dhs_rebuild/figure29_person_time_data_and_curves_4m.png  bands split at 4 months
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("ggplot2", "patchwork"))

model_data <- read.csv(file.path(RESULTS_DIR, "person_time_model_data.csv"),
                       stringsAsFactors = FALSE)
all_curves <- read.csv(file.path(RESULTS_DIR, "person_time_dose_response_curves.csv"),
                       stringsAsFactors = FALSE)

draw <- function(grouping, variable, levels, file) {
curves <- all_curves[all_curves$grouping == grouping, ]
curves$age_group <- factor(curves$age_group, levels = levels)

## ---- one point per survey region, per band ------------------------------------------
model_data$pfpr_pm <- model_data$pfpr * model_data$person_months_w
model_data$band <- model_data[[variable]]
regions <- aggregate(cbind(deaths_w, person_months_w, pfpr_pm) ~ svkey + regkey + band,
                     model_data, sum)
regions$pfpr <- regions$pfpr_pm / regions$person_months_w
regions$child_years <- regions$person_months_w / 12
regions$rate <- 1000 * regions$deaths_w / regions$child_years
regions$age_group <- factor(regions$band, levels = levels)
shares <- aggregate(deaths_w ~ age_group, regions, sum)
shares$share <- shares$deaths_w / sum(shares$deaths_w)
shares$label <- sprintf("%s deaths\n%.0f%% of all under-5 deaths",
                        format(round(shares$deaths_w), big.mark = ","), 100 * shares$share)
n_regions <- length(unique(paste(regions$svkey, regions$regkey)))
n_above_60 <- length(unique(paste(regions$svkey, regions$regkey)[regions$pfpr > 60]))
regions <- regions[regions$pfpr <= 60, ]

top <- ggplot2::ggplot(regions, ggplot2::aes(pfpr, rate)) +
  ggplot2::geom_point(ggplot2::aes(size = child_years), alpha = 0.25, colour = "#1D6F8B") +
  ggplot2::geom_text(data = shares, ggplot2::aes(x = 60, y = Inf, label = label),
                     hjust = 1, vjust = 1.3, size = 2.7, colour = "grey25", inherit.aes = FALSE) +
  ggplot2::facet_wrap(~age_group, nrow = 1, scales = "free_y") +
  ggplot2::scale_x_continuous(limits = c(0, 60)) +
  ggplot2::scale_size_area(max_size = 3, guide = "none") +
  ggplot2::labs(x = NULL, y = "Deaths per 1,000 child-years in the band,\nfive windows combined",
                title = "Observed death rate against prevalence, one point per survey region",
                subtitle = sprintf(paste0("%s survey regions (%d with prevalence above 60%% not shown); ",
                                          "point size is child-years at risk; each panel on its own scale.\n",
                                          "Prevalence is the person-time-weighted MAP PfPR2-10 over the region's windows"),
                                   format(n_regions, big.mark = ","), n_above_60)) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_blank())

bottom <- ggplot2::ggplot(curves, ggplot2::aes(pfpr, hr)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "#1D6F8B", alpha = 0.15) +
  ggplot2::geom_line(colour = "#1D6F8B", linewidth = 0.9) +
  ggplot2::geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.4) +
  ggplot2::facet_wrap(~age_group, nrow = 1) +
  ggplot2::scale_x_continuous(limits = c(0, 60)) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(x = "MAP PfPR2-10 (%)",
                y = sprintf("Mortality hazard ratio\nversus %d%% prevalence (log scale)", PERSON_TIME_AF_REFERENCE),
                title = "Fitted dose-response, one negative-binomial model per band",
                subtitle = "Smooth prevalence effect with segment, window, calendar-year, covariate, country and survey terms; 95% CIs") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(strip.text = ggplot2::element_blank())

combined <- top / bottom + patchwork::plot_layout(heights = c(1, 1.05))
ggplot2::ggsave(file, combined, width = 12, height = 7.2, dpi = 200, bg = "white")
message(grouping, ", deaths by band: ",
        paste(sprintf("%s %.0f%%", shares$age_group, 100 * shares$share), collapse = "; "))
message("Wrote ", basename(file))
}
draw("six bands", "age6", AGE6,
     file.path(RESULTS_DIR, "figure27_person_time_data_and_curves.png"))
draw("six bands (4-month split)", "age6b", AGE6B,
     file.path(RESULTS_DIR, "figure29_person_time_data_and_curves_4m.png"))

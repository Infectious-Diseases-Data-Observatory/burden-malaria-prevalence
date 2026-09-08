#!/usr/bin/env Rscript
# Descriptive country-year HIV pairs; no imputation or model fitting.
# Run from the project root. Uses only the existing UNICEF workbook.
# Add --incidence to plot the reported rates per 1,000 uninfected population.
library(ggplot2)
source("R_cbh/load_pipeline.R")
out <- "results/cbh/hiv_age_pairs"
source_file <- "data/HIV_Epidemiology_Children_Adolescents_2025.xlsx"
incidence <- "--incidence" %in% commandArgs(trailingOnly = TRUE)
metric <- if (incidence) "incidence" else "counts"
indicator <- if (incidence) "Estimated incidence rate (new HIV infection per 1,000 uninfected population)" else
  "Estimated number of people living with HIV"
x <- as.data.frame(suppressMessages(readxl::read_excel(
  source_file, sheet = "Data", skip = 1, guess_max = 100000)))
x <- x[x$Type == "Country" & x$Sex == "Both" &
         x$Indicator == indicator &
         x$Age %in% c("Age 0-14", "Age 15-19"), ]
x$year <- as.integer(x$Year)
x$value_text <- trimws(as.character(x$Value))
x$upper_limit <- grepl("^<", x$value_text)
x$plot_value <- suppressWarnings(as.numeric(gsub("[<,[:space:]]", "", x$value_text)))
stopifnot(!anyDuplicated(x[c("ISO3", "year", "Age")]),
          !any(grepl(">", x$value_text)))
child <- x[x$Age == "Age 0-14", c("ISO3", "year", "Country/Region", "UNICEF Region",
                                    "value_text", "upper_limit", "plot_value")]
names(child) <- c("iso3", "year", "country", "region", "child_source_value",
                  "child_upper_limit", "child_plot_value")
adolescent <- x[x$Age == "Age 15-19", c("ISO3", "year", "value_text", "upper_limit", "plot_value")]
names(adolescent) <- c("iso3", "year", "adolescent_source_value",
                      "adolescent_upper_limit", "adolescent_plot_value")
pairs <- merge(child, adolescent, by = c("iso3", "year"))
valid <- with(pairs, is.finite(child_plot_value) & child_plot_value > 0 &
                is.finite(adolescent_plot_value) & adolescent_plot_value > 0)
# Missing source entries are '.', not zero. Do not add pseudocounts for log axes.
stopifnot(!any(pairs$child_plot_value == 0, na.rm = TRUE),
          !any(pairs$adolescent_plot_value == 0, na.rm = TRUE),
          !any(pairs$iso3 == "NGA"))
excluded <- pairs[!valid, ]
pairs <- pairs[valid, ]
pairs$region_group <- ifelse(pairs$region %in% c("West and Central Africa", "Eastern and Southern Africa"),
                             pairs$region, "Other regions")
pairs$estimate_type <- ifelse(pairs$child_upper_limit | pairs$adolescent_upper_limit,
                              "At least one upper limit", "Point estimates")
export <- pairs
names(export) <- sub("plot_value$", if (incidence) "plot_rate" else "plot_count", names(export))
cbh_atomic_csv(export, file.path(out, paste0("paired_country_year_", metric, ".csv")))
if (incidence) cbh_atomic_csv(excluded, file.path(out, "incidence_missing_pairs.csv"))
latest <- max(pairs$year)
recent <- pairs[pairs$year == latest, ]
panels <- c(sprintf("All years: %s-%s\n%s country-years; %s countries",
                    min(pairs$year), latest, nrow(pairs), length(unique(pairs$iso3))),
            sprintf("%s snapshot\n%s countries", latest, nrow(recent)))
pairs$panel <- panels[1]
recent$panel <- panels[2]
d <- rbind(pairs, recent)
d$panel <- factor(d$panel, levels = panels)
ticks <- if (incidence) c(.01, .1, 1, 10, 100) else c(100, 1000, 10000, 100000, 1000000)
labels <- if (incidence) c("0.01", "0.1", "1", "10", "100") else c("100", "1,000", "10,000", "100,000", "1,000,000")
limits <- range(c(d$child_plot_value, d$adolescent_plot_value)) * c(.8, 1.25)
p <- ggplot(d, aes(adolescent_plot_value, child_plot_value)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#9AA3AA", linewidth = .5) +
  geom_point(aes(colour = region_group, shape = estimate_type), alpha = .55, size = 2) +
  facet_wrap(~panel, ncol = 2) +
  scale_x_log10(breaks = ticks, labels = labels, limits = limits) +
  scale_y_log10(breaks = ticks, labels = labels, limits = limits) +
  scale_colour_manual(values = c("Eastern and Southern Africa" = "#C26736",
                                 "West and Central Africa" = "#176B87", "Other regions" = "#858D96")) +
  scale_shape_manual(values = c("Point estimates" = 16, "At least one upper limit" = 2),
                     breaks = c("Point estimates", "At least one upper limit")) +
  coord_fixed() +
  labs(title = if (incidence) "Adolescent versus child HIV incidence" else "Adolescent versus child HIV: existing country estimates",
       subtitle = if (incidence) "New infections per 1,000 uninfected population, both sexes | Logarithmic axes" else
         "Numbers living with HIV, both sexes | Counts, not prevalence | Logarithmic axes",
       x = if (incidence) "HIV incidence rate, ages 15-19" else "Adolescents aged 15-19 living with HIV",
       y = if (incidence) "HIV incidence rate, ages 0-14" else "Children aged 0-14 living with HIV",
       colour = NULL, shape = NULL,
       caption = if (incidence) paste0(
         "Source: local UNICEF / UNAIDS 2025 estimates workbook; reported incidence rates used directly. Dashed line: equal rates.\n",
         "Open triangles show at least one rate reported as <0.01, plotted at that upper limit. Annual points repeat countries.\n",
         "Nigeria has no paired child series and is absent. ", nrow(excluded), " country-year pairs with missing rates omitted. Incidence is distinct from prevalence.") else
         paste0("Source: local UNICEF / UNAIDS 2025 estimates workbook. Dashed line: equal counts in the two age groups.\n",
                "Open triangles place censored entries at their reported upper limit (e.g. <100 at 100); these are not point estimates.\n",
                "Nigeria has no paired child series and is absent. Country population size contributes to the relationship; annual points repeat countries.")) +
  guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(alpha = 1)),
         shape = guide_legend(order = 2, nrow = 1, override.aes = list(alpha = 1))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#EEF0F2"),
        strip.text = element_text(face = "bold", size = 12, hjust = 0),
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(colour = "#52616B", margin = margin(b = 12)),
        plot.caption = element_text(hjust = 0, size = 9, colour = "#52616B", lineheight = 1.2),
        legend.position = "bottom", legend.box = "vertical",
        plot.margin = margin(15, 20, 12, 15))
path <- file.path(out, paste0("adolescent_vs_child_hiv_", metric, ".png"))
ggsave(path, p, device = ragg::agg_png, width = 12, height = 8.5, units = "in", dpi = 180, bg = "white")
cat("Saved", path, "\n")
print(data.frame(country_years = nrow(pairs), countries = length(unique(pairs$iso3)),
                 exact_pairs = sum(pairs$estimate_type == "Point estimates"),
                 limit_pairs = sum(pairs$estimate_type != "Point estimates"),
                 latest_year = latest, latest_countries = nrow(recent)))

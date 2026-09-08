#!/usr/bin/env Rscript
# Descriptive country-year HIV count pairs; no imputation or model fitting.
# Run from the project root. Uses only the existing UNICEF workbook.
library(ggplot2)
source("R_cbh/load_pipeline.R")
out <- "results/cbh/hiv_age_pairs"
source_file <- "data/HIV_Epidemiology_Children_Adolescents_2025.xlsx"
x <- as.data.frame(suppressMessages(readxl::read_excel(
  source_file, sheet = "Data", skip = 1, guess_max = 100000)))
x <- x[x$Type == "Country" & x$Sex == "Both" &
         x$Indicator == "Estimated number of people living with HIV" &
         x$Age %in% c("Age 0-14", "Age 15-19"), ]
x$year <- as.integer(x$Year)
x$value_text <- trimws(as.character(x$Value))
x$upper_limit <- grepl("^<", x$value_text)
x$plot_count <- suppressWarnings(as.numeric(gsub("[<,[:space:]]", "", x$value_text)))
stopifnot(!anyDuplicated(x[c("ISO3", "year", "Age")]),
          !any(grepl(">", x$value_text)))
child <- x[x$Age == "Age 0-14", c("ISO3", "year", "Country/Region", "UNICEF Region",
                                    "value_text", "upper_limit", "plot_count")]
names(child) <- c("iso3", "year", "country", "region", "child_source_value",
                  "child_upper_limit", "child_plot_count")
adolescent <- x[x$Age == "Age 15-19", c("ISO3", "year", "value_text", "upper_limit", "plot_count")]
names(adolescent) <- c("iso3", "year", "adolescent_source_value",
                      "adolescent_upper_limit", "adolescent_plot_count")
pairs <- merge(child, adolescent, by = c("iso3", "year"))
valid <- with(pairs, is.finite(child_plot_count) & child_plot_count > 0 &
                is.finite(adolescent_plot_count) & adolescent_plot_count > 0)
stopifnot(all(valid), !any(pairs$iso3 == "NGA"))
pairs$region_group <- ifelse(pairs$region %in% c("West and Central Africa", "Eastern and Southern Africa"),
                             pairs$region, "Other regions")
pairs$estimate_type <- ifelse(pairs$child_upper_limit | pairs$adolescent_upper_limit,
                              "At least one upper limit", "Point estimates")
cbh_atomic_csv(pairs, file.path(out, "paired_country_year_counts.csv"))
latest <- max(pairs$year)
recent <- pairs[pairs$year == latest, ]
panels <- c(sprintf("All years: %s-%s\n%s country-years; %s countries",
                    min(pairs$year), latest, nrow(pairs), length(unique(pairs$iso3))),
            sprintf("%s snapshot\n%s countries", latest, nrow(recent)))
pairs$panel <- panels[1]
recent$panel <- panels[2]
d <- rbind(pairs, recent)
d$panel <- factor(d$panel, levels = panels)
counts <- c(100, 1000, 10000, 100000, 1000000)
labels <- c("100", "1,000", "10,000", "100,000", "1,000,000")
limits <- range(c(d$child_plot_count, d$adolescent_plot_count)) * c(.8, 1.25)
p <- ggplot(d, aes(adolescent_plot_count, child_plot_count)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#9AA3AA", linewidth = .5) +
  geom_point(aes(colour = region_group, shape = estimate_type), alpha = .55, size = 2) +
  facet_wrap(~panel, ncol = 2) +
  scale_x_log10(breaks = counts, labels = labels, limits = limits) +
  scale_y_log10(breaks = counts, labels = labels, limits = limits) +
  scale_colour_manual(values = c("Eastern and Southern Africa" = "#C26736",
                                 "West and Central Africa" = "#176B87", "Other regions" = "#858D96")) +
  scale_shape_manual(values = c("Point estimates" = 16, "At least one upper limit" = 2),
                     breaks = c("Point estimates", "At least one upper limit")) +
  coord_fixed() +
  labs(title = "Adolescent versus child HIV: existing country estimates",
       subtitle = "Numbers living with HIV, both sexes | Counts, not prevalence | Logarithmic axes",
       x = "Adolescents aged 15-19 living with HIV",
       y = "Children aged 0-14 living with HIV",
       colour = NULL, shape = NULL,
       caption = paste0("Source: local UNICEF / UNAIDS 2025 estimates workbook. Dashed line: equal counts in the two age groups.\n",
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
path <- file.path(out, "adolescent_vs_child_hiv_counts.png")
ggsave(path, p, device = ragg::agg_png, width = 12, height = 8.5, units = "in", dpi = 180, bg = "white")
cat("Saved", path, "\n")
print(data.frame(country_years = nrow(pairs), countries = length(unique(pairs$iso3)),
                 exact_pairs = sum(pairs$estimate_type == "Point estimates"),
                 limit_pairs = sum(pairs$estimate_type != "Point estimates"),
                 latest_year = latest, latest_countries = nrow(recent)))

#!/usr/bin/env Rscript
# Compare national under-five annual death counts; no refitting or downloads.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
library(ggplot2)
library(patchwork)
spec <- cbh_trial_spec()
out <- file.path("results/cbh", spec$id, "country_burden_2024")
model_path <- file.path(out, "country_totals_2024.csv")
ihme_path <- "data/ihme_malaria_u5_deaths_by_age_country_year.csv"
model <- cbh_read_csv(model_path)
ihme <- cbh_read_csv(ihme_path)
ihme <- ihme[ihme$Year == 2024 & ihme$Sex == "Both" & ihme$Age == "Under 5" &
             ihme$Condition == "Malaria" & ihme$Measure == "Deaths" & ihme$Unit == "Number", ]
ihme$iso3 <- countrycode::countrycode(ihme$Location, "country.name", "iso3c", warn = FALSE)
cbh_unique(ihme, "iso3", "2024 IHME malaria counts")
ihme <- ihme[c("iso3", "Value", "Lower", "Upper")]
names(ihme)[-1] <- c("ihme_malaria_deaths", "ihme_malaria_lower", "ihme_malaria_upper")
z <- merge(model, ihme, by = "iso3", all = TRUE)
z$comparison_status <- ifelse(!is.finite(z$attributable_under5_deaths), "missing_model_MAP",
                              ifelse(!is.finite(z$ihme_malaria_deaths), "missing_IHME", "matched"))
z$model_to_ihme_ratio <- z$attributable_under5_deaths / z$ihme_malaria_deaths
z$difference_deaths <- z$attributable_under5_deaths - z$ihme_malaria_deaths
z$year <- 2024L; z$age <- "Under 5"; z$sex <- "Both"
cbh_atomic_csv(z, file.path(out, "model_vs_ihme_malaria_2024.csv"))
matched <- z[z$comparison_status == "matched", ]
stopifnot(nrow(matched) == sum(model$status == "estimated"),
          all(matched$attributable_under5_deaths > 0), all(matched$ihme_malaria_deaths > 0))
matched$coverage <- factor(ifelse(matched$map_coverage_below_95pct,
  "MAP coverage <95%", "MAP coverage >=95%"), levels = c("MAP coverage >=95%", "MAP coverage <95%"))
palette <- c("MAP coverage >=95%" = "#16747C", "MAP coverage <95%" = "#CC6A30")
theme_set(theme_minimal(base_size = 12))
sty <- theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
             plot.caption = element_text(hjust = 0), legend.position = "bottom")
axis_label <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
limit <- ceiling(max(matched$attributable_under5_deaths, matched$ihme_malaria_deaths)/25000)*25000
linear <- ggplot(matched, aes(ihme_malaria_deaths, attributable_under5_deaths)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(aes(colour = coverage), size = 2.5, alpha = .85) +
  geom_text(data = matched[matched$iso3 %in% c("NGA", "COD", "NER"), ], aes(label = iso3),
            nudge_y = 6500, size = 3.4, colour = "#253746") +
  scale_colour_manual(values = palette, drop = FALSE) +
  scale_x_continuous(labels = function(x) paste0(x/1000, "k"), limits = c(0, limit)) +
  scale_y_continuous(labels = function(x) paste0(x/1000, "k"), limits = c(0, limit)) +
  coord_fixed() + labs(title = "A  Absolute death counts", x = "IHME malaria deaths", y = "Our malaria-attributable deaths", colour = NULL) + sty
logplot <- ggplot(matched, aes(ihme_malaria_deaths, attributable_under5_deaths)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(aes(colour = coverage), size = 2.5, alpha = .85) +
  geom_text(data = matched[matched$iso3 %in% c("NGA", "COD", "NER", "SWZ", "BWA", "COM", "NAM", "DJI", "MRT", "ERI", "RWA", "ZWE", "KEN", "LBR"), ],
    aes(label = iso3), nudge_y = .13, size = 3, check_overlap = TRUE, colour = "#253746") +
  scale_colour_manual(values = palette, drop = FALSE) +
  scale_x_log10(breaks = c(1,10,100,1000,10000,100000), labels = axis_label, limits = c(.3,250000)) +
  scale_y_log10(breaks = c(1,10,100,1000,10000,100000), labels = axis_label, limits = c(.3,250000)) +
  coord_fixed() + labs(title = "B  All countries on logarithmic axes", x = "IHME malaria deaths (log scale)",
                     y = "Our malaria-attributable deaths (log scale)", colour = NULL) + sty
scatter <- (linear | logplot) + plot_layout(guides = "collect") +
  plot_annotation(title = "Under-five malaria mortality by country, 2024",
    subtitle = sprintf("%d matched countries | Dashed line: equal estimates | Above the line: our estimate is higher", nrow(matched)),
    caption = paste("Our estimate: all-cause mortality reduction under national PfPR2-10 -> 0%. IHME: cause-specific malaria deaths.",
      "Point estimates only; joint uncertainty for model country totals is not available. Orange points have limited MAP population coverage.",
      "Ages 2-4 use equal baseline rates and equal death/person-time shares. Three countries lack model estimates because MAP is missing.", sep = "\n"),
    theme = theme(plot.title = element_text(face = "bold", size = 17), plot.caption = element_text(hjust = 0, size = 10)))
scatter <- scatter & theme(legend.position = "bottom")
ggsave(file.path(out, "model_vs_ihme_malaria_2024.png"), scatter, device = ragg::agg_png,
       width = 14, height = 8.5, dpi = 160, bg = "white")

# A second chart labels every matched country, avoiding scatter-label overlap.
matched <- matched[order(matched$ihme_malaria_deaths), ]
matched$country_label <- paste0(matched$country, ifelse(matched$map_coverage_below_95pct, " *", ""))
matched$country_label <- factor(matched$country_label, levels = matched$country_label)
long <- rbind(data.frame(country_label = matched$country_label, deaths = matched$ihme_malaria_deaths, source = "IHME malaria"),
              data.frame(country_label = matched$country_label, deaths = matched$attributable_under5_deaths, source = "Our attributable estimate"))
by_country <- ggplot(matched, aes(y = country_label)) +
  geom_segment(aes(x = ihme_malaria_deaths, xend = attributable_under5_deaths, yend = country_label), colour = "grey65") +
  geom_point(data = long, aes(x = deaths, colour = source, shape = source), size = 2.5) +
  scale_colour_manual(values = c("IHME malaria" = "#253746", "Our attributable estimate" = "#16747C")) +
  scale_shape_manual(values = c("IHME malaria" = 16, "Our attributable estimate" = 17)) +
  scale_x_log10(breaks = c(1,10,100,1000,10000,100000), labels = axis_label) +
  labs(title = "Country-by-country comparison, 2024", subtitle = "Under-five annual death counts; countries ordered by IHME estimate",
    x = "Deaths (log scale)", y = NULL, colour = NULL, shape = NULL,
    caption = "* MAP covers less than 95% of population weight within the available raster footprint.\nModel totals are signed sums, including negative age-band contributions. Point estimates only.") + sty
ggsave(file.path(out, "model_vs_ihme_malaria_by_country_2024.png"), by_country, device = ragg::agg_png,
       width = 11, height = 14, dpi = 160, bg = "white")
special <- z[z$iso3 %in% c("COD", "NGA"), ]
rows <- vapply(seq_len(nrow(special)), function(i) sprintf("| %s | %s | %s | %.2f |",
  special$country[i], format(round(special$attributable_under5_deaths[i]), big.mark = ","),
  format(round(special$ihme_malaria_deaths[i]), big.mark = ","), special$model_to_ihme_ratio[i]), "")
writeLines(c("# Model versus IHME malaria mortality, 2024", "",
  sprintf("The comparison matches %d countries on ISO3, 2024, both sexes and under-five age. The x-axis uses IHME **malaria** deaths, not IHME all-cause deaths. All %d matched counts are positive and appear on both the linear and logarithmic panels.", nrow(matched), nrow(matched)), "",
  "Our estimate is the signed reduction in all-cause mortality predicted under zero national PfPR. IHME reports cause-specific malaria deaths. These are related but different estimands: our fitted association can reflect indirect effects and residual confounding. They should not be interpreted as interchangeable measurements or as independent validation against observed deaths.", "",
  "| Country | Our attributable deaths | IHME malaria deaths | Model / IHME |",
  "|---|---:|---:|---:|", rows, "",
  sprintf("Across the same %d countries, our point estimates sum to %s deaths versus %s for IHME (ratio %.2f).",
    nrow(matched), format(round(sum(matched$attributable_under5_deaths)), big.mark = ","),
    format(round(sum(matched$ihme_malaria_deaths)), big.mark = ","),
    sum(matched$attributable_under5_deaths)/sum(matched$ihme_malaria_deaths)), "",
  "Cape Verde, Lesotho and São Tomé and Príncipe lack model estimates because usable MAP prevalence is missing. Lesotho also has no entry in this malaria export. These countries remain in the comparison CSV with missingness status; missing values are never replaced by zero. Partial MAP coverage is flagged in colour and with an asterisk on the labelled chart.", "",
  "Both charts show point estimates only. IHME's source uncertainty bounds are retained in the CSV. The current pipeline does not provide joint country-total model intervals, and marginal age-band interval endpoints have not been summed. Ages 2-4 retain the authorized equal-rate/equal-person-time assumption; negative age-specific effects remain in model totals.", "",
  "The malaria comparator is the existing `data/ihme_malaria_u5_deaths_by_age_country_year.csv` export. Our denominator uses the newer all-cause export dated 2026-09-09. Their precise IHME release/version alignment has not been verified. See the [burden report](REPORT.md) for exposure and inference limitations.", "",
  "![Scatter comparison](model_vs_ihme_malaria_2024.png)", "",
  "![Every matched country](model_vs_ihme_malaria_by_country_2024.png)", "",
  "[Comparison data](model_vs_ihme_malaria_2024.csv)"), file.path(out, "IHME_COMPARISON.md"))
files <- c(model_path, ihme_path, "R_cbh/burden/03_compare_ihme_2024.R")
cbh_atomic_csv(data.frame(file = files, md5 = vapply(files, cbh_file_hash, "")), file.path(out, "comparison_provenance.csv"))
print(special[c("iso3", "attributable_under5_deaths", "ihme_malaria_deaths", "model_to_ihme_ratio")], row.names = FALSE)
message("Verified and plotted ", nrow(matched), " matched countries.")

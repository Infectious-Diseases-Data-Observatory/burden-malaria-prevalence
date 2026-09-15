# =============================================================================
# 44_deaths_by_age.R — distribution of age at death among under-5 deaths.
#
# Weighted and unweighted counts of deaths by age in completed months (0-59)
# for deaths in the 60 months before interview, across the surveys in the
# panel: the deaths the person-time models (scripts 40-43) are fitted to. Ages
# come from b7 (age at death in months) as recorded. DHS records age at death
# in years from the second birthday, and b7 is then years x 12, so every death
# from 24 months up sits at 24, 36 or 48 exactly; the spike at 12 mixes deaths
# reported as "one year" with those at exactly 12 months. Those are coding
# facts, not heaping to smooth away, and they are why the person-time model's
# age segments are 12-23, 24-35, 36-47 and 48-59 months.
#
# Outputs (results/dhs_rebuild)
#   deaths_by_age_month.csv    weighted and unweighted deaths by month of age
#   deaths_by_age_deciles.csv  weighted deciles of age at death, plus the same by age segment
#   figure25_deaths_by_age.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages("ggplot2")

CACHE <- file.path(DATA_DIR, "deaths_by_age_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
WINDOW_MONTHS <- 60L

registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
map <- read.csv(MAP_REGION_CSV, stringsAsFactors = FALSE)
registry <- registry[registry$svkey %in% map$svkey, , drop = FALSE]

survey_counts <- function(survey) {
  cache <- file.path(CACHE, paste0(survey$svkey, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  br <- readRDS(as.character(survey$local_recode))
  keep <- !is.na(br$b7) & br$b7 < 60 & !is.na(br$v005) & br$v005 > 0 &
    (br$b3 + br$b7) >= (br$v008 - WINDOW_MONTHS) & (br$b3 + br$b7) < br$v008
  age <- as.integer(br$b7[keep])
  weight <- br$v005[keep] / 1e6
  out <- data.frame(svkey = survey$svkey, iso3 = survey$iso3, age_months = 0:59,
                    deaths_w = as.numeric(tapply(weight, factor(age, levels = 0:59), sum)),
                    deaths_n = as.integer(table(factor(age, levels = 0:59))),
                    stringsAsFactors = FALSE)
  out$deaths_w[is.na(out$deaths_w)] <- 0
  rm(br); gc(FALSE)
  saveRDS(out, cache)
  out
}

rows <- list()
for (i in seq_len(nrow(registry))) {
  survey <- registry[i, , drop = FALSE]
  rows[[survey$svkey]] <- tryCatch(survey_counts(survey), error = function(e) {
    message("  skip ", survey$svkey, ": ", conditionMessage(e)); NULL })
  if (i %% 20 == 0) message("  ", i, " surveys")
}
by_survey <- do.call(rbind, rows)
by_age <- aggregate(cbind(deaths_w, deaths_n) ~ age_months, by_survey, sum)
by_age$share_w <- by_age$deaths_w / sum(by_age$deaths_w)
by_age$cumulative_share_w <- cumsum(by_age$share_w)
by_age$segment <- cut(by_age$age_months, c(-1, 0, 2, 5, 11, 23, 35, 47, 59),
                      labels = c("0", "1-2", "3-5", "6-11", "12-23", "24-35", "36-47", "48-59"))
write.csv(by_age, file.path(RESULTS_DIR, "deaths_by_age_month.csv"), row.names = FALSE)

# Weighted deciles of age at death. Ages are whole months, so a decile is the
# first month at which the cumulative share reaches it; the interpolated value
# treats deaths as spread across the month.
deciles <- data.frame(decile = seq(10, 90, by = 10))
deciles$age_months <- vapply(deciles$decile / 100, function(q)
  by_age$age_months[which(by_age$cumulative_share_w >= q)[1]], numeric(1))
deciles$age_months_interpolated <- vapply(deciles$decile / 100, function(q) {
  i <- which(by_age$cumulative_share_w >= q)[1]
  below <- if (i > 1) by_age$cumulative_share_w[i - 1] else 0
  by_age$age_months[i] + (q - below) / by_age$share_w[i]
}, numeric(1))
write.csv(deciles, file.path(RESULTS_DIR, "deaths_by_age_deciles.csv"), row.names = FALSE)

by_segment <- aggregate(cbind(deaths_w, deaths_n) ~ segment, by_age, sum)
by_segment$share_w <- by_segment$deaths_w / sum(by_segment$deaths_w)
message(sprintf("%d surveys; %.0f weighted (%d recorded) under-5 deaths in the 60 months before interview",
                length(rows), sum(by_age$deaths_w), sum(by_age$deaths_n)))
message("\nDeaths by DHS age segment (weighted share):")
print(transform(by_segment, deaths_w = round(deaths_w), share_w = sprintf("%.1f%%", 100 * share_w)),
      row.names = FALSE)
message("\nWeighted deciles of age at death (months):")
print(transform(deciles, age_months_interpolated = round(age_months_interpolated, 1)), row.names = FALSE)
message(sprintf("Weighted mean age at death %.1f months; share in the first month %.1f%%; months 12/24/36/48 (year-coded ages) hold %.1f%% of deaths",
                sum(by_age$age_months * by_age$share_w), 100 * by_age$share_w[1],
                100 * sum(by_age$share_w[by_age$age_months %in% c(12, 24, 36, 48)])))

## ---- figure -------------------------------------------------------------------------------
# Deciles that fall in the same month share one label.
decile_labels <- aggregate(decile ~ age_months, deciles, function(d)
  if (length(d) > 1) sprintf("D%d-D%d", min(d) / 10, max(d) / 10) else sprintf("D%d", d / 10))
decile_labels$x <- vapply(decile_labels$age_months, function(m)
  mean(deciles$age_months_interpolated[deciles$age_months == m]), numeric(1))
plot <- ggplot2::ggplot(by_age, ggplot2::aes(age_months, deaths_w / 1e3)) +
  ggplot2::geom_col(fill = "#1D6F8B", width = 0.9) +
  ggplot2::geom_vline(xintercept = deciles$age_months_interpolated, linetype = "dotted",
                      colour = "grey35") +
  ggplot2::annotate("text", x = decile_labels$x + 0.4, y = max(by_age$deaths_w) / 1e3 * 0.99,
                    label = decile_labels$decile, size = 2.7, colour = "grey30",
                    hjust = 0, vjust = 1) +
  ggplot2::scale_x_continuous(breaks = c(0, 6, 12, 24, 36, 48, 59)) +
  ggplot2::labs(x = "Age at death (completed months, as recorded)",
                y = "Weighted deaths (thousands)",
                title = "Under-5 deaths by age at death, 60 months before interview",
                subtitle = sprintf(paste0(
                  "%d DHS/MIS surveys, %.0f thousand weighted deaths; dotted lines are the weighted deciles.\n",
                  "Ages from the second birthday are recorded in years, so those deaths sit at 24, 36 and 48;\n",
                  "the spike at 12 mixes 'one year' with exactly 12 months"),
                  length(rows), sum(by_age$deaths_w) / 1e3)) +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(RESULTS_DIR, "figure25_deaths_by_age.png"), plot,
                width = 10, height = 4.6, dpi = 200)
message("\nWrote deaths_by_age_month.csv, deaths_by_age_deciles.csv and figure25_deaths_by_age.png")

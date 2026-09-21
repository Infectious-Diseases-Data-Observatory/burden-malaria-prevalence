#!/usr/bin/env Rscript
# Four-source summary of malaria mortality in 2000, 2004, 2015 and 2024, with
# declines to 2024. Reads saved aggregates only; no refitting and no downloads.
# The three under-five sources share one denominator; the World Malaria Report
# is all ages for the WHO African Region and is NOT on the same footing.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
settings <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION", "regional"))
out <- file.path(settings$out, "burden_comparison")
years <- c(2000L, 2004L, 2015L, 2024L)

annual <- cbh_read_csv(file.path(settings$out, "annual_comparison/annual_totals_2000_2024.csv"))
annual <- annual[match(years, annual$year), ]
wmr_af <- cbh_read_csv("results/who_wmr2025_africa_deaths.csv")
wmr_gl <- cbh_read_csv("results/who_wmr2025_global_deaths.csv")
wmr_af <- wmr_af[match(years, wmr_af$year), ]
wmr_gl <- wmr_gl[match(years, wmr_gl$year), ]

# Mortality rates printed in the World Malaria Report 2025 narrative text
# (printed pages 27 and 37), deaths per 100 000 population at risk, all ages.
# The report labels only selected years; 2004 is not published as a number.
wmr_rate_af <- c(140, NA, 62.5, 51.9)
wmr_rate_gl <- c(28.6, NA, 14.9, 13.8)

row <- function(source, scope, unit, v)
  data.frame(source = source, scope = scope, unit = unit,
             y2000 = v[1], y2004 = v[2], y2015 = v[3], y2024 = v[4],
             decline_from_2000 = 100 * (1 - v[4] / v[1]),
             decline_from_2004 = 100 * (1 - v[4] / v[2]),
             decline_from_2015 = 100 * (1 - v[4] / v[3]), row.names = NULL)

u5 <- "under-5, 42 countries"; af <- "all ages, WHO African Region"; gl <- "all ages, global"
tab <- rbind(
  row("PfPR-ACM model", u5, "deaths", annual$model_deaths),
  row("IHME",           u5, "deaths", annual$ihme_malaria_deaths),
  row("UN IGME",        u5, "deaths", annual$who_cacode_deaths),
  row("WHO WMR 2025",   af, "deaths", wmr_af$point),
  row("WHO WMR 2025",   gl, "deaths", wmr_gl$point),
  row("PfPR-ACM model", u5, "deaths per 1000 child-years", annual$model_rate_per100000 / 100),
  row("IHME",           u5, "deaths per 1000 child-years", annual$ihme_malaria_rate_per100000 / 100),
  row("UN IGME",        u5, "deaths per 1000 child-years", annual$who_cacode_rate_per100000 / 100),
  row("WHO WMR 2025",   af, "deaths per 100 000 population at risk", wmr_rate_af),
  row("WHO WMR 2025",   gl, "deaths per 100 000 population at risk", wmr_rate_gl))
cbh_atomic_csv(tab, file.path(out, "source_comparison_2000_2024.csv"))

fmt <- function(v, unit) if (grepl("^deaths$", unit)) format(round(v), big.mark = ",") else
  ifelse(is.na(v), "--", sprintf("%.2f", v))
pc <- function(v) ifelse(is.na(v), "--", sprintf("%+.1f%%", -v))
md <- c("# Malaria mortality by source, 2000-2024", "",
  "Declines are percentage reductions in the stated quantity between that year and 2024; a negative sign is an increase.",
  "The first three sources are under-five deaths in the same 42 countries on a shared denominator.",
  "The World Malaria Report is all ages and a different geography, so its values are context, not a like-for-like comparison.",
  "The report publishes mortality rates only for selected years; 2004 is not among them.",
  "It publishes no death rate per 1000 and no under-five-specific mortality rate; per-1000 death quantities appear only as unquantified model inputs in the Annex 1 methods.",
  "Its denominator, population at risk, is the population of high-endemic areas plus half that of low-endemic areas, with the risk proportion held constant from 2000 to 2024, so it is neither total population nor a child population.", "",
  "| Source | Scope | Unit | 2000 | 2004 | 2015 | 2024 | From 2000 | From 2004 | From 2015 |",
  "|---|---|---|---:|---:|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    r$source, r$scope, r$unit,
    fmt(r$y2000, r$unit), fmt(r$y2004, r$unit), fmt(r$y2015, r$unit), fmt(r$y2024, r$unit),
    sprintf("%.1f%%", r$decline_from_2000),
    ifelse(is.na(r$decline_from_2004), "--", sprintf("%.1f%%", r$decline_from_2004)),
    sprintf("%.1f%%", r$decline_from_2015)))
}
md <- c(md, "", paste("Sources: annual_comparison/annual_totals_2000_2024.csv for the under-five series;",
  "results/who_wmr2025_{africa,global}_deaths.csv for World Malaria Report 2025 death counts;",
  "mortality rates transcribed from the World Malaria Report 2025 narrative text, printed pages 27 and 37."))
writeLines(md, file.path(out, "SOURCE_COMPARISON.md"))
# LaTeX fragment: sources as columns, 2000/2015/2024 values and changes as rows.
# Plain text only; no .tex file is ever written by this pipeline.
pick <- function(src, scope, unit) {
  r <- tab[tab$source == src & tab$scope == scope & tab$unit == unit, ]
  stopifnot(nrow(r) == 1L); r
}
deaths <- list(pick("PfPR-ACM model", u5, "deaths"), pick("IHME", u5, "deaths"),
               pick("UN IGME", u5, "deaths"))
rates  <- list(pick("PfPR-ACM model", u5, "deaths per 1000 child-years"),
               pick("IHME", u5, "deaths per 1000 child-years"),
               pick("UN IGME", u5, "deaths per 1000 child-years"))
num  <- function(x) format(round(x), big.mark = ",", trim = TRUE)
rat  <- function(x) sprintf("%.2f", x)
# The report prints its rates to at most one decimal (140, 62.5); match that.
wmr  <- function(x) sub("\\.0$", "", sprintf("%.1f", x))
# tab stores declines (positive = reduction); print signed change instead.
chg  <- function(x) sprintf("$%s%.1f\\%%$", ifelse(-x < 0, "-", "+"), abs(x))
line <- function(label, vals) paste0("\\quad ", label, " & ", paste(vals, collapse = " & "), " \\\\")
body <- c(
  "\\multicolumn{4}{@{}l}{\\textit{Malaria deaths}} \\\\",
  line("2000", vapply(deaths, function(r) num(r$y2000), "")),
  line("2015", vapply(deaths, function(r) num(r$y2015), "")),
  line("2024", vapply(deaths, function(r) num(r$y2024), "")),
  line("Change 2000--2024", vapply(deaths, function(r) chg(r$decline_from_2000), "")),
  line("Change 2015--2024", vapply(deaths, function(r) chg(r$decline_from_2015), "")),
  "\\addlinespace",
  "\\multicolumn{4}{@{}l}{\\textit{Malaria mortality rate}} \\\\",
  line("2000", vapply(rates, function(r) rat(r$y2000), "")),
  line("2015", vapply(rates, function(r) rat(r$y2015), "")),
  line("2024", vapply(rates, function(r) rat(r$y2024), "")),
  line("Change 2000--2024", vapply(rates, function(r) chg(r$decline_from_2000), "")),
  line("Change 2015--2024", vapply(rates, function(r) chg(r$decline_from_2015), "")))
tex <- c("% Standalone table fragment; uses booktabs (already in the manuscript preamble).",
  "% Generated by R_cbh/reporting/12_source_comparison_table.R from saved aggregates.",
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{Malaria mortality by source in 2000, 2015 and 2024, and change to 2024. Negative values are reductions.}",
  "\\label{tab:source-comparison}", "\\small",
  "\\setlength{\\tabcolsep}{3pt}", "\\renewcommand{\\arraystretch}{1.15}",
  "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lrrr@{}}", "\\toprule",
  paste("&", "\\shortstack{PfPR-ACM\\\\model}", "& IHME &", "\\shortstack{UN\\\\IGME}", "\\\\"),
  "\\midrule", body, "\\bottomrule", "\\end{tabular*}",
  "\\par\\vspace{0.5em}", "\\begin{minipage}{\\linewidth}", "\\footnotesize",
  paste("All three columns are deaths before age five in the same 42 countries. Rates use the same annual",
        "under-five person-years, implied by the IHME all-cause death count and rate, so within a year the",
        "three sources share one denominator."),
  paste("Model deaths are IHME all-cause deaths in each age band multiplied by the estimated malaria-attributable",
        "fraction; IHME and UN IGME are cause-specific malaria death estimates. These are different estimands."),
  paste("Counts and rates move in opposite directions after 2015 because under-five person-time in these",
        "countries grew by 12% between 2015 and 2024."),
  "\\end{minipage}", "\\end{table}")
writeLines(tex, file.path(out, "source_comparison.latex.txt"))

inputs <- c(file.path(settings$out, "annual_comparison/annual_totals_2000_2024.csv"),
  "results/who_wmr2025_africa_deaths.csv", "results/who_wmr2025_global_deaths.csv",
  "R_cbh/reporting/12_source_comparison_table.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")),
  file.path(out, "source_comparison_provenance.csv"))
cat("Source comparison table written:", file.path(out, "SOURCE_COMPARISON.md"), "\n")

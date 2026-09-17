#!/usr/bin/env Rscript
# Present the plan's national/annual/state estimates; no new models or TeX files.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
settings <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional"))
root <- settings$out;out <- file.path(root,"tables")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
paths <- file.path(root,c("annual_comparison/country_estimates_2004_2024.csv",
  "annual_comparison/annual_totals_2004_2024.csv","nigeria_states/state_totals_2024.csv"))
x <- cbh_read_csv(paths[1]);x <- x[x$year==2024,]
cbh_unique(x,"iso3","2024 national estimates")
stopifnot(nrow(x)==42,all(is.finite(x$model_deaths)),all(is.finite(x$ihme_malaria_deaths)),
  all(is.finite(x$who_cacode_deaths)))
x$difference_from_ihme <- x$model_deaths-x$ihme_malaria_deaths
x <- x[order(-abs(x$difference_from_ihme),x$iso3),]
x$rank_absolute_difference_from_ihme <- seq_len(nrow(x))
cols <- c("rank_absolute_difference_from_ihme","iso3","country","year","model_deaths",
  "ihme_malaria_deaths","who_cacode_deaths","difference_from_ihme")
cbh_atomic_csv(x[cols],file.path(out,"country_comparison_2024.csv"))
fmt <- function(z)format(round(z),big.mark=",",scientific=FALSE,trim=TRUE)
total <- colSums(x[c("model_deaths","ihme_malaria_deaths","who_cacode_deaths")])
header <- paste0("| Country | ",cbh_paper_model_label()," | IHME | UN IGME |")
body <- function(d)sprintf("| %s | %s | %s | %s |",d$country,fmt(d$model_deaths),
  fmt(d$ihme_malaria_deaths),fmt(d$who_cacode_deaths))
total_row <- sprintf("| **All 42 countries** | **%s** | **%s** | **%s** |",fmt(total[1]),fmt(total[2]),fmt(total[3]))
note <- "Deaths before age five in 2024. Countries are ordered by the absolute difference between the PfPR-ACM model and IHME; the total always includes all 42 estimable countries. UN IGME is the direct under-five CA-CODE 2026 series, not a WHO all-age proxy. Model values use the revised 18-variable regional-adjustment MAP gamma=2 fits. Point estimates are rounded to the nearest death; totals are calculated before rounding. These sources estimate different quantities; their agreement is descriptive."
writeLines(c("# Country comparison, 2024","",header,"|---|---:|---:|---:|",body(x),total_row,"",note),
  file.path(out,"country_comparison_2024.md"))
writeLines(c("# Ten largest country differences, 2024","",header,"|---|---:|---:|---:|",body(x[1:10,]),total_row,"",note),
  file.path(out,"country_comparison_2024_top10.md"))
# Match the existing manuscript table's selection, with current source labels.
# Plain-text LaTeX source lets the author update embedded tables without editing TeX here.
tex_country <- function(x)gsub("&","\\\\&",x,fixed=TRUE)
tex <- c("% Replacement table source; saved as plain text. No manuscript file is edited.",
  "\\begin{table}[htbp]","\\centering",
  "\\caption{Ten countries with the largest absolute differences between the PfPR-ACM model and IHME estimates of malaria-attributable deaths before age five in 2024.}",
  "\\label{tab:country}","\\begin{tabular}{lrrr}","\\toprule",
  "Country & PfPR-ACM model & IHME & UN IGME \\\\","\\midrule",
  sprintf("%s & %s & %s & %s \\\\",tex_country(ifelse(x$iso3[1:10]=="COD","DR Congo",x$country[1:10])),fmt(x$model_deaths[1:10]),
    fmt(x$ihme_malaria_deaths[1:10]),fmt(x$who_cacode_deaths[1:10])),"\\midrule",
  sprintf("All 42 countries & %s & %s & %s \\\\",fmt(total[1]),fmt(total[2]),fmt(total[3])),
  "\\bottomrule","\\end{tabular}","\\par\\smallskip",
  "\\begin{minipage}{\\linewidth}\\footnotesize",
  "UN IGME denotes the direct under-five CA-CODE 2026 series. Totals include all 42 estimable countries, not only the ten shown. Death counts are rounded after estimation. The model uses the revised regional-adjustment MAP fits with $\\gamma=2$.",
  "\\end{minipage}","\\end{table}")
writeLines(tex,file.path(out,"country_comparison_2024.latex.txt"))
annual <- cbh_read_csv(paths[2]);states <- cbh_read_csv(paths[3])
stopifnot(nrow(annual)==21,nrow(states)==37,
  max(abs(unlist(annual[annual$year==2024,c("model_deaths","ihme_malaria_deaths","who_cacode_deaths")])-total))<1e-6)
writeLines(c("# Current primary tables","",
  "These tables report the age-band effects and national/annual/state burden specified in Sections 3–4 of the analysis plan.","",
  "- [Age-band results](age_band_results.md): observed deaths and percentages summing to 100.0%, EDF, 40%→20% and 20%→0% hazard ratios. [CSV](age_band_results.csv); [LaTeX source as text](age_band_results.latex.txt).",
  "- [2024 country comparison](country_comparison_2024.md), all 42 countries. [Ten-country manuscript presentation](country_comparison_2024_top10.md); [CSV](country_comparison_2024.csv); [LaTeX source as text](country_comparison_2024.latex.txt).",
  "- [Annual deaths and rates, 2004–2024](../annual_comparison/README.md). [CSV](../annual_comparison/annual_totals_2004_2024.csv).",
  "- [Nigerian state estimates, 2024](../nigeria_states/README.md). [CSV](../nigeria_states/state_totals_2024.csv).", "",
  "LaTeX source is provided in `.latex.txt` files so the author can paste it into the existing manuscript. No `.tex` files are created or edited."),file.path(out,"README.md"))
inputs <- c(paths,"R_cbh/reporting/09_burden_tables.R","R_cbh/reporting/labels.R")
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"burden_table_provenance.csv"))
message("Country tables and current table index written; no TeX files changed")

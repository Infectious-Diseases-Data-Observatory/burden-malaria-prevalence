# Shared aggregate-only table generator; no model fitting or manuscript writes.
cbh_primary_age_table <- function(root, ages) {
  paths <- file.path(root,c("fit_diagnostics.csv","pfpr_40_to_20_contrasts.csv",
    "pfpr_20_to_zero_contrasts.csv","pfpr_edf.csv","primary_sample.csv"))
  read_age <- function(path) {
    x <- cbh_read_csv(path)
    cbh_unique(x,"age_band",basename(path))
    stopifnot(nrow(x)==length(ages),setequal(x$age_band,ages),all(x$series=="map_full"))
    x[match(ages,x$age_band),]
  }
  d <- read_age(paths[1]);h <- read_age(paths[2]);z <- read_age(paths[3]);e <- read_age(paths[4])
  sample <- cbh_read_csv(paths[5])
  total <- sum(d$deaths)
  stopifnot(total==sample$deaths,sum(d$rows)==sample$records)
  # Largest-remainder rounding at 0.1 percentage point: the displayed age
  # shares sum to 100.0%, without altering raw counts or unrounded proportions.
  units <- 1000*d$deaths/total
  apportioned <- floor(units)
  remainder <- as.integer(round(1000-sum(apportioned)))
  if(remainder>0L) {
    ix <- order(units-apportioned,decreasing=TRUE)[seq_len(remainder)]
    apportioned[ix] <- apportioned[ix]+1L
  }
  pct <- apportioned/10
  stopifnot(sum(apportioned)==1000,all(abs(pct-100*d$deaths/total)<.1))
  # In the source zero-contrast table, lower_95/upper_95 are AF limits.
  # Compute HR limits from the log-HR covariance, then cross-check 1-AF.
  lo0 <- exp(z$log_hazard_ratio-1.96*z$standard_error)
  hi0 <- exp(z$log_hazard_ratio+1.96*z$standard_error)
  stopifnot(max(abs(exp(z$log_hazard_ratio)-z$hazard_ratio_20_to_zero))<1e-12,
    max(abs(lo0-(1-z$upper_95)))<1e-12,max(abs(hi0-(1-z$lower_95)))<1e-12)
  data <- data.frame(age_band=ages,records=d$rows,observed_deaths=d$deaths,
    observed_u5_death_share_pct=100*d$deaths/total,displayed_death_share_pct=pct,pfpr_edf=e$edf,
    hr_40_to_20=h$hazard_ratio_40_to_20,hr_40_to_20_lower_95=h$lower_95,hr_40_to_20_upper_95=h$upper_95,
    hr_20_to_zero=z$hazard_ratio_20_to_zero,hr_20_to_zero_lower_95=lo0,hr_20_to_zero_upper_95=hi0,
    zero_below_observed_support=z$zero_below_observed_support)
  fmt <- function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE)
  hr <- function(x,lo,hi,dash="–") sprintf("%.3f (%.3f%s%.3f)",x,lo,dash,hi)
  note <- paste0("Deaths are observed deaths in the primary analysis sample; percentages use all ",fmt(total),
    " under-five deaths as the denominator. Displayed percentages use largest-remainder rounding to one decimal place and sum to 100.0%. ",
    "Records are child–age-band observations, so a child can contribute multiple records. EDF: effective degrees of freedom of the PfPR spline. ",
    "Both contrast columns are adjusted mortality hazard ratios for reducing PfPR from the first value to the second. ",
    "Intervals are conditional on fitted smoothing parameters, exposure and the fixed HIV imputation; zero PfPR is below observed exposure support in every band.")
  md <- c("| Completed months | Records | Deaths (% of U5 deaths) | PfPR EDF | HR: 40% to 20% (95% interval) | HR: 20% to 0% (95% interval) |",
    "|---|---:|---:|---:|---:|---:|",
    vapply(seq_along(ages),function(i) sprintf("| %s | %s | %s (%.1f%%) | %.2f | %s | %s |",ages[i],fmt(d$rows[i]),
      fmt(d$deaths[i]),pct[i],e$edf[i],hr(data$hr_40_to_20[i],h$lower_95[i],h$upper_95[i]),
      hr(data$hr_20_to_zero[i],lo0[i],hi0[i])),""),
    sprintf("| **Total** | **%s** | **%s (100.0%%)** | — | — | — |",fmt(sum(d$rows)),fmt(total)),"",note)
  out <- file.path(root,"tables")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  cbh_atomic_csv(data,file.path(out,"age_band_results.csv"))
  writeLines(c("# Primary age-band results","",md),file.path(out,"age_band_results.md"))
  tex_ages <- ifelse(ages=="<1","$<1$",gsub("-","--",ages,fixed=TRUE))
  tex <- c(
    "% Standalone table fragment; uses booktabs (already present in the manuscript).",
    "% Generated from the saved primary MAP gamma=2 aggregate results.",
    "\\begin{table}[tbp]",
    "\\centering",
    "\\caption{Age-specific associations between malaria prevalence and all-cause mortality.}",
    "\\label{tab:primary-age-band-effects}",
    "\\small",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\renewcommand{\\arraystretch}{1.15}",
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lrrrcc@{}}",
    "\\toprule",
    "\\shortstack[l]{Age\\\\(months)} & \\shortstack{Child-band\\\\records} & \\shortstack{Deaths\\\\(\\% of U5 deaths)} & \\shortstack{PfPR\\\\EDF} & \\shortstack{HR (95\\% interval)\\\\40\\% $\\to$ 20\\%} & \\shortstack{HR (95\\% interval)\\\\20\\% $\\to$ 0\\%} \\\\",
    "\\midrule",
    vapply(seq_along(ages),function(i) sprintf("%s & %s & %s (%.1f\\%%) & %.2f & %s & %s \\\\",tex_ages[i],fmt(d$rows[i]),
      fmt(d$deaths[i]),pct[i],e$edf[i],hr(data$hr_40_to_20[i],h$lower_95[i],h$upper_95[i],"--"),
      hr(data$hr_20_to_zero[i],lo0[i],hi0[i],"--")),""),
    "\\midrule",
    sprintf("Total & %s & %s (100.0\\%%) & -- & -- & -- \\\\",fmt(sum(d$rows)),fmt(total)),
    "\\bottomrule",
    "\\end{tabular*}",
    "\\par\\vspace{0.5em}",
    "\\begin{minipage}{\\linewidth}",
    "\\footnotesize",
    paste0("Deaths are observed deaths in the primary analysis sample; percentages use all ",fmt(total)," under-five deaths as the denominator. "),
    "Displayed percentages use largest-remainder rounding to one decimal place and sum to 100.0\\%.",
    "Records are child--age-band observations; a child can contribute multiple records.",
    "EDF: effective degrees of freedom of the PfPR spline. HR: adjusted mortality hazard ratio for reducing PfPR from the first value to the second.",
    "Models were fitted separately by age band using MAP prevalence and $\\gamma=2$.",
    "Intervals condition on fitted smoothing parameters, exposure and the fixed HIV imputation; zero PfPR is below observed exposure support in every band.",
    "\\end{minipage}",
    "\\end{table}")
  writeLines(trimws(tex,which="right"),file.path(out,"age_band_results.tex"))
  inputs <- c(paths,"R_cbh/reporting/age_band_table.R")
  cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"provenance.csv"))
  invisible(md)
}

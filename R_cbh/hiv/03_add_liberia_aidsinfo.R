#!/usr/bin/env Rscript
# Child HIV incidence for Liberia from UNAIDS counts (added 24 September 2026).
# The UNAIDS/UNICEF workbook has no incidence series for Liberia, but AIDSinfo publishes
# its annual new HIV infections among children aged 0-14 (UNAIDS epidemiological
# estimates 2026; saved in data/unaids_aidsinfo/). The workbook's child rate ("new HIV
# infections per 1,000 uninfected population", age 0-14) is reproduced by
#   1000 x new infections aged 0-14 / (under-5 population - children aged 0-4 living with HIV),
# checked here against every country-year with a published rate. Liberia's rate is built
# the same way, with the IHME-implied under-5 person-years used by the burden step and
# the workbook's Liberia count of children aged 0-4 living with HIV.
# Writes two panels without touching the frozen ones:
#   child_incidence_country_year_lbr.csv           base panel + Liberia (primary)
#   child_incidence_country_year_extended_lbr.csv  extended panel, Liberia replaced
source("R_cbh/load_pipeline.R")
library(data.table)
hiv_dir <- "data/derived_cbh/hiv_incidence"; out <- "results/cbh/hiv_incidence/liberia_aidsinfo"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
paths <- c(base = file.path(hiv_dir, "child_incidence_country_year.csv"),
  extended = file.path(hiv_dir, "child_incidence_country_year_extended.csv"),
  aidsinfo = "data/unaids_aidsinfo/liberia_new_hiv_infections_children_0_14.csv",
  workbook = "data/HIV_Epidemiology_Children_Adolescents_2025.xlsx",
  ihme = list.files("data", pattern = "2026-09-09 10-58-22[.]csv$", full.names = TRUE))
stopifnot(all(file.exists(paths)), length(paths) == 5L)
years <- 2000:2024
num <- function(x) suppressWarnings(as.numeric(gsub(",", "", x)))

## ---- inputs ---------------------------------------------------------------------------------
wb <- as.data.table(suppressMessages(readxl::read_excel(paths[["workbook"]], sheet = "Data", skip = 1, guess_max = 100000)))
wb <- wb[Type == "Country" & Sex == "Both"][, `:=`(year = as.integer(Year), value = num(Value))]
ni <- wb[Age == "Age 0-14" & Indicator == "Estimated number of annual new HIV infections", .(iso3 = ISO3, year, ni = value)]
pl04 <- wb[Age == "Age 0-4" & Indicator == "Estimated number of people living with HIV", .(iso3 = ISO3, year, plhiv_0_4 = value)]
rate <- wb[Age == "Age 0-14" & Indicator == "Estimated incidence rate (new HIV infection per 1,000 uninfected population)",
  .(iso3 = ISO3, year, published_rate = value)]
ihme <- fread(paths[["ihme"]])[Sex == "Both" & Condition == "All causes" & Measure == "Deaths" & Age == "Under 5" & Year %in% years]
ihme[, iso3 := countrycode::countrycode(Location, "country.name", "iso3c", warn = FALSE)]
ihme <- dcast(ihme[!is.na(iso3)], iso3 + Year ~ Unit, value.var = "Value")
u5 <- ihme[, .(iso3, year = Year, under5_person_years = Number / `Rate (per 100,000)` * 1e5)]
stopifnot(all(is.finite(u5$under5_person_years)), all(u5$under5_person_years > 0))

## ---- check the definition against published rates ---------------------------------------------
v <- Reduce(function(a, b) merge(a, b, by = c("iso3", "year")), list(ni, rate, u5))
v <- merge(v, pl04, by = c("iso3", "year"), all.x = TRUE)
v <- v[is.finite(ni) & is.finite(published_rate) & published_rate > 0 & ni >= 200]
v[, derived_rate := 1000 * ni / (under5_person_years - fifelse(is.finite(plhiv_0_4), plhiv_0_4, 0))]
v[, ratio := derived_rate / published_rate]
q <- quantile(v$ratio, c(.05, .25, .5, .75, .95))
stopifnot(nrow(v) > 300, abs(q[["50%"]] - 1) < .05)
cbh_atomic_csv(as.data.frame(v[order(iso3, year)]), file.path(out, "definition_check_country_years.csv"))
by_country <- v[, .(country_years = .N, median_ratio = median(ratio), min_ratio = min(ratio), max_ratio = max(ratio)), by = iso3][order(iso3)]
cbh_atomic_csv(as.data.frame(by_country), file.path(out, "definition_check_by_country.csv"))

## ---- Liberia ------------------------------------------------------------------------------------
a <- fread(paths[["aidsinfo"]])[year %in% years]
stopifnot(all(a$iso3 == "LBR"), setequal(a$year, years), all(is.finite(a$estimate)))
l <- merge(merge(a[, .(year, ni = estimate, ni_lower = lower, ni_upper = upper)], u5[iso3 == "LBR"], by = "year"),
  pl04[iso3 == "LBR", .(year, plhiv_0_4)], by = "year", all.x = TRUE)
stopifnot(nrow(l) == length(years), all(is.finite(l$under5_person_years)), all(is.finite(l$plhiv_0_4)))
l[, uninfected := under5_person_years - plhiv_0_4]
l[, `:=`(rate = 1000 * ni / uninfected, lower = 1000 * ni_lower / uninfected, upper = 1000 * ni_upper / uninfected)]
cbh_atomic_csv(as.data.frame(l), file.path(out, "liberia_child_incidence.csv"))
make_rows <- function(template) {
  sig <- unique(template$imputation_signature); stopifnot(length(sig) == 1L)
  data.table(iso3 = "LBR", country_name = "Liberia", region = "West and Central Africa", year = l$year,
    child_rate = l$rate, child_source_value = sprintf("%.4f", l$rate), adolescent_source_value = NA_character_,
    hiv_incidence_per1000 = l$rate, lower_95 = l$lower, upper_95 = l$upper,
    hiv_incidence_status = "derived_aidsinfo_counts", imputation_signature = sig)
}
base <- fread(paths[["base"]]); ext <- fread(paths[["extended"]])
stopifnot(!"LBR" %in% base$iso3, "LBR" %in% ext$iso3)
panel <- rbind(base, make_rows(base)[, names(base), with = FALSE]); setorder(panel, iso3, year)
panel_ext <- rbind(ext[iso3 != "LBR"], make_rows(ext)[, names(ext), with = FALSE]); setorder(panel_ext, iso3, year)
cbh_unique(as.data.frame(panel), c("iso3", "year"), "Panel with Liberia"); cbh_unique(as.data.frame(panel_ext), c("iso3", "year"), "Extended panel with Liberia")
fwrite(panel, file.path(hiv_dir, "child_incidence_country_year_lbr.csv"))
fwrite(panel_ext, file.path(hiv_dir, "child_incidence_country_year_extended_lbr.csv"))
old <- ext[iso3 == "LBR" & year %in% years, .(year, extended_model_rate = hiv_incidence_per1000)]
cmp <- merge(l[, .(year, derived_rate = rate)], old, by = "year")
cbh_atomic_csv(as.data.frame(cmp), file.path(out, "liberia_derived_vs_extended_model.csv"))
nb <- rate[iso3 %in% c("SLE", "GIN", "CIV") & year %in% c(2005, 2010, 2015, 2020, 2024)]
fmt <- function(x) sprintf("%.2f", x)
writeLines(c("# Child HIV incidence for Liberia from UNAIDS counts", "",
  "The UNAIDS/UNICEF workbook (`data/HIV_Epidemiology_Children_Adolescents_2025.xlsx`) has no incidence series for Liberia, so Liberia was excluded from the complete-case primary. AIDSinfo publishes Liberia's annual new HIV infections among children aged 0–14 (UNAIDS epidemiological estimates 2026; retrieved 24 September 2026 and saved in `data/unaids_aidsinfo/`).", "",
  sprintf("**Definition check.** The workbook's child rate (new infections per 1,000 uninfected population, age 0–14) is reproduced by 1,000 × new infections aged 0–14 / (under-5 person-years − children aged 0–4 living with HIV): across %d country-years in %d countries with at least 200 new infections, derived/published has median %.3f (IQR %.3f–%.3f; 5–95%% %.3f–%.3f). The under-5 person-years are the IHME-implied denominators of the burden step. Dividing by the whole population aged 0–14 instead gives about 0.38 of the published rate, so the published child rate is effectively per uninfected child under 5.",
    nrow(v), uniqueN(v$iso3), q[["50%"]], q[["25%"]], q[["75%"]], q[["5%"]], q[["95%"]]), "",
  "**Liberia.** The same formula with Liberia's AIDSinfo counts, IHME-implied under-5 person-years and the workbook's children aged 0–4 living with HIV:", "",
  "| Year | New infections 0–14 | Rate per 1,000 (range from the count bounds) |", "|---|---:|---:|",
  sprintf("| %d | %s | %s (%s–%s) |", l$year, format(round(l$ni), big.mark = ","), fmt(l$rate), fmt(l$lower), fmt(l$upper))[l$year %in% c(2000, 2005, 2010, 2015, 2020, 2024)], "",
  paste0("Published rates in neighbouring countries for comparison: ", paste(sprintf("%s %d: %s", nb$iso3, nb$year, fmt(nb$published_rate)), collapse = "; "), "."), "",
  "**Panels.** `child_incidence_country_year_lbr.csv` is the frozen base panel plus Liberia 2000–2024 (status `derived_aidsinfo_counts`), used by the primary from v7; `child_incidence_country_year_extended_lbr.csv` is the extended panel with Liberia's latent-adolescent imputation replaced by these values, used by the imputed-covariate sensitivity. The frozen panels are unchanged. Caveats: the counts come from the 2026 estimates round while the other countries' rates come from the 2025 workbook; Liberia's rate is fixed (its range is shown here but not propagated).", "",
  "[Liberia series](liberia_child_incidence.csv) · [definition check by country](definition_check_by_country.csv) · [comparison with the extended model's imputation](liberia_derived_vs_extended_model.csv)"),
  file.path(out, "REPORT.md"))
prov <- c(unname(paths), "R_cbh/hiv/03_add_liberia_aidsinfo.R")
cbh_atomic_csv(data.frame(file = prov, md5 = vapply(prov, cbh_file_hash, "")), file.path(out, "provenance.csv"))
print(round(q, 3)); print(l[year %in% c(2000, 2005, 2010, 2015, 2020, 2024), .(year, ni = round(ni), rate = round(rate, 3))]); print(cmp[year %in% c(2005, 2015, 2024)])
message("Liberia child HIV incidence written")

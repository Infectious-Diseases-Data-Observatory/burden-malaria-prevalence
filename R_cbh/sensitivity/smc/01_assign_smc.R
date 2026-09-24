#!/usr/bin/env Rscript
# Admin-1 SMC switch-on years for the SMC before/after analysis.
# Aggregates the cluster file to admin-1 units (analysis regions; Nigerian states),
# pooled over rounds since a cluster's switch-on year is fixed per area, and builds
# the Nigerian child-to-state crosswalk from the survey files. Writes unit tables
# and the crosswalk to the ignored data folder and aggregate copies (no cluster
# identifiers) to results.
source("R_cbh/load_pipeline.R")
source("R_cbh/sensitivity/smc/settings.R")
library(data.table)
st <- cbh_smc_settings()
for (p in c(st$out, st$private)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
alias <- function(iso3, key) { a <- match(paste(iso3, key), paste(st$aliases$iso3, st$aliases$source_key)); key[!is.na(a)] <- st$aliases$regkey[a[!is.na(a)]]; key }
x <- fread(st$source)
stopifnot(all(x$country %in% names(st$countries)), !anyDuplicated(x[, .(country, survey, cluster_id)]),
          all(x$smc_grain %in% c("admin2", "admin1", "national", "unmatched", "none")),
          all(is.na(x$smc_first_year) == x$smc_grain %in% c("unmatched", "none")))
blank <- x[!nzchar(trimws(adm1_name)), .N, by = .(country, survey, smc_grain)]
x <- x[nzchar(trimws(adm1_name))]
x[, iso3 := st$countries[country]]
x[, regkey := alias(iso3, cbh_rkey(adm1_name))]
x[, unit_type := fifelse(iso3 == "NGA", "state", "analysis_region")]
# Countries that have introduced SMC: any recorded switch-on year.
start <- x[!is.na(smc_first_year), .(national_start = min(smc_first_year)), by = iso3]
countries <- sort(start$iso3)
x <- x[iso3 %in% countries & smc_grain != "unmatched"]
x[, first := fifelse(is.na(smc_first_year), Inf, as.numeric(smc_first_year))]
x[, confirmed := smc_grain %in% c("admin2", "admin1", "none")]
years <- 2000:2025
majority_year <- function(f) { if (!length(f)) return(NA_integer_)
  y <- years[vapply(years, function(t) mean(f <= t) >= st$majority, TRUE)]; if (length(y)) min(y) else NA_integer_ }
units <- x[, .(unit_type = unit_type[1], clusters = .N, clusters_admin2 = sum(smc_grain == "admin2"), clusters_admin1 = sum(smc_grain == "admin1"),
  clusters_national = sum(smc_grain == "national"), clusters_none = sum(smc_grain == "none"), rounds = uniqueN(survey),
  first_year_main = majority_year(first), first_year_confirmed = majority_year(first[confirmed]),
  earliest_cluster_year = suppressWarnings(as.integer(min(first[is.finite(first)])))), by = .(iso3, regkey)]
units <- merge(units, start, by = "iso3")
units[, status_basis := fifelse((clusters_admin2 + clusters_admin1 + clusters_none) / clusters >= st$majority, "confirmed", "national")]
units[, earliest_cluster_year := fifelse(is.finite(earliest_cluster_year), earliest_cluster_year, NA_integer_)]
setorder(units, iso3, regkey)
coverage <- x[, .(year = years, coverage = vapply(years, function(t) mean(first <= t), 0),
  coverage_confirmed = vapply(years, function(t) if (any(confirmed)) mean(first[confirmed] <= t) else NA_real_, 0)), by = .(iso3, regkey)]
stopifnot(!anyDuplicated(units[, .(iso3, regkey)]), all(nzchar(units$regkey)), all(coverage$coverage >= 0 & coverage$coverage <= 1))

## ---- Nigerian child-to-state crosswalk -----------------------------------------------------------
# The builder numbers clusters by first appearance (psu = match(v021, unique(v021))) and child i is
# row i of the recode it read, so the state is taken row by row and the psu codes are re-derived
# and checked against the shard before use.
registry <- fread("data/derived_dhs/survey_registry.csv")
cfg <- cbh_config(); dmeta <- readRDS(file.path(cfg$output_dir, "manifest.rds"))$manifest
mmeta <- readRDS("data/derived_mics/cbh_build/manifest.rds")$manifest
src <- st$nigeria_state_sources
cross <- rbindlist(lapply(seq_len(nrow(src)), function(i) {
  s <- src[i, ]
  if (s$file == "registry") {
    r <- readRDS(registry$local_recode[registry$svkey == s$survey]); state <- as.character(r[[s$state]])
    shard <- readRDS(file.path(cfg$output_dir, dmeta$file[dmeta$survey == s$survey]))$data
  } else {
    r <- readRDS(sprintf("data/derived_mics/recodes/%s.rds", s$survey))
    h <- haven::read_sav(s$file, col_select = c(s$cluster, s$state))
    h <- unique(data.table(cluster = as.integer(h[[s$cluster]]), state = as.character(haven::as_factor(h[[s$state]]))))
    stopifnot(!anyDuplicated(h$cluster)); state <- h$state[match(as.integer(r$v001), h$cluster)]
    shard <- readRDS(file.path("data/derived_mics/cbh_build", mmeta$file[mmeta$survey == s$survey]))$data
  }
  psu <- paste(s$survey, "p", match(r$v021, unique(r$v021)), sep = ":")
  d <- data.table(child_id = paste(s$survey, "c", seq_len(nrow(r)), sep = ":"), psu = psu, regkey = alias("NGA", cbh_rkey(state)))
  k <- match(as.character(shard$child_id), d$child_id)
  stopifnot(!anyNA(k), identical(as.character(shard$psu), d$psu[k]), !anyNA(d$regkey[k]), all(nzchar(d$regkey[k])))
  d <- unique(d[k, .(child_id, regkey)]); d[, survey := s$survey]
  message(s$survey, ": ", nrow(d), " children mapped to ", uniqueN(d$regkey), " states (psu codes reproduced)")
  d
}))
miss <- setdiff(unique(cross$regkey), units[iso3 == "NGA", regkey])
if (length(miss)) stop("Nigerian states without an SMC record: ", paste(miss, collapse = ", "))

fwrite(units, file.path(st$private, "unit_smc.csv")); fwrite(coverage, file.path(st$private, "unit_smc_coverage.csv"))
fwrite(cross, file.path(st$private, "nigeria_child_state.csv"))
cbh_atomic_csv(as.data.frame(units), file.path(st$out, "admin1_smc_years.csv"))
cbh_atomic_csv(as.data.frame(coverage), file.path(st$out, "admin1_smc_coverage.csv"))
cbh_atomic_csv(as.data.frame(cross[, .(children = .N), by = .(survey, state = regkey)][order(survey, state)]), file.path(st$out, "nigeria_state_children.csv"))
if (nrow(blank)) cbh_atomic_csv(as.data.frame(blank), file.path(st$out, "clusters_without_admin1_label.csv"))
print(units[, .(units = .N, confirmed = sum(status_basis == "confirmed"), first_main = paste(sort(unique(first_year_main)), collapse = "/")), by = iso3])
message("SMC countries: ", paste(countries, collapse = " "), "; ", nrow(units), " admin-1 units; Nigerian crosswalk ", nrow(cross), " children")

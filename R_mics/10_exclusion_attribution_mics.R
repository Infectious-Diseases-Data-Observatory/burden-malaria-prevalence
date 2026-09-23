#!/usr/bin/env Rscript
# Disjoint attribution of MICS complete-case exclusions to the FIRST missing covariate,
# in the same order as R_cbh/covariates/14_exclusion_attribution.R, so MICS and DHS
# tallies can be combined in the study flow. Aggregate counts only.
source("R_cbh/load_pipeline.R"); source("R_cbh/analysis/model.R")
source("R_cbh/primary/settings.R"); source("R_cbh/primary/specification.R")
library(data.table)
st <- cbh_primary_settings("regional_mics"); spec <- cbh_primary_regional_spec(st)
order_vars <- c("log_hiv_incidence","log_health_expenditure_pc","political_stability","log_gdp_pc",
  "wasting_pct","stunting_pct","facility_delivery_pct","electricity_pct","mean_wealth_quintile",
  "dtp3_pct","measles_pct","short_birth_interval_pct","improved_water_pct","improved_sanitation_pct",
  "urban_pct","mean_maternal_education_years","mean_maternal_age_first_birth")
stopifnot(setequal(order_vars, spec$covariates))
wide <- cbh_read_csv(st$mics_overlay); hiv <- cbh_read_csv(st$hiv_panel)
meta <- readRDS(file.path(st$mics_output_dir, "manifest.rds")); stopifnot(meta$complete)
m <- meta$manifest[meta$manifest$status %in% c("built","cached"), ]
tallies <- list()
for (i in seq_len(nrow(m))) {
  o <- readRDS(file.path(st$mics_output_dir, m$file[i])); stopifnot(identical(o$signature, m$signature[i]))
  d <- cbh_attach_incidence(o$data, hiv); d <- d[d$model_ready, , drop = FALSE]; if (!nrow(d)) next
  j <- match(paste(d$survey, d$regkey), paste(wide$survey, wide$regkey)); stopifnot(!anyNA(j))
  for (v in spec$regional) d[[v]] <- wide[[v]][j]
  valid <- vapply(order_vars, function(v) is.finite(d[[v]]), logical(nrow(d)))
  first <- apply(!valid, 1, function(r) if (any(r)) order_vars[which(r)[1]] else "retained")
  tallies[[i]] <- data.table(survey = m$survey[i], country = m$country[i], reason = first, death = d$death)[
    , .(records = .N, deaths = sum(death)), by = .(survey, country, reason)]
}
t <- rbindlist(tallies)
out <- "results/mics_inventory"
fwrite(t[reason != "retained"], file.path(out, "exclusion_attribution_mics_by_survey.csv"))
print(t[reason != "retained", .(records = sum(records), deaths = sum(deaths), surveys = uniqueN(survey)), by = reason][order(-records)])
cat("MICS model-ready records:", sum(t$records), "| retained:", sum(t[reason == "retained"]$records), "\n")

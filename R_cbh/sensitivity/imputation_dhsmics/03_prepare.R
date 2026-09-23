#!/usr/bin/env Rscript
# Analysis dataset for the DHS+MICS imputed-covariate sensitivity (v6). Mirrors the
# imputed branch of R_cbh/primary/00_prepare_regional.R, over both survey manifests:
# regional covariates from the jointly imputed overlay, the three national series
# (political stability, log GDP, log health expenditure) from the imputed national
# panel by country and band entry year, and child HIV incidence from the extended
# panel. Every MAP-eligible record must be retained. Covariates are re-standardised
# on this sample. No TeX writes; the prepared data stay under the ignored data folder.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/primary/specification.R")
source("R_cbh/sensitivity/imputation_dhsmics/settings.R")
library(data.table)
st <- cbh_imputed_dhsmics_settings(); cfg <- cbh_config()
spec <- cbh_primary_regional_spec(st); stopifnot(length(spec$covariates) == 17L, "urban_pct" %in% spec$covariates)
for (p in c(st$out, st$private)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
provenance <- rbind(cbh_read_csv(file.path(st$imputation_results, "regional_imputation_provenance.csv")),
                    cbh_read_csv(file.path(st$imputation_results, "national_imputation_provenance.csv")))
stopifnot(identical(unname(vapply(provenance$file, cbh_file_hash, "")), provenance$md5), file.exists(st$hiv_panel))
dmeta <- readRDS(file.path(cfg$output_dir, "manifest.rds")); mmeta <- readRDS(file.path(st$mics_output_dir, "manifest.rds"))
stopifnot(dmeta$complete, mmeta$complete)
source_files <- c(file.path(cfg$output_dir, "manifest.rds"), file.path(st$mics_output_dir, "manifest.rds"),
  file.path(st$imputation_results, c("regional_imputation_provenance.csv", "national_imputation_provenance.csv")),
  st$imputed_overlay, st$imputed_national, st$hiv_panel,
  "R_cbh/sensitivity/imputation_dhsmics/03_prepare.R", "R_cbh/sensitivity/imputation_dhsmics/settings.R",
  "R_cbh/primary/specification.R", "R_cbh/covariates/regional.R", "R_cbh/analysis/model.R")
signature <- cbh_hash(list(files = source_files, hashes = vapply(source_files, cbh_file_hash, ""), spec = spec))
signature_path <- file.path(st$private, "preparation_signature.rds")
if (file.exists(st$data) && file.exists(signature_path)) {
  cache <- readRDS(signature_path)
  if (identical(cache$signature, signature) && identical(cache$data_md5, cbh_file_hash(st$data))) {
    message("Verified current prepared dataset; retaining it."); quit(save = "no", status = 0)
  }
}
wide <- cbh_read_csv(st$imputed_overlay); cbh_unique(wide, c("survey", "regkey"), "Imputed overlay")
national <- cbh_read_csv(st$imputed_national); cbh_unique(national, c("iso3", "year"), "Imputed national panel")
hiv <- cbh_read_csv(st$hiv_panel)
flag_cols <- grep("_model_imputed$", names(wide), value = TRUE); stopifnot(length(flag_cols) == 13L)
required <- c("death", "age_band", "pfpr_pct", "calendar_year", "band_years", "survey", "country", "region", spec$covariates)
keep_cols <- unique(c(required, "child_id", "entry_year", "survey_year", "regkey"))
national_vars <- c(political_stability = "political_stability", log_gdp_pc = "log_gdp_pc",
                   log_health_expenditure_pc = "log_health_expenditure_pc")
stopifnot(all(national_vars %in% spec$covariates))
shards <- rbind(
  data.frame(programme = "DHS", dir = cfg$output_dir, dmeta$manifest[dmeta$manifest$status %in% c("built", "cached"), c("survey", "country", "file", "signature")]),
  data.frame(programme = "MICS", dir = st$mics_output_dir, mmeta$manifest[mmeta$manifest$status %in% c("built", "cached"), c("survey", "country", "file", "signature")]))
pieces <- selections <- counts_by_survey <- list()
max_diff <- setNames(numeric(length(national_vars)), national_vars)
for (i in seq_len(nrow(shards))) {
  s <- shards[i, ]
  object <- readRDS(file.path(s$dir, s$file)); stopifnot(identical(object$signature, s$signature))
  d <- cbh_attach_incidence(object$data, hiv)
  ready <- d$model_ready
  j <- match(paste(d$survey, d$regkey), paste(wide$survey, wide$regkey))
  if (any(ready & is.na(j))) stop("Imputed overlay has no row for ", paste(unique(paste(d$survey, d$regkey)[ready & is.na(j)]), collapse = ", "))
  for (v in spec$regional) d[[v]] <- wide[[v]][j]
  k <- match(paste(d$country, d$entry_year), paste(national$iso3, national$year))
  stopifnot(!anyNA(k[ready]))
  for (v in national_vars) {
    # Observed shard values must equal the panel's; only gaps change.
    both <- ready & is.finite(d[[v]])
    if (any(both)) max_diff[v] <- max(max_diff[v], abs(d[[v]][both] - national[[v]][k[both]]))
    d[[v]] <- national[[v]][k]
  }
  counts_by_survey[[i]] <- data.frame(programme = s$programme, survey = s$survey, country = s$country, records = sum(ready),
    regional_any_model_imputed = sum(rowSums(as.matrix(wide[j[ready], flag_cols, drop = FALSE])) > 0),
    political_stability_imputed = sum(national$political_stability_source[k[ready]] != "wgi_observed"),
    gdp_imputed = sum(national$log_gdp_pc_source[k[ready]] != "wdi_observed"),
    health_expenditure_imputed = sum(national$log_health_expenditure_pc_source[k[ready]] != "ghed_observed"),
    hiv_no_adolescent_series = sum(d$hiv_incidence_status[ready] == "imputed_no_adolescent_series"))
  complete <- ready & complete.cases(d[required])
  for (v in required) if (is.numeric(d[[v]])) complete <- complete & is.finite(d[[v]])
  selections[[i]] <- data.frame(programme = s$programme, survey = s$survey, country = s$country,
    model_ready_records = sum(ready), records = sum(complete), deaths = sum(d$death[complete]))
  pieces[[i]] <- d[complete, keep_cols, drop = FALSE]
  if (i %% 20L == 0L) message("Assembled ", i, "/", nrow(shards), " survey shards")
}
stopifnot(all(max_diff < 1e-8))
selection <- cbh_bind(selections); counts_by_survey <- cbh_bind(counts_by_survey)
# No model-ready record may be lost to covariate availability.
lost <- selection[selection$records != selection$model_ready_records, ]
if (nrow(lost)) stop("Records lost to covariates: ", paste(lost$survey, collapse = ", "))
d <- as.data.frame(rbindlist(pieces)); rm(pieces, object); gc(FALSE)
d$age_band <- factor(d$age_band, levels = cfg$age_bands$age_band)
for (v in c("survey", "country", "region")) d[[v]] <- factor(d[[v]])
stopifnot(!anyNA(d[required]), all(d$death %in% 0:1), all(d$band_years > 0),
  uniqueN(d, by = c("child_id", "age_band")) == nrow(d), all(table(d$age_band) > 0))
counts <- data.frame(records = nrow(d), children = uniqueN(d$child_id), deaths = sum(d$death),
  regions = nlevels(d$region), surveys = nlevels(d$survey), countries = nlevels(d$country))
dhs_rows <- !startsWith(as.character(d$survey), "MC_")
stopifnot(sum(dhs_rows) == st$expected_dhs_records, sum(d$death[dhs_rows]) == st$expected_dhs_deaths,
  uniqueN(d$region[dhs_rows]) == st$expected_dhs_regions,
  counts$records == st$expected_records, counts$deaths == st$expected_deaths, counts$regions == st$expected_regions,
  counts$surveys == st$expected_surveys, counts$countries == st$expected_countries)
scaling <- data.frame(variable = spec$covariates, mean = vapply(d[spec$covariates], mean, 0),
  sd = vapply(d[spec$covariates], sd, 0), row.names = NULL)
stopifnot(all(is.finite(scaling$sd) & scaling$sd > 0))
for (i in seq_len(nrow(scaling))) { v <- scaling$variable[i]; d[[paste0("z_", v)]] <- (d[[v]] - scaling$mean[i]) / scaling$sd[i] }
stopifnot(all(c("z_urban_pct", "z_mean_maternal_age_first_birth", "z_wasting_pct", "z_stunting_pct") %in% names(d)))
prepared <- list(signature = signature, data = d, scaling = scaling, specification = spec, counts = counts, selection = selection,
  source_hashes = data.frame(file = source_files, md5 = vapply(source_files, cbh_file_hash, "")))
cbh_atomic_rds(prepared, st$data)
cbh_atomic_rds(list(signature = signature, data_md5 = cbh_file_hash(st$data)), signature_path)
cbh_atomic_csv(counts, file.path(st$out, "prepared_sample.csv"))
cbh_atomic_csv(scaling, file.path(st$out, "covariate_scaling.csv"))
cbh_atomic_csv(selection, file.path(st$out, "prepared_selection_by_survey.csv"))
cbh_atomic_csv(prepared$source_hashes, file.path(st$out, "preparation_provenance.csv"))
cbh_atomic_csv(counts_by_survey, file.path(st$out, "imputation_record_counts_by_survey.csv"))
total <- aggregate(counts_by_survey[setdiff(names(counts_by_survey), c("programme", "survey", "country"))],
  list(programme = counts_by_survey$programme), sum)
total <- rbind(total, cbind(programme = "All", as.data.frame(lapply(total[-1], sum))))
cbh_atomic_csv(total, file.path(st$out, "imputation_record_counts.csv"))
print(counts); print(total)
message("DHS+MICS imputed-covariate dataset prepared: every MAP-eligible record retained")

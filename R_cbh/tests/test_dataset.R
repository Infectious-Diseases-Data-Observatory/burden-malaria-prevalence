# Synthetic records only; run from the project root with Rscript.
source("R_cbh/load_pipeline.R")
cfg <- cbh_config()
cmc <- function(year, month = 1) (year - 1900) * 12 + month
survey <- data.frame(svkey = "TESTFL", iso3 = "TST", year = 2020)
rule <- data.frame(svkey = "TESTFL", region_var = "v024", cmc_offset_months = 0L,
  calendar = "gregorian", group_donor = NA_character_, group_fine = NA_character_,
  group_coarse = NA_character_, strata_var = NA_character_)
boundary <- data.frame(svkey = "TESTFL", region = "North", regkey = "north")
empty_overrides <- data.frame(svkey = character(), region_var = character(),
  source_label = character(), regkey = character(), region_id = character())
fixture <- function(ages, death_month = rep(NA_real_, length(ages)), interview = cmc(2020, 6)) {
  n <- length(ages)
  data.frame(caseid = paste0("synthetic", seq_len(n)), bidx = 1L, bord = 2L,
    v005 = seq_len(n) * 1e6, v008 = interview, v011 = interview - ages - 300,
    v021 = "synthetic_psu", v022 = "synthetic_stratum", v024 = "North",
    v025 = "rural", v133 = 10L, v190 = "middle", b0 = "single birth",
    b3 = interview - ages, b4 = "female", b5 = ifelse(is.na(death_month), "yes", "no"),
    b6 = ifelse(is.na(death_month), NA, 200 + death_month), b7 = death_month, b11 = 30)
}
make <- function(br, rules = rule) {
  geo <- list(regkey = rep("north", nrow(br)), region = rep("TST:TESTFL:north", nrow(br)))
  cbh_make_bands(br, survey, rules, geo, cfg)
}
expect_error <- function(expr, pattern = NULL) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(e)))
}

# Boundary deaths enter exactly the intended band and never later bands.
deaths <- c(0, 1, 5, 6, 11, 12, 23, 24, 35, 36, 47, 48, 59, 60)
d <- make(fixture(rep(60, length(deaths)), deaths))$data
event <- d[d$death == 1, ]
expected <- findInterval(deaths, c(0, 1, 6, 12, 24, 36, 48, 60))
stopifnot(nrow(event) == 13, all(event$age_band_index == expected[as.integer(sub(".*:", "", event$child_id))]))
for (i in seq_along(deaths)) {
  z <- d[d$child_id == paste0("TESTFL:c:", i), ]
  stopifnot(all(z$age_band_index <= expected[i]))
}
# Complete potential bands use identical eligibility for dead and living children.
z <- make(fixture(c(24, 18, 18), c(NA, 14, NA)))$data
stopifnot(sum(z$child_id == "TESTFL:c:1") == 4, sum(z$child_id == "TESTFL:c:2") == 3,
          sum(z$child_id == "TESTFL:c:3") == 3, sum(z$death) == 0)
# Entry at 60 months before interview is included; birth before the window is allowed.
z <- make(fixture(c(60, 61, 72, 120)))$data
stopifnot(any(z$child_id == "TESTFL:c:1" & z$age_band_index == 1),
          !any(z$child_id == "TESTFL:c:2" & z$age_band_index == 1),
          any(z$child_id == "TESTFL:c:2" & z$age_band_index == 2),
          any(z$child_id == "TESTFL:c:3" & z$age_band_index == 4),
          !any(z$child_id == "TESTFL:c:4"), all(z$band_years == exp(z$log_band_years)))
empty <- make(fixture(0))$data
stopifnot(nrow(empty) == 0, all(c("death", "age_band", "entry_year") %in% names(empty)))

# Reported unit precision is preserved, including whole-year reports and B7 fallback.
a <- cbh_death_age(c("died on day of birth", "day: 30", "days: 31", "months: 6", "303", "999", NA),
                   c(0, 0, 1, 6, 37, 4, NA))
stopifnot(identical(a$lo, c(0, 0, 1, 6, 36, 4, NA_real_)), a$hi[5] == 48,
          a$source[6] == "b7_imputed_months", is.na(a$lo[7]))
labelled <- structure(c(100, 206, 303), labels = c("died on day of birth" = 100), class = c("haven_labelled", "vctrs_vctr", "double"))
stopifnot(identical(cbh_death_age(labelled, c(0, 6, 40))$lo, c(0, 6, 36)))
stopifnot(identical(cbh_death_age(factor(c("months: 6", "303")), c(6, 40))$lo, c(6, 36)))

# Invalid children are excluded with a counted reason; duplicates exclude both copies.
br <- fixture(rep(60, 7))
br$b5[1] <- "unknown"; br$b5[2] <- "no"
br$v005[3] <- 0; br$b3[4] <- br$v008[4] + 1
br$caseid[6:7] <- "same_synthetic_mother"
z <- make(br)
stopifnot(nrow(z$data) == 7, sum(z$child_checks$children[z$child_checks$status != "valid"]) == 6)
br <- fixture(60, 6); br$b7 <- NA; br$b6 <- 300
stopifnot(make(br)$child_checks$status == "ambiguous_death_band")

# Explicit calendar offset: equal dates after converting older Ethiopian CMCs.
br <- fixture(60)
gregorian <- make(br)$data
br[c("v008", "v011", "b3")] <- br[c("v008", "v011", "b3")] - 92
et_rule <- rule; et_rule$cmc_offset_months <- 92; et_rule$calendar <- "ethiopian_plus_92_approx"
converted <- make(br, et_rule)$data
stopifnot(identical(gregorian$band_entry_cmc, converted$band_entry_cmc),
          identical(gregorian$maternal_age_birth, converted$maternal_age_birth))
expect_error(make(br), "calendar")
actual_rules <- cbh_read_csv(cfg$survey_rules)
stopifnot(actual_rules$cmc_offset_months[actual_rules$svkey == "ET81FL"] == 92,
          actual_rules$cmc_offset_months[actual_rules$svkey == "ET8AFL"] == 0)

# Exact/canonical geographic matches; ambiguity and close spelling do not auto-match.
g <- cbh_match_regions(c("Province du Nord", "Nort", "North"), boundary)
stopifnot(g$regkey[1] == "north", is.na(g$regkey[2]), g$regkey[3] == "north")
ambiguous <- data.frame(region = c("North", "North Province"), regkey = c("a", "b"))
stopifnot(is.na(cbh_match_regions("Nord", ambiguous)$regkey))
stopifnot(is.na(cbh_match_regions("North excluding City", boundary)$regkey))
g <- cbh_match_regions(character(), boundary)
stopifnot(nrow(g) == 0)

# Isolated end-to-end fixture including cached input joins and output reading.
temp_root <- tempfile("cbh-tests-"); dir.create(temp_root)
dir.create(file.path(temp_root, "data"))
test_cfg <- cfg; test_cfg$root <- temp_root; test_cfg$output_dir <- file.path(temp_root, "data", "derived_cbh")
write_input <- function(d, name) {
  path <- file.path(temp_root, "data", paste0(name, ".csv"))
  write.csv(d, path, row.names = FALSE, na = "")
  path
}
raw_path <- file.path(temp_root, "data", "test.rds")
saveRDS(fixture(c(60, 24)), raw_path)
registry <- cbind(survey, local_recode = raw_path)
test_cfg$registry <- write_input(registry, "registry")
test_cfg$survey_rules <- write_input(rule, "rules")
test_cfg$region_overrides <- write_input(empty_overrides, "overrides")
test_cfg$boundary_regions <- write_input(boundary, "boundaries")
map <- data.frame(svkey = "TESTFL", regkey = "north", year = 2015:2020, pfpr2_10 = c(0, 10, 20, 30, 40, 50))
test_cfg$annual_map <- write_input(map, "map")
for (name in names(test_cfg$annual_panels)) {
  spec <- test_cfg$annual_panels[[name]]
  p <- data.frame(iso3 = "TST", year = 2015:2020)
  for (v in unname(spec$columns)) p[[v]] <- 50
  for (v in unname(spec$source_status)) p[[v]] <- c(rep("assumed_pre_series_zero", 3), rep("wuenic_estimate", 3))
  if (name == "hiv") p <- p[p$year != 2015, ]
  test_cfg$annual_panels[[name]]$path <- write_input(p, name)
}
test_cfg$statcompiler <- file.path(temp_root, "data", "absent.csv")
geo <- cbh_geography(fixture(60), survey, rule, boundary, empty_overrides, registry)
stopifnot(geo$regkey == "north")
ov <- data.frame(svkey = "TESTFL", region_var = "v024", source_label = "North", regkey = "north", region_id = "stable_north")
stopifnot(cbh_geography(fixture(60), survey, rule, boundary, ov, registry)$region == "TST:stable_north")
# A declared donor can group fine recode regions into a coarser MAP geography.
donor_path <- file.path(temp_root, "data", "donor.rds")
donor_data <- fixture(c(60, 60)); donor_data$v024 <- c("North A", "North B"); donor_data$szone <- "North"
saveRDS(donor_data, donor_path)
donor_registry <- rbind(registry, data.frame(svkey = "DONORFL", iso3 = "TST", year = 2020, local_recode = donor_path))
donor_rule <- rule; donor_rule$group_donor <- "DONORFL"; donor_rule$group_fine <- "v024"; donor_rule$group_coarse <- "szone"
donor_target <- fixture(60); donor_target$v024 <- "North A"
stopifnot(cbh_geography(donor_target, survey, donor_rule, boundary, empty_overrides, donor_registry)$regkey == "north")
inputs <- cbh_external_inputs(test_cfg)
z <- cbh_attach_external(gregorian, survey, boundary, inputs, test_cfg)
stopifnot(z$pfpr_pct[1] == 0, z$model_ready[1], is.na(z$hiv_prev_pct[1]),
          is.na(z$hib3_pct[1]), all(z$exposure_year == z$entry_year))
inputs$map <- inputs$map[inputs$map$year != 2015, ]
stopifnot(is.na(cbh_attach_external(gregorian, survey, boundary, inputs, test_cfg)$pfpr_pct[1]))
inputs$map <- map
# Duplicate survey-region measurements are flagged, never arbitrarily selected.
inputs$statcompiler <- data.frame(svkey = "TESTFL", variable = "imp_water", level = "subnational",
                                  CharacteristicLabel = c("North", "North"), Value = c(40, 50))
dup <- cbh_attach_external(gregorian, survey, boundary, inputs, test_cfg)
stopifnot(all(is.na(dup$imp_water_survey_region_pct)), all(dup$imp_water_survey_region_pct_status == "ambiguous_source"))
stopifnot(nrow(cbh_attach_external(empty, survey, boundary, inputs, test_cfg)) == 0)
cbh_validate_bands(z, test_cfg)
map_bad <- rbind(map, map[1, ])
bad_cfg <- test_cfg; bad_cfg$annual_map <- write_input(map_bad, "duplicate_map")
expect_error(cbh_external_inputs(bad_cfg), "duplicate keys")
bad_cfg <- test_cfg; bad_cfg$output_dir <- file.path(temp_root, "public")
expect_error(cbh_private_output(bad_cfg), "ignored data")
m1 <- cbh_build(test_cfg); stopifnot(m1$status == "built", m1$rows == 11)
m2 <- cbh_build(test_cfg); stopifnot(m2$status == "cached", identical(m1$signature, m2$signature))
d <- cbh_load_analysis(test_cfg$output_dir, required_covariates = "hiv_prev_pct")
stopifnot(all(is.finite(d$hiv_prev_pct)), abs(mean(d$analysis_weight) - 1) < 1e-12,
          length(levels(d$age_band)) == 7, !is.ordered(d$age_band),
          attr(d, "selection_report")$missing_requested_covariates > 0)
map$pfpr2_10[1] <- 5; write.csv(map, test_cfg$annual_map, row.names = FALSE)
m3 <- cbh_build(test_cfg); stopifnot(m3$status == "built", m3$signature != m1$signature)
test_cfg$entry_lookback_months <- 59L
m4 <- cbh_build(test_cfg); stopifnot(m4$status == "built", m4$rows < m3$rows)
expect_error(cbh_load_analysis(test_cfg$output_dir, surveys = "UNKNOWN"), "absent")
meta <- readRDS(file.path(test_cfg$output_dir, "manifest.rds")); meta$complete <- FALSE
saveRDS(meta, file.path(test_cfg$output_dir, "manifest.rds"))
expect_error(cbh_load_analysis(test_cfg$output_dir), "did not complete")
unlink(temp_root, recursive = TRUE)
cat("All child-band dataset checks passed (synthetic data only).\n")

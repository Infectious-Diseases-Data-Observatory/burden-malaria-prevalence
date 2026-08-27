# =============================================================================
# 03_build_analysis_dataset.R
# Build one aggregate survey-region analysis dataset containing:
#   * MAP PfPR2-10
#   * DHS all-under-5, neonatal and post-neonatal mortality + exposure
#   * DHS survey-region covariates
#   * nearest-year national World Bank covariates
#   * transparent missingness, eligibility and single-imputation fields
#
# Standard run (from authorised DHS recodes and outputs of scripts 01-02):
#   Rscript R_dhs/03_build_analysis_dataset.R
#
# Migration/verification run starting from the existing aggregate panel:
#   Rscript R_dhs/03_build_analysis_dataset.R --from-legacy-aggregate
#
# When the survey registry is available, migration mode reads the locally
# authorised recodes only to derive the three new vaccine aggregates. No
# respondent-level records are written to the analysis outputs.
# =============================================================================

source("R_dhs/00_config.R")

args <- commandArgs(trailingOnly = TRUE)
legacy_mode <- "--from-legacy-aggregate" %in% args

COVARIATE_CATALOG <- data.frame(
  variable = c(
    "pct_urban", "dtp3_reg", "measles", "pentavalent3_reg",
    "pcv3_reg", "rotavirus_complete_reg", "facility", "educ_yrs", "wealth_q",
    "excl_bf", "stunting", "underweight", "wasting", "birth_int",
    "mage1", "imp_water", "imp_sanit", "elec_dhs",
    "hib3_wuenic", "pcv3_wuenic", "rotac_wuenic", "log_hiv_prev",
    "log_gdp", "log_hexp_pc", "polstab"
  ),
  level = c(
    rep("survey-region", 18),
    rep("national-exact-year", 4),
    rep("national-nearest-year", 3)
  ),
  type = c(
    "proportion", "proportion", "proportion", "proportion", "proportion",
    "proportion", "proportion", "continuous", "continuous", "proportion",
    "proportion", "proportion", "proportion", "proportion", "continuous",
    "proportion", "proportion", "proportion",
    "proportion", "proportion", "proportion", "continuous",
    "continuous", "continuous", "continuous"
  )
)

yes_received <- function(x) {
  grepl(
    "yes|received|reported|card|marked|mother",
    tolower(as.character(x))
  )
}

region_covariates <- function(br, region_var) {
  if (!region_var %in% names(br)) return(NULL)
  weight <- suppressWarnings(as.numeric(br$v005)) / 1e6
  region <- rkey(as.character(br[[region_var]]))
  aggregate_value <- function(value, eligible) {
    weighted_mean_by_region(value, weight, region, eligible)
  }

  age_months <- suppressWarnings(
    as.numeric(br$v008) - as.numeric(br$b3)
  )
  alive <- if ("b5" %in% names(br)) {
    tolower(as.character(br$b5)) == "yes"
  } else {
    rep(TRUE, nrow(br))
  }

  urban <- if ("v025" %in% names(br)) {
    100 * aggregate_value(
      as.numeric(tolower(as.character(br$v025)) == "urban"),
      !is.na(br$v025)
    )
  } else numeric(0)

  immunisation <- function(variable) {
    if (!variable %in% names(br)) return(numeric(0))
    labels <- tolower(as.character(br[[variable]]))
    eligible <- alive & is.finite(age_months) &
      age_months >= 12 & age_months <= 23 &
      !is.na(labels) & labels != "missing"
    100 * aggregate_value(as.numeric(yes_received(labels)), eligible)
  }
  dtp3_reg <- immunisation("h7")
  measles <- immunisation("h9")

  immunisation_series <- function(variables, required_doses) {
    available <- intersect(variables, names(br))
    if (!length(available)) return(numeric(0))
    labels <- lapply(available, function(variable) {
      tolower(as.character(br[[variable]]))
    })
    observed <- lapply(labels, function(x) {
      !is.na(x) & x != "missing" & nzchar(x)
    })
    received <- lapply(labels, yes_received)
    observed_any <- Reduce(`|`, observed)
    dose_count <- Reduce(`+`, lapply(received, as.integer))
    eligible <- alive & is.finite(age_months) &
      age_months >= 12 & age_months <= 23 & observed_any
    100 * aggregate_value(
      as.numeric(dose_count >= required_doses),
      eligible
    )
  }
  pentavalent3_reg <- immunisation_series(c("h51", "h52", "h53"), 3)
  pcv3_reg <- immunisation_series(c("h54", "h55", "h56"), 3)

  # DHS schedules contain either two or three rotavirus doses. Treat H59 as
  # evidence of a three-dose schedule only when it has observed values for an
  # eligible child; otherwise use completion of H57-H58.
  rotavirus_doses <- 2L
  if ("h59" %in% names(br)) {
    h59_labels <- tolower(as.character(br$h59))
    h59_observed <- alive & is.finite(age_months) &
      age_months >= 12 & age_months <= 23 &
      !is.na(h59_labels) & h59_labels != "missing" & nzchar(h59_labels)
    if (any(h59_observed)) rotavirus_doses <- 3L
  }
  rotavirus_complete_reg <- immunisation_series(
    c("h57", "h58", "h59"),
    rotavirus_doses
  )

  facility <- numeric(0)
  if ("m15" %in% names(br)) {
    labels <- tolower(as.character(br$m15))
    usable <- !is.na(labels) & labels != "missing" & nzchar(labels)
    facility_pattern <- paste(
      "hospital|clinic|health|dispensary|maternit|cms|pmi|hut|post|doctor",
      "nursing|centre|center|sector|infirmary|polyclin|hopital|hôpital",
      "clinique|sante|santé|cabinet",
      sep = "|"
    )
    home_pattern <- "home|house|domicile|maison|parent"
    in_facility <- !grepl(home_pattern, labels) &
      grepl(facility_pattern, labels)
    facility <- 100 * aggregate_value(as.numeric(in_facility), usable)
  }

  mother_id <- if ("caseid" %in% names(br)) {
    as.character(br$caseid)
  } else {
    paste(br$v001, br$v002, br$v003)
  }
  first_mother <- !duplicated(mother_id)
  household_id <- paste(br$v001, br$v002)
  first_household <- !duplicated(household_id)

  education <- suppressWarnings(as.numeric(as.character(br$v133)))
  education[education < 0 | education > 25] <- NA_real_
  educ_yrs <- aggregate_value(
    education,
    first_mother & is.finite(education)
  )

  wealth <- match(
    tolower(as.character(br$v190)),
    c("poorest", "poorer", "middle", "richer", "richest")
  )
  wealth_q <- aggregate_value(
    as.numeric(wealth),
    first_mother & is.finite(wealth)
  )

  excl_bf <- numeric(0)
  if ("v404" %in% names(br)) {
    under_six_months <- alive & is.finite(age_months) & age_months < 6
    breastfeeding <- grepl(
      "yes|breast", tolower(as.character(br$v404))
    )
    other_food_vars <- intersect(
      c(
        "v409", "v410", "v411", "v411a", "v412", "v412a", "v413",
        "v414a", "v414b", "v414c", "v414e", "v414f", "v414g",
        "v414h", "v414i", "v414j", "v414k", "v414l", "v414m",
        "v414n", "v414o", "v414p", "v414v"
      ),
      names(br)
    )
    other_food <- if (length(other_food_vars)) {
      Reduce(
        `|`,
        lapply(other_food_vars, function(v) yes_received(br[[v]]))
      )
    } else {
      rep(FALSE, nrow(br))
    }
    excl_bf <- 100 * aggregate_value(
      as.numeric(breastfeeding & !other_food),
      under_six_months & !is.na(br$v404)
    )
  }

  anthropometry <- function(variable) {
    if (!variable %in% names(br)) return(numeric(0))
    z <- suppressWarnings(as.numeric(br[[variable]]))
    z[z >= 9990 | z <= -600] <- NA_real_
    100 * aggregate_value(as.numeric(z < -200), alive & is.finite(z))
  }
  stunting <- anthropometry("hw5")
  underweight <- anthropometry("hw8")
  wasting <- anthropometry("hw11")

  birth_int <- if ("b11" %in% names(br)) {
    interval <- suppressWarnings(as.numeric(br$b11))
    100 * aggregate_value(
      as.numeric(interval < 24),
      is.finite(interval)
    )
  } else numeric(0)

  maternal_age <- suppressWarnings(as.numeric(br$v212))
  maternal_age[maternal_age < 8 | maternal_age > 45] <- NA_real_
  mage1 <- if ("v212" %in% names(br)) {
    aggregate_value(
      maternal_age,
      first_mother & is.finite(maternal_age)
    )
  } else numeric(0)

  household_fraction <- function(variable, include_pattern, exclude_pattern = NULL) {
    if (!variable %in% names(br)) return(numeric(0))
    labels <- tolower(as.character(br[[variable]]))
    usable <- first_household & !is.na(labels) &
      labels != "missing" & nzchar(labels)
    included <- grepl(include_pattern, labels)
    if (!is.null(exclude_pattern)) {
      included <- included & !grepl(exclude_pattern, labels)
    }
    100 * aggregate_value(as.numeric(included), usable)
  }
  imp_water <- household_fraction(
    "v113",
    "pipe|tap|standpipe|borehole|tube ?well|protected|rain|bottled|sachet",
    "unprotected"
  )
  imp_sanit <- household_fraction(
    "v116",
    "flush|septic|sewer|ventilated|vip|slab|composting",
    "without slab|open pit|no facil|bush|field|hanging|bucket|somewhere"
  )
  elec_dhs <- household_fraction("v119", "yes")

  values <- list(
    pct_urban = urban,
    dtp3_reg = dtp3_reg,
    measles = measles,
    pentavalent3_reg = pentavalent3_reg,
    pcv3_reg = pcv3_reg,
    rotavirus_complete_reg = rotavirus_complete_reg,
    facility = facility,
    educ_yrs = educ_yrs,
    wealth_q = wealth_q,
    excl_bf = excl_bf,
    stunting = stunting,
    underweight = underweight,
    wasting = wasting,
    birth_int = birth_int,
    mage1 = mage1,
    imp_water = imp_water,
    imp_sanit = imp_sanit,
    elec_dhs = elec_dhs
  )
  region_keys <- Reduce(union, lapply(values, names))
  if (!length(region_keys)) return(NULL)

  out <- data.frame(regkey = region_keys)
  for (variable in names(values)) {
    out[[variable]] <- pick_named(values[[variable]], region_keys)
  }
  out
}

build_from_raw <- function() {
  if (!file.exists(SURVEY_REGISTRY_CSV)) {
    stop("Survey registry missing. Run script 01.")
  }
  if (!file.exists(MAP_REGION_CSV)) {
    stop("MAP survey-region estimates missing. Run script 02.")
  }
  registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
  map <- read.csv(MAP_REGION_CSV, stringsAsFactors = FALSE)

  rows <- list()
  for (i in seq_len(nrow(registry))) {
    survey <- registry[i, , drop = FALSE]
    survey_map <- map[map$svkey == survey$svkey, , drop = FALSE]
    recode <- survey$local_recode
    if (!is.character(recode) || is.na(recode) || !file.exists(recode) ||
        !nrow(survey_map)) next

    br <- tryCatch(readRDS(recode), error = function(e) NULL)
    if (is.null(br)) next
    region_var <- best_region_var(br, survey_map$regkey)
    mortality <- mortality_by_region(br, region_var)
    covariates <- region_covariates(br, region_var)
    if (is.null(mortality) || is.null(covariates)) next

    region <- merge(mortality, covariates, by = "regkey", all = TRUE)
    region <- merge(
      region,
      unique(survey_map[, c(
        "svkey", "iso3", "year", "regkey", "region",
        "pfpr2_10", "population_weight"
      )]),
      by = "regkey"
    )
    if (!nrow(region)) next
    rows[[survey$svkey]] <- region
    message(
      "Aggregated ", survey$svkey, ": ", nrow(region),
      " matched survey-regions."
    )
  }
  if (!length(rows)) stop("No survey-regions could be assembled.")
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

attach_regional_vaccine_covariates <- function(data) {
  vaccine_variables <- c(
    "pentavalent3_reg", "pcv3_reg", "rotavirus_complete_reg"
  )
  for (variable in vaccine_variables) data[[variable]] <- NA_real_

  if (!file.exists(SURVEY_REGISTRY_CSV)) {
    warning(
      "Survey registry unavailable; new regional vaccine covariates remain ",
      "missing. Run script 01 to inventory local DHS recodes."
    )
    return(data)
  }

  registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
  rows <- list()
  for (svkey in unique(data$svkey)) {
    survey <- registry[
      registry$svkey == svkey & registry$recode_available,
      ,
      drop = FALSE
    ]
    if (!nrow(survey)) next
    recode <- survey$local_recode[1]
    if (!is.character(recode) || is.na(recode) || !file.exists(recode)) next

    br <- tryCatch(readRDS(recode), error = function(e) NULL)
    if (is.null(br)) next
    target_keys <- unique(data$regkey[data$svkey == svkey])
    region_var <- best_region_var(br, target_keys)
    covariates <- region_covariates(br, region_var)
    if (is.null(covariates)) next

    covariates$svkey <- svkey
    rows[[svkey]] <- covariates[, c(
      "svkey", "regkey", vaccine_variables
    )]
  }
  if (!length(rows)) return(data)

  vaccine_panel <- do.call(rbind, rows)
  key <- paste(data$svkey, data$regkey)
  vaccine_key <- paste(vaccine_panel$svkey, vaccine_panel$regkey)
  matched <- match(key, vaccine_key)
  for (variable in vaccine_variables) {
    data[[variable]] <- vaccine_panel[[variable]][matched]
  }
  message(
    "Attached direct DHS vaccine covariates from ",
    length(unique(vaccine_panel$svkey)), " surveys."
  )
  data
}

build_from_legacy_aggregate <- function() {
  path <- file.path(REPO_ROOT, "results", "component2_region_data_full.csv")
  if (!file.exists(path)) stop("Legacy aggregate panel not found: ", path)
  out <- read.csv(path, stringsAsFactors = FALSE)
  if ("m1mo5y" %in% names(out)) {
    out$postneonatal_mortality <- out$m1mo5y
  }
  out$nnmr <- out$u5mr - out$postneonatal_mortality

  # Restore original missing values using the legacy imputation flags so the
  # new <=5% eligibility rule is evaluated on pre-imputation data.
  for (variable in COVARIATE_CATALOG$variable) {
    flag <- paste0(variable, "_imp")
    if (variable %in% names(out) && flag %in% names(out)) {
      out[[variable]][as.logical(out[[flag]])] <- NA_real_
    }
  }
  attach_regional_vaccine_covariates(out)
}

attach_national_covariates <- function(data) {
  panels <- list(
    dtp3 = read_or_fetch_wb("wb_dtp3.csv", "SH.IMM.IDPT", "dtp3"),
    gdp_pc = read_or_fetch_wb(
      "wb_gdp_pc.csv", "NY.GDP.PCAP.CD", "gdp_pc"
    ),
    hexp_gdp = read_or_fetch_wb(
      "wb_hexp_gdp.csv", "SH.XPD.CHEX.GD.ZS", "hexp_gdp"
    ),
    hexp_pc = read_or_fetch_wb(
      "wb_hexp_pc.csv", "SH.XPD.CHEX.PC.CD", "hexp_pc"
    ),
    # The Bank archived PV.EST; the Worldwide Governance Indicators are now
    # served under GOV_WGI_ codes.
    polstab = read_or_fetch_wb("wb_polstab.csv", "GOV_WGI_PV_EST", "polstab"),
    elec = read_or_fetch_wb("wb_elec.csv", "EG.ELC.ACCS.ZS", "elec")
  )
  for (name in names(panels)) {
    data[[name]] <- mapply(
      function(iso3, year) {
        nearest_panel_value(panels[[name]], iso3, year, name)
      },
      data$iso3,
      data$year
    )
  }
  data$log_gdp <- log(data$gdp_pc)
  data$log_hexp_pc <- log(data$hexp_pc)
  data
}

attach_unicef_immunisation <- function(data) {
  vaccine_variables <- c(
    "hib3_wuenic", "pcv3_wuenic", "rotac_wuenic"
  )
  status_variables <- paste0(vaccine_variables, "_status")
  if (!file.exists(UNICEF_IMMUNISATION_CSV)) {
    for (variable in vaccine_variables) data[[variable]] <- NA_real_
    for (variable in status_variables) data[[variable]] <- "panel_unavailable"
    warning(
      "Compact UNICEF immunisation panel unavailable; run ",
      "R_dhs/02b_extract_unicef_immunisation.R."
    )
    return(data)
  }

  panel <- read.csv(
    UNICEF_IMMUNISATION_CSV,
    stringsAsFactors = FALSE
  )
  required <- c("iso3", "year", vaccine_variables, status_variables)
  missing_columns <- setdiff(required, names(panel))
  if (length(missing_columns)) {
    stop(
      "Compact UNICEF panel is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  panel_key <- paste(panel$iso3, panel$year, sep = "|")
  if (anyDuplicated(panel_key)) {
    stop("Compact UNICEF panel has duplicate country-year rows.")
  }

  data_key <- paste(data$iso3, data$year, sep = "|")
  matched <- match(data_key, panel_key)
  for (variable in c(vaccine_variables, status_variables)) {
    data[[variable]] <- panel[[variable]][matched]
  }
  message(
    "Matched UNICEF immunisation estimates to ",
    sum(is.finite(matched)), " of ", nrow(data), " survey-region-years."
  )
  data
}

# Derived child (0-14) HIV prevalence from UNAIDS 2025 estimates (via UNICEF).
# The workbook reports no prevalence rate, so it is derived as the estimated
# number of children 0-14 living with HIV divided by the World Bank 0-14
# population (SP.POP.0014.TO), joined by ISO3 and exact survey year. The
# prevalence spans ~3 orders of magnitude and is entered on the log scale, so
# the analysis covariate is log_hiv_prev. Countries whose UNAIDS child series
# is absent (only 15-19 published) stay missing and fall under the standard
# <=5% imputation rule with an explicit status flag.
attach_hiv_prevalence <- function(data) {
  data$hiv_prev <- NA_real_
  data$log_hiv_prev <- NA_real_
  data$hiv_prev_status <- "no_under15_series"
  if (!file.exists(HIV_XLSX)) {
    warning(
      "HIV epidemiology workbook not found: ", HIV_XLSX,
      "; child HIV prevalence remains missing."
    )
    return(data)
  }
  required_packages("readxl")

  raw <- as.data.frame(suppressMessages(readxl::read_excel(
    HIV_XLSX, sheet = "Data", skip = 1, guess_max = 100000
  )))
  parse_estimate <- function(x) {
    x <- gsub(",", "", trimws(as.character(x)))
    below <- grepl("^<", x)
    x <- gsub("[<>]", "", x)
    value <- suppressWarnings(as.numeric(x))
    value[below] <- value[below] / 2   # e.g. "<500" -> 250 midpoint
    value
  }
  plhiv <- raw[
    raw$Sex == "Both" & raw$Age == "Age 0-14" &
      raw$Indicator == "Estimated number of people living with HIV",
    c("ISO3", "Year", "Value"),
    drop = FALSE
  ]
  plhiv$iso3 <- toupper(plhiv$ISO3)
  plhiv$year <- suppressWarnings(as.integer(plhiv$Year))
  plhiv$plhiv <- parse_estimate(plhiv$Value)
  plhiv <- plhiv[
    nchar(plhiv$iso3) == 3L & is.finite(plhiv$year) &
      is.finite(plhiv$plhiv) & plhiv$plhiv > 0,
    ,
    drop = FALSE
  ]
  plhiv <- plhiv[!duplicated(paste(plhiv$iso3, plhiv$year, sep = "|")), , drop = FALSE]

  population <- read_or_fetch_wb("wb_pop_0_14.csv", "SP.POP.0014.TO", "pop_0_14")
  matched_pop <- match(
    paste(plhiv$iso3, plhiv$year, sep = "|"),
    paste(population$iso3, population$year, sep = "|")
  )
  plhiv$pop_0_14 <- population$pop_0_14[matched_pop]
  plhiv$hiv_prev <- 100 * plhiv$plhiv / plhiv$pop_0_14
  panel <- plhiv[
    is.finite(plhiv$hiv_prev) & plhiv$hiv_prev > 0,
    c("iso3", "year", "plhiv", "pop_0_14", "hiv_prev"),
    drop = FALSE
  ]
  write.csv(panel, HIV_PANEL_CSV, row.names = FALSE)

  matched <- match(
    paste(data$iso3, data$year, sep = "|"),
    paste(panel$iso3, panel$year, sep = "|")
  )
  data$hiv_prev <- panel$hiv_prev[matched]
  data$log_hiv_prev <- log(data$hiv_prev)
  data$hiv_prev_status[is.finite(matched)] <- "unaids_2025_estimate"

  covered_countries <- sort(unique(panel$iso3))
  summary_panel <- data.frame(
    variable = "log_hiv_prev",
    source = "UNAIDS 2025 estimates (0-14 PLHIV) / World Bank SP.POP.0014.TO",
    matched_region_years = sum(is.finite(matched)),
    total_region_years = nrow(data),
    countries_with_series = length(covered_countries),
    countries_without_series = paste(
      sort(setdiff(unique(data$iso3), covered_countries)), collapse = " "
    ),
    min_prevalence_pct = min(panel$hiv_prev),
    max_prevalence_pct = max(panel$hiv_prev)
  )
  write.csv(summary_panel, HIV_SUMMARY_CSV, row.names = FALSE)

  message(
    "Matched child HIV prevalence to ", sum(is.finite(matched)), " of ",
    nrow(data), " survey-region-years (no under-15 UNAIDS series: ",
    paste(sort(setdiff(unique(data$iso3), covered_countries)), collapse = " "),
    ")."
  )
  data
}

if (legacy_mode) {
  message("Building from the existing aggregate panel for migration validation.")
  analysis <- build_from_legacy_aggregate()
} else {
  required_packages(c("DHS.rates", "countrycode"))
  analysis <- build_from_raw()
  analysis <- attach_national_covariates(analysis)
}
analysis <- attach_unicef_immunisation(analysis)
analysis <- attach_hiv_prevalence(analysis)

analysis$pfpr10 <- analysis$pfpr2_10 / 10
analysis$nnmr <- if ("nnmr" %in% names(analysis)) {
  analysis$nnmr
} else {
  analysis$u5mr - analysis$postneonatal_mortality
}

missingness <- apply_missingness_rule(
  analysis,
  COVARIATE_CATALOG,
  threshold = MAX_MISSING
)
analysis <- missingness$data
catalog <- missingness$catalog

country_mean <- tapply(
  analysis$pfpr2_10,
  analysis$iso3,
  mean,
  na.rm = TRUE
)
analysis$country_mean_pfpr <- unname(country_mean[analysis$iso3])
analysis$country_mean_pfpr_gt_1 <- analysis$country_mean_pfpr > PFPR_FLOOR
analysis$pfpr_region_ge_1 <- is.finite(analysis$pfpr2_10) &
  analysis$pfpr2_10 >= PFPR_FLOOR

included_covariates <- catalog$variable[catalog$included_in_main]
analysis$complete_case_eligible <- if (length(included_covariates)) {
  rowSums(vapply(
    included_covariates,
    function(variable) analysis[[paste0(variable, "_imputed")]],
    logical(nrow(analysis))
  )) == 0
} else {
  TRUE
}

analysis$main_sample <- analysis$country_mean_pfpr_gt_1 &
  analysis$pfpr_region_ge_1 &
  is.finite(analysis$exposure) & analysis$exposure > 0 &
  is.finite(analysis$u5mr) & analysis$u5mr > 0 &
  is.finite(analysis$nnmr) & analysis$nnmr > 0 &
  is.finite(analysis$postneonatal_mortality) &
  analysis$postneonatal_mortality > 0

year_center <- round(mean(analysis$year[analysis$main_sample], na.rm = TRUE))
analysis$year_c <- analysis$year - year_center

write.csv(catalog, COVARIATE_CSV, row.names = FALSE)
write.csv(analysis, ANALYSIS_CSV, row.names = FALSE)
saveRDS(analysis, ANALYSIS_RDS)

inclusion <- data.frame(
  stage = c(
    "assembled survey-region panel",
    "country mean PfPR >1%",
    "region PfPR >=1%",
    "main shared-outcome sample",
    "complete-case eligible-covariate sensitivity"
  ),
  region_years = c(
    nrow(analysis),
    sum(analysis$country_mean_pfpr_gt_1),
    sum(analysis$pfpr_region_ge_1),
    sum(analysis$main_sample),
    sum(analysis$main_sample & analysis$complete_case_eligible)
  )
)
write.csv(
  inclusion,
  file.path(RESULTS_DIR, "analysis_inclusion_counts.csv"),
  row.names = FALSE
)

message(
  "Analysis dataset: ", nrow(analysis), " survey-region-years; ",
  length(unique(analysis$svkey)), " surveys; ",
  length(unique(analysis$iso3)), " countries."
)
message(
  "Main sample: ", sum(analysis$main_sample), " region-years; ",
  length(unique(analysis$iso3[analysis$main_sample])),
  " countries; year centered at ", year_center, "."
)
message(
  sum(catalog$included_in_main), " of ", nrow(catalog),
  " candidate covariates passed the <=5% missingness rule."
)

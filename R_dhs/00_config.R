# =============================================================================
# 00_config.R
# Shared configuration and functions for the rebuilt DHS/MIS-only analysis.
#
# Run all scripts from the repository root. Raw/access-controlled DHS records
# remain under data/ and are never written to results/. Only aggregate
# survey-region outputs and non-disclosive summaries are logged.
# =============================================================================

options(stringsAsFactors = FALSE)

REPO_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(REPO_ROOT, "R_dhs", "00_config.R"))) {
  stop("Run the rebuilt pipeline from the repository root.")
}

DATA_DIR       <- file.path(REPO_ROOT, "data")
DERIVED_DIR    <- file.path(DATA_DIR, "derived_dhs")
RESULTS_DIR    <- file.path(REPO_ROOT, "results", "dhs_rebuild")
DHS_LOCAL_DIR  <- file.path(DATA_DIR, "dhs")
DHS_CACHE_DIR  <- path.expand("~/.rdhs_cache/datasets_reformatted")
BOUNDARY_DIR   <- file.path(DATA_DIR, "dhs_boundaries")
MAP_RASTER_DIR <- file.path(DATA_DIR, "map_annual")
MAP_REGION_DIR <- file.path(DATA_DIR, "map_region_cache")
UNICEF_GLOBAL_CSV <- file.path(
  DATA_DIR, "fusion_GLOBAL_DATAFLOW_UNICEF_1.0_all.csv"
)
UNICEF_IMMUNISATION_CSV <- file.path(
  DERIVED_DIR, "unicef_immunisation_country_year.csv"
)
UNICEF_IMMUNISATION_SUMMARY_CSV <- file.path(
  RESULTS_DIR, "unicef_immunisation_extraction_summary.csv"
)
HIV_XLSX <- file.path(
  DATA_DIR, "HIV_Epidemiology_Children_Adolescents_2025.xlsx"
)
HIV_PANEL_CSV <- file.path(
  DERIVED_DIR, "hiv_prevalence_country_year.csv"
)
HIV_SUMMARY_CSV <- file.path(
  RESULTS_DIR, "hiv_prevalence_extraction_summary.csv"
)

for (d in c(DERIVED_DIR, RESULTS_DIR, BOUNDARY_DIR, MAP_RASTER_DIR, MAP_REGION_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

DHS_START_YEAR <- 2000L
DHS_END_YEAR   <- 2025L
PFPR_FLOOR     <- 1
MAX_MISSING    <- 0.05
# Cut year for the era subgroup splits. 2013 rather than 2010: it is the median
# region-year of the panel and halves the rows the models are fitted on almost
# exactly (553 before, 558 after), where a 2010 cut gives 367 against 744 and
# leaves the earlier era both less precise and with fewer distinct years for the
# calendar smooth. It is also the year the panel is centred on.
ERA_CUT <- 2013L
ERA_EARLY <- paste("before", ERA_CUT)
ERA_LATE <- paste(ERA_CUT, "onwards")

AF_REFERENCE   <- 1

MAP_DATASET_ID <- "Malaria__202508_Global_Pf_Parasite_Rate"
MAP_AFRICA_EXTENT <- matrix(c(-18, -35, 52, 38), nrow = 2)
GPW_TIF <- file.path(
  DATA_DIR, "pop",
  "gpw_v4_population_density_rev11_2020_2.5m.tif"
)

SURVEY_REGISTRY_CSV <- file.path(DERIVED_DIR, "survey_registry.csv")
MAP_REGION_CSV      <- file.path(DERIVED_DIR, "map_pfpr_by_survey_region.csv")
MAP_STATUS_CSV      <- file.path(RESULTS_DIR, "map_extraction_status.csv")
ANALYSIS_CSV        <- file.path(DERIVED_DIR, "dhs_analysis_dataset.csv")
ANALYSIS_RDS        <- file.path(DERIVED_DIR, "dhs_analysis_dataset.rds")
COVARIATE_CSV       <- file.path(RESULTS_DIR, "covariate_missingness.csv")
MODEL_BUNDLE_RDS    <- file.path(DERIVED_DIR, "main_model_bundle.rds")

required_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "Missing required packages: ", paste(missing, collapse = ", "),
      ". Install only after confirming they are on the approved package list."
    )
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# DHS Stata recodes are encoded in Latin-1, but some of them reach us decoded as
# Windows-1250, so a Latin-1 accented letter arrives as the Central-European
# letter sharing its byte: Kasai's i-diaeresis (byte 0xEF) becomes d-caron and
# Thies's e-grave (byte 0xE8) becomes c-caron. The region name then fails to
# join to the boundary's DHSREGEN and the region is silently dropped. Repair it
# by mapping each letter that exists in CP1250 but NOT in Latin-1 back to the
# Latin-1 letter at the same byte. Correctly encoded names are left untouched,
# so the substitution is a no-op for all but the affected surveys.
.CP1250_TO_LATIN1 <- local({
  bytes <- as.raw(128:255)
  decode <- function(from) {
    vapply(bytes, function(b) {
      tryCatch(iconv(rawToChar(as.raw(b)), from, "UTF-8"),
               error = function(e) NA_character_)
    }, character(1))
  }
  cp1250 <- decode("CP1250")
  latin1 <- decode("ISO-8859-1")
  keep <- !is.na(cp1250) & !is.na(latin1) &
    grepl("^[[:alpha:]]$", cp1250) & grepl("^[[:alpha:]]$", latin1) &
    is.na(suppressWarnings(iconv(cp1250, "UTF-8", "ISO-8859-1")))
  list(from = paste(cp1250[keep], collapse = ""),
       to = paste(latin1[keep], collapse = ""))
})

rkey <- function(x) {
  s <- as.character(x)
  if (nzchar(.CP1250_TO_LATIN1$from)) {
    s <- chartr(.CP1250_TO_LATIN1$from, .CP1250_TO_LATIN1$to, s)
  }
  # sub = "" drops any character iconv cannot transliterate instead of
  # returning NA for the whole string, which would silently void the key and
  # drop the region. Transliterable accents are unaffected.
  gsub(
    "[^a-z0-9]", "",
    tolower(iconv(s, "", "ASCII//TRANSLIT", sub = ""))
  )
}

survey_key <- function(filename) {
  x <- toupper(sub("\\..*$", "", basename(filename)))
  paste0(substr(x, 1, 2), substr(x, 5, 8))
}

# --- reconciling region names between a recode and its boundary --------------
# A recode and its boundary routinely name the same region differently, and
# exact key matching drops every such region silently. Five kinds of difference
# occur in this corpus:
#
#   language    an English boundary against a French or Portuguese recode
#               ("nord" for North, "regiao norte" for North)
#   word order  "Cuanza Norte" against "North Cuanza", "Kigali Ville" against
#               "Kigali City"
#   qualifier   "Centre (sans Ouagadougou)" against "Centre"; "Maritime
#               (sans agglomeration de Lome)" against "Maritime excluding
#               Greater Lome area"
#   spelling    "Matebeleland" for "Matabeleland", "Liberville" for
#               "Libreville", "extreme nor" truncated from "Extreme-Nord"
#   renaming    Namibia's Caprivi became Zambezi; Mauritanian boundaries carry
#               the English translation of the Arabic name
#
# Six stages run in order, each removing what it matched before the next runs.
# Nothing is ever guessed: every stage requires its pairing to be unique in both
# directions, and the three heuristic stages additionally require the two sides
# to describe the same number of units. That equal-count condition is what stops
# a mis-assignment where one boundary polygon spans several recode regions -
# Mali 2001 folds Kidal, Gao and Timbouctou into one polygon, Tanzania 2004
# reports 26 regions against 8 zones - and in both a pairing would attach one
# region's mortality to a polygon covering several.

# Words carrying no geographic content. "et"/"and" matter: they let
# "Tiris Zemour et Inchiri" agree with "Tiris Zemour and Inchiri".
REGION_STOPWORDS <- c(
  "de", "du", "des", "da", "do", "dos", "das", "d", "la", "le", "les", "el",
  "of", "the", "and", "et", "e", "y", "a", "o", "au", "aux",
  "region", "regiao", "regioes", "regions", "province", "provincia",
  "provincias", "provinces", "state", "states", "zone", "zones", "pvk"
)

# French and Portuguese forms folded onto one English token. Deliberately
# minimal: only words that actually occur in these region names, and NOT the
# English inflections (northern, southern), which the prefix stage already
# reconciles and which can be genuine distinct regions.
REGION_TOKEN_LEXICON <- c(
  nord = "north", norte = "north",
  sud = "south", sul = "south",
  est = "east", leste = "east", este = "east",
  ouest = "west", oeste = "west",
  centre = "central", center = "central", centro = "central",
  extreme = "far",
  ville = "city", cidade = "city",
  grande = "greater", agglomeration = "area",
  hiperendemica = "hyperendemic", hiperendemico = "hyperendemic",
  mesoendemica = "mesoendemic", mesoendemico = "mesoendemic",
  estavel = "stable", instavel = "unstable"
)

# A word introducing an exclusion clause; it and everything after it is dropped.
REGION_EXCLUSION_WORDS <- c(
  "sans", "without", "excluding", "excl", "exclu", "exceto", "sem",
  "notsurveyed"
)

# Pairs no rule could derive. Each is an official rename or a boundary carrying
# the English translation of a local name, so they are applied unconditionally
# rather than being gated like the heuristic stages.
REGION_NAME_SYNONYMS <- as.data.frame(rbind(
  # Namibia renamed Caprivi to Zambezi in 2013; the recode predates the boundary.
  c("caprivi", "zambezi"),
  # Comoros: recodes use the French names, boundaries the Comorian ones.
  c("moheli", "mwali"),
  c("anjouan", "ndzuwani"),
  c("ndzouani", "ndzuwani"),
  c("grandecomore", "ngazidja"),
  # Mauritania: boundaries translate the Arabic names.
  c("hodhechargui", "easternbasinregion"),
  c("hodhelchargui", "easternbasinregion"),
  c("hodhgharbi", "westernbasinregion"),
  c("hodhelgharbi", "westernbasinregion"),
  # Ethiopian recodes abbreviate Benishangul-Gumuz.
  c("bengumz", "benishangulgumuz")
), stringsAsFactors = FALSE)
names(REGION_NAME_SYNONYMS) <- c("recode", "boundary")

# Accented letters are folded to their base letter directly rather than through
# iconv's TRANSLIT, which on macOS writes the diacritic as a SEPARATE ASCII
# character ("regiao" comes back as "regi~ao"). That is harmless for rkey(),
# which strips punctuation, but here it would split a word in two and destroy
# the token.
REGION_ACCENTED <- paste0(
  "àáâãäåèéêë",
  "ìíîïòóôõö",
  "ùúûüçñýÿ"
)
REGION_UNACCENTED <- "aaaaaaeeeeiiiiooooouuuucnyy"

# Reduce a raw region label to a sorted bag of meaningful, language-normalised
# tokens, so that word order and language stop mattering.
region_tokens <- function(label) {
  text <- as.character(label)
  if (nzchar(.CP1250_TO_LATIN1$from)) {
    text <- chartr(.CP1250_TO_LATIN1$from, .CP1250_TO_LATIN1$to, text)
  }
  text <- tolower(text)
  text <- chartr(REGION_ACCENTED, REGION_UNACCENTED, text)
  text <- gsub("[^a-z0-9]+", " ", text)
  tokens <- strsplit(trimws(text), " +")[[1]]
  tokens <- tokens[nzchar(tokens)]

  clause <- which(tokens %in% REGION_EXCLUSION_WORDS)
  if (length(clause)) {
    tokens <- if (clause[1] == 1L) character(0) else tokens[seq_len(clause[1] - 1L)]
  }
  tokens <- tokens[!tokens %in% REGION_STOPWORDS]

  translated <- REGION_TOKEN_LEXICON[tokens]
  tokens <- ifelse(is.na(translated), tokens, translated)

  # Split written-together compass compounds so "Northeast" agrees with
  # "nord-est" once both halves are translated.
  tokens <- unlist(lapply(tokens, function(token) {
    parts <- regmatches(token, regexec("^(north|south)(east|west)$", token))[[1]]
    if (length(parts) == 3L) parts[2:3] else token
  }))
  sort(unique(tokens))
}

region_signature <- function(labels) {
  vapply(labels, function(label) paste(region_tokens(label), collapse = ""),
         character(1), USE.NAMES = FALSE)
}

match_region_keys <- function(recode_labels, boundary_labels,
                              min_prefix = 3L, max_edits = 2L) {
  frame <- function(labels) {
    labels <- as.character(labels)
    labels <- labels[!is.na(labels) & nzchar(trimws(labels))]
    out <- data.frame(label = labels, key = rkey(labels),
                      stringsAsFactors = FALSE)
    out <- out[nzchar(out$key) & !is.na(out$key), , drop = FALSE]
    out[!duplicated(out$key), , drop = FALSE]
  }
  recode <- frame(recode_labels)
  boundary <- frame(boundary_labels)
  pairs <- data.frame(from = character(0), to = character(0),
                      how = character(0), stringsAsFactors = FALSE)
  if (!nrow(recode) || !nrow(boundary)) return(pairs)

  equal_counts <- nrow(recode) == nrow(boundary)
  take <- function(from_keys, to_keys, how) {
    if (!length(from_keys)) return(invisible(NULL))
    pairs <<- rbind(pairs, data.frame(from = from_keys, to = to_keys,
                                      how = rep(how, length(from_keys)),
                                      stringsAsFactors = FALSE))
    recode <<- recode[!recode$key %in% from_keys, , drop = FALSE]
    boundary <<- boundary[!boundary$key %in% to_keys, , drop = FALSE]
  }

  # 1. exact
  exact <- intersect(recode$key, boundary$key)
  take(exact, exact, "exact")

  # 2. curated synonyms
  if (nrow(recode) && nrow(boundary)) {
    hit <- REGION_NAME_SYNONYMS[
      REGION_NAME_SYNONYMS$recode %in% recode$key &
        REGION_NAME_SYNONYMS$boundary %in% boundary$key, , drop = FALSE
    ]
    hit <- hit[!duplicated(hit$recode) & !duplicated(hit$boundary), , drop = FALSE]
    take(hit$recode, hit$boundary, "synonym")
  }

  # 3. language- and order-insensitive token signature. Ungated by unit count,
  #    because equal signatures are evidence rather than a heuristic, but the
  #    signature must be unique on both sides so the pairing stays one to one.
  if (nrow(recode) && nrow(boundary)) {
    recode$signature <- region_signature(recode$label)
    boundary$signature <- region_signature(boundary$label)
    shared <- intersect(
      recode$signature[!duplicated(recode$signature) &
                         !(duplicated(recode$signature, fromLast = TRUE))],
      boundary$signature[!duplicated(boundary$signature) &
                           !(duplicated(boundary$signature, fromLast = TRUE))]
    )
    shared <- shared[nzchar(shared)]
    if (length(shared)) {
      take(recode$key[match(shared, recode$signature)],
           boundary$key[match(shared, boundary$signature)], "canonical")
    }
  }

  # 4. prefix, 5. edit distance, 6. elimination - heuristics, so equal counts
  #    are required.
  if (equal_counts && nrow(recode) && nrow(boundary)) {
    unique_pairing <- function(candidates, how) {
      for (key in names(candidates)) {
        if (!key %in% recode$key) next
        hit <- intersect(candidates[[key]], boundary$key)
        if (length(hit) != 1L) next
        suitors <- sum(vapply(candidates, function(x) hit %in% x, logical(1)))
        if (suitors != 1L) next
        take(key, hit, how)
      }
    }

    prefix_candidates <- lapply(recode$key, function(key) {
      boundary$key[
        (nchar(key) >= min_prefix & startsWith(boundary$key, key)) |
          (nchar(boundary$key) >= min_prefix & startsWith(key, boundary$key))
      ]
    })
    names(prefix_candidates) <- recode$key
    unique_pairing(prefix_candidates, "prefix")

    if (nrow(recode) && nrow(boundary)) {
      fuzzy_candidates <- lapply(recode$key, function(key) {
        distance <- as.integer(utils::adist(key, boundary$key))
        allowed <- pmin(max_edits,
                        ceiling(0.25 * pmin(nchar(key), nchar(boundary$key))))
        boundary$key[distance <= allowed & distance > 0L]
      })
      names(fuzzy_candidates) <- recode$key
      unique_pairing(fuzzy_candidates, "fuzzy")
    }

    # A single leftover on each side is settled by elimination, not guessed.
    if (nrow(recode) == 1L && nrow(boundary) == 1L) {
      take(recode$key, boundary$key, "elimination")
    }
  }

  pairs
}

# Roll a recode's regions up to a coarser boundary's units, using a grouping
# read from a sibling survey of the same country that carries both levels. The
# target's own region names are reconciled against the donor's first, so a
# spelling change between rounds does not break the mapping. Returns NULL unless
# EVERY region the target reports can be placed, because a partial grouping
# would quietly drop respondents.
group_regions_to_boundary <- function(br, region_var, boundary_labels,
                                      survey, registry) {
  target_labels <- unique(as.character(br[[region_var]]))
  target_labels <- target_labels[!is.na(target_labels) & nzchar(target_labels)]
  target_keys <- unique(rkey(target_labels))

  siblings <- registry[registry$iso3 == survey$iso3 &
                         registry$svkey != survey$svkey, , drop = FALSE]
  if (!nrow(siblings)) return(NULL)
  siblings <- siblings[order(abs(siblings$year - survey$year)), , drop = FALSE]

  for (i in seq_len(nrow(siblings))) {
    path <- as.character(siblings$local_recode[i])
    if (!is.character(path) || is.na(path) || !file.exists(path)) next
    donor <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(donor)) next
    grouping <- derive_region_grouping(donor, boundary_labels, target_keys)
    rm(donor)
    gc(FALSE)
    if (is.null(grouping)) next

    alignment <- match_region_keys(target_labels, grouping$fine_label)
    if (nrow(alignment) < length(target_keys)) next
    boundary_for_target <- grouping$boundary[match(alignment$to, grouping$fine)]
    values <- boundary_for_target[
      match(rkey(as.character(br[[region_var]])), alignment$from)
    ]
    observed <- !is.na(br[[region_var]]) & nzchar(as.character(br[[region_var]]))
    if (any(is.na(values[observed]))) next

    return(list(
      values = values,
      source = sprintf("%s %s->%s", siblings$svkey[i],
                       grouping$donor_fine[1], grouping$donor_coarse[1])
    ))
  }
  NULL
}

best_region_var <- function(br, boundary_labels, prefer = "v024") {
  boundary_keys <- unique(rkey(as.character(boundary_labels)))
  boundary_keys <- boundary_keys[nzchar(boundary_keys) & !is.na(boundary_keys)]
  if (!length(boundary_keys)) return(prefer)

  levels_of <- function(v) {
    values <- unique(as.character(br[[v]]))
    values[!is.na(values) & nzchar(values)]
  }
  # Score with the full reconciliation rather than a bare key intersection, so a
  # variable that names the boundary's units in another language or word order
  # is not passed over. Senegal 2015 keeps its four zones in szone while v024
  # holds fourteen regions; Mali 2001 keeps the boundary's seven units in v023.
  reconciled <- function(v) {
    if (!v %in% names(br)) return(0L)
    values <- levels_of(v)
    if (length(values) < 2L || length(values) > 60L) return(0L)
    nrow(match_region_keys(values, boundary_labels))
  }

  if (prefer %in% names(br) && reconciled(prefer) == length(boundary_keys)) {
    return(prefer)
  }

  candidates <- names(br)[vapply(
    br, function(x) is.character(x) || is.factor(x), logical(1)
  )]
  if (!length(candidates)) return(prefer)
  scores <- vapply(candidates, reconciled, integer(1))
  if (max(scores) < 1L) return(prefer)

  best <- candidates[scores == max(scores)]
  # Among equally reconciling variables take the one with fewest levels: it is
  # the one actually at the boundary's granularity rather than a finer variable
  # that happens to share some names.
  sizes <- vapply(best, function(v) length(levels_of(v)), integer(1))
  best <- best[which.min(sizes)]

  preferred_score <- reconciled(prefer)
  if (preferred_score >= max(scores)) return(prefer)
  # Only displace the preferred variable when the alternative reconciles with at
  # least three fifths of the boundary, so a marginal gain cannot swap it out.
  if (max(scores) / length(boundary_keys) >= 0.6) best else prefer
}

# --- boundaries published coarser than the recode ----------------------------
# Some boundaries carry fewer units than the recode reports regions: Tanzania
# 2004 reports 26 regions against 8 zones, the Senegal continuous surveys of
# 2012 and 2014 report 14 against 4. Because MAP prevalence is computed per
# boundary polygon, the recode's regions must be rolled UP to the boundary's
# unit before anything can be joined. Pushing the zone's prevalence down onto
# each region instead would repeat one prevalence across many rows and
# pseudo-replicate it.
#
# The grouping is never invented. It is read out of a sibling survey of the same
# country whose recode carries BOTH levels at once - Senegal 2015 holds v024
# (14 regions) beside szone (4 zones), Tanzania 2010 holds v024 beside v023 - so
# the grouping is DHS's own definition rather than our guess.
derive_region_grouping <- function(donor, boundary_labels, target_keys) {
  boundary_keys <- unique(rkey(as.character(boundary_labels)))
  boundary_keys <- boundary_keys[nzchar(boundary_keys) & !is.na(boundary_keys)]
  candidates <- names(donor)[vapply(
    donor, function(x) is.character(x) || is.factor(x), logical(1)
  )]
  levels_of <- function(v) {
    values <- unique(as.character(donor[[v]]))
    values[!is.na(values) & nzchar(values)]
  }

  # A variable of the donor that reconciles one-to-one with every boundary unit.
  coarse <- NULL
  for (v in candidates) {
    values <- levels_of(v)
    if (length(values) < 2L || length(values) > 3L * length(boundary_keys)) next
    crosswalk <- match_region_keys(values, boundary_labels)
    if (nrow(crosswalk) == length(boundary_keys) &&
          length(unique(crosswalk$to)) == length(boundary_keys)) {
      coarse <- list(variable = v, crosswalk = crosswalk)
      break
    }
  }
  if (is.null(coarse)) return(NULL)

  # A finer variable of the donor that nests cleanly inside it AND actually
  # names the regions the target survey reports.
  for (v in candidates) {
    if (identical(v, coarse$variable)) next
    values <- levels_of(v)
    if (length(values) <= length(boundary_keys) || length(values) > 60L) next
    pairs <- unique(data.frame(
      fine_label = as.character(donor[[v]]),
      fine = rkey(as.character(donor[[v]])),
      coarse = rkey(as.character(donor[[coarse$variable]])),
      stringsAsFactors = FALSE
    ))
    pairs <- pairs[!duplicated(pairs$fine) | !duplicated(pairs$coarse), , drop = FALSE]
    pairs <- pairs[nzchar(pairs$fine) & nzchar(pairs$coarse), , drop = FALSE]
    if (!nrow(pairs)) next
    # The grouping must be many-to-one: a region spanning two boundary units
    # cannot be rolled up without splitting its mortality.
    if (any(table(pairs$fine) > 1L)) next
    if (mean(target_keys %in% pairs$fine) < 0.6) next
    pairs$boundary <- coarse$crosswalk$to[
      match(pairs$coarse, coarse$crosswalk$from)
    ]
    pairs <- pairs[!is.na(pairs$boundary), , drop = FALSE]
    if (!nrow(pairs)) next
    return(data.frame(fine = pairs$fine, fine_label = pairs$fine_label,
                      boundary = pairs$boundary,
                      donor_fine = v, donor_coarse = coarse$variable,
                      stringsAsFactors = FALSE))
  }
  NULL
}

weighted_mean_by_region <- function(value, weight, region, eligible) {
  keep <- eligible & is.finite(value) & is.finite(weight) &
    !is.na(region) & nzchar(region)
  if (!any(keep)) return(numeric(0))
  numerator <- tapply(weight[keep] * value[keep], region[keep], sum)
  denominator <- tapply(weight[keep], region[keep], sum)
  numerator / denominator
}

pick_named <- function(x, keys) {
  if (!length(x)) return(rep(NA_real_, length(keys)))
  unname(x[keys])
}

# Reference window for the mortality outcome, in months, passed to
# DHS.rates::chmort as Period. 12 rather than the DHS default of 60.
#
# The exposure is MAP prevalence in the SURVEY year, so a 60-month mortality
# window centres the outcome about 2.5 years before the exposure it is regressed
# on. A 12-month window nearly removes that mismatch. 27_horizon_lag_selection.R
# shows the dose-response is stable across horizons once the model structure is
# held fixed, so this costs little in the estimate; what it costs is precision,
# roughly doubling the interval on a regional rate (median 95% width 69.5 per
# 1000 against 35.7 at 60 months). The neonatal negative control is cleanest
# here: +0.09% per 10 PfPR points, p = 0.94, against +1.33%, p = 0.21 at 60.
#
# Sensitivity across 12, 24, 36, 48 and 60 months lives in script 27.
CHMORT_PERIOD <- 12L

# DHS.rates takes the sample strata from v022, and two different things go wrong
# with it. Some recodes ship v022 entirely empty (DR Congo 2007, Senegal 2008),
# and chmort aborts with "missing values in `strata'". Others ship it populated
# but malformed: Malawi 2004's clusters are not nested within its strata, and
# chmort aborts with "clusters not nested in strata at top level". Either way the
# whole survey used to vanish with no diagnostic.
#
# So rather than picking one stratum up front, try each candidate in turn and
# keep the first that chmort accepts. v022 is always tried first, so every survey
# whose v022 works is unaffected. Only the rate and the weighted exposure are
# read out of chmort, and neither depends on the stratification - it enters the
# standard errors alone - so the fallback cannot move the point estimates.
strata_candidates <- function(br) {
  usable <- function(v) v %in% names(br) && !all(is.na(br[[v]]))
  candidates <- character(0)
  if (usable("v022")) candidates <- c(candidates, "v022")
  if (usable("v023")) candidates <- c(candidates, "v023")
  candidates
}

mortality_by_region <- function(br, region_var) {
  required_packages("DHS.rates")
  if (!region_var %in% names(br)) return(NULL)

  candidates <- strata_candidates(br)
  # Last resort: the standard DHS design stratification of region by residence.
  if (all(c("v024", "v025") %in% names(br))) {
    br$strata_design <- as.integer(factor(paste(br$v024, br$v025)))
    candidates <- c(candidates, "strata_design")
  }
  if (!length(candidates)) return(NULL)

  rates <- NULL
  for (stratum in candidates) {
    rates <- tryCatch(
      suppressMessages(
        DHS.rates::chmort(br, Class = region_var, Strata = stratum,
                        Period = CHMORT_PERIOD)
      ),
      error = function(e) NULL
    )
    if (!is.null(rates)) break
  }
  if (is.null(rates)) return(NULL)

  u5 <- rates[
    grepl("^U5MR", rownames(rates)),
    c("Class", "R", "WN"),
    drop = FALSE
  ]
  nn <- rates[
    grepl("^NNMR", rownames(rates)),
    c("Class", "R"),
    drop = FALSE
  ]
  if (!nrow(u5) || !nrow(nn)) return(NULL)
  names(u5)[2:3] <- c("u5mr", "exposure")
  names(nn)[2] <- "nnmr"
  out <- merge(u5, nn, by = "Class")
  data.frame(
    regkey = rkey(out$Class),
    u5mr = as.numeric(out$u5mr),
    nnmr = as.numeric(out$nnmr),
    postneonatal_mortality = as.numeric(out$u5mr - out$nnmr),
    exposure = as.numeric(out$exposure)
  )
}

local_recode_path <- function(no_extension) {
  candidates <- c(
    file.path(DHS_LOCAL_DIR, paste0(no_extension, ".rds")),
    file.path(DHS_CACHE_DIR, paste0(no_extension, ".rds"))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) normalizePath(hit[1], winslash = "/") else NA_character_
}

wb_fetch <- function(indicator, date = "2000:2025") {
  required_packages(c("httr", "jsonlite"))
  url <- sprintf(
    paste0(
      "https://api.worldbank.org/v2/country/all/indicator/%s",
      "?date=%s&format=json&per_page=20000"
    ),
    indicator, date
  )
  for (attempt in seq_len(4)) {
    response <- tryCatch(
      httr::GET(url, httr::timeout(180)),
      error = function(e) NULL
    )
    if (!is.null(response) && httr::status_code(response) == 200) {
      parsed <- jsonlite::fromJSON(
        httr::content(response, "text", encoding = "UTF-8"),
        simplifyDataFrame = TRUE
      )
      # A retired or renamed indicator still returns HTTP 200, with the reason
      # in a one-element "message" payload rather than the usual two-element
      # [metadata, data] list. Report that instead of failing on a subscript.
      if (is.data.frame(parsed) && "message" %in% names(parsed)) {
        stop(
          "World Bank rejected indicator ", indicator, ": ",
          tryCatch(parsed$message[[1]]$value[1],
                   error = function(e) "unknown reason")
        )
      }
      if (!is.list(parsed) || length(parsed) < 2 || is.null(parsed[[2]])) {
        stop("World Bank returned no data for indicator ", indicator)
      }
      payload <- parsed[[2]]
      out <- data.frame(
        iso3 = payload$countryiso3code,
        year = as.integer(payload$date),
        value = as.numeric(payload$value)
      )
      return(out[is.finite(out$value) & nchar(out$iso3) == 3, ])
    }
    message("World Bank retry ", attempt, " for ", indicator)
  }
  stop("World Bank fetch failed: ", indicator)
}

read_or_fetch_wb <- function(filename, indicator, value_name) {
  path <- file.path(DATA_DIR, filename)
  if (file.exists(path)) {
    out <- read.csv(path, stringsAsFactors = FALSE)
  } else {
    message("Fetching World Bank indicator ", indicator, " -> ", filename)
    # An indicator the Bank has archived must not abort the whole rebuild. The
    # covariate simply stays missing and the existing missingness rule drops it,
    # which is what already happens to any series we cannot populate.
    out <- tryCatch(wb_fetch(indicator), error = function(e) {
      warning("World Bank indicator ", indicator, " unavailable (",
              conditionMessage(e), "); ", value_name,
              " will be left missing.", call. = FALSE)
      NULL
    })
    if (is.null(out)) {
      empty <- data.frame(iso3 = character(0), year = integer(0),
                          stringsAsFactors = FALSE)
      empty[[value_name]] <- numeric(0)
      return(empty)
    }
    names(out)[names(out) == "value"] <- value_name
    write.csv(out, path, row.names = FALSE)
  }
  if (!value_name %in% names(out) && "value" %in% names(out)) {
    names(out)[names(out) == "value"] <- value_name
  }
  out
}

nearest_panel_value <- function(panel, iso3, year, value_name) {
  rows <- panel[
    panel$iso3 == iso3 & is.finite(panel[[value_name]]),
    ,
    drop = FALSE
  ]
  if (!nrow(rows)) return(NA_real_)
  rows[[value_name]][which.min(abs(rows$year - year))]
}

single_impute <- function(x, country) {
  x[!is.finite(x)] <- NA_real_
  country_median <- ave(
    x, country,
    FUN = function(z) {
      value <- suppressWarnings(median(z, na.rm = TRUE))
      if (is.finite(value)) value else NA_real_
    }
  )
  overall <- suppressWarnings(median(x, na.rm = TRUE))
  out <- x
  missing <- !is.finite(out)
  out[missing] <- country_median[missing]
  out[!is.finite(out)] <- overall
  out
}

apply_missingness_rule <- function(data, catalog, threshold = MAX_MISSING) {
  stopifnot(all(c("variable", "level", "type") %in% names(catalog)))
  rows <- vector("list", nrow(catalog))

  for (i in seq_len(nrow(catalog))) {
    variable <- catalog$variable[i]
    if (!variable %in% names(data)) data[[variable]] <- NA_real_
    raw <- suppressWarnings(as.numeric(data[[variable]]))
    raw[!is.finite(raw)] <- NA_real_
    missing <- !is.finite(raw)
    missing_prop <- mean(missing)
    include <- is.finite(missing_prop) && missing_prop <= threshold &&
      any(is.finite(raw))

    analysis_name <- paste0(variable, "_analysis")
    imputed_name <- paste0(variable, "_imputed")
    data[[imputed_name]] <- missing
    data[[analysis_name]] <- if (include) {
      single_impute(raw, data$iso3)
    } else {
      rep(NA_real_, nrow(data))
    }

    rows[[i]] <- data.frame(
      variable = variable,
      level = catalog$level[i],
      type = catalog$type[i],
      missing_n = sum(missing),
      total_n = length(raw),
      missing_prop = missing_prop,
      included_in_main = include,
      imputation = if (include && any(missing)) {
        "country median; overall median fallback"
      } else if (include) {
        "none required"
      } else {
        "not imputed; excluded"
      }
    )
  }
  list(data = data, catalog = do.call(rbind, rows))
}

logit_percent <- function(x) {
  p <- pmin(pmax(x / 100, 0.005), 0.995)
  log(p / (1 - p))
}

make_ridge_matrix <- function(data, catalog, preprocessing = NULL) {
  included <- catalog$variable[as.logical(catalog$included_in_main)]
  if (!length(included)) stop("No covariates passed the missingness rule.")

  transformed <- lapply(included, function(variable) {
    x <- data[[paste0(variable, "_analysis")]]
    type <- catalog$type[match(variable, catalog$variable)]
    if (identical(type, "proportion")) logit_percent(x) else as.numeric(x)
  })
  names(transformed) <- included

  if (is.null(preprocessing)) {
    means <- vapply(transformed, mean, numeric(1), na.rm = TRUE)
    sds <- vapply(transformed, stats::sd, numeric(1), na.rm = TRUE)
    keep <- is.finite(means) & is.finite(sds) & sds > 0
    means <- means[keep]
    sds <- sds[keep]
    included <- names(means)
    preprocessing <- list(
      variables = included,
      means = means,
      sds = sds,
      types = setNames(
        catalog$type[match(included, catalog$variable)],
        included
      )
    )
  } else {
    included <- preprocessing$variables
  }

  matrix_out <- vapply(included, function(variable) {
    x <- data[[paste0(variable, "_analysis")]]
    if (identical(preprocessing$types[[variable]], "proportion")) {
      x <- logit_percent(x)
    }
    (x - preprocessing$means[[variable]]) / preprocessing$sds[[variable]]
  }, numeric(nrow(data)))
  matrix_out <- as.matrix(matrix_out)
  colnames(matrix_out) <- included

  list(matrix = matrix_out, preprocessing = preprocessing)
}

MODEL_SPECS <- data.frame(
  specification = c(
    "linear_no_interaction",
    "spline_no_interaction",
    "linear_time_interaction",
    "spline_time_interaction",
    "full_te_surface"
  ),
  prevalence_form = c("linear", "spline", "linear", "spline", "tensor"),
  time_interaction = c(FALSE, FALSE, TRUE, TRUE, TRUE)
)

# Whether the prevalence slope is allowed to vary by country. Set to FALSE:
# a country-specific random slope is not interpretable alongside the
# population-average effect the burden extrapolation needs, and the fixed slope
# is the quantity actually reported. The cost is that inference is no longer
# protected against slope heterogeneity between countries - regions within a
# country are treated as independent evidence about one common slope, which
# narrows the interval on that slope. `07_sensitivity_model_structure.R` and
# `15_specification_forest.R` both carry the with-slope fit as a sensitivity so
# the size of that effect stays visible.
INCLUDE_COUNTRY_PFPR_SLOPE <- FALSE

model_formula <- function(
    specification,
    include_country_slope = INCLUDE_COUNTRY_PFPR_SLOPE,
    spline_k = 6,
    year_k = 8,
    extra_terms = NULL) {
  prevalence <- switch(
    specification,
    linear_no_interaction = "pfpr10",
    spline_no_interaction = sprintf("s(pfpr10, k=%d)", spline_k),
    linear_time_interaction = "pfpr10 + pfpr10:year_c",
    spline_time_interaction = sprintf(
      "s(pfpr10, k=%d) + ti(pfpr10, year_c, k=c(%d,%d))",
      spline_k, spline_k, min(year_k, 6)
    ),
    full_te_surface = sprintf(
      "te(pfpr10, year_c, k=c(%d,%d))",
      spline_k, min(year_k, 6)
    ),
    stop("Unknown model specification: ", specification)
  )
  random_terms <- c(
    "s(country, bs='re')",
    if (include_country_slope) "s(country, pfpr10, bs='re')"
  )
  terms <- c(
    prevalence,
    "G",
    if (!identical(specification, "full_te_surface")) {
      sprintf("s(year_c, k=%d)", year_k)
    },
    # Additional unpenalised terms, e.g. a standardised auxiliary outcome used
    # as a covariate (36_neonatal_as_covariate.R). Standardised so that zero is
    # the mean, which is where newdata_at_mean() evaluates them.
    extra_terms,
    random_terms,
    "offset(log(exposure))"
  )
  as.formula(paste("deaths ~", paste(terms, collapse = " + ")))
}

fit_ridge_gam <- function(
    data,
    outcome,
    catalog,
    specification,
    method = "ML",
    preprocessing = NULL,
    include_country_slope = INCLUDE_COUNTRY_PFPR_SLOPE,
    spline_k = 6,
    year_k = 8,
    extra_terms = NULL) {
  required_packages("mgcv")
  keep <- is.finite(data[[outcome]]) & data[[outcome]] > 0 &
    is.finite(data$exposure) & data$exposure > 0 &
    is.finite(data$pfpr10) & is.finite(data$year_c) &
    !is.na(data$iso3)
  dd <- data[keep, , drop = FALSE]
  dd$country <- factor(dd$iso3)
  dd$deaths <- round(dd[[outcome]] / 1000 * dd$exposure)

  ridge <- make_ridge_matrix(dd, catalog, preprocessing)
  dd$G <- ridge$matrix
  penalty <- list(G = list(diag(ncol(ridge$matrix))))

  model <- mgcv::gam(
    model_formula(
      specification,
      include_country_slope = include_country_slope,
      spline_k = spline_k,
      year_k = year_k,
      extra_terms = extra_terms
    ),
    family = mgcv::nb(),
    method = method,
    paraPen = penalty,
    data = dd
  )
  list(
    model = model,
    data = dd,
    preprocessing = ridge$preprocessing,
    specification = specification,
    outcome = outcome,
    include_country_slope = include_country_slope,
    spline_k = spline_k,
    year_k = year_k,
    extra_terms = extra_terms
  )
}

model_summary_row <- function(fit) {
  model <- fit$model
  summary_model <- summary(model)
  linear <- "pfpr10" %in% rownames(summary_model$p.table)
  spline_row <- grep(
    "^s\\(pfpr10\\)",
    rownames(summary_model$s.table),
    value = TRUE
  )
  tensor_row <- grep(
    "^ti\\(pfpr10,year_c\\)",
    rownames(summary_model$s.table),
    value = TRUE
  )
  full_tensor_row <- grep(
    "^te\\(pfpr10,year_c\\)",
    rownames(summary_model$s.table),
    value = TRUE
  )
  interaction_linear <- "pfpr10:year_c" %in% rownames(summary_model$p.table)
  beta <- se <- p <- edf <- interaction_beta <- interaction_se <-
    interaction_edf <- interaction_p <- full_surface_edf <-
    full_surface_p <- NA_real_
  if (linear) {
    beta <- summary_model$p.table["pfpr10", "Estimate"]
    se <- summary_model$p.table["pfpr10", "Std. Error"]
    p <- summary_model$p.table[
      "pfpr10",
      grep("^Pr\\(", colnames(summary_model$p.table), value = TRUE)[1]
    ]
  } else if (length(spline_row)) {
    edf <- summary_model$s.table[spline_row[1], "edf"]
    p <- summary_model$s.table[spline_row[1], "p-value"]
  }
  if (interaction_linear) {
    interaction_beta <- summary_model$p.table["pfpr10:year_c", "Estimate"]
    interaction_se <- summary_model$p.table["pfpr10:year_c", "Std. Error"]
    interaction_p <- summary_model$p.table[
      "pfpr10:year_c",
      grep("^Pr\\(", colnames(summary_model$p.table), value = TRUE)[1]
    ]
  } else if (length(tensor_row)) {
    interaction_edf <- summary_model$s.table[tensor_row[1], "edf"]
    interaction_p <- summary_model$s.table[tensor_row[1], "p-value"]
  }
  if (length(full_tensor_row)) {
    full_surface_edf <- summary_model$s.table[full_tensor_row[1], "edf"]
    full_surface_p <- summary_model$s.table[full_tensor_row[1], "p-value"]
  }
  data.frame(
    specification = fit$specification,
    outcome = fit$outcome,
    n = nrow(model$model),
    countries = nlevels(model$model$country),
    edf_total = sum(model$edf),
    AIC = AIC(model),
    pfpr_beta = beta,
    pfpr_se = se,
    pct_change_per_10 = if (is.finite(beta)) 100 * (exp(beta) - 1) else NA,
    pct_change_lo = if (is.finite(beta)) 100 * (exp(beta - 1.96 * se) - 1) else NA,
    pct_change_hi = if (is.finite(beta)) 100 * (exp(beta + 1.96 * se) - 1) else NA,
    pfpr_smooth_edf = edf,
    pfpr_p = p,
    time_interaction_beta = interaction_beta,
    time_interaction_se = interaction_se,
    time_interaction_edf = interaction_edf,
    time_interaction_p = interaction_p,
    full_surface_edf = full_surface_edf,
    full_surface_p = full_surface_p
  )
}

newdata_at_mean <- function(model, pfpr10, year_c = 0) {
  n <- length(pfpr10)
  country_level <- levels(model$model$country)[1]
  g_names <- colnames(model$model$G)
  out <- data.frame(
    pfpr10 = pfpr10,
    year_c = rep(year_c, length.out = n),
    exposure = 1,
    country = factor(rep(country_level, n), levels = levels(model$model$country))
  )
  out$G <- matrix(0, nrow = n, ncol = length(g_names))
  colnames(out$G) <- g_names
  # Any further model variable (an extra_terms covariate) is held at zero, its
  # standardised mean, so the prediction stays population-average.
  extra <- setdiff(all.vars(stats::formula(model)),
                   c("deaths", "pfpr10", "year_c", "G", "country", "exposure"))
  for (v in extra) if (!v %in% names(out)) out[[v]] <- 0
  out
}

population_lpmatrix <- function(model, newdata) {
  X <- predict(model, newdata, type = "lpmatrix")
  random_columns <- grep(
    "^s\\(country\\)|^s\\(country,pfpr10\\)",
    colnames(X)
  )
  if (length(random_columns)) X[, random_columns] <- 0
  X
}

link_prediction <- function(model, newdata) {
  X <- population_lpmatrix(model, newdata)
  beta <- coef(model)
  fit <- as.numeric(X %*% beta)
  se <- sqrt(rowSums((X %*% vcov(model)) * X))
  data.frame(fit = fit, se = se)
}

link_difference <- function(model, high, reference) {
  Xh <- population_lpmatrix(model, high)
  Xr <- population_lpmatrix(model, reference)
  dX <- Xh - Xr
  beta <- coef(model)
  fit <- as.numeric(dX %*% beta)
  se <- sqrt(rowSums((dX %*% vcov(model)) * dX))
  data.frame(fit = fit, se = se)
}

af_from_model <- function(model, prevalence, year_c = 0, reference = AF_REFERENCE) {
  high <- newdata_at_mean(model, prevalence / 10, year_c)
  low <- newdata_at_mean(model, rep(reference / 10, length(prevalence)), year_c)
  delta <- link_difference(model, high, low)
  data.frame(
    prevalence = prevalence,
    af = pmax(1 - exp(-delta$fit), 0),
    lo = pmax(1 - exp(-(delta$fit - 1.96 * delta$se)), 0),
    hi = pmax(1 - exp(-(delta$fit + 1.96 * delta$se)), 0)
  )
}

read_analysis_data <- function() {
  if (file.exists(ANALYSIS_RDS)) {
    readRDS(ANALYSIS_RDS)
  } else if (file.exists(ANALYSIS_CSV)) {
    out <- read.csv(ANALYSIS_CSV, stringsAsFactors = FALSE)
    for (v in intersect(
      c("main_sample", "pfpr_region_ge_1", "country_mean_pfpr_gt_1",
        "complete_case_eligible"),
      names(out)
    )) out[[v]] <- as.logical(out[[v]])
    out
  } else {
    stop("Analysis dataset not found. Run R_dhs/03_build_analysis_dataset.R.")
  }
}

## ---- which outcome the brms scripts (30-31) model ---------------------------
# Post-neonatal mortality is the primary outcome. Setting BRMS_OUTCOME=neonatal
# in the environment refits the identical Stan programs to neonatal mortality,
# the negative control, and every output file gains a "_neonatal" suffix so the
# primary results are never overwritten. The neonatal attributable fraction is
# deliberately NOT clipped at zero: a control should be allowed to sit on zero
# with its interval straddling it, and clipping would push the summary above
# zero by construction.
brms_outcome <- function() {
  choice <- tolower(Sys.getenv("BRMS_OUTCOME", "postneonatal"))
  table <- list(
    postneonatal = list(name = "postneonatal",
                        column = "postneonatal_mortality",
                        label = "post-neonatal mortality", suffix = "",
                        title_suffix = "", clip = TRUE),
    neonatal = list(name = "neonatal", column = "nnmr",
                    label = "neonatal mortality", suffix = "_neonatal",
                    title_suffix = " (neonatal negative control)",
                    clip = FALSE)
  )
  if (!choice %in% names(table)) {
    stop("BRMS_OUTCOME must be one of: ", paste(names(table), collapse = ", "))
  }
  table[[choice]]
}

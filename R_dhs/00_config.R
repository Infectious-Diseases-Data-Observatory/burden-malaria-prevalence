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

best_region_var <- function(br, target_keys, prefer = "v024") {
  target_keys <- unique(target_keys[nzchar(target_keys)])
  coverage <- function(v) {
    if (!v %in% names(br)) return(0)
    mean(target_keys %in% rkey(unique(as.character(br[[v]]))))
  }
  if (prefer %in% names(br) && coverage(prefer) >= 0.5) return(prefer)

  candidates <- names(br)[vapply(
    br,
    function(x) is.character(x) || is.factor(x),
    logical(1)
  )]
  if (!length(candidates)) return(prefer)
  scores <- vapply(
    candidates,
    function(v) sum(target_keys %in% rkey(unique(as.character(br[[v]])))),
    integer(1)
  )
  best <- candidates[which.max(scores)]
  if (max(scores) > 0 && coverage(best) >= 0.6) best else prefer
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

# DHS.rates takes the sample strata from v022. A few recodes ship it empty
# (DR Congo 2007 and Senegal 2008 among them), and chmort then aborts with
# "missing values in `strata'", which silently cost us the whole survey. Fall
# back to the next usable stratum definition: v023 where the recode populates
# it, otherwise the standard DHS design stratification of region by urban/rural.
# Only the rate and the weighted exposure are read out of chmort, and neither
# depends on the stratification - it enters the standard errors alone - so the
# fallback cannot move the point estimates.
strata_variable <- function(br) {
  usable <- function(v) v %in% names(br) && !all(is.na(br[[v]]))
  if (usable("v022")) return("v022")
  if (usable("v023")) return("v023")
  if (all(c("v024", "v025") %in% names(br))) return(NA_character_)
  NULL
}

mortality_by_region <- function(br, region_var) {
  required_packages("DHS.rates")
  if (!region_var %in% names(br)) return(NULL)
  strata <- strata_variable(br)
  if (is.null(strata)) return(NULL)
  if (is.na(strata)) {
    br$strata_fallback <- as.integer(factor(paste(br$v024, br$v025)))
    strata <- "strata_fallback"
  }
  rates <- tryCatch(
    suppressMessages(
      DHS.rates::chmort(br, Class = region_var, Strata = strata)
    ),
    error = function(e) NULL
  )
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

model_formula <- function(
    specification,
    include_country_slope = TRUE,
    spline_k = 6,
    year_k = 8) {
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
    include_country_slope = TRUE,
    spline_k = 6,
    year_k = 8) {
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
      year_k = year_k
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
    year_k = year_k
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

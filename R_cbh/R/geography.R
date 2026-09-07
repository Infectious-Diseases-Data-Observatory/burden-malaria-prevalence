# Deterministic geographic normalization migrated from R_dhs/00_config.R.
# No fuzzy, prefix or elimination matching is accepted in this pipeline.

cbh_match_regions <- function(labels, boundaries) {
  cbh_require(boundaries, c("region", "regkey"), "Boundary region table")
  labels <- unique(as.character(labels))
  labels <- labels[!is.na(labels) & nzchar(trimws(labels))]
  target <- unique(boundaries[c("region", "regkey")])
  cbh_unique(target, "regkey", "Boundary regions")
  out <- data.frame(source_label = labels, source_key = cbh_rkey(labels),
                    regkey = rep(NA_character_, length(labels)),
                    match_method = rep("unmatched", length(labels)))
  # Require unique normalized target names; the cached regkey is the join key.
  tkey <- cbh_rkey(target$region)
  unambiguous <- function(x) !duplicated(x) & !duplicated(x, fromLast = TRUE)
  hit <- match(out$source_key, tkey)
  good <- !is.na(hit) & unambiguous(tkey)[hit]
  out$regkey[good] <- target$regkey[hit[good]]
  out$match_method[good] <- "exact"
  syn <- REGION_NAME_SYNONYMS
  hit <- match(out$source_key, syn$recode)
  dest <- match(syn$boundary[hit], tkey)
  good <- is.na(out$regkey) & !is.na(dest) & unambiguous(tkey)[dest]
  out$regkey[good] <- target$regkey[dest[good]]
  out$match_method[good] <- "curated_synonym"
  tsig <- cbh_region_signature(target$region)
  hit <- match(cbh_region_signature(out$source_label), tsig)
  good <- is.na(out$regkey) & !is.na(hit) & unambiguous(tsig)[hit] & nzchar(tsig[hit])
  out$regkey[good] <- target$regkey[hit[good]]
  out$match_method[good] <- "canonical"
  out
}

cbh_geography <- function(br, survey, rules, boundaries, overrides, registry) {
  variable <- rules$region_var
  if (!variable %in% names(br)) stop("Configured region variable is unavailable.")
  labels <- cbh_column(br, variable)
  crosswalk <- cbh_match_regions(labels, boundaries)
  if (!is.na(rules$group_donor) && nzchar(rules$group_donor)) {
    donor <- registry[registry$svkey == rules$group_donor, , drop = FALSE]
    if (nrow(donor) != 1L || !file.exists(donor$local_recode)) stop("Configured geographic donor is unavailable.")
    d <- cbh_read_recode(donor$local_recode, c(rules$group_fine, rules$group_coarse))
    cbh_require(d, c(rules$group_fine, rules$group_coarse), "Geographic donor")
    pairs <- unique(data.frame(fine = cbh_column(d, rules$group_fine),
                               coarse = cbh_column(d, rules$group_coarse)))
    pairs <- pairs[complete.cases(pairs) & nzchar(pairs$fine) & nzchar(pairs$coarse), ]
    if (anyDuplicated(pairs$fine)) stop("Donor geography is not a many-to-one grouping.")
    coarse <- cbh_match_regions(pairs$coarse, boundaries)
    pairs$regkey <- coarse$regkey[match(pairs$coarse, coarse$source_label)]
    # Exact/canonical matching also reconciles target fine labels to the donor.
    fine <- cbh_match_regions(labels, data.frame(region = pairs$fine, regkey = pairs$fine))
    crosswalk$regkey <- pairs$regkey[match(fine$regkey, pairs$fine)]
    crosswalk$match_method <- ifelse(is.na(crosswalk$regkey), "unmatched", "explicit_donor_grouping")
  }
  ov <- overrides[overrides$svkey == survey$svkey & overrides$region_var == variable, , drop = FALSE]
  if (nrow(ov)) {
    cbh_unique(ov, "source_label", "Region overrides")
    if (any(!ov$regkey %in% boundaries$regkey)) stop("Region override targets an unknown boundary key.")
    j <- match(crosswalk$source_label, ov$source_label)
    crosswalk$regkey[!is.na(j)] <- ov$regkey[j[!is.na(j)]]
    crosswalk$match_method[!is.na(j)] <- "explicit_override"
  }
  # Conservatively distinguish survey boundary versions unless explicitly harmonized.
  crosswalk$region <- ifelse(is.na(crosswalk$regkey), NA_character_,
                            paste(survey$iso3, survey$svkey, crosswalk$regkey, sep = ":"))
  if (nrow(ov)) {
    j <- match(crosswalk$source_label, ov$source_label)
    stable <- !is.na(j) & !is.na(ov$region_id[j]) & nzchar(ov$region_id[j])
    crosswalk$region[stable] <- paste(survey$iso3, ov$region_id[j[stable]], sep = ":")
  }
  j <- match(labels, crosswalk$source_label)
  crosswalk$svkey <- rep(survey$svkey, nrow(crosswalk))
  crosswalk$region_var <- rep(variable, nrow(crosswalk))
  list(regkey = crosswalk$regkey[j], region = crosswalk$region[j], crosswalk = crosswalk)
}
.CP1250_TO_LATIN1 <- local({
    bytes <- as.raw(128:255)
    decode <- function(from) {
        vapply(bytes, function(b) {
            tryCatch(iconv(rawToChar(as.raw(b)), from, "UTF-8"), error = function(e) NA_character_)
        }, character(1))
    }
    cp1250 <- decode("CP1250")
    latin1 <- decode("ISO-8859-1")
    keep <- !is.na(cp1250) & !is.na(latin1) & grepl("^[[:alpha:]]$", cp1250) & grepl("^[[:alpha:]]$", 
        latin1) & is.na(suppressWarnings(iconv(cp1250, "UTF-8", "ISO-8859-1")))
    list(from = paste(cp1250[keep], collapse = ""), to = paste(latin1[keep], collapse = ""))
})

cbh_rkey <- function(x) {
    s <- as.character(x)
    if (nzchar(.CP1250_TO_LATIN1$from)) {
        s <- chartr(.CP1250_TO_LATIN1$from, .CP1250_TO_LATIN1$to, s)
    }
    gsub("[^a-z0-9]", "", tolower(iconv(s, "", "ASCII//TRANSLIT", sub = "")))
}

REGION_STOPWORDS <- c("de", "du", "des", "da", "do", "dos", "das", "d", "la", "le", "les", "el", 
    "of", "the", "and", "et", "e", "y", "a", "o", "au", "aux", "region", "regiao", "regioes", "regions", 
    "province", "provincia", "provincias", "provinces", "state", "states", "zone", "zones", "pvk")

REGION_TOKEN_LEXICON <- c(nord = "north", norte = "north", sud = "south", sul = "south", est = "east", 
    leste = "east", este = "east", ouest = "west", oeste = "west", centre = "central", center = "central", 
    centro = "central", extreme = "far", ville = "city", cidade = "city", grande = "greater", agglomeration = "area", 
    hiperendemica = "hyperendemic", hiperendemico = "hyperendemic", mesoendemica = "mesoendemic", 
    mesoendemico = "mesoendemic", estavel = "stable", instavel = "unstable")

REGION_NAME_SYNONYMS <- as.data.frame(rbind(c("caprivi", "zambezi"), c("moheli", "mwali"), c("anjouan", 
    "ndzuwani"), c("ndzouani", "ndzuwani"), c("grandecomore", "ngazidja"), c("hodhechargui", "easternbasinregion"), 
    c("hodhelchargui", "easternbasinregion"), c("hodhgharbi", "westernbasinregion"), c("hodhelgharbi", 
        "westernbasinregion"), c("bengumz", "benishangulgumuz")), stringsAsFactors = FALSE)

names(REGION_NAME_SYNONYMS) <- c("recode", "boundary")

REGION_ACCENTED <- paste0("àáâãäåèéêë", "ìíîïòóôõö", "ùúûüçñýÿ")

REGION_UNACCENTED <- "aaaaaaeeeeiiiiooooouuuucnyy"

cbh_region_tokens <- function(label) {
    text <- as.character(label)
    if (nzchar(.CP1250_TO_LATIN1$from)) {
        text <- chartr(.CP1250_TO_LATIN1$from, .CP1250_TO_LATIN1$to, text)
    }
    text <- tolower(text)
    text <- chartr(REGION_ACCENTED, REGION_UNACCENTED, text)
    text <- gsub("[^a-z0-9]+", " ", text)
    tokens <- strsplit(trimws(text), " +")[[1]]
    tokens <- tokens[nzchar(tokens)]
    # Preserve exclusion clauses: a region excluding a city is a different area.
    tokens <- tokens[!tokens %in% REGION_STOPWORDS]
    translated <- REGION_TOKEN_LEXICON[tokens]
    tokens <- ifelse(is.na(translated), tokens, translated)
    tokens <- unlist(lapply(tokens, function(token) {
        parts <- regmatches(token, regexec("^(north|south)(east|west)$", token))[[1]]
        if (length(parts) == 3L) 
            parts[2:3]
        else token
    }))
    sort(unique(tokens))
}

cbh_region_signature <- function(labels) {
    vapply(labels, function(label) paste(cbh_region_tokens(label), collapse = ""), character(1), USE.NAMES = FALSE)
}

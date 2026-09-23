#!/usr/bin/env Rscript
# Inventory of MICS microdata: birth-history content and availability of the
# primary model's covariates, per survey. Reads the extracted SPSS files under
# data/MICS_extracted (git-ignored); writes aggregate tables only.
suppressPackageStartupMessages({ library(haven); library(data.table) })
root <- "data/MICS_extracted"
out <- "results/mics_inventory"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
meta <- fread("data/derived_mics/variable_metadata.csv")
surveys <- sort(list.dirs(root, recursive = FALSE, full.names = FALSE))

first_of <- function(d, cands) { nm <- names(d); j <- match(tolower(cands), tolower(nm)); j <- j[!is.na(j)]
  if (length(j)) nm[j[1]] else NA_character_ }
lab_hits <- function(s, file, pattern) {
  z <- meta[survey == s & get("file") == file & grepl(pattern, label, ignore.case = TRUE)]
  z$variable
}
pct_nonmissing <- function(path, var) {
  if (is.na(var) || !file.exists(path)) return(NA_real_)
  x <- tryCatch(read_sav(path, col_select = all_of(var))[[1]], error = function(e) NULL)
  if (is.null(x)) return(NA_real_)
  x <- zap_missing(x); v <- suppressWarnings(as.numeric(x))
  round(100 * mean(!is.na(v) & !(v %in% c(9, 99, 999, 9999, 99.9, 99.99, 999.9))), 1)
}
safe_read <- function(path, cols) {
  if (!file.exists(path)) return(NULL)
  have <- names(read_sav(path, n_max = 0)); cols <- intersect(cols, have)
  if (!length(cols)) return(NULL)
  as.data.table(lapply(read_sav(path, col_select = all_of(cols)), function(x) as.numeric(zap_missing(x))))
}

rows <- list()
for (s in surveys) {
  f <- function(x) file.path(root, s, x)
  files <- basename(list.files(file.path(root, s), pattern = "\\.sav$"))
  rnd <- as.integer(sub(".*MICS(\\d).*", "\\1", s))
  iso3 <- sub("^([A-Z]{3}).*", "\\1", s)
  year <- as.integer(sub("^[A-Z]{3}(\\([^)]*\\))?_(\\d{4})_.*", "\\2", s))
  subnational <- grepl("\\(", s) & !grepl("including current", s)
  r <- list(survey = s, iso3 = iso3, year = year, round = rnd, subnational = subnational,
            n_files = length(files), has_bh = "bh.sav" %in% files)
  # ---- birth history ------------------------------------------------------------
  # MICS naming varies by round and country: CMC birth date may be BH4C or CCDOB or
  # only month/year; age at death may be imputed months (BH9C) or units + number
  # (BH9U/BH9N, BH9UC/BH9NC, BH9A/BH9B, Portuguese HN10_U/HN10_I); interview date may
  # be WDOI, WM6C, CMCDOIW or month/year. Names are matched case-insensitively.
  if (r$has_bh) {
    raw <- read_sav(f("bh.sav"))
    nm <- names(raw); low <- tolower(nm)
    col <- function(...) { for (a in tolower(c(...))) { j <- match(a, low); if (!is.na(j)) return(as.numeric(zap_missing(raw[[j]]))) }; NULL }
    cmc_from <- function(mo, yr) { if (is.null(mo) || is.null(yr)) return(NULL)
      yr <- ifelse(yr < 100, yr + 1900, yr); ifelse(mo %in% 1:12 & yr > 1900 & yr < 2030, (yr - 1900) * 12 + mo, NA) }
    dob <- col("BH4C", "CCDOB"); dob_src <- if (!is.null(dob)) "cmc" else "month_year"
    if (is.null(dob)) dob <- cmc_from(col("BH4MC", "BH4M", "HN5_MES"), col("BH4YC", "BH4Y", "HN5_ANO"))
    doi <- col("WDOI", "WM6C", "CMCDOIW")
    if (is.null(doi)) doi <- cmc_from(col("WM6M"), col("WM6Y"))
    alive <- col("BH5", "HN6")
    agem <- col("BH9C")
    if (is.null(agem)) {
      u <- col("BH9UC", "BH9U", "BH9A", "HN10_U"); n <- col("BH9NC", "BH9N", "BH9B", "HN10_I")
      if (!is.null(u) && !is.null(n)) agem <- ifelse(u == 1 & n < 98, n / 30.4375,
        ifelse(u == 2 & n < 98, n, ifelse(u == 3 & n < 98, 12 * n, NA)))
    }
    r$bh_dob_source <- dob_src
    r$bh_has_core <- !is.null(dob) && !is.null(doi) && !is.null(alive) && !is.null(agem)
    if (isTRUE(r$bh_has_core)) {
      # Survival coded 1 = alive, 2 = dead in every MICS round inspected.
      died <- !is.na(alive) & alive == 2
      recent <- !is.na(dob) & !is.na(doi) & (doi - dob) >= 0 & (doi - dob) < 60
      r$births <- nrow(raw); r$deaths <- sum(died)
      r$deaths_age_known <- sum(died & is.finite(agem))
      r$deaths_under5 <- sum(died & is.finite(agem) & agem < 60)
      r$births_last5y <- sum(recent)
      r$deaths_last5y_under5 <- sum(recent & died & is.finite(agem) & agem < 60)
      r$interview_year <- as.integer(1900 + floor((median(doi, na.rm = TRUE) - 1) / 12))
      r$alive_codes <- paste(sort(unique(na.omit(alive))), collapse = ",")
    }
    rm(raw); gc(FALSE)
  }
  # ---- region -------------------------------------------------------------------
  src <- if ("wm.sav" %in% files) f("wm.sav") else if ("hh.sav" %in% files) f("hh.sav") else NA
  if (!is.na(src)) {
    h <- read_sav(src, n_max = 0); rv <- first_of(h, c("HH7", "HHREG", "WMREG", "HH7A"))
    r$region_var <- rv
    if (!is.na(rv)) {
      x <- read_sav(src, col_select = all_of(rv))[[1]]; labs <- attr(x, "labels")
      r$n_regions <- length(unique(na.omit(as.numeric(zap_missing(x)))))
      r$region_labels <- paste(head(names(labs)[labs %in% unique(as.numeric(x))], 60), collapse = " | ")
    }
  }
  # ---- covariate sources (candidate variables and completeness) -------------------
  if ("wm.sav" %in% files) {
    w <- read_sav(f("wm.sav"), n_max = 0)
    r$edu_var <- first_of(w, c("welevel", "WB4", "WM11", "melevel"))
    r$wealth_var <- first_of(w, c("windex5", "WLTHIND5", "wlthind5"))
    r$urban_var <- first_of(w, c("HH6"))
    r$first_birth_var <- first_of(w, c("WDOBFC"))
    # Place of delivery for the last birth: MN8 in MICS3, MN18 in MICS4-5, MN20 in MICS6. Match by
    # name first; the label fallback excludes antenatal, attendant and sibling-mortality questions.
    std_del <- if (rnd <= 3) c("MN8", "MN7_A") else if (rnd <= 5) c("MN18", "MN17_PLACE") else c("MN20", "MN19")
    del <- first_of(w, std_del)
    if (is.na(del)) {
      z <- meta[survey == s & file == "wm.sav" & grepl(paste0("place of delivery|where did you give birth|",
        "lieu.*accouch|o[uù].*accouch|local do parto|onde.*(deu [àa] luz|parto)"), label, ignore.case = TRUE) &
        !grepl("antenatal|prenatal|pr[ée]natal|pr[ée]-natal|assist|attend|sister|brother|s[oœ]ur|fr[eè]re|irm[aã]", label,
        ignore.case = TRUE)]$variable
      del <- if (length(z)) z[1] else NA
    }
    r$delivery_var <- del
  }
  if ("hh.sav" %in% files) {
    hh <- read_sav(f("hh.sav"), n_max = 0)
    r$water_var <- first_of(hh, c("WS1"))
    # Toilet type by label: the variable number varies (WS7 in MICS3 and early MICS4, WS8 in later MICS4-5,
    # WS11 in MICS6), and neighbouring items record sharing or the number of households.
    tl <- meta[survey == s & file == "hh.sav" & grepl("kind of toilet|type of toilet|toilet facility|type de toilette|lieux d.?aisance|type de latrine|tipo de (casa de banho|sanit|retrete)|instala..o sanit", label, ignore.case = TRUE) &
      !grepl("shared|partag|compartilh|households using|number of|nombre de|n.mero de|location|emplacement|localiza", label, ignore.case = TRUE)]$variable
    tl <- c(tl[grepl("^WS", tl, ignore.case = TRUE)], tl[!grepl("^WS", tl, ignore.case = TRUE)])   # reported before observed
    r$toilet_var <- if (length(tl)) tl[1] else if (rnd <= 3) first_of(hh, "WS7") else if (rnd <= 5) first_of(hh, "WS8") else first_of(hh, "WS11")
    el <- lab_hits(s, "hh.sav", "^electricity$|has electricity|electricity$|[ée]lectricit[ée]|electricidade")
    r$electricity_var <- if (length(el)) el[1] else first_of(hh, c("HC8", "HC8A", "HC9A"))
  }
  if ("ch.sav" %in% files) {
    dtp <- lab_hits(s, "ch.sav", "dpt|dtp|dtc|penta|diphth|dipht[ée]rie|difteria|t[ée]tracoq")
    r$dtp_vars <- length(dtp)
    # Third dose: the dose digit follows the vaccine code (IM3D3D, IM3PENTA3D, IM6PENTA3D);
    # a bare "3" also matches the module number, so match the vaccine code and dose together.
    # MICS3 codes doses by letter (IM4A/B/C, or IM5A/B/C where DPT is combined with HepB), and some
    # MICS3 labels are wrong (Zimbabwe 2009), so the variable name decides dose 3 there.
    r$dtp3_var <- { z <- dtp[grepl("(penta|dpt|dtp|dtc|[^a-z]d)3", dtp, ignore.case = TRUE) |
      (rnd <= 3 & grepl("^im[0-9]c", dtp, ignore.case = TRUE))]
      if (!length(z)) z <- dtp[grepl("third|(dpt|dtp|dtc|penta)[^0-9]{0,8}3|dose ?3", meta[survey == s & file == "ch.sav" &
        variable %in% dtp]$label[match(dtp, meta[survey == s & file == "ch.sav" & variable %in% dtp]$variable)],
        ignore.case = TRUE)]
      if (length(z)) z[1] else NA }
    ms <- lab_hits(s, "ch.sav", "measles|mcv|mmr|rougeole|sarampo|rub[ée]ol|rubella|\\bmr ?1\\b|\\brr ?1\\b|\\bvar\\b")
    r$measles_var <- if (length(ms)) ms[1] else NA
  }
  zfile <- intersect(c("ch.sav", "who_z.sav"), files)
  haz <- whz <- character()
  for (zf in zfile) {
    haz <- c(haz, meta[survey == s & file == zf & (variable %in% c("HAZ2", "HAZ", "zhfa", "haz") |
      grepl("height.for.age|length.for.age|taille.pour.[âa]ge|altura.para.idade", label, ignore.case = TRUE))]$variable)
    whz <- c(whz, meta[survey == s & file == zf & (variable %in% c("WHZ2", "WHZ", "zwfl", "whz") |
      grepl("weight.for.height|weight.for.length|poids.pour.taille|peso.para.altura", label, ignore.case = TRUE))]$variable)
  }
  r$haz_var <- if (length(haz)) haz[1] else NA
  r$whz_var <- if (length(whz)) whz[1] else NA
  # Raw measurements, from which z-scores could be computed if none are supplied.
  if ("ch.sav" %in% files) {
    wt <- lab_hits(s, "ch.sav", "^weight|child's weight|weight \\(kilograms\\)|weight in kilo")
    ht <- lab_hits(s, "ch.sav", "^height|^length|height \\(cm|length \\(cm|height or length")
    r$raw_anthro <- length(wt) > 0 && length(ht) > 0
  }
  rows[[s]] <- as.data.table(lapply(r, function(x) if (length(x) == 0) NA else x))
}
inv <- rbindlist(rows, fill = TRUE)
fwrite(inv, file.path(out, "survey_inventory.csv"))
cat("surveys:", nrow(inv), "| with birth history:", sum(inv$has_bh, na.rm = TRUE),
    "| with usable core BH variables:", sum(inv$bh_has_core %in% TRUE), "\n")

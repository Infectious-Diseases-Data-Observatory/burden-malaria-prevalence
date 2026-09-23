#!/usr/bin/env Rscript
# The 13 regional covariates for each MICS survey x analysis region, computed from the
# microdata to match the DHS definitions used by the primary analysis, then the same two
# fallbacks: national WUENIC DTP3/measles for the survey year, then the within-survey mean of
# available regions. Writes data/derived_mics/regional_covariates_wide_mics.csv plus a long
# table with denominators, and an audit in results/mics_inventory/covariates/.
suppressPackageStartupMessages({ library(haven); library(data.table) })
source("R_cbh/R/geography.R")
root <- "data/MICS_extracted"; D <- "data/derived_mics"
aud <- "results/mics_inventory/covariates"; dir.create(aud, recursive = TRUE, showWarnings = FALSE)
setup <- fread(file.path(D, "survey_setup.csv")); reg <- fread(file.path(D, "survey_registry_mics.csv"))
rmap <- fread(file.path(D, "region_map.csv")); inv <- fread("results/mics_inventory/survey_inventory.csv")
built <- fread(file.path(D, "cbh_build/survey_manifest.csv"))[status %in% c("built", "cached")]$survey   # cached = reused shard
setup <- setup[svkey %in% built]
setup[, year := reg$year[match(svkey, reg$svkey)]]; setup[, iso3 := substr(folder, 1, 3)]
card_only <- c("MC_GIN2016", "MC_COM2022", "MC_TCD2019")   # recall doses unusable: WUENIC fallback (user decision)
# Primary and total secondary durations (years) by country, used to turn level + grade into years.
dur <- fread(text = "iso3,P,S
BEN,6,7\nCAF,6,7\nCIV,6,7\nCMR,6,7\nCOD,6,6\nCOG,6,7\nCOM,6,7\nGHA,6,6\nGIN,6,7\nGMB,6,6\nGNB,6,6\nKEN,8,4
LSO,7,5\nMDG,5,7\nMLI,6,6\nMOZ,7,5\nMRT,6,7\nMWI,8,4\nNGA,6,6\nSEN,6,7\nSLE,6,6\nSOM,8,4\nSSD,8,4\nSTP,6,6
SWZ,7,5\nTCD,6,7\nTGO,6,7\nZWE,7,6")

num <- function(x) as.numeric(zap_missing(zap_labels(x)))
getv <- function(d, ...) { low <- tolower(names(d)); for (a in tolower(c(...))) { j <- match(a, low); if (!is.na(j)) return(d[[j]]) }; NULL }
labs_of <- function(x) { l <- attr(x, "labels"); if (is.null(l)) return(setNames(character(), character())); setNames(tolower(names(l)), l) }
wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (!any(ok)) return(c(NA, 0, 0)); c(sum(x[ok] * w[ok]) / sum(w[ok]), sum(ok), sum(w[ok])) }
region_of <- function(d, s) {
  if (startsWith(s$region_var, "CONST:")) return(rep(sub("CONST:", "", s$region_var), nrow(d)))
  v <- sub("^WM:", "", s$region_var); x <- getv(d, v, "HH7", "hh7", "HHB", "HHREG", "WMREG", "region", "zone")
  if (is.null(x)) return(rep(NA_character_, nrow(d)))
  # Labels can differ in case or accents between files of the same survey, so match on normalised keys.
  m <- rmap[svkey == s$svkey]; m$analysis_region[match(cbh_rkey(as.character(as_factor(x))), cbh_rkey(m$mics_label))]
}
# Education: classify the level label, then years = offset + grade (capped), with fallbacks.
edu_years <- function(level, grade, iso) {
  P <- dur[iso3 == iso]$P; S <- dur[iso3 == iso]$S; LS <- ceiling(S / 2)
  L <- labs_of(level); lv <- num(level); g <- num(grade); g[!is.finite(g) | g > 20] <- NA
  txt <- ifelse(is.na(lv), NA, L[as.character(lv)])
  cat <- rep(NA_character_, length(lv))
  pat <- list(none = "none|aucun|nenhum|n[ãa]o|never|pre.?school|pre.?primary|pr[ée].?scolaire|pr[ée].?escolar|maternelle|kindergarten|ece|eccde|early childhood|alfabetiz|koranic|coranique|non.?formal|non.?standard|adult educ",
    higher = "higher|tertiary|university|sup[ée]rieur|superior|bacharelato|post.?secondary|diploma|college|polytechnic",
    voc = "voc|tech|profission|commercial|nursing|teacher|vei|iei|foundation certificate",
    upper = "upper sec|senior sec|sss|shs|secondaire 2|lyc[ée]e|esg2|2nd cycle|second cycle|high$",
    lower = "lower sec|junior sec|jss|jhs|middle|intermediate|secondaire 1|coll[èe]ge|esg1|1er cycle|post.?primary",
    secondary = "second|secund|secondaire|high",
    primary = "primary|primaire|prim[áa]rio|basic|b[áa]sico|ep1|ep2|fondamental")
  for (k in names(pat)) { hit <- is.na(cat) & !is.na(txt) & grepl(pat[[k]], txt); cat[hit] <- k }
  cat[!is.na(txt) & grepl("missing|manquant|dk|nsp|ns$|no response|non r[ée]ponse|incoherent|sem re", txt) & is.na(cat)] <- "missing"
  yrs <- rep(NA_real_, length(lv))
  cap <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  yrs[cat == "none"] <- 0
  i <- cat == "primary" & !is.na(cat); yrs[i] <- ifelse(is.na(g[i]), P / 2, cap(g[i], 0, P))
  i <- cat == "lower" & !is.na(cat); yrs[i] <- P + ifelse(is.na(g[i]), LS / 2, ifelse(g[i] > LS, cap(g[i] - P, 0, LS), cap(g[i], 0, LS)))
  i <- cat == "upper" & !is.na(cat); yrs[i] <- P + LS + ifelse(is.na(g[i]), (S - LS) / 2, ifelse(g[i] > S - LS, cap(g[i] - P - LS, 0, S - LS), cap(g[i], 0, S - LS)))
  i <- cat == "secondary" & !is.na(cat); yrs[i] <- P + ifelse(is.na(g[i]), S / 2, ifelse(g[i] > S, cap(g[i] - P, 0, S), cap(g[i], 0, S)))
  i <- cat == "voc" & !is.na(cat); yrs[i] <- P + S / 2
  i <- cat == "higher" & !is.na(cat); yrs[i] <- P + S + ifelse(is.na(g[i]), 2, cap(g[i], 0, 6))
  list(years = yrs, category = cat)
}
card_dose <- function(v) { x <- num(v); ifelse(is.na(x), NA, ifelse((x >= 1 & x <= 31) | x %in% c(44, 66), 1, ifelse(x == 0, 0, NA))) }
classify_labels <- function(x, improved, unimproved) {
  L <- labs_of(x); v <- num(x); t <- L[as.character(v)]
  out <- rep(NA_real_, length(v)); out[grepl(unimproved, t)] <- 0; out[is.na(out) & grepl(improved, t)] <- 1; out
}
water_code <- function(x) {  # JMP 2017: piped, borehole, protected well/spring, rain, tanker, cart, kiosk, bottled, sachet
  v <- num(x); imp <- v %in% c(11:14, 21, 31, 41, 51, 61, 71, 72, 91, 92); unimp <- v %in% c(32, 42, 81, 96)
  byl <- classify_labels(x, "pip|tap|robinet|torneira|borehole|forage|furo|tube|protected|prot[ée]g|protegid|rain|pluie|chuva|bottle|bouteille|garraf|sachet|tank|citerne|kiosk|fontaine|standpipe|chafariz",
                         "unprotected|non prot|n[ãa]o protegid|surface|river|rivi|rio|lake|lac|pond|dam|stream|canal|marigot|ruisseau|other|autre|outro")
  ifelse(imp, 1, ifelse(unimp, 0, byl))
}
sanit_code <- function(x) {  # JMP 2017: flush to sewer/septic/pit/unknown, VIP, pit with slab, composting
  v <- num(x); imp <- v %in% c(11, 12, 13, 15, 18, 21, 22, 31); unimp <- v %in% c(14, 23, 41, 51, 95, 96)
  byl <- classify_labels(x, "sewer|[ée]gout|septic|septique|fossa|pit latrine with slab|avec dalle|com laje|ventilated|vip|am[ée]lior|composting|compost|flush to pit",
                         "without slab|sans dalle|sem laje|open pit|traditional|tradicional|bucket|seau|balde|hanging|suspend|no facility|bush|field|nature|brousse|mato|elsewhere|ailleurs|other|autre")
  ifelse(imp, 1, ifelse(unimp, 0, byl))
}
facility_code <- function(x) {
  classify_labels(x, "hospital|h[ôo]pital|hospit|health|sant[ée]|sa[úu]de|clinic|clinique|cl[íi]nica|cent|post|poste|posto|dispens|matern|facility|cs |csps|cm |pmi|private med|m[ée]dical",
                  "home|domicile|maison|casa|resid|tbA|traditional|on the way|en route|route|other|autre|outro")
}

long <- list(); nat <- list()
for (i in seq_len(nrow(setup))) {
  s <- setup[i]; f <- function(x) file.path(root, s$folder, x); iso <- s$iso3
  bh <- read_sav(f("bh.sav")); wm <- read_sav(f("wm.sav")); hh <- read_sav(f("hh.sav")); ch <- read_sav(f("ch.sav"))
  add <- function(var, reg, x, w, pop) {
    d <- data.table(reg, x, w)[!is.na(reg)]
    r <- d[, { m <- wmean(x, w); .(value = m[1], observed_n = m[2], weighted_n = m[3], eligible_n = .N) }, by = reg]
    long[[length(long) + 1]] <<- data.table(survey = s$svkey, country = iso, regkey = cbh_rkey(r$reg), region = r$reg,
      variable = var, value = r$value, observed_n = r$observed_n, weighted_n = r$weighted_n, eligible_n = r$eligible_n,
      population = pop, source = "MICS_microdata_weighted_summary")
    tot <- wmean(d$x, d$w); nat[[length(nat) + 1]] <<- data.table(survey = s$svkey, country = iso, variable = var, national = tot[1], n = tot[2])
  }
  # ---- mothers with a birth in the 60 months before interview (one row per woman) -----
  use_ref <- is.null(getv(bh, "HH1", "WM1"))
  key <- function(d) if (use_ref) paste(num(getv(d, "XX1")), 0, num(getv(d, "MEMID", "LN"))) else
    paste(num(getv(d, "HH1", "WM1")), num(if (is.null(getv(d, "HH2", "WM2"))) 0 else getv(d, "HH2", "WM2")), num(getv(d, "LN", "WM4", "WM3")))
  bdob <- getv(bh, "BH4C", "CCDOB"); bdob <- if (!is.null(bdob)) num(bdob) else {
    mo <- num(getv(bh, "BH4MC", "BH4M", "HN5_MES")); yr <- num(getv(bh, "BH4YC", "BH4Y", "HN5_ANO")); yr <- ifelse(yr < 100, yr + 1900, yr)
    ifelse(mo %in% 1:12 & yr > 1900 & yr < 2030, (yr - 1900) * 12 + mo, NA) }
  bkey <- key(bh); wkey <- key(wm)
  wdoi <- num(getv(wm, "WDOI", "WM6C", "CMCDOIW")); if (!length(wdoi) || all(is.na(wdoi))) { mo <- num(getv(wm, "WM6M")); yr <- num(getv(wm, "WM6Y")); wdoi <- (yr - 1900) * 12 + mo }
  wdob <- num(getv(wm, "WDOB", "WM8C")); ww <- num(getv(wm, "wmweight", "WMWEIGHT"))
  first <- tapply(bdob, bkey, min, na.rm = TRUE); first[!is.finite(first)] <- NA
  recent <- tapply(seq_along(bkey), bkey, function(j) any(is.finite(bdob[j]))) & TRUE
  lastb <- tapply(bdob, bkey, max, na.rm = TRUE)
  mother <- wkey %in% names(lastb) & is.finite(wdoi) & (wdoi - lastb[wkey]) >= 0 & (wdoi - lastb[wkey]) < 60
  mother[is.na(mother)] <- FALSE
  wreg <- region_of(wm, s)
  afb <- floor((first[wkey] - wdob) / 12); afb[!(afb >= 8 & afb <= 49)] <- NA
  add("mean_maternal_age_first_birth", wreg[mother], afb[mother], ww[mother], "mothers_birth_last_5y")
  lvl <- getv(wm, "WB6A", "WB4", "WB4X", "WM11", "wm11", "MELEVEL"); grd <- getv(wm, "WB6B", "WB5", "WB5X", "WM12", "wm12")
  if (!is.null(lvl)) { e <- edu_years(lvl, if (is.null(grd)) rep(NA, nrow(wm)) else grd, iso)
    # Women who never attended school skip the level question in MICS; DHS counts them as 0 years.
    wmm <- data.table(variable = names(wm), label = vapply(wm, function(x) { l <- attr(x, "label"); if (is.null(l)) "" else l }, ""))
    ev <- wmm[grepl("ever attended|jamais fr[ée]quent|d[ée]j[àa] fr[ée]quent|alguma vez frequent|ever been to school|attended school", label, ignore.case = TRUE)]$variable
    if (length(ev)) { a <- num(getv(wm, ev[1])); e$years[a == 2 & is.na(e$years)] <- 0; e$category[a == 2 & is.na(e$category)] <- "none" }
    wl <- getv(wm, "welevel", "melevel")
    if (!is.null(wl)) { tl <- labs_of(wl)[as.character(num(wl))]
      nz <- grepl("^none$|pre.?primary or none|^aucun|^nenhum|^sans|no education|none or pre", tl) & is.na(e$years)
      e$years[nz] <- 0; e$category[nz] <- "none" }
    add("mean_maternal_education_years", wreg[mother], e$years[mother], ww[mother], "mothers_birth_last_5y")
    ed_cat <- table(e$category[mother], useNA = "ifany") }
  wq <- num(getv(wm, "windex5", "WLTHIND5")); add("mean_wealth_quintile", wreg[mother], ifelse(wq %in% 1:5, wq, NA)[mother], ww[mother], "mothers_birth_last_5y")
  hh6 <- getv(wm, "HH6")
  if (is.null(hh6) && !is.null(getv(hh, "HH6")) && !use_ref) {
    hk <- paste(num(getv(hh, "HH1")), num(getv(hh, "HH2"))); wk2 <- paste(num(getv(wm, "HH1", "WM1")), num(getv(wm, "HH2", "WM2")))
    hh6 <- getv(hh, "HH6")[match(wk2, hk)] }
  # Urban by label: some surveys add categories such as "slum (informal settlement)", which are urban.
  urb <- if (is.null(hh6)) rep(if (grepl("Dakar", s$folder)) 100 else NA, nrow(wm)) else {
    L6 <- labs_of(hh6); t6 <- L6[as.character(num(hh6))]
    ifelse(grepl("urban|urbain|urbano|slum|informal|bidonville|city|ville|cidade", t6), 100, ifelse(grepl("rural", t6), 0,
      ifelse(num(hh6) == 1, 100, ifelse(num(hh6) == 2, 0, NA)))) }
  add("urban_pct", wreg[mother], urb[mother], ww[mother], "mothers_birth_last_5y")
  # ---- facility delivery: last birth in the 2 years before interview (MICS module) ----
  dv <- inv[survey == s$folder]$delivery_var
  if (length(dv) && !is.na(dv) && nzchar(dv) && !is.null(getv(wm, dv))) {
    fx <- facility_code(getv(wm, dv)); add("facility_delivery_pct", wreg, 100 * fx, ww, "last_birth_2y") }
  # ---- short birth interval: non-first births in last 60 months, preceding interval 7-23 m --
  ord <- order(bkey, bdob); pb <- rep(NA_real_, length(bdob))
  pb[ord] <- ave(bdob[ord], bkey[ord], FUN = function(z) c(NA, head(z, -1)))
  bdoi <- wdoi[match(bkey, wkey)]; bw <- ww[match(bkey, wkey)]; breg <- wreg[match(bkey, wkey)]
  intv <- bdob - pb; recentb <- is.finite(bdoi - bdob) & (bdoi - bdob) >= 0 & (bdoi - bdob) < 60 & is.finite(intv) & intv >= 7
  add("short_birth_interval_pct", breg[recentb], 100 * (intv[recentb] <= 23), bw[recentb], "nonfirst_births_last_5y")
  # ---- households ------------------------------------------------------------------
  hw <- num(getv(hh, "hhweight", "HHWEIGHT")); hreg <- region_of(hh, s)
  if (!is.null(getv(hh, "HH46", "HH9")) && "survey" %in% tolower(names(hh))) { }   # (NGA 2021 NICS rows handled below)
  keep <- rep(TRUE, nrow(hh)); sv <- getv(hh, "survey"); if (!is.null(sv)) { lab <- tolower(as.character(as_factor(sv))); if (any(grepl("mics", lab))) keep <- grepl("mics", lab) }
  hwm <- getv(hh, "hhweightMICS"); if (!is.null(hwm)) hw <- num(hwm)
  wv <- getv(hh, "WS1"); if (!is.null(wv)) add("improved_water_pct", hreg[keep], 100 * water_code(wv)[keep], hw[keep], "households")
  hm <- data.table(variable = names(hh), label = vapply(hh, function(x) { l <- attr(x, "label"); if (is.null(l)) "" else l }, ""))
  tcand <- hm[grepl("kind of toilet|type of toilet|toilet facility|type de toilette|lieux d.?aisance|type de latrine|tipo de (casa de banho|sanit|retrete)|instala..o sanit", label, ignore.case = TRUE) &
              !grepl("shared|partag|compartilh|households using|number of|nombre de|n.mero de|location|emplacement|localiza", label, ignore.case = TRUE)]$variable
  tcand <- c(tcand[grepl("^WS", tcand, ignore.case = TRUE)], tcand[!grepl("^WS", tcand, ignore.case = TRUE)])   # reported (WS) before observed items
  tv <- if (length(tcand)) getv(hh, tcand[1]) else getv(hh, inv[survey == s$folder]$toilet_var); if (!is.null(tv)) add("improved_sanitation_pct", hreg[keep], 100 * sanit_code(tv)[keep], hw[keep], "households")
  ev <- getv(hh, inv[survey == s$folder]$electricity_var, "HC8", "HC8A", "HC9A")
  if (!is.null(ev)) { e <- num(ev); L <- labs_of(ev)
    yes <- if (any(grepl("off.?grid|generator|g[ée]n[ée]rat|solar|solaire", L))) e %in% as.numeric(names(L)[grepl("yes|oui|sim|grid|r[ée]seau|off|generator|solar|solaire", L) & !grepl("^no|non$|n[ãa]o", L)]) else e == 1
    elec <- ifelse(is.na(e) | e >= 7, NA, 100 * yes); add("electricity_pct", hreg[keep], elec[keep], hw[keep], "households") }
  # ---- children: vaccination 12-23 months and anthropometry 0-59 months -------------
  cw <- num(getv(ch, "chweight", "CHWEIGHT")); cage <- num(getv(ch, "CAGE", "cage")); creg <- region_of(ch, s)
  ckeep <- rep(TRUE, nrow(ch)); sv <- getv(ch, "survey"); if (!is.null(sv)) { lab <- tolower(as.character(as_factor(sv))); if (any(grepl("mics", lab))) ckeep <- grepl("mics", lab) }
  cwm <- getv(ch, "chweightMICS"); if (!is.null(cwm)) cw <- num(cwm)
  a1223 <- is.finite(cage) & cage >= 12 & cage <= 23 & ckeep
  if (!s$svkey %in% card_only) {
    chm <- data.table(variable = names(ch), label = vapply(ch, function(x) { l <- attr(x, "label"); if (is.null(l)) "" else l }, ""))
    d3 <- inv[survey == s$folder]$dtp3_var
    rc <- chm[grepl("(times|fois|vezes|number of|nombre de|quantas).*(dpt|dtp|dtc|penta|diph)|(dpt|dtp|dtc|penta|diph).*(times|fois|vezes)", label, ignore.case = TRUE)]$variable
    card <- if (length(d3) && !is.na(d3) && !is.null(getv(ch, d3))) card_dose(getv(ch, d3)) else rep(NA, nrow(ch))
    rec3 <- if (length(rc)) { r <- num(getv(ch, rc[length(rc)])); ifelse(r >= 3 & r <= 7, 1, ifelse(r >= 0 & r < 3, 0, NA)) } else rep(NA, nrow(ch))
    dtp <- ifelse(card %in% 1, 1, ifelse(rec3 %in% 1, 1, ifelse(card %in% 0 | rec3 %in% 0, 0, NA)))
    add("dtp3_pct", creg[a1223], 100 * dtp[a1223], cw[a1223], "children_12_23m")
    mday <- chm[grepl("(measles|rougeole|sarampo|\\bmr\\b|\\brr\\b|\\bvar\\b|mcv)", label, ignore.case = TRUE) & grepl("(^|[^a-z])(day|jour|dia)([^a-z]|$)", label, ignore.case = TRUE)]$variable
    if (!length(mday)) mday <- chm[grepl("^IM[0-9]*(M|MEAS|MR|N|VAR|M1|N1)D$", variable, ignore.case = TRUE)]$variable
    mrec <- chm[grepl("(measles|rougeole|sarampo|\\bmr\\b|\\brr\\b|\\bvar\\b)", label, ignore.case = TRUE) & grepl("ever|d[ée]j[àa]|j[áa]|received|given|reçu|recebeu", label, ignore.case = TRUE) & !grepl("times|fois|vezes", label, ignore.case = TRUE)]$variable
    mc <- if (length(mday)) card_dose(getv(ch, mday[1])) else rep(NA, nrow(ch))
    mr <- if (length(mrec)) { r <- num(getv(ch, mrec[1])); ifelse(r == 1, 1, ifelse(r == 2, 0, NA)) } else rep(NA, nrow(ch))
    ms <- ifelse(mc %in% 1 | mr %in% 1, 1, ifelse(mc %in% 0 | mr %in% 0, 0, NA))
    add("measles_pct", creg[a1223], 100 * ms[a1223], cw[a1223], "children_12_23m")
  }
  zs <- ch; zkey <- NULL
  if (is.null(getv(ch, "HAZ2", "haz2", "HAZ21", "zhaz")) && file.exists(f("who_z.sav"))) {       # MICS3: WHO z-scores in who_z.sav
    z <- read_sav(f("who_z.sav")); zk <- paste(num(getv(z, "HH1")), num(getv(z, "HH2")), num(getv(z, "LN")))
    ck <- paste(num(getv(ch, "HH1", "UF1")), num(getv(ch, "HH2", "UF2")), num(getv(ch, "LN", "UF4"))); j <- match(ck, zk)
    haz <- num(getv(z, "HAZ2", "haz2"))[j]; whz <- num(getv(z, "WHZ2", "whz2"))[j]
    hf <- num(getv(z, "hazflag", "HAZFLAG"))[j]; wf <- num(getv(z, "whzflag", "WHZFLAG"))[j]
  } else { haz <- num(getv(ch, "HAZ2", "haz2", "HAZ21", "zhaz")); whz <- num(getv(ch, "WHZ2", "whz2", "WHZ21", "zwhz"))
    hf <- num(getv(ch, "HAZFLAG", "hazflag", "HAZFLAG1")); wf <- num(getv(ch, "WHZFLAG", "whzflag", "WHZFLAG1"))
    if (!length(hf)) hf <- rep(NA, nrow(ch)); if (!length(wf)) wf <- rep(NA, nrow(ch)) }
  if (length(haz) && any(is.finite(haz))) {
    a059 <- is.finite(cage) & cage >= 0 & cage <= 59 & ckeep
    okh <- a059 & is.finite(haz) & abs(haz) <= 6 & !(hf %in% 1); okw <- a059 & is.finite(whz) & abs(whz) <= 5 & !(wf %in% 1)
    add("stunting_pct", creg[okh], 100 * (haz[okh] < -2), cw[okh], "children_0_59m_valid")
    add("wasting_pct", creg[okw], 100 * (whz[okw] < -2), cw[okw], "children_0_59m_valid")
  }
  message(sprintf("%-15s done", s$svkey))
}
L <- rbindlist(long); N <- rbindlist(nat)
vars <- c("mean_maternal_age_first_birth","mean_maternal_education_years","mean_wealth_quintile","urban_pct","dtp3_pct",
          "measles_pct","facility_delivery_pct","short_birth_interval_pct","improved_water_pct","improved_sanitation_pct",
          "electricity_pct","wasting_pct","stunting_pct")
# Every analysis region of every built survey gets a row, including regions with no observations.
frame <- unique(rmap[svkey %in% setup$svkey, .(survey = svkey, region = analysis_region)])[, regkey := cbh_rkey(region)]
W <- dcast(L, survey + regkey ~ variable, value.var = "value")
W <- merge(frame[, .(survey, regkey)], W, by = c("survey", "regkey"), all.x = TRUE)
for (v in vars) if (!v %in% names(W)) W[, (v) := NA_real_]
W[, country := setup$iso3[match(survey, setup$svkey)]]; W[, survey_year := setup$year[match(survey, setup$svkey)]]
# Fallback 1: WUENIC national DTP3 / MCV1 for the survey year (from the local UNICEF snapshot).
u <- fread("data/fusion_GLOBAL_DATAFLOW_UNICEF_1.0_all.csv", select = c("REF_AREA","INDICATOR","SEX","TIME_PERIOD","OBS_VALUE","UNIT_MEASURE","OBS_STATUS","DATA_SOURCE","AGE"))
u <- u[SEX == "_T" & UNIT_MEASURE == "PCNT" & AGE == "M12T23" & OBS_STATUS == "E" & grepl("WHO/UNICEF estimates of national immunization coverage", DATA_SOURCE, fixed = TRUE) &
       grepl("IM_DTP3|IM_MCV1", INDICATOR)]
u[, ind := ifelse(grepl("DTP3", INDICATOR), "dtp3_pct", "measles_pct")]
u <- u[, .(value = as.numeric(OBS_VALUE[1])), by = .(iso3 = sub(":.*", "", REF_AREA), year = as.integer(TIME_PERIOD), ind)]
for (v in c("dtp3_pct", "measles_pct")) {
  W[, paste0(v, "_before_imputation") := get(v)]
  if (v %in% c("dtp3_pct", "measles_pct")) W[survey %in% card_only, (v) := NA_real_]
  nv <- u[ind == v][match(paste(W$country, W$survey_year), paste(iso3, year))]$value
  miss <- !is.finite(W[[v]])
  W[, paste0(v, "_imputation_source") := ifelse(miss & is.finite(nv), "UNICEF_WUENIC_country_year", ifelse(miss, "missing", "retained_regional_estimate"))]
  W[miss, (v) := nv[miss]]
}
# Fallback 2: within-survey arithmetic mean of available regions (each region counted once).
for (v in setdiff(vars, c())) {
  W[, paste0(v, "_regional_mean_imputed") := FALSE]
  W[, (paste0(v, "_regional_mean_imputed")) := !is.finite(get(v)) & any(is.finite(get(v))), by = survey]
  W[, (v) := ifelse(is.finite(get(v)), get(v), mean(get(v)[is.finite(get(v))])), by = survey]
  W[!is.finite(get(v)), (v) := NA_real_]
}
fwrite(W, file.path(D, "regional_covariates_wide_mics.csv"))
fwrite(L, file.path(D, "regional_covariates_long_mics.csv"))
fwrite(N, file.path(aud, "national_values_by_survey.csv"))
comp <- W[, lapply(.SD, function(x) round(100 * mean(is.finite(x)))), by = survey, .SDcols = vars]
fwrite(comp, file.path(aud, "completeness_by_survey.csv"))
cat("regions:", nrow(W), "| surveys:", uniqueN(W$survey), "\n")
print(W[, lapply(.SD, function(x) sum(!is.finite(x))), .SDcols = vars])

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
# Vaccine columns the label search misses, or where one antigen is split over two columns (audit 2026-09-28).
vacc_extra <- list(
  MC_MRT2011 = list(dtp_card = "IM3PE3D", dtp_times = "IM10B", meas_rec = "IM16"), # Penta columns; IM16 label mis-decoded (Cyrillic)
  MC_COG2014 = list(dtp_card = "IM3PT3D"),                                         # pentavalent 3 card column (recall: DTC IM12 and Penta IM15D)
  MC_BEN2021 = list(meas_card = "IM6RRD", meas_rec = "IM26B"),                     # measles-rubella (RR)
  MC_GHA2017 = list(meas_rec = "IM26A"),                                           # "MR1" defeats \bmr\b
  MC_GNB2018 = list(meas_rec = "IM26"),                                            # "VAS"
  MC_MDG2018 = list(meas_rec = "IM26"),                                            # "ROR"
  MC_TGO2017 = list(meas_rec = "IM26"),                                            # "ROR"
  MC_MOZ2008 = list(dtp_times = "IM16", meas_card = "IM6D"),                       # labels name no vaccine ('Quantas vezes recebeu?', 'Sarampo')
  MC_ZWE2009 = list(dtp_card = "im5cd"))                                           # DPT-HepB3 row of the MICS3 card (as MWI2006 IM5CD); im4cd is 99 when the dose is recorded there
# Surveys outside v7 exempt from the 98% explicit-evidence check (their children with no information still count as 0).
vacc_evidence_exempt <- c(MC_SSD2010 = "after excluding children outside the AG2Y 0-1 module filter, about 6% of the rest have a blank module (non-response, counted as 0)",
                          MC_SOM2011NE = "item nonresponse: 'ever vaccinated' / 'ever DPT' coded 9 (missing) for about 5% of children 12-23 months")
# Primary and total secondary durations (years) by country, used to turn level + grade into years.
dur <- fread(text = "iso3,P,S
BEN,6,7\nCAF,6,7\nCIV,6,7\nCMR,6,7\nCOD,6,6\nCOG,6,7\nCOM,6,7\nGHA,6,6\nGIN,6,7\nGMB,6,6\nGNB,6,6\nKEN,8,4
LSO,7,5\nMDG,5,7\nMLI,6,6\nMOZ,7,5\nMRT,6,7\nMWI,8,4\nNGA,6,6\nSEN,6,7\nSLE,6,6\nSOM,8,4\nSSD,8,4\nSTP,6,6
SWZ,7,5\nTCD,6,7\nTGO,6,7\nZWE,7,6")
emap <- fread("R_mics/education_level_map.csv", encoding = "UTF-8")   # per-survey level -> years map (every built MICS survey)
stopifnot(!anyDuplicated(emap[, .(svkey, code)]), all(xor(is.na(emap$flat), is.na(emap$offset))), all(is.na(emap$flat) == !is.na(emap$cap)),
          all(is.na(emap$flat) == !is.na(emap$base)), all(is.na(emap$flat) == !is.na(emap$na_add)))

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
# Education (generic fallback for surveys not in education_level_map.csv): classify the level label, then years = offset + grade (capped).
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
# Education from the reviewed per-survey map: graded level = offset + clamp(grade - base, 0, cap); level without a grade =
# offset + na_add; ungraded level (none, pre-primary, Koranic, Mahadra, non-formal, literacy, out-of-sequence vocational) = flat.
edu_years_map <- function(level, grade, sv) {
  M <- emap[svkey == sv]; lv <- num(level); g <- num(grade); g[!is.finite(g) | g >= 90] <- NA
  j <- match(lv, M$code); txt <- labs_of(level)[as.character(lv)]
  miss <- is.na(j) & !is.na(lv) & (lv %in% c(7, 8, 9, 97, 98, 99) | grepl("missing|manqu|^dk$|nsp|^ns$|no response|non r[ée]ponse|incoh|sem re", txt))
  bad <- is.na(j) & !is.na(lv) & !miss
  if (any(bad)) stop(sv, ": level codes not in education_level_map.csv: ", paste(sort(unique(lv[bad])), collapse = ", "))
  flat <- M$flat[j]; rel <- g - M$base[j]
  # A grade beyond the level's length but inside its cumulative span is counted from school entry (MOZ2008 superior '12').
  cum <- !is.na(j) & is.na(flat) & !is.na(g) & rel > M$cap[j] & g >= M$offset[j] & g <= M$offset[j] + M$cap[j]
  yrs <- ifelse(is.na(j), NA_real_, ifelse(!is.na(flat), flat, ifelse(cum, g,
           M$offset[j] + ifelse(is.na(g), M$na_add[j], pmin(pmax(rel, 0), M$cap[j])))))
  list(years = yrs, category = ifelse(is.na(j), ifelse(miss, "missing", NA_character_), M$label[j]), row = j,
       over = !is.na(j) & is.na(flat) & !is.na(g) & !cum & ((rel < 0 & g != 0) | rel > M$cap[j]),   # grade 0 = no grade completed at that level
       graded = !is.na(j) & is.na(flat) & !is.na(rel) & rel >= 1 & rel <= M$cap[j])
}
card_dose <- function(v) { x <- num(v); ifelse(is.na(x), NA, ifelse((x >= 1 & x <= 31) | x %in% c(44, 66), 1, ifelse(x == 0, 0, NA))) }
rec3_code <- function(v) { r <- num(v); ifelse(r >= 3 & r <= 7, 1, ifelse(r >= 0 & r < 3, 0, NA)) }
yesno_code <- function(v) { r <- num(v); ifelse(r == 1, 1, ifelse(r == 2, 0, NA)) }
cols_of <- function(d, vs, f) { vs <- vs[!is.na(vs)]; if (!length(vs)) return(matrix(NA_real_, nrow(d), 1))
  matrix(vapply(vs, function(v) f(getv(d, v)), numeric(nrow(d))), nrow = nrow(d)) }
any_dose <- function(M) { y <- rowSums(matrix(M %in% 1, nrow(M))) > 0; n <- rowSums(matrix(M %in% 0, nrow(M))) > 0; ifelse(y, 1, ifelse(n, 0, NA_real_)) }
classify_labels <- function(x, improved, unimproved) {
  L <- labs_of(x); v <- num(x); t <- L[as.character(v)]
  out <- rep(NA_real_, length(v)); out[grepl(unimproved, t)] <- 0; out[is.na(out) & grepl(improved, t)] <- 1; out
}
water_code <- function(x, sv = "") {  # JMP 2017: piped, borehole, protected well/spring, rain, tanker, cart, kiosk, bottled, sachet
  v <- num(x); imp <- v %in% c(11:14, 21, 31, 41, 51, 61, 71, 72, 91, 92); unimp <- v %in% c(32, 42, 81, 96)
  byl <- classify_labels(x, "pip|tap|robinet|torneira|borehole|forage|furo|tube|protected|prot[ée]g|protegid|rain|pluie|chuva|bottle|bouteille|garraf|sachet|tank|citerne|kiosk|fontaine|standpipe|chafariz",
                         "unprotected|non prot|n[ãa]o protegid|surface|river|rivi|rio|lake|lac|pond|dam|stream|canal|marigot|ruisseau|other|autre|outro")
  out <- ifelse(imp, 1, ifelse(unimp, 0, byl))
  if (sv == "MC_MRT2011") out[v %in% c(33, 34)] <- 0   # puits traditionnel couvert / non couvert: unlined wells (JMP needs lining + cover; 33 = user decision)
  if (sv == "MC_COD2017") out[v == 62] <- 1            # bidon/bassin/seau livré à domicile = delivered water (JMP 2017 improved)
  if (sv == "MC_MOZ2008") out[v == 32] <- 1            # 'Sem bomba manual' continues 31 'poço ou furo protegido': protected, no pump
  out
}
sanit_code <- function(x, sv = "") {  # JMP 2017: flush to sewer/septic/pit/unknown, VIP, pit with slab, composting
  v <- num(x); imp <- v %in% c(11, 12, 13, 15, 18, 21, 22, 31); unimp <- v %in% c(14, 23, 41, 51, 95, 96)
  byl <- classify_labels(x, "sewer|[ée]gout|septic|septique|fossa|pit latrine with slab|avec dalle|com laje|ventilated|vip|am[ée]lior|composting|compost|flush to pit|w\\. slab",
                         "without slab|sans dalle|sem laje|open pit|traditional|tradicional|bucket|seau|balde|hanging|suspend|no facility|bush|field|nature|brousse|mato|elsewhere|ailleurs|other|autre")
  out <- ifelse(imp, 1, ifelse(unimp, 0, byl)); t <- labs_of(x)[as.character(v)]
  out[grepl("endroit inconnu|lieu inconnu|unknown place|dk where|don.t know where|ne sait pas o|desconhecido", t)] <- 1   # flush to unknown place = improved (JMP 2017); MRT2011 code 14
  if (sv == "MC_GHA2017") out[v == 24] <- NA_real_   # 'pit latrine with seat' (slab unknown): out of the denominator, as is 61 'mobile toilet' (user decision: status quo)
  if (sv == "MC_MDG2018") out[v == 24] <- 1          # 'latrine à fosse avec dalle non lavable': a slab, improved (user decision: status quo)
  out
}
# Health facility as in DHS RH_DELP_C_DHF: public (21-29), private medical (31-39), NGO/mission (41-49), sector unknown (76) = 1;
# home (11-19), en route, other (81-96) = 0; 97-99 missing. Codes are used only when 11/12 are homes (standard MICS layout);
# otherwise (MC_MOZ2008, MC_SSD2010) labels decide, as before; 'phc' adds SSD2010 code 2 'PHCF' (primary health care facility).
# MC_COD2017/MC_TCD2019 labels look shifted by one from code 33 (no 33 or 76; 36 'MATERNITE PRIVEE', 96 'AUTRE PRIVE MEDICAL',
# no plain 'Autre'): codes are trusted, so 96 = other = 0, as in v7 (read as a facility, COD2017 would gain about 1 point).
fac_home <- "home|domicile|maison|casa|resid"
facility_code <- function(x) {
  v <- num(x); L <- labs_of(x); t <- L[as.character(v)]
  std <- any(names(L) %in% c("11", "12")) && all(grepl(fac_home, L[names(L) %in% c("11", "12")]))
  out <- if (std) ifelse(v >= 11 & v <= 19, 0, ifelse((v >= 21 & v <= 49) | v == 76, 1, ifelse(v >= 81 & v <= 96, 0, NA_real_))) else
    classify_labels(x, "hospital|h[ôo]pital|hospit|health|sant[ée]|sa[úu]de|clinic|clinique|cl[íi]nica|cent|post|poste|posto|dispens|matern|facility|cs |csps|cm |pmi|private med|m[ée]dical|phc",
                    "home|domicile|maison|casa|resid|tba|traditional|on the way|en route|route|other|autre|outro")
  out[grepl("on (the )?way|en route|roadside|vehicle|in the open|a caminho", t)] <- 0
  out
}
# Stop when more than tol% of answered codes are unclassified, or (facility) a label contradicts its class.
check_coding <- function(sv, what, x, cls, w, allow_na = numeric(), tol = 0.5) {
  v <- num(x); t <- labs_of(x)[as.character(v)]; ok <- is.finite(v) & !(v %in% 97:99) & !grepl("^(missing|manquant|sem info|em falta|no response|non r[ée]ponse)", t) & is.finite(w) & w > 0
  un <- ok & is.na(cls) & !(v %in% allow_na); unc <- 100 * sum(w[un]) / sum(w[ok])
  if (unc > tol) stop(sprintf("%s %s: %.2f%% of answers unclassified (codes %s)", sv, what, unc, paste(sort(unique(v[un])), collapse = ",")))
  if (what == "facility") {
    bad <- ok & ((cls %in% 1 & grepl(fac_home, t) & !grepl("matern|nursing", t)) |
                 (cls %in% 0 & grepl("hospital|h[ôoe]pital|clinic|clinique|cl[íi]nica|health|sant[ée]|sa[úu]de|dispens|matern|mission|cham", t) &
                  !grepl("on (the )?way|en route|^autre|^other|^outro|autre priv", t)))
    if (any(bad)) stop(sprintf("%s facility: label/class conflict: %s", sv, paste(unique(t[bad]), collapse = " | ")))
  }
  invisible(unc)
}
# Weighted share of interviewed households with a coded value.
obs_share <- function(sv, what, cls, w, done, min_share = 0.95) { ok <- done & is.finite(w) & w > 0; sh <- sum(w[ok & is.finite(cls)]) / sum(w[ok])
  if (sh < min_share) stop(sprintf("%s: %s observed for %.1f%% of interviewed households", sv, what, 100 * sh)); sh }

long <- list(); nat <- list(); vacc_diag <- list(); wb7_done <- character()
for (i in seq_len(nrow(setup))) {
  s <- setup[i]; f <- function(x) file.path(root, s$folder, x); iso <- s$iso3; rnd <- inv[survey == s$folder]$round
  stopifnot(length(rnd) == 1, rnd %in% 2:6)
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
  lvn <- c("WB6A", "WB4", "WB4X", "WM11", "wm11", "MELEVEL"); lvn <- lvn[tolower(lvn) %in% tolower(names(wm))][1]
  lvl <- if (is.na(lvn)) NULL else getv(wm, lvn); grd <- getv(wm, "WB6B", "WB5", "WB5X", "WM12", "wm12")
  if (!is.null(lvl)) { g0 <- if (is.null(grd)) rep(NA, nrow(wm)) else grd; mapped <- s$svkey %in% emap$svkey
    e <- if (mapped) edu_years_map(lvl, g0, s$svkey) else { warning(s$svkey, ": not in education_level_map.csv, generic conversion"); edu_years(lvl, g0, iso) }
    # Never attended: the item immediately before the level question (MICS6 WB5, MICS4/5 WB3, MICS3/4 WM10), chosen by
    # position, not label text (labels differ by language/encoding: GNB2018 'Alguma vez já frequentou', COD2017 'Déjé fréquenté').
    evn <- unname(c(WB6A = "WB5", WB4 = "WB3", WB4X = "WB3", WM11 = "WM10")[toupper(lvn)]); evx <- if (is.na(evn)) NULL else getv(wm, evn)
    a <- rep(NA_real_, nrow(wm))
    if (!is.null(evx)) { if (!grepl("^(no|non|n[ãa]o|nã£o)$", trimws(labs_of(evx)["2"]))) stop(s$svkey, ": ", evn, " code 2 is not 'no'"); a <- num(evx) }
    e$years[a %in% 2 & is.na(e$years)] <- 0; e$category[a %in% 2 & is.na(e$category)] <- "none"
    wl <- getv(wm, "welevel", "melevel")
    if (!is.null(wl)) { tl <- labs_of(wl)[as.character(num(wl))]
      nz <- grepl("(^|[^a-z])(none|aucun|nenhum|nunca)([^a-z]|$)|no education|sans instruction", tl) & is.na(e$years) & !(a %in% 1)
      e$years[nz] <- 0; e$category[nz] <- "none" }
    # MICS6 WB6B is the highest grade ATTENDED; WB7 asks whether it was completed. DHS v133 counts completed years, so a
    # graded level not completed loses one year (user decision). MICS3-5 WB7 is a literacy item and must not match.
    w7 <- getv(wm, "WB7"); w7c <- !is.null(w7) && isTRUE(grepl("compl[eéè]t|ach[eéè]v|conclu", attr(w7, "label"), ignore.case = TRUE))
    if (w7c && !(rnd >= 6 && identical(toupper(lvn), "WB6A"))) stop(s$svkey, ": WB7 reads as a completion item outside the MICS6 layout")
    if (rnd >= 6 && identical(toupper(lvn), "WB6A") && !is.null(w7) && !w7c) stop(s$svkey, ": MICS6 WB7 is not a grade-completion item: ", attr(w7, "label"))
    if (w7c) { nc <- num(w7) %in% 2 & (if (mapped) e$graded else { gg <- if (is.null(grd)) rep(NA_real_, nrow(wm)) else num(grd); is.finite(gg) & gg >= 1 & gg <= 20 &
        e$category %in% c("primary", "lower", "upper", "secondary", "higher") & is.finite(e$years) & e$years > 0 })
      e$years[nc] <- e$years[nc] - 1; wb7_done <- c(wb7_done, s$svkey) }
    mw <- ww[mother]; na_sh <- sum(mw[is.na(e$years[mother])], na.rm = TRUE) / sum(mw, na.rm = TRUE)
    if (mapped) {
      if (na_sh > 0.02) stop(s$svkey, ": ", round(100 * na_sh, 1), "% of mothers lack education years")
      ov <- data.table(j = e$row[mother], o = e$over[mother])[!is.na(j), .(n = .N, sh = mean(o)), by = j][n >= 20 & sh > 0.10]
      if (nrow(ov)) stop(s$svkey, ": grades outside mapped range for level codes ", paste(emap[svkey == s$svkey]$code[ov$j], collapse = ", "))
    } else if (na_sh > 0.02) warning(s$svkey, ": ", round(100 * na_sh, 1), "% of mothers lack education years (generic conversion)")
    if (any(mother & a %in% 2) && mean(e$years[mother & a %in% 2] == 0, na.rm = TRUE) < 0.99) stop(s$svkey, ": never-attended mothers with non-zero years")
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
  if (s$svkey == "MC_MOZ2008") {
    # MICS3 Mozambique splits the place of delivery: MN7_A intended place, MN7_B delivered there (1) or elsewhere (6), MN8 the
    # actual place, asked only when MN7_B = 6. Place = MN7_A if MN7_B = 1, MN8 if MN7_B = 6; denominator = women answering MN7_B.
    pa <- getv(wm, "MN7_A"); p8 <- getv(wm, "MN8"); b <- num(getv(wm, "MN7_B"))
    stopifnot(!is.null(pa), !is.null(p8), identical(unname(attr(pa, "labels")), unname(attr(p8, "labels"))))
    pl <- labelled(ifelse(b %in% 1, num(pa), ifelse(b %in% 6, num(p8), NA_real_)), attr(pa, "labels")); mn <- !is.na(b)
    fx <- facility_code(pl); check_coding(s$svkey, "facility", pl[mn], fx[mn], ww[mn], allow_na = c(15, 23))   # 15 brigadas móveis, 23 farmácia: excluded, as before
    add("facility_delivery_pct", wreg[mn], 100 * fx[mn], ww[mn], "last_birth_2y")
  } else if (length(dv) && !is.na(dv) && nzchar(dv) && !is.null(getv(wm, dv))) {
    fx <- facility_code(getv(wm, dv)); check_coding(s$svkey, "facility", getv(wm, dv), fx, ww); add("facility_delivery_pct", wreg, 100 * fx, ww, "last_birth_2y") }
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
  # Interviewed households: result HH46 (MICS6) / HH9 (MICS3-5) = 1 when so labelled (KEN2009MOM weights every household); else weight > 0.
  hr <- getv(hh, if (rnd >= 6) "HH46" else "HH9")
  hdone <- if (!is.null(hr) && isTRUE(grepl("complet|rempli|accept|preench", labs_of(hr)["1"]))) num(hr) %in% 1 else is.finite(hw) & hw > 0
  wv <- getv(hh, "WS1")
  if (!is.null(wv)) { wc <- water_code(wv, s$svkey); check_coding(s$svkey, "water", wv[keep], wc[keep], hw[keep], allow_na = if (s$svkey == "MC_SOM2006") 52:54 else numeric())   # SOM2006 roof top / berkad / balli: unclassified, as in v7
    obs_share(s$svkey, "water", wc[keep], hw[keep], hdone[keep], min_share = if (s$svkey == "MC_SOM2006") 0.8 else 0.95)
    add("improved_water_pct", hreg[keep], 100 * wc[keep], hw[keep], "households") }
  hm <- data.table(variable = names(hh), label = vapply(hh, function(x) { l <- attr(x, "label"); if (is.null(l)) "" else l }, ""))
  tcand <- hm[grepl("kind of toilet|type of toilet|toilet facility|type de toilette|lieux d.?aisance|type de latrine|tipo de (casa de banho|sanit|retrete)|instala..o sanit", label, ignore.case = TRUE) &
              !grepl("shared|partag|compartilh|households using|number of|nombre de|n.mero de|location|emplacement|localiza", label, ignore.case = TRUE)]$variable
  tcand <- c(tcand[grepl("^WS", tcand, ignore.case = TRUE)], tcand[!grepl("^WS", tcand, ignore.case = TRUE)])   # reported (WS) before observed items
  tv <- if (length(tcand)) getv(hh, tcand[1]) else getv(hh, inv[survey == s$folder]$toilet_var)
  if (!is.null(tv)) { sc <- sanit_code(tv, s$svkey)
    # MC_MOZ2008 asks WS7 only if the household has (WS6A) or can use (WS6B) a toilet; the others report where they defecate
    # (WS9A: beach, bush, other) = no facility = unimproved (JMP; the DHS denominator is all households).
    if (s$svkey == "MC_MOZ2008") sc[is.na(num(tv)) & ((num(getv(hh, "WS6A")) %in% 2 & num(getv(hh, "WS6B")) %in% 2) | is.finite(num(getv(hh, "WS9A"))))] <- 0
    check_coding(s$svkey, "sanitation", tv[keep], sc[keep], hw[keep], allow_na = if (s$svkey == "MC_GHA2017") c(24, 61) else if (s$svkey == "MC_MDG2012S") 24:25 else numeric())   # MDG2012S 'dales' slab latrines: unclassified, as in v7
    obs_share(s$svkey, "sanitation", sc[keep], hw[keep], hdone[keep])
    add("improved_sanitation_pct", hreg[keep], 100 * sc[keep], hw[keep], "households") }
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
    # DHS/MICS tabulation: denominator = every living child 12-23 months with a completed under-5 interview; vaccinated only on
    # positive evidence (card date, 44 marked on card, 66 mother reported, recall DTP count 3+, measles 'yes'); anything else = 0.
    chm <- data.table(variable = names(ch), label = vapply(ch, function(x) { l <- attr(x, "label"); if (is.null(l)) "" else l }, ""))
    rx <- getv(ch, if (rnd >= 6) "UF17" else "UF9"); if (is.null(rx)) stop(s$svkey, ": no under-5 interview result item")
    if (!isTRUE(grepl("complet|rempli|accept|preench", labs_of(rx)["1"]))) stop(s$svkey, ": under-5 interview result code 1 is not 'completed'")
    cdone <- a1223 & num(rx) %in% 1
    # Card seen by round, not by label ('viu', 'vi'; GNB2014 IM5 is a campaign item): MICS6 IM5 1-3 (card and/or other document), MICS3-5 IM1 = 1.
    seen_var <- if (rnd >= 6) "IM5" else "IM1"; sx <- getv(ch, seen_var); if (is.null(sx)) stop(s$svkey, ": no card-seen item ", seen_var)
    seen <- num(sx) %in% (if (rnd >= 6) 1:3 else 1)
    if (rnd >= 6 && s$svkey != "MC_COD2017" && any(cdone & num(sx) %in% 5)) stop(s$svkey, ": IM5 = 5 outside MC_COD2017 (only COD2017's 'card kept at the health centre' is excluded)")
    x <- vacc_extra[[s$svkey]]
    bad <- setdiff(toupper(unlist(x)), toupper(names(ch))); if (length(bad)) stop(s$svkey, ": vacc_extra column(s) not in ch.sav: ", paste(bad, collapse = ", "))
    d3 <- unique(c(inv[survey == s$folder]$dtp3_var, x$dtp_card)); d3 <- d3[!is.na(d3) & vapply(d3, function(v) !is.null(getv(ch, v)), TRUE)]
    rc <- chm[grepl("(times|fois|vezes|number of|nombre de|quantas).*(dpt|dtp|dtc|penta|diph)|(dpt|dtp|dtc|penta|diph).*(times|fois|vezes)", label, ignore.case = TRUE)]$variable
    card <- any_dose(cols_of(ch, d3, card_dose))
    if (length(d3)) { R3 <- cols_of(ch, d3, num); dated <- a1223 & rowSums(matrix((R3 >= 1 & R3 <= 31) | R3 %in% 44, nrow(ch)), na.rm = TRUE) > 0
      if (sum(dated) >= 20 && mean(seen[dated]) < 0.99) stop(sprintf("%s: %.1f%% of card-dated DTP3 doses not flagged seen by %s", s$svkey, 100 * (1 - mean(seen[dated])), seen_var))
      card[seen & !(card %in% 1)] <- 0 }   # dose not recorded on a seen card (blank, 97/98/99) = not given, as in the MICS and DHS tabulations
    rcu <- unique(c(rc, x$dtp_times))   # every DTP/DTCoq and pentavalent count (user decision); only COG2014 gains children (IM12 DTC); MOZ2008 IM16 by name
    rec3 <- any_dose(cols_of(ch, rcu, rec3_code))
    # Every DTP-containing dose-3 card column must be used: IM3D3D, IM6PENTA3D, ... and the MICS3 rows IM4CD/IM5CD (ZWE2009 im5cd).
    d3_all <- chm[grepl("^IM[0-9]+([A-Z]*3|C)D$", variable, ignore.case = TRUE) & grepl("dpt|dtp|dtc|penta|diph|t[ée]tracoq", paste(variable, label), ignore.case = TRUE)]$variable
    if (length(setdiff(toupper(d3_all), toupper(d3)))) stop(s$svkey, ": DTP-containing dose-3 card column(s) not used: ", paste(setdiff(d3_all, d3), collapse = ", "))
    if (!length(rcu)) stop(s$svkey, ": no DTP 'number of times' recall column")
    mday <- chm[grepl("(measles|rougeole|sarampo|\\bmr\\b|\\brr\\b|\\bvar\\b|mcv)", label, ignore.case = TRUE) & grepl("(^|[^a-z])(day|jour|dia)([^a-z]|$)", label, ignore.case = TRUE)]$variable
    if (!length(mday)) mday <- chm[grepl("^IM[0-9]*(M|MEAS|MR|N|VAR|M1|N1)D$", variable, ignore.case = TRUE)]$variable
    mrec <- chm[grepl("(measles|rougeole|sarampo|\\bmr\\b|\\brr\\b|\\bvar\\b)", label, ignore.case = TRUE) & grepl("ever|d[ée]j[àa]|j[áa]|received|given|reçu|recebeu", label, ignore.case = TRUE) & !grepl("times|fois|vezes", label, ignore.case = TRUE)]$variable
    mcard <- unique(c(if (length(mday)) mday[1], x$meas_card))   # first measles-containing dose only; never the M2D/ND/RR2 columns
    mrecu <- unique(c(if (length(mrec)) mrec[1], x$meas_rec))
    mc <- any_dose(cols_of(ch, mcard, card_dose)); mr <- any_dose(cols_of(ch, mrecu, yesno_code))
    if (length(mcard)) mc[seen & !(mc %in% 1)] <- 0   # only when a measles card column exists
    if (!length(mrecu)) stop(s$svkey, ": no measles 'ever received' recall column")
    if (!is.null(getv(ch, "IM6MEASD")) && any(a1223 & card_dose(getv(ch, "IM6MEASD")) %in% 1 & !mc %in% 1)) stop(s$svkey, ": IM6MEASD shows measles doses the card columns miss")
    # 'Ever' items, by questionnaire structure: DTP/Penta ever = the yes/no item paired with each DTP count by name (IM12A <- IM11A,
    # IM21 <- IM20, IM16 <- IM15, IM10B <- IM10A, IM15D <- IM15C; names interleave, e.g. BEN2014 IM11, IM11A, IM12, IM12A, and
    # SOM2011 has a DPT 'free ORS packet' item IM11A before IM12); any vaccination ever = IM11 (MICS6), IM10 (MICS3 layout, DPT3
    # card column IM4CD/IM5CD), IM6 (MICS4-5, including STP2014) or IM0A (NGA2016); every non-card-only survey has one.
    yn <- function(v) { L <- labs_of(getv(ch, v)); isTRUE(grepl("^(yes|oui|sim)", L["1"]) & grepl("^n", L["2"])) }
    ever_of <- function(v) { m <- regmatches(toupper(v), regexec("^IM([0-9]+)([A-Z]?)$", toupper(v)))[[1]]; if (!length(m)) return(NA_character_)
      n <- as.integer(m[2]); cand <- c(paste0("IM", n - 1, m[3]), if (m[3] %in% LETTERS[-1]) paste0("IM", n, LETTERS[match(m[3], LETTERS) - 1]))
      cand <- cand[vapply(cand, function(v) !is.null(getv(ch, v)) && yn(v) && grepl("dpt|dtp|dtc|penta|diph|t[ée]tracoq|injec", chm$label[match(v, toupper(chm$variable))], ignore.case = TRUE), TRUE)]
      if (length(cand)) cand[1] else NA_character_ }
    dtp_ever <- unname(vapply(rcu, ever_of, "")); if (anyNA(dtp_ever)) stop(s$svkey, ": no DTP 'ever received' item paired with ", paste(rcu[is.na(dtp_ever)], collapse = ", "))
    av <- c(if (rnd >= 6) "IM11" else if (grepl("^IM[45]CD$", toupper(d3[1]))) "IM10" else "IM6", "IM0A"); av <- av[vapply(av, function(v) !is.null(getv(ch, v)), TRUE)][1]
    if (!is.na(av) && !(yn(av) && grepl("vac|immun", chm$label[match(toupper(av), toupper(chm$variable))], ignore.case = TRUE))) stop(s$svkey, ": ", av, " is not an 'ever vaccinated' item")
    anyx <- if (is.na(av)) rep(NA_real_, nrow(ch)) else num(getv(ch, av))
    E <- cols_of(ch, dtp_ever, num); ever_no <- rowSums(E == 1, na.rm = TRUE) == 0 & rowSums(E == 2, na.rm = TRUE) > 0   # MRT2011: IM11 asked only if IM10A is not yes
    dk <- function(vs) rowSums(matrix(cols_of(ch, vs, num) %in% c(8, 98), nrow(ch))) > 0
    # Explicit evidence first (card or recall 0/1, 'never' answers, don't know = 0); the survey must explain at least 98% of children.
    dtp <- ifelse(card %in% 1 | rec3 %in% 1, 1, ifelse(card %in% 0 | rec3 %in% 0 | anyx %in% c(2, 8) | ever_no | dk(c(dtp_ever, rcu)), 0, NA))
    ms <- ifelse(mc %in% 1 | mr %in% 1, 1, ifelse(mc %in% 0 | mr %in% 0 | anyx %in% c(2, 8) | dk(mrecu), 0, NA))
    if (length(d3) && any(cdone & seen & is.na(dtp))) stop(s$svkey, ": children with a seen card but no DTP3 value")
    # Children the questionnaire did not ask are excluded, not counted as unvaccinated (user decisions): MC_COD2017 card kept at the
    # health centre (IM5 = 5, recall skipped); MC_SSD2010 children outside the module's age filter (reported age AG2Y not 0-1)
    # with an entirely blank module (28 September 2026; filter-failing children who were asked keep their answers).
    k5 <- if (s$svkey == "MC_COD2017") { if (!isTRUE(grepl("centre de sant", labs_of(sx)["5"]))) stop("MC_COD2017: IM5 code 5 is not 'card kept at the health centre'")
      num(sx) %in% 5 } else if (s$svkey == "MC_SSD2010") { imv <- grep("^IM", names(ch), value = TRUE, ignore.case = TRUE)
      !(num(getv(ch, "AG2Y")) %in% 0:1) & rowSums(vapply(imv, function(v) !is.na(num(getv(ch, v))), logical(nrow(ch)))) == 0 } else rep(FALSE, nrow(ch))
    okc <- cdone & !is.na(creg) & !k5 & is.finite(cw) & cw > 0
    shd <- sum(cw[okc & is.finite(dtp)]) / sum(cw[okc]); shm <- sum(cw[okc & is.finite(ms)]) / sum(cw[okc])
    if (min(shd, shm) < 0.98 && !s$svkey %in% names(vacc_evidence_exempt)) stop(sprintf("%s: explicit vaccination evidence for only %.1f%% (DTP3) / %.1f%% (measles) of children 12-23 months", s$svkey, 100 * shd, 100 * shm))
    vacc_diag[[length(vacc_diag) + 1]] <- data.table(survey = s$svkey, n_12_23 = sum(okc), dtp_evidence_share = shd, measles_evidence_share = shm,
      k5_excluded = sum(k5 & cdone & !is.na(creg)), exempt = s$svkey %in% names(vacc_evidence_exempt), dtp_card = paste(d3, collapse = "+"), dtp_recall = paste(rcu, collapse = "+"), dtp_ever = paste(dtp_ever, collapse = "+"),
      any_ever = ifelse(is.na(av), "", av), meas_card = paste(mcard, collapse = "+"), meas_recall = paste(mrecu, collapse = "+"))
    dtp[cdone & !k5 & is.na(dtp)] <- 0; ms[cdone & !k5 & is.na(ms)] <- 0   # no information at all = not vaccinated
    dtp[!cdone | k5] <- NA; ms[!cdone | k5] <- NA
    add("dtp3_pct", creg[a1223], 100 * dtp[a1223], cw[a1223], "children_12_23m")
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
L <- rbindlist(long); N <- rbindlist(nat); VD <- rbindlist(vacc_diag)
# Regression lock on the v9 national values (28 September 2026) of the surveys whose vaccine columns were changed; not an
# independent validation. The audit's values excluded don't-know and no-information children (57.1, 71.9, 67.2 Penta recall
# only, 54.3, 55.1, 73.4, 87.2, 77.1); here those children count as 0 and COG2014 also uses DTC recall. ZWE2009: WUENIC 73.
exp <- data.table(survey = c("MC_MRT2011", "MC_MRT2011", "MC_COG2014", "MC_BEN2021", "MC_MDG2018", "MC_TGO2017", "MC_GHA2017", "MC_GNB2018", "MC_MOZ2008", "MC_MOZ2008", "MC_ZWE2009"),
  variable = c("dtp3_pct", "measles_pct", "dtp3_pct", "measles_pct", "measles_pct", "measles_pct", "measles_pct", "measles_pct", "dtp3_pct", "measles_pct", "dtp3_pct"),
  value = c(55.65, 70.45, 65.27, 53.55, 54.67, 73.27, 86.49, 75.54, 73.31, 73.73, 67.28))
chk <- merge(exp, N, by = c("survey", "variable")); stopifnot(nrow(chk) == nrow(exp), all(abs(chk$national - chk$value) < 0.15))
# The WB7 completed-grade adjustment ran exactly in the MICS6 surveys with the WB6A/WB6B/WB7 layout.
w7x <- setup$svkey[inv$round[match(setup$folder, inv$survey)] >= 6]
if (!setequal(wb7_done, w7x)) stop("WB7 adjustment ran in ", paste(setdiff(wb7_done, w7x), collapse = ", "), " / not in ", paste(setdiff(w7x, wb7_done), collapse = ", "))
stopifnot(all(setdiff(setup$svkey, card_only) %in% VD$survey), VD[!(exempt)]$dtp_evidence_share >= 0.98, VD[!(exempt)]$measles_evidence_share >= 0.98)
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
fwrite(VD, file.path(aud, "vaccination_evidence_by_survey.csv"))
comp <- W[, lapply(.SD, function(x) round(100 * mean(is.finite(x)))), by = survey, .SDcols = vars]
fwrite(comp, file.path(aud, "completeness_by_survey.csv"))
cat("regions:", nrow(W), "| surveys:", uniqueN(W$survey), "| WB7 completed-grade adjustment:", length(wb7_done), "surveys\n")
print(W[, lapply(.SD, function(x) sum(!is.finite(x))), .SDcols = vars])

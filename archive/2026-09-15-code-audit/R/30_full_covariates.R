# =============================================================================
# 30_full_covariates.R — assemble a FULL predictor set for neonatal & post-
# neonatal mortality, region-level (DHS, admin-1) where possible + national
# (World Bank) otherwise, then impute missing. Extends the region-health panel
# (R/23) with the user's requested predictors.
#
# DHS BR (design-weighted by v005, region = best_region_var vs panel):
#   measles      % children 12-23mo who received measles (h9)
#   excl_bf      % children <6mo exclusively breastfed (v404 & no other food/liquid 24h)
#   stunting     % children <5y HAZ<-2 (hw5)     underweight WAZ<-2 (hw8)   wasting WHZ<-2 (hw11)
#   birth_int    % births with preceding interval <24mo (b11)
#   mage1        mean maternal age at first birth (v212)     [+ carried: dtp3_reg,facility,educ_yrs,wealth_q]
# World Bank (national, nearest survey year):
#   hexp_gdp %GDP (SH.XPD.CHEX.GD.ZS), hexp_pc per-capita USD (SH.XPD.CHEX.PC.CD),
#   polstab political stability (PV.EST), electricity access (EG.ELC.ACCS.ZS)
# Missing -> imputed (country-median, then overall-median) with *_imp flags kept.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(rdhs); library(mgcv)}); options(rappdir_permission = TRUE)
CACHE <- path.expand("~/.rdhs_cache/datasets_reformatted")
XDIR  <- file.path(DATA, "panel_cache_extra"); dir.create(XDIR, showWarnings = FALSE)

## ---- World Bank pulls (cached) ---------------------------------------------
wb_get <- function(ind, nm) { f <- file.path(DATA, paste0("wb_", nm, ".csv"))
  if (file.exists(f)) return(read.csv(f, stringsAsFactors = FALSE))
  x <- tryCatch(wb_fetch(ind), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) { cat("  WB FAIL:", ind, "\n"); return(NULL) }
  x <- x[, c("iso3","year","value")]; names(x)[3] <- nm; write.csv(x, f, row.names = FALSE); x }
cat("[WB] fetching health expenditure / political stability / electricity ...\n")
wb <- list(hexp_gdp = wb_get("SH.XPD.CHEX.GD.ZS","hexp_gdp"), hexp_pc = wb_get("SH.XPD.CHEX.PC.CD","hexp_pc"),
           polstab  = wb_get("PV.EST","polstab"),             elec    = wb_get("EG.ELC.ACCS.ZS","elec"))
for (nm in names(wb)) cat(sprintf("    %-9s %s\n", nm, if (is.null(wb[[nm]])) "-- unavailable --" else sprintf("%d rows", nrow(wb[[nm]]))))

## ---- enumerate surveys + panel region keys ---------------------------------
ds <- dhs_datasets(fileFormat = "FL"); brs <- ds[ds$FileType == "Births Recode", ]
brs$iso3 <- countrycode::countrycode(brs$CountryName, "country.name", "iso3c", warn = FALSE)
brs$reg  <- countrycode::countrycode(brs$iso3, "iso3c", "region", warn = FALSE)
brs$yr   <- as.integer(brs$SurveyYear)
sv <- brs[!is.na(brs$reg) & brs$reg == "Sub-Saharan Africa" & brs$yr >= 2000 &
            brs$SurveyType %in% c("DHS","MIS","AIS"), ]
sv$svkey <- paste0(substr(sv$FileName,1,2), substr(sv$FileName,5,8)); sv$noext <- sub("\\..*$","",sv$FileName)
pan <- read.csv(file.path(RESULTS, "component2_region_data_expanded_health.csv"), stringsAsFactors = FALSE)
pk <- split(pan$regkey, pan$svkey)

reg_wmean <- function(val, w, reg, ok) { k <- ok & is.finite(val) & is.finite(w) & !is.na(reg) & nzchar(reg)
  if (!any(k)) return(numeric(0)); tapply(w[k]*val[k], reg[k], sum) / tapply(w[k], reg[k], sum) }
yn <- function(x) grepl("yes|received|reported|card|marked", tolower(as.character(x)))     # DHS "given/received"

extra_rows <- function(i) {
  s <- sv[i,]; cf <- file.path(XDIR, paste0(s$svkey, ".rds")); tk <- pk[[s$svkey]]
  if (file.exists(cf)) { cc <- readRDS(cf); if (is.null(tk) || mean(cc$regkey %in% tk) >= 0.6) return(cc) }
  f <- file.path(CACHE, paste0(s$noext, ".rds")); if (!file.exists(f)) return(NULL)
  b <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(b) || !"v024" %in% names(b)) return(NULL)
  rv <- if (length(tk)) best_region_var(b, tk) else "v024"
  w <- as.numeric(b$v005)/1e6; reg <- rkey(as.character(b[[rv]]))
  agem <- suppressWarnings(as.numeric(b$v008) - as.numeric(b$b3)); alive <- tolower(as.character(b$b5)) == "yes"
  R <- function(val, ok) reg_wmean(val, w, reg, ok)
  # measles (12-23mo)
  measles <- if ("h9" %in% names(b)) 100*R(as.numeric(yn(b$h9)), alive & is.finite(agem) & agem>=12 & agem<=23 & !is.na(b$h9) & tolower(as.character(b$h9))!="missing") else numeric(0)
  # exclusive breastfeeding (<6mo): currently BF and nothing else in last 24h
  excl_bf <- numeric(0)
  if ("v404" %in% names(b)) { u6 <- alive & is.finite(agem) & agem < 6
    bf <- grepl("yes|breast", tolower(as.character(b$v404)))
    other_vars <- intersect(c("v409","v410","v411","v411a","v412","v412a","v413","v414a","v414b","v414c","v414e","v414f","v414g","v414h","v414i","v414j","v414k","v414l","v414m","v414n","v414o","v414p","v414v"), names(b))
    gave <- Reduce(`|`, lapply(other_vars, function(v) yn(b[[v]])), FALSE)
    excl_bf <- 100*R(as.numeric(bf & !gave), u6 & !is.na(b$v404)) }
  # anthropometry (<5y): z-scores x100, flag >=9990 missing
  anth <- function(col, thr = -200) { if (!col %in% names(b)) return(numeric(0))
    z <- suppressWarnings(as.numeric(b[[col]])); z[z >= 9990 | z <= -600] <- NA
    100*R(as.numeric(z < thr), alive & is.finite(z)) }
  stunting <- anth("hw5"); underweight <- anth("hw8"); wasting <- anth("hw11")
  # short preceding birth interval (<24mo) and maternal age at first birth
  bint <- if ("b11" %in% names(b)) { bi <- suppressWarnings(as.numeric(b$b11)); 100*R(as.numeric(bi < 24), is.finite(bi)) } else numeric(0)
  mid <- if ("caseid" %in% names(b)) as.character(b$caseid) else paste(b$v001,b$v002,b$v003)
  fm  <- !duplicated(mid); fh <- !duplicated(paste(b$v001, b$v002))      # mother- / household-level dedup
  v212 <- suppressWarnings(as.numeric(b$v212)); v212[v212 > 45 | v212 < 8] <- NA
  mage1 <- if ("v212" %in% names(b)) R(v212, fm & is.finite(v212)) else numeric(0)
  # household water / sanitation / electricity (region %, household-deduped)
  hh_frac <- function(col, pat, excl=NULL) { if (!col %in% names(b)) return(numeric(0))
    lab <- tolower(as.character(b[[col]])); ok <- fh & !is.na(b[[col]]) & nzchar(lab) & lab != "missing"
    imp <- grepl(pat, lab) & (if (is.null(excl)) TRUE else !grepl(excl, lab)); 100*R(as.numeric(imp), ok) }
  imp_water <- hh_frac("v113", "pipe|tap|standpipe|borehole|tube ?well|protected|rain|bottled|sachet", "unprotected")
  imp_sanit <- hh_frac("v116", "flush|septic|sewer|ventilated|vip|slab|composting", "without slab|open pit|no facil|bush|field|hanging|bucket|somewhere")
  elec_dhs  <- hh_frac("v119", "yes")
  regs <- Reduce(union, list(names(measles),names(excl_bf),names(stunting),names(underweight),names(wasting),
                             names(bint),names(mage1),names(imp_water),names(imp_sanit),names(elec_dhs)))
  if (!length(regs)) return(NULL)
  pk2 <- function(v) if (length(v)) unname(v[regs]) else NA_real_
  r <- data.frame(svkey=s$svkey, regkey=regs, measles=pk2(measles), excl_bf=pk2(excl_bf),
                  stunting=pk2(stunting), underweight=pk2(underweight), wasting=pk2(wasting),
                  birth_int=pk2(bint), mage1=pk2(mage1), imp_water=pk2(imp_water),
                  imp_sanit=pk2(imp_sanit), elec_dhs=pk2(elec_dhs), stringsAsFactors = FALSE)
  saveRDS(r, cf); r
}
cat(sprintf("[DHS] building extra region covariates for %d surveys ...\n", nrow(sv)))
X <- do.call(rbind, lapply(seq_len(nrow(sv)), function(i)
  tryCatch(extra_rows(i), error = function(e) { message(sv$svkey[i], ": ", conditionMessage(e)); NULL })))
cat(sprintf("[DHS] extra region rows: %d from %d surveys\n", nrow(X), length(unique(X$svkey))))

## ---- merge everything onto the panel ---------------------------------------
d <- merge(pan, X, by = c("svkey","regkey"), all.x = TRUE)
nrst <- function(tab, iso, yr) { if (is.null(tab)) return(NA_real_); s <- tab[tab$iso3==iso & is.finite(tab[[3]]),]
  if (!nrow(s)) NA_real_ else s[[3]][which.min(abs(s$year - yr))] }
for (nm in names(wb)) d[[nm]] <- mapply(function(is,yr) nrst(wb[[nm]], is, yr), d$iso3, d$year)
d$log_hexp_pc <- log(d$hexp_pc)

## ---- coverage report + imputation (country-median -> overall-median) -------
vars <- c("dtp3_reg","measles","facility","educ_yrs","wealth_q","excl_bf","stunting","underweight",
          "wasting","birth_int","mage1","imp_water","imp_sanit","elec_dhs","hexp_gdp","log_hexp_pc","elec")
cat("\n=== covariate completeness (of", nrow(d), "region-years) ===\n")
for (v in vars) cat(sprintf("  %-12s %3.0f%%  mean %8.2f\n", v, 100*mean(is.finite(d[[v]])), mean(d[[v]], na.rm=TRUE)))
impute <- function(x, grp) { m <- ave(x, grp, FUN = function(z) median(z, na.rm=TRUE))
  x2 <- ifelse(is.finite(x), x, m); ifelse(is.finite(x2), x2, median(x, na.rm=TRUE)) }
for (v in vars) { d[[paste0(v,"_imp")]] <- !is.finite(d[[v]]); d[[v]] <- impute(d[[v]], d$iso3) }
write.csv(d, file.path(RESULTS, "component2_region_data_full.csv"), row.names = FALSE)
cat("\nsaved: results/component2_region_data_full.csv\n")

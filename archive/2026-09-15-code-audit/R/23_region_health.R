# =============================================================================
# 23_region_health.R — REGION-LEVEL health-access covariates from DHS microdata,
# to replace the national WB DTP3 (which is constant within each survey and so
# cannot adjust for SUBNATIONAL health-system variation). For every BR recode in
# the expanded panel we compute, per admin-1 region (v024, design-weighted v005):
#   dtp3_reg  = % of living children 12-23 mo who received DPT3 (card or recall)
#   facility  = % of recent births delivered in a health facility (m15)
#   educ_yrs  = mean single-year maternal education (v133), one row per mother
#   wealth_q  = mean household wealth quintile (v190 -> 1..5), one row per mother
# Then merge into results/component2_region_data_expanded.csv, build a composite
# health-access index (PC1), and refit the primary model with vs without these
# region covariates to see how the malaria (PfPR) effect moves.
#
# Reads the rdhs reformatted cache directly (no re-download). Resumable via
# data/panel_cache_health/.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(mgcv)})
CACHE <- path.expand("~/.rdhs_cache/datasets_reformatted")
HDIR  <- file.path(DATA, "panel_cache_health"); dir.create(HDIR, showWarnings = FALSE)

## ---- enumerate the same surveys as R/18 (BR, SSA, 2000+) --------------------
suppressMessages(library(rdhs)); options(rappdir_permission = TRUE)
ds <- dhs_datasets(fileFormat = "FL"); br <- ds[ds$FileType == "Births Recode", ]
br$iso3 <- countrycode::countrycode(br$CountryName, "country.name", "iso3c", warn = FALSE)
br$reg  <- countrycode::countrycode(br$iso3, "iso3c", "region", warn = FALSE)
br$yr   <- as.integer(br$SurveyYear)
sv <- br[!is.na(br$reg) & br$reg == "Sub-Saharan Africa" & br$yr >= 2000 &
           br$SurveyType %in% c("DHS","MIS","AIS"), ]
sv$svkey <- paste0(substr(sv$FileName, 1, 2), substr(sv$FileName, 5, 8))
sv$noext <- sub("\\..*$", "", sv$FileName)

# region keys per survey from the (recovered) panel — used to detect the recode
# region variable that matches the panel's granularity (e.g. sreg1 for Tanzania).
.pan <- read.csv(file.path(RESULTS, "component2_region_data_expanded.csv"), stringsAsFactors = FALSE)
pk <- split(.pan$regkey, .pan$svkey)

## ---- region-weighted mean helper (returns named vector keyed by regkey) -----
reg_wmean <- function(val, w, reg, ok) {
  k <- ok & is.finite(val) & is.finite(w) & !is.na(reg) & nzchar(reg)
  if (!any(k)) return(numeric(0))
  tapply(w[k] * val[k], reg[k], sum) / tapply(w[k], reg[k], sum)
}
FAC <- "hospital|clinic|health|dispensary|maternit|matern|cms|pmi|hut|post|doctor|nursing|centre|center|sector|infirmary|polyclin|hopital|hôpital|clinique|sante|santé|cabinet"
HOME <- "home|house|domicile|maison|parent"

health_rows <- function(i) {
  s <- sv[i, ]; cf <- file.path(HDIR, paste0(s$svkey, ".rds")); tk <- pk[[s$svkey]]
  if (file.exists(cf)) { cc <- readRDS(cf)                         # reuse cache unless stale vs panel keys
    if (is.null(tk) || mean(cc$regkey %in% tk) >= 0.6) return(cc)  # (old v024 Tanzania rows -> recompute)
  }
  f <- file.path(CACHE, paste0(s$noext, ".rds")); if (!file.exists(f)) return(NULL)
  b <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(b) || !"v024" %in% names(b)) return(NULL)
  rv  <- if (length(tk)) best_region_var(b, tk) else "v024"        # region var matching the panel granularity
  w   <- as.numeric(b$v005) / 1e6
  reg <- rkey(as.character(b[[rv]]))
  # DPT3 among living children 12-23 months
  dpt3 <- numeric(0)
  if ("h7" %in% names(b) && "b3" %in% names(b)) {
    agem  <- as.numeric(b$v008) - as.numeric(b$b3)
    alive <- tolower(as.character(b$b5)) == "yes"
    h7    <- tolower(as.character(b$h7))
    elig  <- alive & is.finite(agem) & agem >= 12 & agem <= 23 & !is.na(h7) & h7 != "missing"
    recv  <- grepl("card|mother|report|marked", h7)
    dpt3  <- 100 * reg_wmean(as.numeric(recv), w, reg, elig)
  }
  # facility delivery (births with a usable m15)
  facility <- numeric(0)
  if ("m15" %in% names(b)) {
    m15 <- tolower(as.character(b$m15)); mok <- !is.na(m15) & m15 != "missing" & nzchar(m15)
    fac <- as.numeric(!grepl(HOME, m15) & grepl(FAC, m15))
    facility <- 100 * reg_wmean(fac, w, reg, mok)
  }
  # mother-level covariates: dedupe to one row per mother
  mid <- if ("caseid" %in% names(b)) as.character(b$caseid) else paste(b$v001, b$v002, b$v003)
  fm  <- !duplicated(mid)
  edu <- suppressWarnings(as.numeric(as.character(b$v133))); edu[edu > 25 | edu < 0] <- NA
  educ_yrs <- reg_wmean(edu, w, reg, fm & is.finite(edu))
  wq  <- match(tolower(as.character(b$v190)), c("poorest","poorer","middle","richer","richest"))
  wealth_q <- reg_wmean(as.numeric(wq), w, reg, fm & is.finite(wq))
  regs <- Reduce(union, list(names(dpt3), names(facility), names(educ_yrs), names(wealth_q)))
  if (!length(regs)) return(NULL)
  pick <- function(v) if (length(v)) unname(v[regs]) else NA_real_
  r <- data.frame(svkey = s$svkey, regkey = regs, dtp3_reg = pick(dpt3), facility = pick(facility),
                  educ_yrs = pick(educ_yrs), wealth_q = pick(wealth_q), stringsAsFactors = FALSE)
  saveRDS(r, cf); r
}

cat(sprintf("[1] building region health covariates for %d surveys ...\n", nrow(sv)))
rows <- list()
for (i in seq_len(nrow(sv))) {
  r <- tryCatch(health_rows(i), error = function(e) { message(sv$svkey[i], ": ", conditionMessage(e)); NULL })
  rows[[sv$svkey[i]]] <- r
}
H <- do.call(rbind, rows)
cat(sprintf("[1] region health rows: %d from %d surveys\n", nrow(H), length(unique(H$svkey))))
for (v in c("dtp3_reg","facility","educ_yrs","wealth_q"))
  cat(sprintf("    %-9s non-missing %.0f%%  mean %.1f\n", v, 100*mean(is.finite(H[[v]])), mean(H[[v]], na.rm=TRUE)))

## ---- merge into the expanded panel ------------------------------------------
d <- read.csv(file.path(RESULTS, "component2_region_data_expanded.csv"), stringsAsFactors = FALSE)
d <- merge(d, H, by = c("svkey","regkey"), all.x = TRUE)
cat(sprintf("\n[2] panel rows with all 4 region health covariates: %d of %d\n",
            sum(complete.cases(d[, c("dtp3_reg","facility","educ_yrs","wealth_q")])), nrow(d)))
cat("    correlations among region health covariates:\n")
print(round(cor(d[, c("dtp3_reg","facility","educ_yrs","wealth_q")], use = "pairwise.complete.obs"), 2))
write.csv(d, file.path(RESULTS, "component2_region_data_expanded_health.csv"), row.names = FALSE)

## ---- composite health-access index (PC1, higher = better access) ------------
cc  <- complete.cases(d[, c("dtp3_reg","facility","educ_yrs","wealth_q")])
pc  <- prcomp(scale(d[cc, c("dtp3_reg","facility","educ_yrs","wealth_q")]))
ld  <- pc$rotation[, 1]; if (mean(ld) < 0) ld <- -ld            # orient so higher = better
d$hai <- NA_real_
d$hai[cc] <- as.numeric(scale(d[cc, c("dtp3_reg","facility","educ_yrs","wealth_q")]) %*% ld)
cat(sprintf("\n[3] health-access index PC1: %.0f%% of variance; loadings %s\n",
            100 * summary(pc)$importance[2, 1],
            paste(sprintf("%s=%.2f", names(ld), ld), collapse = "  ")))

## ---- refit primary model: baseline vs +region-health, SAME sample -----------
need <- c("m1mo5y","pfpr10","log_gdp","pct_urban","year_c","iso3","svkey",
          "dtp3","dtp3_reg","facility","educ_yrs","wealth_q")
fitd <- d[complete.cases(d[, need]) & is.finite(d$exposure) & d$exposure > 0 & d$pfpr2_10 >= 1, ]
fitd$deaths  <- round(fitd$m1mo5y / 1000 * fitd$exposure)
fitd$country <- factor(fitd$iso3); fitd$svkey <- factor(fitd$svkey)
fitd$hai <- as.numeric(scale(fitd[, c("dtp3_reg","facility","educ_yrs","wealth_q")]) %*% ld)
cat(sprintf("\n[4] refit sample (all covariates present): %d region-years, %d countries, %d surveys\n",
            nrow(fitd), nlevels(fitd$country), nlevels(fitd$svkey)))

RE  <- "s(country,bs=\"re\") + s(country,pfpr10,bs=\"re\") + offset(log(exposure))"
fit <- function(core) gam(as.formula(paste("deaths ~", core, "+ s(year_c) +", RE)),
                          family = nb(), method = "REML", data = fitd)
eff <- function(m) { b <- summary(m)$p.table["pfpr10", ]
  sprintf("%+.1f%% (%+.1f to %+.1f)", (exp(b[1])-1)*100, (exp(b[1]-1.96*b[2])-1)*100, (exp(b[1]+1.96*b[2])-1)*100) }
m_base <- fit("pfpr10 + dtp3 + log_gdp + pct_urban")                                         # national DTP3 only
m_full <- fit("pfpr10 + dtp3_reg + facility + educ_yrs + wealth_q + log_gdp + pct_urban")     # + region health
m_hai  <- fit("pfpr10 + hai + log_gdp + pct_urban")                                          # composite index

cat("\n=== PfPR effect on post-neonatal mortality (per +10 PfPR2-10 pts) ===\n")
cat(sprintf("  baseline (national DTP3)          : %s   AIC %.1f\n", eff(m_base), AIC(m_base)))
cat(sprintf("  + region health (4 covars)        : %s   AIC %.1f\n", eff(m_full), AIC(m_full)))
cat(sprintf("  + composite health-access index   : %s   AIC %.1f\n", eff(m_hai),  AIC(m_hai)))
cat("\n[full model] covariate coefficients (per unit, on log-mortality scale):\n")
print(round(summary(m_full)$p.table[c("dtp3_reg","facility","educ_yrs","wealth_q","log_gdp","pct_urban"), c("Estimate","Std. Error","Pr(>|z|)")], 4))
saveRDS(list(m_base = m_base, m_full = m_full, m_hai = m_hai, loadings = ld,
             center = round(mean(d$year, na.rm = TRUE))), file.path(RESULTS, "region_health_models.rds"))
cat("\nsaved: results/component2_region_data_expanded_health.csv + region_health_models.rds\n")

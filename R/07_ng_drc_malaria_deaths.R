# =============================================================================
# 07_ng_drc_malaria_deaths.R — malaria-attributable U5 deaths for Nigeria & DR
# Congo, NATIONAL and SUBNATIONAL (sum over admin-1), by the two methods:
#
#   Component 1 (share):  malaria deaths = share_hat(PfPR, GDP, DTP3)/100 x D
#                         share_hat from lm(share_u5 ~ pfpr_pct + log_gdp + dtp3)
#   Component 2 (AF):     malaria deaths = AF(PfPR) x D,  AF = 1 - exp(-beta*p/10)
#                         beta = pfpr10 fixed effect of the DHS U5MR mixed model
#
# D = all-cause U5 deaths (national U5MR x births, IGME/WB). Subnational splits D
# across admin-1 by population weight (GPW), so the subnational denominator sums
# EXACTLY to the national one; only the malaria FRACTION varies with local PfPR.
# National vs subnational therefore differ only through the shape of share()/AF()
# applied to disaggregated vs aggregate prevalence (share is linear -> ~equal;
# AF is concave -> subnational sum <= national point estimate).
#
# Uncertainty: ONLY the malaria coefficient (C1 pfpr_pct slope; C2 pfpr10 beta).
# Prevalence, GDP, DTP3, births, U5MR and population weights are treated as fixed.
# Because the coefficient is shared across admin-1 units, the CI on a total is the
# total re-evaluated at the coefficient's CI bounds (fully correlated), not a sum
# of independent variances.
# =============================================================================
source("R/00_utils.R")
suppressMessages({library(terra); library(sf); library(malariaAtlas); library(lme4)})

ISOS <- c(NGA = "Nigeria", COD = "DR Congo")

## ---- national inputs --------------------------------------------------------
c1 <- read.csv(file.path(RESULTS, "component1_country_data.csv"), stringsAsFactors = FALSE)
c1$log_gdp <- log(c1$gdp_pc)
nat <- c1[c1$iso3 %in% names(ISOS), ]
nat$allcause_u5 <- nat$u5mr / 1000 * nat$births            # D: all-cause U5 deaths

## ---- Component 1 coefficients (share_u5 model; shared spec 00_utils.R) -------
m1  <- fit_c1_share(c1, "share_u5")
b1  <- coef(m1); se_p1 <- summary(m1)$coefficients["pfpr_pct", "Std. Error"]
share_hat <- function(p, lg, dt, bpf) (b1["(Intercept)"] + bpf * p + b1["log_gdp"] * lg + b1["dtp3"] * dt)

## ---- Component 2 coefficient (pfpr10 fixed effect, U5MR LMM; shared spec) ----
d2 <- read.csv(file.path(RESULTS, "component2_region_data.csv"), stringsAsFactors = FALSE)
m2 <- fit_c2_lmm(d2, "u5mr")$m                       # no-stunting default (full country coverage)
beta <- fixef(m2)["pfpr10"]; se_b2 <- sqrt(vcov(m2)["pfpr10", "pfpr10"])
AF <- function(p, b) 1 - exp(-b * p / 10)

cat(sprintf("C1 pfpr_pct slope = %.4f pp/PfPR-pt (SE %.4f)\n", b1["pfpr_pct"], se_p1))
cat(sprintf("C2 pfpr10 beta    = %.4f (SE %.4f) => +%.1f%% U5MR per +10 PfPR pts\n\n",
            beta, se_b2, (exp(beta) - 1) * 100))

## ---- admin-1 PfPR + population weight (cached; getShp + GPW otherwise) -------
ADM <- file.path(DATA, "pfpr_pop_admin1_ng_cd_2024.csv")
if (file.exists(ADM)) {
  adm <- read.csv(ADM, stringsAsFactors = FALSE)
} else {
  pf_full <- read_pfpr(); den_full <- terra::rast(GPW_TIF)
  extract_admin1 <- function(iso, label) {
    v   <- terra::makeValid(terra::vect(sf::st_make_valid(getShp(ISO = iso, admin_level = "admin1"))))
    pf  <- terra::crop(pf_full, v)
    den <- terra::resample(terra::crop(den_full, pf), pf, method = "bilinear")
    w   <- terra::mask(den, pf)
    an  <- terra::extract(pf * w, v, fun = sum, na.rm = TRUE, exact = TRUE, ID = FALSE)
    aw  <- terra::extract(w,      v, fun = sum, na.rm = TRUE, exact = TRUE, ID = FALSE)
    data.frame(country = label, iso3 = iso, area = terra::values(v)$name_1,
               pfpr_pct = 100 * an[[1]] / aw[[1]], pop_w = aw[[1]], stringsAsFactors = FALSE)
  }
  adm <- do.call(rbind, lapply(names(ISOS), function(i) extract_admin1(i, ISOS[[i]])))
  adm <- adm[is.finite(adm$pfpr_pct) & is.finite(adm$pop_w), ]
  write.csv(adm, ADM, row.names = FALSE)
}

## ---- estimate deaths, national and subnational-sum --------------------------
estimate <- function(iso) {
  n  <- nat[nat$iso3 == iso, ]; D <- n$allcause_u5
  a  <- adm[adm$iso3 == iso, ]; w <- a$pop_w / sum(a$pop_w); Du <- D * w   # allocate D by pop
  pbar <- sum(w * a$pfpr_pct)                                              # pop-wtd mean admin1 PfPR
  # Component 1 (linear share) — totals at national and subnational prevalence
  c1_tot <- function(bpf, p, Dvec) sum(share_hat(p, n$log_gdp, n$dtp3, bpf) / 100 * Dvec)
  # Component 2 (concave AF)
  c2_tot <- function(b, p, Dvec) sum(AF(p, b) * Dvec)
  row <- function(level, p, Dvec) data.frame(
    country = ISOS[[iso]], level = level, allcause_u5 = round(D),
    c1 = c1_tot(b1["pfpr_pct"], p, Dvec),
    c1_lo = c1_tot(b1["pfpr_pct"] - 1.96 * se_p1, p, Dvec),
    c1_hi = c1_tot(b1["pfpr_pct"] + 1.96 * se_p1, p, Dvec),
    c2 = c2_tot(beta, p, Dvec),
    c2_lo = c2_tot(beta - 1.96 * se_b2, p, Dvec),
    c2_hi = c2_tot(beta + 1.96 * se_b2, p, Dvec),
    ihme = n$deaths, pfpr = ifelse(level == "National", n$pfpr_pct, pbar),
    n_units = ifelse(level == "National", 1L, nrow(a)), stringsAsFactors = FALSE)
  rbind(row("National", n$pfpr_pct, D),
        row(sprintf("Subnational (sum of %d admin-1)", nrow(a)), a$pfpr_pct, Du))
}
res <- do.call(rbind, lapply(names(ISOS), estimate))
write.csv(res, file.path(RESULTS, "ng_drc_malaria_deaths.csv"), row.names = FALSE)

## ---- print table ------------------------------------------------------------
fmt <- function(x, lo, hi) sprintf("%s [%s-%s]", format(round(x), big.mark=","),
                                   format(round(lo), big.mark=","), format(round(hi), big.mark=","))
out <- with(res, data.frame(Country = country, Level = level, PfPR = sprintf("%.1f%%", pfpr),
  `AllCause_U5` = format(allcause_u5, big.mark = ","),
  `Component1` = fmt(c1, c1_lo, c1_hi),
  `Component2` = fmt(c2, c2_lo, c2_hi),
  `IHME_ref` = format(ihme, big.mark = ","), check.names = FALSE))
cat("Malaria-attributable UNDER-5 deaths (point [95% CI from malaria coefficient only]):\n\n")
print(out, row.names = FALSE)
cat("\nsaved: results/ng_drc_malaria_deaths.csv\n")

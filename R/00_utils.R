# =============================================================================
# 00_utils.R — shared configuration and helper functions.
# Source this first (all other R/ scripts start with source("R/00_utils.R")).
# Run scripts from the repository root.
# =============================================================================
suppressMessages({
  library(terra); library(sf); library(jsonlite); library(httr); library(countrycode)
})

## ---- paths & constants ------------------------------------------------------
DATA    <- "data"
RESULTS <- "results"
dir.create(RESULTS, showWarnings = FALSE)

PFPR_TIF <- file.path(DATA, "pfpr_2to10_africa_2024.tif")                 # MAP PfPR2-10 (proportion), year 2024
GPW_TIF  <- file.path(DATA, "pop", "gpw_v4_population_density_rev11_2020_2.5m.tif")  # GPW density weight
MAP_YEAR <- 2024                                                          # latest MAP release year (used as ~2025 proxy)

## ---- small helpers ----------------------------------------------------------
# normalise a place name to a join key (lowercase, ASCII, alphanumeric only)
rkey <- function(x) gsub("[^a-z0-9]", "", tolower(iconv(as.character(x), "", "ASCII//TRANSLIT")))

## ---- World Bank API (long format, with retries) ----------------------------
wb_fetch <- function(indicator, date = "2000:2024") {
  u <- sprintf("https://api.worldbank.org/v2/country/all/indicator/%s?date=%s&format=json&per_page=20000",
               indicator, date)
  for (attempt in 1:4) {
    r <- tryCatch(GET(u, timeout(180)), error = function(e) NULL)
    if (!is.null(r) && status_code(r) == 200) {
      j <- fromJSON(content(r, "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)[[2]]
      d <- data.frame(iso3 = j$countryiso3code, year = as.integer(j$date),
                      value = as.numeric(j$value), stringsAsFactors = FALSE)
      return(d[!is.na(d$value) & nchar(d$iso3) == 3, ])
    }
    message("  World Bank retry ", attempt, " for ", indicator)
  }
  stop("World Bank fetch failed: ", indicator)
}

# most recent non-NA year per country
wb_latest <- function(indicator, date = "2010:2024") {
  d <- wb_fetch(indicator, date)
  d <- d[order(d$iso3, -d$year), ]
  d[!duplicated(d$iso3), c("iso3", "year", "value")]
}

# Sub-Saharan Africa country list (World Bank region classification)
ssa_countries <- function() {
  j <- fromJSON("https://api.worldbank.org/v2/country?format=json&per_page=400")[[2]]
  d <- data.frame(iso3 = j$id, country = trimws(j$name),
                  region = trimws(j$region$value), stringsAsFactors = FALSE)
  d[d$region == "Sub-Saharan Africa", c("iso3", "country")]
}

## ---- MAP PfPR raster --------------------------------------------------------
read_pfpr <- function() { stopifnot(file.exists(PFPR_TIF)); terra::rast(PFPR_TIF) }

# population-weighted national/subnational PfPR2-10 (%) for polygons `v` (SpatVector)
pfpr_by_admin <- function(pf, v) {
  stopifnot(file.exists(GPW_TIF))
  den <- terra::resample(terra::crop(terra::rast(GPW_TIF), pf), pf, method = "bilinear")
  w   <- terra::mask(den, pf); num <- pf * w
  an  <- terra::extract(num, v, fun = sum, na.rm = TRUE, exact = TRUE, ID = FALSE)
  aw  <- terra::extract(w,   v, fun = sum, na.rm = TRUE, exact = TRUE, ID = FALSE)
  100 * an[[1]] / aw[[1]]
}

## ---- age-standardise parasite prevalence to the PfPR2-10 reference ----------
# DHS measures children ~6-59 months (0.5-5y); MAP PfPR is standardised to ages 2-10.
# Uses malariaAtlas::convertPrevalence (Smith et al. 2007 age-prevalence model).
# Input & output are PROPORTIONS (0-1).
to_pfpr210 <- function(p_prop) {
  out <- rep(NA_real_, length(p_prop)); ok <- is.finite(p_prop)
  if (any(ok)) {                                        # age args must match prevalence length
    n <- sum(ok)
    out[ok] <- suppressMessages(as.numeric(malariaAtlas::convertPrevalence(
      p_prop[ok], rep(0.5, n), rep(5, n), rep(2, n), rep(10, n))))
  }
  out
}

## ---- IHME / GBD export parser (classic GBD or Data Explorer schema) ---------
# Returns iso3, deaths (Number), rate (per 100k, if present), year.
parse_ihme <- function(csv) {
  d <- read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
  names(d) <- tolower(trimws(names(d)))
  alias <- c(location = "location_name", condition = "cause_name", measure = "measure_name",
             unit = "metric_name", sex = "sex_name", age = "age_name", value = "val")
  for (s in names(alias)) if (s %in% names(d) && !(alias[s] %in% names(d))) names(d)[names(d) == s] <- alias[s]
  norm <- function(x) tolower(trimws(as.character(x)))
  f <- d
  if ("measure_name" %in% names(f)) f <- f[grepl("death",   norm(f$measure_name)), ]
  if ("cause_name"   %in% names(f)) f <- f[grepl("malaria", norm(f$cause_name)),   ]
  if ("sex_name"     %in% names(f)) f <- f[grepl("both",    norm(f$sex_name)),     ]
  ac <- grep("age.*name|^age_name$|^age$", names(f), value = TRUE)[1]
  if (!is.na(ac)) f <- f[grepl("under 5|under-5|<5|0-4|0 to 4", norm(f[[ac]])), ]
  if ("year" %in% names(f)) { yr <- max(f$year); f <- f[f$year == yr, ] } else yr <- NA
  iso_of <- function(x) countrycode(x, "country.name", "iso3c", warn = FALSE)
  num  <- f[grepl("number", norm(f$metric_name)), ]
  rate <- f[grepl("rate",   norm(f$metric_name)), ]
  out <- data.frame(iso3 = iso_of(num$location_name), deaths = as.numeric(num$val),
                    year = yr, stringsAsFactors = FALSE)
  out$rate <- as.numeric(rate$val)[match(out$iso3, iso_of(rate$location_name))]
  out[!is.na(out$iso3), ]
}

## ---- DHS.rates: all-cause mortality by region -------------------------------
# Returns regkey, u5mr (5q0) and m1mo5y (U5MR - neonatal), per 1,000 live births.
mort_by_region <- function(br, regvar) {
  if (!regvar %in% names(br)) return(NULL)
  r <- tryCatch(suppressMessages(DHS.rates::chmort(br, Class = regvar)), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  u <- r[grepl("^U5MR", rownames(r)), c("Class", "R")]; names(u)[2] <- "u5mr"
  n <- r[grepl("^NNMR", rownames(r)), c("Class", "R")]; names(n)[2] <- "nnmr"
  m <- merge(u, n, by = "Class")
  data.frame(regkey = rkey(m$Class), u5mr = m$u5mr, m1mo5y = m$u5mr - m$nnmr, stringsAsFactors = FALSE)
}

## ---- DHS PR recode: prevalence + region covariates --------------------------
# Design-weighted (hv005) by region: RDT & microscopy positivity in tested 6-59mo
# children (the tested set is 6-59mo de facto by design), % urban (hv025),
# % stunted (hc70 HAZ < -2). Returns one row per region.
prev_cov_by_region <- function(pr, regvar) {
  if (!regvar %in% names(pr)) return(NULL)
  w <- as.numeric(pr$hv005) / 1e6; reg <- as.character(pr[[regvar]])
  posfrac <- function(col) { v <- tolower(as.character(col)); t <- v %in% c("negative", "positive") & !is.na(reg)
    if (!any(t)) return(NULL); 100 * tapply(w[t] * (v[t] == "positive"), reg[t], sum) / tapply(w[t], reg[t], sum) }
  rdt <- posfrac(pr$hml35); if (is.null(rdt)) return(NULL)
  mic <- if ("hml32" %in% names(pr)) posfrac(pr$hml32) else NULL
  frac_if <- function(flag, ok) { k <- ok & !is.na(reg); if (!any(k)) return(NULL)
    100 * tapply(w[k] * flag[k], reg[k], sum) / tapply(w[k], reg[k], sum) }
  urb <- if ("hv025" %in% names(pr)) frac_if(tolower(as.character(pr$hv025)) == "urban", rep(TRUE, length(reg))) else NULL
  st  <- if ("hc70"  %in% names(pr)) { h <- as.numeric(pr$hc70); h[h >= 9990] <- NA; frac_if(h < -200, !is.na(h)) } else NULL
  regs <- names(rdt)
  pick <- function(v) if (is.null(v)) NA_real_ else as.numeric(v[regs])
  data.frame(regkey = rkey(regs), region = regs, rdt = as.numeric(rdt),
             mic = pick(mic), pct_urban = pick(urb), stunting = pick(st), stringsAsFactors = FALSE)
}

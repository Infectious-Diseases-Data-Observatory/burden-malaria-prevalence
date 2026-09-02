# =============================================================================
# 40_build_person_time.R — deaths and person-time by region, window and age.
#
# The person-time design. For every child in a Births Recode, the 60 months
# before that child's own interview are cut into five 12-month windows (window 1
# ends at the interview) and the child's life into DHS.rates' eight age
# segments. A child contributes person-months to each (window, segment) cell its
# age and calendar time pass through, censored at the interview month or the
# month of death, and one death to the cell in which it died. Sums are weighted
# by the sampling weight and taken to region x window x segment, with cluster
# sums kept for a delete-one-cluster jackknife.
#
# Unlike the synthetic-cohort probabilities the panel carries, these are
# observed counts and observed time at risk, so a Poisson or negative-binomial
# regression on them has the right likelihood, zero-death cells are data, and
# each window can be paired with the MAP prevalence of its own calendar year.
#
# Conventions
#   * months are DHS century-month codes; a child's age in month t is t - b3;
#   * the interview month itself is incomplete and excluded (as chmort does);
#   * a death at age b7 falls in month b3 + b7, which counts as a full month at
#     risk; children with no death are at risk to the month before interview;
#   * windows are defined per child from its own interview month, so a survey's
#     windows spread over its fieldwork period; the calendar-year composition
#     of each region-window's person-time is recorded for pairing with MAP.
#
# Outputs (data/derived_dhs unless noted)
#   person_time_region_window_segment.csv   weighted deaths and person-months
#   person_time_window_year_shares.csv      person-time share by calendar year
#   data/person_time_cluster_sums.rds       cluster-level sums for the jackknife
#   results/dhs_rebuild/person_time_validation.csv   against chmort, 60 months
# =============================================================================

source("R_dhs/00_config.R")
required_packages("DHS.rates")

N_WINDOWS <- 5L
WINDOW_MONTHS <- 12L
CACHE <- file.path(DATA_DIR, "person_time_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
PERSON_TIME_CSV <- file.path(DERIVED_DIR, "person_time_region_window_segment.csv")
YEAR_SHARE_CSV <- file.path(DERIVED_DIR, "person_time_window_year_shares.csv")
CLUSTER_SUMS_RDS <- file.path(DATA_DIR, "person_time_cluster_sums.rds")
VALIDATION_SURVEYS <- 4L

cmc_year <- function(cmc) 1900L + (cmc - 1L) %/% 12L
cmc_january <- function(year) (year - 1900L) * 12L + 1L

# Ethiopian DHS recodes give every date in the Ethiopian calendar, whose
# century-month codes run 92 months (7 years 8 months) behind the Gregorian
# ones; the five surveys ET41FL-ET81FL all show interview years eight years
# before their registry year. Windows must be in Gregorian months to pair with
# MAP's calendar years, so the offset is detected from the modal interview year
# and added to every date. Rates are unaffected (both dates shift together).
ETHIOPIAN_OFFSET_MONTHS <- 92L
calendar_offset_months <- function(v008, survey_year) {
  modal_year <- cmc_year(as.integer(stats::median(v008, na.rm = TRUE)))
  gap <- survey_year - modal_year
  if (gap >= 7L && gap <= 9L) ETHIOPIAN_OFFSET_MONTHS else 0L
}

## ---- one survey ----------------------------------------------------------------------
person_time_survey <- function(br, region, cmc_offset = 0L) {
  keep <- !is.na(br$v005) & br$v005 > 0 & !is.na(region) & nzchar(region) &
    is.finite(br$b3) & is.finite(br$v008)
  br <- br[keep, , drop = FALSE]
  region <- region[keep]
  weight <- br$v005 / 1e6
  cluster <- as.integer(br$v021)
  b3 <- as.integer(br$b3) + cmc_offset
  v008 <- as.integer(br$v008) + cmc_offset
  b7 <- as.integer(br$b7)
  death_month <- b3 + b7                       # NA for children alive
  # exclusive end of time at risk
  end <- ifelse(is.na(b7), v008, pmin(death_month + 1L, v008))

  cell_rows <- list()
  share_rows <- list()
  for (w in seq_len(N_WINDOWS)) {
    hi <- v008 - WINDOW_MONTHS * (w - 1L)      # exclusive
    lo <- hi - WINDOW_MONTHS
    window_deaths_any <- !is.na(death_month) & death_month >= lo & death_month < hi &
      death_month < v008
    for (s in seq_along(UNDER5_SEGMENTS)) {
      seg_lo <- UNDER5_SEGMENTS[[s]][1]
      seg_hi <- UNDER5_SEGMENTS[[s]][2]
      start <- pmax(b3 + seg_lo, lo)
      stop <- pmin(b3 + seg_hi, hi, end)
      months <- pmax(0L, stop - start)
      died <- window_deaths_any & !is.na(b7) & b7 >= seg_lo & b7 < seg_hi
      use <- months > 0L | died
      if (!any(use)) next
      key <- paste(region[use], cluster[use], sep = "|")
      pm_w <- rowsum(weight[use] * months[use], key)
      d_w <- rowsum(weight[use] * died[use], key)
      pm_n <- rowsum(as.numeric(months[use]), key)
      d_n <- rowsum(as.numeric(died[use]), key)
      parts <- do.call(rbind, strsplit(rownames(pm_w), "|", fixed = TRUE))
      cell_rows[[length(cell_rows) + 1L]] <- data.frame(
        regkey = parts[, 1], cluster = as.integer(parts[, 2]), window = w,
        segment = s, seg_lo = seg_lo, seg_hi = seg_hi,
        person_months_w = as.numeric(pm_w), deaths_w = as.numeric(d_w),
        person_months_n = as.numeric(pm_n), deaths_n = as.numeric(d_n),
        stringsAsFactors = FALSE)
    }
    # calendar-year composition of the window's under-5 person-time
    e_lo <- pmax(lo, b3)
    e_hi <- pmin(hi, end, b3 + 60L)
    exposed <- e_hi > e_lo
    if (any(exposed)) {
      y1 <- cmc_year(lo[exposed])
      jan_next <- cmc_january(y1 + 1L)
      m1 <- pmax(0L, pmin(e_hi[exposed], jan_next) - e_lo[exposed])
      m2 <- pmax(0L, e_hi[exposed] - pmax(e_lo[exposed], jan_next))
      key1 <- paste(region[exposed], y1, sep = "|")
      key2 <- paste(region[exposed], y1 + 1L, sep = "|")
      shares <- rbind(
        data.frame(key = key1, pm = weight[exposed] * m1, stringsAsFactors = FALSE),
        data.frame(key = key2, pm = weight[exposed] * m2, stringsAsFactors = FALSE))
      shares <- shares[shares$pm > 0, , drop = FALSE]
      tot <- rowsum(shares$pm, shares$key)
      parts <- do.call(rbind, strsplit(rownames(tot), "|", fixed = TRUE))
      share_rows[[length(share_rows) + 1L]] <- data.frame(
        regkey = parts[, 1], window = w, calendar_year = as.integer(parts[, 2]),
        person_months_w = as.numeric(tot), stringsAsFactors = FALSE)
    }
  }
  cells <- do.call(rbind, cell_rows)
  shares <- do.call(rbind, share_rows)
  shares <- shares[order(shares$regkey, shares$window, shares$calendar_year), ]
  list(cells = cells, shares = shares)
}

## ---- validation against chmort at 60 months --------------------------------------
# From the person-time table, a segment's death probability is
# 1 - exp(-rate x width) and U5MR = 1000 x (1 - prod(1 - q)). chmort uses the
# synthetic-cohort ratio of deaths to children exposed instead, so agreement
# should be close but not exact; a difference beyond a few per 1000 would mean
# a construction error.
piecewise_rates <- function(cells) {
  agg <- aggregate(cbind(deaths_w, person_months_w) ~ regkey + segment + seg_lo + seg_hi,
                   cells, sum)
  out <- do.call(rbind, lapply(split(agg, agg$regkey), function(d) {
    d <- d[order(d$segment), ]
    rate <- d$deaths_w / d$person_months_w
    q <- 1 - exp(-rate * (d$seg_hi - d$seg_lo))
    q[!is.finite(q)] <- 0
    data.frame(regkey = d$regkey[1],
               u5mr = 1000 * (1 - prod(1 - q)),
               nnmr = 1000 * q[1],
               postneonatal = 1000 * (1 - prod(1 - q)) - 1000 * q[1],
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}
validate_against_chmort <- function(br, region, cells, svkey) {
  br$region_final <- region
  reference <- tryCatch(
    suppressMessages(DHS.rates::chmort(
      br[!is.na(br$region_final), , drop = FALSE], Class = "region_final",
      Period = 60, Strata = strata_candidates(br)[1])),
    error = function(e) NULL)
  if (is.null(reference)) return(NULL)
  mine <- piecewise_rates(cells)
  rows <- lapply(c("U5MR", "NNMR"), function(label) {
    block <- reference[grepl(paste0("^", label), rownames(reference)), c("Class", "R")]
    m <- merge(data.frame(regkey = rkey(block$Class), chmort = as.numeric(block$R)),
               data.frame(regkey = mine$regkey, person_time = mine[[tolower(label)]]),
               by = "regkey")
    if (!nrow(m)) return(NULL)
    data.frame(svkey = svkey, outcome = label, regions = nrow(m),
               correlation = stats::cor(m$chmort, m$person_time),
               median_ratio = stats::median(m$person_time / m$chmort),
               max_abs_difference = max(abs(m$person_time - m$chmort)),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## ---- all surveys -------------------------------------------------------------------------
registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
map <- read.csv(MAP_REGION_CSV, stringsAsFactors = FALSE)
registry <- registry[registry$svkey %in% map$svkey, , drop = FALSE]
message("Building person-time for ", nrow(registry), " surveys")

all_cells <- list(); all_shares <- list(); validation <- list()
validated <- 0L
for (i in seq_len(nrow(registry))) {
  survey <- registry[i, , drop = FALSE]
  cache <- file.path(CACHE, paste0(survey$svkey, ".rds"))
  if (file.exists(cache)) {
    got <- readRDS(cache)
  } else {
    survey_map <- map[map$svkey == survey$svkey, , drop = FALSE]
    br <- tryCatch(readRDS(as.character(survey$local_recode)), error = function(e) NULL)
    if (is.null(br)) { message("  ", survey$svkey, ": recode unreadable"); next }
    region <- tryCatch(survey_region_vector(br, survey_map, survey, registry),
                       error = function(e) NULL)
    if (is.null(region) || all(is.na(region))) {
      message("  ", survey$svkey, ": no region assignment"); rm(br); gc(FALSE); next
    }
    cmc_offset <- calendar_offset_months(br$v008, survey$year)
    if (cmc_offset != 0L) message("  ", survey$svkey, ": dates shifted by ",
                                  cmc_offset, " months (Ethiopian calendar)")
    got <- person_time_survey(br, region, cmc_offset)
    got$cmc_offset <- cmc_offset
    got$cells$svkey <- survey$svkey; got$cells$iso3 <- survey$iso3
    got$cells$survey_year <- survey$year
    got$shares$svkey <- survey$svkey
    if (validated < VALIDATION_SURVEYS) {
      got$validation <- tryCatch(validate_against_chmort(br, region, got$cells, survey$svkey),
                                 error = function(e) NULL)
      if (!is.null(got$validation)) validated <- validated + 1L
    }
    saveRDS(got, cache)
    rm(br); gc(FALSE)
    message(sprintf("  %s: %d regions, %.0f weighted deaths, %.0fk person-months",
                    survey$svkey, length(unique(got$cells$regkey)),
                    sum(got$cells$deaths_w), sum(got$cells$person_months_w) / 1e3))
  }
  all_cells[[survey$svkey]] <- got$cells
  all_shares[[survey$svkey]] <- got$shares
  if (!is.null(got$validation)) validation[[survey$svkey]] <- got$validation
}

cluster_sums <- do.call(rbind, all_cells)
rownames(cluster_sums) <- NULL
saveRDS(cluster_sums, CLUSTER_SUMS_RDS)

region_cells <- aggregate(
  cbind(deaths_w, person_months_w, deaths_n, person_months_n) ~
    svkey + iso3 + survey_year + regkey + window + segment + seg_lo + seg_hi,
  cluster_sums, sum)
region_cells$age_group <- segment_age_group(region_cells$seg_lo)
region_cells <- region_cells[order(region_cells$svkey, region_cells$regkey,
                                   region_cells$window, region_cells$segment), ]
write.csv(region_cells, PERSON_TIME_CSV, row.names = FALSE)

shares <- do.call(rbind, all_shares)
shares <- aggregate(person_months_w ~ svkey + regkey + window + calendar_year, shares, sum)
totals <- aggregate(person_months_w ~ svkey + regkey + window, shares, sum)
names(totals)[4] <- "window_total"
shares <- merge(shares, totals, by = c("svkey", "regkey", "window"))
shares$share <- shares$person_months_w / shares$window_total
shares <- shares[order(shares$svkey, shares$regkey, shares$window, shares$calendar_year), ]
write.csv(shares[, c("svkey", "regkey", "window", "calendar_year", "share",
                     "person_months_w")],
          YEAR_SHARE_CSV, row.names = FALSE)

if (length(validation)) {
  validation <- do.call(rbind, validation)
  rownames(validation) <- NULL
  write.csv(validation, file.path(RESULTS_DIR, "person_time_validation.csv"),
            row.names = FALSE)
  message("\nPerson-time life table against chmort (60 months):")
  print(transform(validation, correlation = round(correlation, 4),
                  median_ratio = round(median_ratio, 3),
                  max_abs_difference = round(max_abs_difference, 2)), row.names = FALSE)
}
message(sprintf(paste0(
  "\nPerson-time table: %d region-window-segment cells, %d regions, %d surveys; ",
  "%.0f weighted deaths over %.1f million weighted person-months"),
  nrow(region_cells), length(unique(paste(region_cells$svkey, region_cells$regkey))),
  length(unique(region_cells$svkey)), sum(region_cells$deaths_w),
  sum(region_cells$person_months_w) / 1e6))
message("Wrote ", basename(PERSON_TIME_CSV), ", ", basename(YEAR_SHARE_CSV),
        ", ", basename(CLUSTER_SUMS_RDS))

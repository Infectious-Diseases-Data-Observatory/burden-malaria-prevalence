# Audit of the improved-water covariate: for every survey in the analysis
# dataset, tabulate the household-level labels of v113 (source of drinking
# water) as the births recode carries them, flag which labels the covariate's
# regular expression counts as improved, and recompute the region shares.
setwd("/Users/jameswatson/Documents/Claude Projects/MIS:DHS malaria prevalence")
source("R_dhs/00_config.R")
OUT <- "/private/tmp/claude-501/-Users-jameswatson-Documents-Claude-Projects-MIS-DHS-malaria-prevalence/9c576ede-4a2a-47fd-a386-9f921fb6a27c/scratchpad"
INCLUDE <- "pipe|tap|standpipe|borehole|tube ?well|protected|rain|bottled|sachet"
EXCLUDE <- "unprotected"

reg <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
map <- read.csv(MAP_REGION_CSV, stringsAsFactors = FALSE)
OTHER <- list(v116 = list(include = "flush|septic|sewer|ventilated|vip|slab|composting",
                          exclude = "without slab|open pit|no facil|bush|field|hanging|bucket|somewhere"),
              v119 = list(include = "yes", exclude = NULL))
d <- read.csv(file.path(DERIVED_DIR, "dhs_analysis_dataset.csv"), stringsAsFactors = FALSE)
surveys <- unique(d$svkey)

one <- function(k) {
  f <- reg$local_recode[reg$svkey == k][1]
  if (!file.exists(f)) return(list(labels = NULL, note = data.frame(svkey = k, note = "recode missing")))
  b <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(b)) return(list(labels = NULL, note = data.frame(svkey = k, note = "readRDS failed")))
  if (!"v113" %in% names(b)) return(list(labels = NULL, note = data.frame(svkey = k, note = "no v113")))
  w <- suppressWarnings(as.numeric(b$v005)) / 1e6
  fh <- !duplicated(paste(b$v001, b$v002))
  raw <- b$v113
  lab <- tolower(as.character(raw))
  code <- if (is.factor(raw)) as.character(as.integer(raw)) else if (!is.null(attr(raw, "labels"))) as.character(unclass(raw)) else NA_character_
  usable <- fh & !is.na(lab) & lab != "missing" & nzchar(lab)
  matched <- grepl(INCLUDE, lab) & !grepl(EXCLUDE, lab)
  survey <- reg[reg$svkey == k, , drop = FALSE][1, ]
  survey_map <- map[map$svkey == k, , drop = FALSE]
  region <- survey_region_vector(b, survey_map, survey, reg)
  hh <- data.frame(label = lab[usable], n_hh = 1, w_hh = w[usable], stringsAsFactors = FALSE)
  tab <- aggregate(cbind(n_hh, w_hh) ~ label, data = hh, FUN = sum)
  tab$matched <- grepl(INCLUDE, tab$label) & !grepl(EXCLUDE, tab$label)
  tab$numeric_code <- grepl("^[0-9]+$", tab$label)
  # DHS-IV style recodes carry only group headings for the codes: 1x piped, 2x well water (open),
  # 3x covered well/borehole, 4x surface water, 51 rainwater, 61 tanker truck, 71 bottled, 96 other, 97 not de jure
  group_rule <- function(code) { c2 <- suppressWarnings(as.integer(code)); tens <- c2 %/% 10
    ifelse(is.na(c2), NA, ifelse(c2 %in% c(97, 99), NA, ifelse(tens %in% c(1, 3) | c2 %in% c(51, 71), TRUE, FALSE))) }
  tab$group_rule_improved <- ifelse(tab$numeric_code, group_rule(tab$label), NA)
  tab$svkey <- k
  tab$v113_class <- paste(class(raw), collapse = "/")
  # region shares as the pipeline computes them, plus the share of households
  # whose label matches nothing in either pattern and is not obviously unimproved
  region_share <- tapply(w[usable] * matched[usable], region[usable], sum) / tapply(w[usable], region[usable], sum)
  is_code <- grepl("^[0-9]+$", lab)
  corrected_flag <- ifelse(is_code, group_rule(lab), matched)
  dejure <- grepl("dejure|de jure", lab)
  keep <- usable & !dejure & !is.na(corrected_flag)
  corrected <- tapply(w[keep] * corrected_flag[keep], region[keep], sum) / tapply(w[keep], region[keep], sum)
  rs <- data.frame(svkey = k, regkey = names(region_share), recomputed = 100 * as.numeric(region_share),
                   corrected = 100 * as.numeric(corrected[names(region_share)]),
                   share_numeric_codes = 100 * as.numeric(tapply(w[usable] * is_code[usable], region[usable], sum) / tapply(w[usable], region[usable], sum)))
  dropped <- fh & (is.na(lab) | lab == "missing" | !nzchar(lab))
  note <- data.frame(svkey = k, note = sprintf("%d households, %.1f%% dropped as missing/blank, %d distinct labels, %.1f%% matched",
                                               sum(fh), 100 * mean(dropped[fh]), nrow(tab), 100 * sum(tab$w_hh[tab$matched]) / sum(tab$w_hh)))
  # the two other household covariates: label tables only
  other <- do.call(rbind, lapply(names(OTHER), function(v) {
    if (!v %in% names(b)) return(NULL)
    l <- tolower(as.character(b[[v]])); u <- fh & !is.na(l) & l != "missing" & nzchar(l)
    t <- aggregate(cbind(n_hh, w_hh) ~ label, data = data.frame(label = l[u], n_hh = 1, w_hh = w[u], stringsAsFactors = FALSE), FUN = sum)
    t$matched <- grepl(OTHER[[v]]$include, t$label) & (if (is.null(OTHER[[v]]$exclude)) TRUE else !grepl(OTHER[[v]]$exclude, t$label))
    t$numeric_code <- grepl("^[0-9]+$", t$label); t$svkey <- k; t$variable <- v; t }))
  list(labels = tab, regions = rs, note = note, other = other)
}
res <- parallel::mclapply(surveys, function(k) tryCatch(one(k), error = function(e) list(labels = NULL, note = data.frame(svkey = k, note = conditionMessage(e)))),
                          mc.cores = 4)
labels <- do.call(rbind, lapply(res, `[[`, "labels"))
regions <- do.call(rbind, lapply(res, function(r) r$regions))
notes <- do.call(rbind, lapply(res, `[[`, "note"))
other <- do.call(rbind, lapply(res, function(r) r$other))
write.csv(other, file.path(OUT, "v116_v119_labels_by_survey.csv"), row.names = FALSE)
write.csv(labels, file.path(OUT, "v113_labels_by_survey.csv"), row.names = FALSE)
write.csv(regions, file.path(OUT, "v113_region_shares.csv"), row.names = FALSE)
write.csv(notes, file.path(OUT, "v113_notes.csv"), row.names = FALSE)
# distinct labels across surveys with total weighted households and surveys using them
distinct <- aggregate(cbind(w_hh, n_hh) ~ label + matched + numeric_code, data = labels, FUN = sum)
distinct$surveys <- as.integer(table(labels$label)[distinct$label])
distinct <- distinct[order(-distinct$w_hh), ]
write.csv(distinct, file.path(OUT, "v113_distinct_labels.csv"), row.names = FALSE)
cat("surveys:", length(surveys), " with labels:", length(unique(labels$svkey)), " distinct labels:", nrow(distinct), "\n")
cat("DONE\n")

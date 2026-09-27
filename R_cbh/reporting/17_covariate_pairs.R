#!/usr/bin/env Rscript
# Supplementary figure (added 27 September 2026, at the user's request): bivariate scatter plots of the
# 17 adjustment variables in the primary models, one point per survey-region, coloured by survey
# programme (DHS or MICS). The 13 regional variables are constant within a survey-region; the four
# national annual series vary with band-entry year, so each survey-region shows their mean over its
# child-band records (the values the models see, averaged). Variables are on the modelling scale
# before standardisation (natural log for HIV incidence, GDP and health expenditure). Lower panels:
# scatter (x = column variable, y = row variable); diagonal: density by programme, scaled to the panel
# height; upper panels: Pearson r across survey-regions (all, DHS, MICS), shaded by |r|. Writes
# programme-level summaries only (correlations, means, SDs and observed ranges), no survey-region
# values. Versions without MICS surveys (the DHS-only history) are skipped.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
st <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION", "regional_mics"))
coverage <- fread(file.path(st$out, "survey_map/survey_coverage.csv"))
if (!"MICS" %in% coverage$type) { message("covariate_pairs: skipped, ", st$id, " has no MICS surveys"); quit(save = "no", status = 0) }
stopifnot(!anyDuplicated(coverage$survey), all(coverage$type %in% c("DHS", "MICS")))
spec <- cbh_primary_regional_spec(st); covs <- spec$covariates
stopifnot(length(covs) == 17L, identical(fread(file.path(st$out, "covariate_scaling.csv"))$variable, covs))

prepared <- readRDS(st$data)
stopifnot(identical(cbh_file_hash(st$data), unique(fread(file.path(st$out, "fit_manifest.csv"))$prepared_data_md5)))
d <- as.data.table(prepared$data[c("survey", "region", covs)]); rm(prepared); invisible(gc(FALSE))
stopifnot(nrow(d) == st$expected_records, !anyNA(d), all(unique(d$survey) %in% coverage$survey))
const <- d[, lapply(.SD, uniqueN), by = region, .SDcols = spec$regional]
stopifnot(all(as.matrix(const[, -1]) == 1L))
sr <- d[, c(list(records = .N), lapply(.SD, mean)), by = .(survey, region), .SDcols = covs]
stopifnot(!anyDuplicated(sr$region), sum(sr$records) == nrow(d))
sr[, source := factor(coverage$type[match(survey, coverage$survey)], levels = c("DHS", "MICS"))]
rm(d); invisible(gc(FALSE))
out <- file.path(st$out, "covariate_pairs"); dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Grouping as in the covariate forest plot (14_covariate_forest.R), so related variables sit together.
ord <- c("mean_maternal_age_first_birth", "mean_maternal_education_years", "mean_wealth_quintile", "urban_pct", "electricity_pct",
  "improved_water_pct", "improved_sanitation_pct", "dtp3_pct", "measles_pct", "facility_delivery_pct", "short_birth_interval_pct",
  "wasting_pct", "stunting_pct", "log_hiv_incidence", "log_gdp_pc", "log_health_expenditure_pc", "political_stability")
# Labels wrapped to fit one panel width (they are used for both the column and the row strips).
lab <- c(mean_maternal_age_first_birth = "Age at\nfirst birth\n(years)", mean_maternal_education_years = "Maternal\neducation\n(years)",
  mean_wealth_quintile = "Wealth-\nquintile\nscore", urban_pct = "Urban\n(%)", electricity_pct = "Electricity\n(%)",
  improved_water_pct = "Improved\nwater (%)", improved_sanitation_pct = "Improved\nsanitation\n(%)", dtp3_pct = "DTP3\n(%)",
  measles_pct = "Measles\n(%)", facility_delivery_pct = "Facility\ndelivery\n(%)", short_birth_interval_pct = "Birth\ninterval\n<24 mo (%)",
  wasting_pct = "Wasting\n(%)", stunting_pct = "Stunting\n(%)", log_hiv_incidence = "Child HIV\nincidence\n(/1,000, log)",
  log_gdp_pc = "GDP per\ncapita\n(US$, log)", log_health_expenditure_pc = "Health exp.\nper capita\n(US$, log)", political_stability = "Political\nstability\n(WGI)")
stopifnot(setequal(ord, covs), setequal(names(lab), covs))

# Correlations: unweighted across survey-regions (as plotted) by programme, weighted by records (closer to
# what the model's design sees), and across surveys (one record-weighted mean per survey), which matters
# for the national series because they are nearly constant within a survey.
wcor <- function(x, y, w) { w <- w / sum(w); mx <- sum(w * x); my <- sum(w * y); sum(w * (x - mx) * (y - my)) / sqrt(sum(w * (x - mx)^2) * sum(w * (y - my)^2)) }
sv <- sr[, c(list(records = sum(records)), lapply(.SD, function(v) sum(v * records) / sum(records))), by = .(survey, source), .SDcols = covs]
pairs <- CJ(i = seq_along(ord), j = seq_along(ord))[i < j]
cors <- pairs[, { x <- sr[[ord[i]]]; y <- sr[[ord[j]]]; s <- sr$source
  .(variable_1 = ord[i], variable_2 = ord[j], r_all = cor(x, y), r_dhs = cor(x[s == "DHS"], y[s == "DHS"]), r_mics = cor(x[s == "MICS"], y[s == "MICS"]),
    r_all_record_weighted = wcor(x, y, sr$records), r_surveys = cor(sv[[ord[i]]], sv[[ord[j]]])) }, by = .(i, j)][, !c("i", "j")]
cbh_atomic_csv(as.data.frame(cors), file.path(out, "covariate_correlations.csv"))
summ <- rbindlist(lapply(ord, function(v) sr[, .(variable = v, survey_regions = .N, surveys = uniqueN(survey), mean = mean(get(v)), sd = sd(get(v)),
  min = min(get(v)), max = max(get(v))), by = source]))
total <- rbindlist(lapply(ord, function(v) data.table(variable = v, total_sd = sd(sr[[v]]))))
summ <- merge(summ, total, by = "variable")[order(match(variable, ord), source)]
diffs <- dcast(summ, variable ~ source, value.var = "mean")[, .(variable, mean_dhs = DHS, mean_mics = MICS)]
diffs <- merge(diffs, total, by = "variable")[, std_difference_mics_minus_dhs := (mean_mics - mean_dhs) / total_sd][order(match(variable, ord))]
cbh_atomic_csv(as.data.frame(summ), file.path(out, "covariate_summary_by_source.csv"))
cbh_atomic_csv(as.data.frame(diffs), file.path(out, "covariate_difference_by_source.csv"))

# Panel data for facet_grid(row ~ column, free scales): row variable on y, column variable on x.
lv <- unname(lab[ord]); rng <- lapply(setNames(ord, ord), function(v) range(sr[[v]]))
lower <- rbindlist(lapply(seq_along(ord), function(i) rbindlist(lapply(seq_len(i - 1L), function(j)
  sr[, .(row = lv[i], col = lv[j], x = get(ord[j]), y = get(ord[i]), source)]))))
setorder(lower, source)   # factor order: MICS rows last, so drawn over DHS
dg <- rbindlist(lapply(seq_along(ord), function(i) { v <- ord[i]; r <- rng[[v]]
  dens <- rbindlist(lapply(levels(sr$source), function(s) { k <- density(sr[source == s][[v]], from = r[1], to = r[2], n = 256)
    data.table(source = s, x = k$x, dens = k$y) }))
  dens[, `:=`(row = lv[i], col = lv[i], y = r[1] + .9 * diff(r) * dens / max(dens))] }))
upper <- rbindlist(lapply(seq_len(nrow(cors)), function(k) { a <- match(cors$variable_1[k], ord); b <- match(cors$variable_2[k], ord)
  rx <- rng[[ord[b]]]; ry <- rng[[ord[a]]]
  data.table(row = lv[a], col = lv[b], x = mean(rx), y_all = ry[1] + .68 * diff(ry), y_dhs = ry[1] + .40 * diff(ry), y_mics = ry[1] + .16 * diff(ry),
    r_all = cors$r_all[k], r_dhs = cors$r_dhs[k], r_mics = cors$r_mics[k]) }))
for (z in list(lower, dg, upper)) z[, `:=`(row = factor(row, levels = lv), col = factor(col, levels = lv))]
lower[, source := factor(source, levels = c("DHS", "MICS"))]; dg[, source := factor(source, levels = c("DHS", "MICS"))]

n_src <- sr[, .(regions = .N, surveys = uniqueN(survey)), by = source][order(source)]
src_lab <- setNames(sprintf("%s: %d surveys, %s survey-regions", n_src$source, n_src$surveys, format(n_src$regions, big.mark = ",")), as.character(n_src$source))
cols <- c(DHS = "#2F6DB0", MICS = "#E07B1F")        # points and densities
text_cols <- c(DHS = "#1F4E8C", MICS = "#A04A0C")   # darker variants for the r values (contrast on white and grey)
f2 <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.2f", round(x, 2) + 0))
# Interior breaks only, so tick labels of neighbouring panels do not run together.
inner <- function(l) { b <- scales::breaks_extended(5)(l); w <- diff(l); b[b >= l[1] + .1 * w & b <= l[2] - .1 * w] }
p <- ggplot() +
  geom_rect(data = upper, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = abs(r_all))) +
  geom_text(data = upper, aes(x, y_all, label = f2(r_all)), size = 3, fontface = "bold") +
  geom_text(data = upper, aes(x, y_dhs, label = f2(r_dhs)), size = 2.5, fontface = "bold", colour = text_cols[["DHS"]]) +
  geom_text(data = upper, aes(x, y_mics, label = f2(r_mics)), size = 2.5, fontface = "bold", colour = text_cols[["MICS"]]) +
  geom_point(data = lower, aes(x, y, colour = source), size = .55, alpha = .7, stroke = 0) +
  geom_line(data = dg, aes(x, y, colour = source), linewidth = .5) +
  scale_colour_manual(values = cols, labels = src_lab, name = NULL, guide = guide_legend(override.aes = list(size = 2.5, alpha = 1, linewidth = 1))) +
  scale_fill_gradient(low = "white", high = "grey75", limits = c(0, 1), breaks = c(0, .5, 1),
    name = "Upper panels: Pearson r across all survey-regions (black; grey shading by absolute value), DHS (blue) and MICS (orange) survey-regions",
    guide = guide_colourbar(title.position = "top", barwidth = unit(40, "mm"), barheight = unit(3, "mm"))) +
  scale_x_continuous(breaks = inner) + scale_y_continuous(breaks = inner) +
  facet_grid(row ~ col, scales = "free", switch = "both") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside", strip.background = element_blank(), strip.text = element_text(size = 7, face = "bold", lineheight = .9),
    strip.text.y.left = element_text(angle = 0, hjust = 1), panel.grid = element_blank(), panel.spacing = unit(1.5, "mm"),
    axis.text = element_text(size = 6), axis.ticks = element_line(linewidth = .25), legend.position = "top", legend.text = element_text(size = 8.5),
    legend.title = element_text(size = 8.5), legend.box.spacing = unit(2, "mm"))
ggsave(file.path(out, "sfig_covariate_pairs.png"), p, width = 13, height = 13.4, dpi = 300, device = ragg::agg_png, bg = "white")

# Surveys whose vaccination values are one national value in every region (WUENIC fallbacks and others).
flat <- sr[, .(regions = .N, dtp3 = uniqueN(dtp3_pct) == 1, measles = uniqueN(measles_pct) == 1), by = .(survey, source)][regions > 1 & (dtp3 | measles)]
mics_wide <- "data/derived_mics/regional_covariates_wide_mics.csv"
wsrc <- fread(mics_wide)[, .(all_wuenic = all(dtp3_pct_imputation_source == "UNICEF_WUENIC_country_year" & measles_pct_imputation_source == "UNICEF_WUENIC_country_year")), by = survey][all_wuenic == TRUE]
wuenic <- sr[survey %in% wsrc$survey, .(regions = .N), by = survey]
stopifnot(all(wuenic$survey %in% flat[dtp3 & measles]$survey))
flat_dtp3_other <- flat[dtp3 & !survey %in% wuenic$survey]   # flat DTP3 that is not a documented fallback
dtp3_note <- if (nrow(flat_dtp3_other)) sprintf(" DTP3 is %s%% in every region of %s.", paste(sr[survey %in% flat_dtp3_other$survey, unique(round(dtp3_pct))], collapse = "/"),
  paste(flat_dtp3_other$survey, collapse = ", ")) else ""
rmax <- cors[which.max(abs(r_all))]
writeLines(c("# Supplementary figure: adjustment variables by survey programme", "",
  sprintf(paste("Bivariate distributions of the 17 adjustment variables in the primary models across %s survey-regions: %s from %d DHS surveys (blue) and %s from %d MICS surveys (orange).",
    "Lower panels: one point per survey-region, with the column variable on the horizontal axis and the row variable on the vertical axis.",
    "Diagonal: density by programme, scaled to the panel height (no density axis).",
    "Upper panels: Pearson correlation across all survey-regions (black; shading shows |r|), and within DHS (blue) and MICS (orange) survey-regions; correlations are unweighted.",
    "The 13 regional variables are survey-region summaries; in %d MICS surveys (%d survey-regions) DTP3 and measles coverage are single national WUENIC values for every region, which appear as lines.%s",
    "The four national annual series (child HIV incidence per 1,000 from UNAIDS/UNICEF estimates, model-imputed for Nigeria and Comoros and derived from UNAIDS counts for Liberia; GDP and current health expenditure per capita in current US$; all natural log; and WGI political stability)",
    "are averaged over each survey-region's child-band records. They are nearly constant across the regions of a survey, so their panels show one cluster per survey and their correlations mainly reflect differences between surveys, weighted by each survey's number of regions.",
    "Values are on the modelling scale before standardisation."),
    format(nrow(sr), big.mark = ","), format(n_src[source == "DHS"]$regions, big.mark = ","), n_src[source == "DHS"]$surveys,
    format(n_src[source == "MICS"]$regions, big.mark = ","), n_src[source == "MICS"]$surveys, nrow(wuenic), sum(wuenic$regions), dtp3_note)), file.path(out, "CAPTION.md"))
inputs <- c(st$data, file.path(st$out, c("fit_manifest.csv", "covariate_scaling.csv", "survey_map/survey_coverage.csv")), mics_wide, "R_cbh/reporting/17_covariate_pairs.R", "R_cbh/primary/specification.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")), file.path(out, "provenance.csv"))
message(sprintf("Covariate pairs written: %s (%d survey-regions; largest |r| %.2f, %s vs %s); flat vaccination surveys: %s", out, nrow(sr), rmax$r_all,
  rmax$variable_1, rmax$variable_2, paste(sprintf("%s (%s)", flat$survey, ifelse(flat$dtp3 & flat$measles, "both", ifelse(flat$dtp3, "DTP3", "measles"))), collapse = ", ")))

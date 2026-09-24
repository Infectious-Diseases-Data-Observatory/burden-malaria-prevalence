#!/usr/bin/env Rscript
# Supplementary figures (added 24 September 2026): country and survey-region random intercepts of
# the seven primary age-band models, as multiplicative effects on the band hazard, exp(u), with 95%
# intervals from the coefficient covariance conditional on the smoothing parameters (Vp).
# Countries: forest plot grouped by UN M49 sub-region. Survey-regions (about 1,226 per band):
# ranked caterpillar plot per band, coloured by sub-region. Each panel reports the term's effective
# degrees of freedom and the estimated random-effect SD. Reads the saved fits; no refitting.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(mgcv)})
st <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION", "regional_mics"))
out <- file.path(st$out, "random_effects"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
ages <- cbh_config()$age_bands$age_band
band_label <- setNames(ifelse(ages == "<1", "<1 month", paste(ages, "months")), ages)
manifest <- fread(file.path(st$out, "fit_manifest.csv"))[series == "map_full"][match(ages, age_band)]
stopifnot(nrow(manifest) == 7L)

re_table <- function(f, term) {
  s <- f$smooth[[match(sprintf("s(%s)", term), vapply(f$smooth, `[[`, "", "label"))]]
  ix <- s$first.para:s$last.para; lev <- levels(f$model[[term]])
  stopifnot(inherits(s, "random.effect"), length(ix) == length(lev))
  data.table(level = lev, u = unname(coef(f)[ix]), se = sqrt(diag(f$Vp)[ix]))
}
pieces <- lapply(seq_len(nrow(manifest)), function(i) {
  path <- manifest$model_file[i]; stopifnot(identical(cbh_file_hash(path), manifest$md5[i]))
  f <- readRDS(path)$fit; a <- manifest$age_band[i]
  vc <- gam.vcomp(f)$vc
  sdv <- setNames(vc[, "std.dev"], sub("^s\\((.*)\\)$", "\\1", rownames(vc)))
  edf <- vapply(c("country", "region"), function(v) sum(f$edf[{s <- f$smooth[[match(sprintf("s(%s)", v), vapply(f$smooth, `[[`, "", "label"))]]; s$first.para:s$last.para}]), 0)
  res <- list(country = cbind(age_band = a, re_table(f, "country")), region = cbind(age_band = a, re_table(f, "region")),
    summary = data.table(age_band = a, term = c("country", "region"), levels = c(nlevels(f$model$country), nlevels(f$model$region)),
      edf = unname(edf), re_sd = unname(sdv[c("country", "region")])))
  rm(f); gc(FALSE); res
})
country <- rbindlist(lapply(pieces, `[[`, "country")); region <- rbindlist(lapply(pieces, `[[`, "region"))
summ <- rbindlist(lapply(pieces, `[[`, "summary"))
country[, `:=`(effect = exp(u), lower_95 = exp(u - 1.96 * se), upper_95 = exp(u + 1.96 * se))]
region[, `:=`(effect = exp(u), lower_95 = exp(u - 1.96 * se), upper_95 = exp(u + 1.96 * se))]
country[, excludes_1 := lower_95 > 1 | upper_95 < 1]; region[, excludes_1 := lower_95 > 1 | upper_95 < 1]

# Country names as in Figure 1 (DRC, CAF and RC abbreviated) and UN M49 sub-regions.
names_map <- c(AGO="Angola",BDI="Burundi",BEN="Benin",BFA="Burkina Faso",BWA="Botswana",CAF="CAF",CIV="Côte d'Ivoire",
  CMR="Cameroon",COD="DRC",COG="RC",COM="Comoros",CPV="Cabo Verde",DJI="Djibouti",ERI="Eritrea",ETH="Ethiopia",GAB="Gabon",
  GHA="Ghana",GIN="Guinea",GMB="The Gambia",GNB="Guinea-Bissau",GNQ="Equatorial Guinea",KEN="Kenya",LBR="Liberia",LSO="Lesotho",
  MDG="Madagascar",MLI="Mali",MOZ="Mozambique",MRT="Mauritania",MUS="Mauritius",MWI="Malawi",NAM="Namibia",NER="Niger",
  NGA="Nigeria",RWA="Rwanda",SEN="Senegal",SLE="Sierra Leone",SOM="Somalia",SSD="South Sudan",STP="São Tomé and Príncipe",
  SWZ="Eswatini",TCD="Chad",TGO="Togo",TZA="Tanzania",UGA="Uganda",ZAF="South Africa",ZMB="Zambia",ZWE="Zimbabwe")
m49 <- list(`West Africa`=c("BEN","BFA","CIV","CPV","GHA","GIN","GMB","GNB","LBR","MLI","MRT","NER","NGA","SEN","SLE","TGO"),
  `Central Africa`=c("AGO","CAF","CMR","COD","COG","GAB","GNQ","STP","TCD"),
  `East Africa`=c("BDI","COM","DJI","ERI","ETH","KEN","MDG","MOZ","MUS","MWI","RWA","SOM","SSD","TZA","UGA","ZMB","ZWE"),
  `Southern Africa`=c("BWA","LSO","NAM","SWZ","ZAF"))
sub_of <- setNames(rep(names(m49), lengths(m49)), unlist(m49))
country[, `:=`(iso3 = level, name = names_map[level], subregion = factor(sub_of[level], levels = names(m49)))]
region[, iso3 := sub(":.*$", "", level)]
region[, subregion := factor(sub_of[iso3], levels = names(m49))]
stopifnot(!anyNA(country$name), !anyNA(country$subregion), !anyNA(region$subregion))
cbh_atomic_csv(as.data.frame(country[, .(age_band, iso3, name, subregion, u, se, effect, lower_95, upper_95, excludes_1)]), file.path(out, "country_random_effects.csv"))
cbh_atomic_csv(as.data.frame(region[, .(age_band, region = level, iso3, subregion, u, se, effect, lower_95, upper_95, excludes_1)]), file.path(out, "region_random_effects.csv"))
cbh_atomic_csv(as.data.frame(summ), file.path(out, "random_effect_summary.csv"))

facet_lab <- function(tm) {
  x <- summ[term == tm]
  setNames(sprintf("%s\nEDF %.1f/%d, SD %.2f", band_label[x$age_band], x$edf, x$levels, x$re_sd), x$age_band)
}
theme_re <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(), strip.text.x = element_text(face = "bold"),
  panel.border = element_rect(fill = NA, colour = "grey85"), panel.spacing.x = grid::unit(8, "pt"), legend.position = "bottom")

# Country forest plot.
order_names <- unique(country[order(subregion, name)]$name)
country[, `:=`(label = factor(name, levels = rev(order_names)), band = factor(age_band, levels = ages))]
cl <- facet_lab("country")
p1 <- ggplot(country, aes(effect, label)) +
  geom_vline(xintercept = 1, colour = "grey55", linewidth = .4) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), width = 0, orientation = "y", colour = "#1F4E79", linewidth = .45) +
  geom_point(aes(shape = excludes_1), colour = "#1F4E79", fill = "#1F4E79", size = 1.8, stroke = .6) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), labels = c(`TRUE` = "95% interval excludes 1", `FALSE` = "95% interval includes 1"), name = NULL) +
  scale_x_log10(labels = function(x) sub("\\.?0+$", "", sprintf("%.2f", x))) +
  facet_grid(subregion ~ band, scales = "free_y", space = "free_y", switch = "y", labeller = labeller(band = cl)) +
  labs(x = "Country random intercept, exp(u): multiplicative effect on the band hazard (log scale)", y = NULL) +
  theme_re + theme(panel.grid.major.y = element_blank(), strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"), axis.text.x = element_text(size = 8))
ggsave(file.path(out, "sfig_country_random_effects.png"), p1, width = 16, height = 10, dpi = 300, device = ragg::agg_png, bg = "white")

# Survey-region caterpillar plot.
region[, band := factor(age_band, levels = ages)]
region[, rank := frank(effect, ties.method = "first"), by = band]
rl <- facet_lab("region")
cols <- c(`West Africa` = "#C34D26", `Central Africa` = "#7B4FA3", `East Africa` = "#009E73", `Southern Africa` = "#E69F00")
p2 <- ggplot(region, aes(rank, effect, colour = subregion)) +
  geom_hline(yintercept = 1, colour = "grey45", linewidth = .4) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95), colour = "grey70", linewidth = .2, alpha = .3) +
  geom_point(size = .55) +
  scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
  scale_y_log10(labels = function(x) sub("\\.?0+$", "", sprintf("%.2f", x))) +
  facet_wrap(~band, nrow = 2, labeller = labeller(band = rl), scales = "free_x") +
  labs(x = "Survey-regions ranked by estimated random intercept", y = "Survey-region random intercept, exp(u) (log scale)") +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  theme_re + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
ggsave(file.path(out, "sfig_region_random_effects.png"), p2, width = 15, height = 8.5, dpi = 300, device = ragg::agg_png, bg = "white")

common <- paste0("Estimates are exp(u), the multiplicative effect of the random intercept on the band hazard, conditional on the fixed effects (PfPR[2–10] and calendar-year splines and 17 covariates) and the other random intercepts, from the seven separately fitted age-band models of the primary analysis (",
  st$id, "). Intervals are 95% intervals from the coefficient covariance conditional on the estimated smoothing parameters. Panel titles give the effective degrees of freedom of the random-effect term out of its number of levels and the estimated random-effect standard deviation (gam.vcomp); an EDF near zero means the term is shrunk almost entirely to zero, which the gamma = 2 penalty encourages.")
writeLines(c("# Supplementary figure: country random intercepts by age band", "",
  paste("Country random intercepts for the 37 countries in each age-band model, grouped by UN M49 sub-region (country names as in Figure 1).", common)),
  file.path(out, "CAPTION_country.md"))
writeLines(c("# Supplementary figure: survey-region random intercepts by age band", "",
  paste(sprintf("Survey-region random intercepts (%s–%s survey-regions per age band; region identifiers are specific to each survey), ranked within each band and coloured by UN M49 sub-region; grey vertical lines are 95%% intervals.",
    format(min(summ[term == "region"]$levels), big.mark = ","), format(max(summ[term == "region"]$levels), big.mark = ",")), common)),
  file.path(out, "CAPTION_region.md"))
inputs <- c(manifest$model_file, file.path(st$out, "fit_manifest.csv"), "R_cbh/reporting/15_random_effects_forest.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")), file.path(out, "provenance.csv"))
print(summ); message("Random-effect figures written: ", out)

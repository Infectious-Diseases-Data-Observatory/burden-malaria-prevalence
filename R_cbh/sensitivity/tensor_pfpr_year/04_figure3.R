#!/usr/bin/env Rscript
# Figure 3 analogue for the tensor-product sensitivity (added 25 September 2026): the main-text
# rate figure (R_cbh/reporting/11_burden_comparison_figure.R) with the PfPR-ACM values from the
# te(PfPR, calendar year) models (no region random intercept) instead of the primary.
# A = 2024 country rates versus IHME; B = Nigerian states versus IHME, 2024; C = annual under-five
# malaria mortality 2000-2024 versus IHME and UN IGME, with the primary as a reference line.
# Hazard ratios are year-specific, evaluated at mid-year and held at December 2023 for 2024.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/reporting/labels.R")
source("R_cbh/sensitivity/tensor_pfpr_year/components.R")
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
base <- cbh_primary_settings("regional_mics"); root <- base$out
out <- "results/cbh/tensor_pfpr_year_dhsmics_map_gamma2_v1/figure3"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
ages <- cbh_config()$age_bands$age_band
comp <- cbh_tensor_components(); z_of <- function(a) comp[[sprintf("tensor_age_%d", match(a, ages))]]
model_label <- paste(cbh_paper_model_label(), "(tensor product)"); primary_label <- paste(cbh_paper_model_label(), "(primary)")
paths <- c(age = file.path(root, "annual_comparison/country_age_estimates_2000_2024.csv"), country = file.path(root, "annual_comparison/country_estimates_2000_2024.csv"),
  fig5 = file.path(root, "annual_comparison/figure5_data.csv"), state_age = file.path(root, "nigeria_states/state_age_estimates_2024.csv"),
  state_tot = file.path(root, "nigeria_states/state_totals_2024.csv"))
stopifnot(all(file.exists(paths)))

# National deaths, 2000-2024.
age_in <- fread(paths[["age"]]); cty <- fread(paths[["country"]])
nat <- rbindlist(lapply(ages, function(a) { z <- age_in[age_band == a]; k <- cbh_tensor_lhr(z_of(a), z$pfpr_pct, 0, cbh_tensor_eval_year(z$year))
  z[, .(iso3, year, deaths = allcause_deaths * (1 - exp(k$est)))] }))[, .(model_deaths = sum(deaths)), by = .(iso3, year)]
nat <- merge(nat, cty[, .(iso3, year, country, under5_person_years, ihme_malaria_deaths, ihme_malaria_rate_per100000, who_cacode_deaths)], by = c("iso3", "year"))
stopifnot(nrow(nat) == 42L * 25L)
cbh_atomic_csv(as.data.frame(nat[order(iso3, year)]), file.path(out, "country_year_tensor.csv"))

# Nigerian states, 2024 (state IHME all-cause deaths by band and state PfPR).
sa <- fread(paths[["state_age"]]); st <- fread(paths[["state_tot"]]); stopifnot(nrow(st) == 37L, all(sa$year == 2024))
sa[, deaths := {k <- cbh_tensor_lhr(z_of(age_band[1]), pfpr_pct, 0, cbh_tensor_eval_year(2024)); ihme_deaths * (1 - exp(k$est))}, by = age_band]
states <- merge(sa[, .(model_deaths = sum(deaths)), by = state], st[, .(state, under5_person_years, ihme_malaria_deaths, ihme_malaria_rate_per100000, primary_deaths = attributable_under5_deaths)], by = "state")
states[, `:=`(model_rate = 1000 * model_deaths / under5_person_years, ihme_rate = ihme_malaria_rate_per100000 / 100)]
cbh_atomic_csv(as.data.frame(states[order(state)]), file.path(out, "nigeria_states_2024_tensor.csv"))

# Annual series (common IHME-implied denominator, as in the main figure).
f5 <- fread(paths[["fig5"]]); stopifnot(setequal(unique(f5$source), c(cbh_paper_model_label(), "IHME", "UN IGME")))
tens <- nat[, .(deaths = sum(model_deaths), under5_person_years = sum(under5_person_years)), by = year][, `:=`(source = model_label, rate_per100000 = 1e5 * deaths / under5_person_years)]
long <- rbind(f5[source != cbh_paper_model_label(), .(year, source, rate_per100000)], tens[, .(year, source, rate_per100000)],
  f5[source == cbh_paper_model_label(), .(year, source = primary_label, rate_per100000)])
sources <- c(model_label, "IHME", "UN IGME", primary_label)
long[, `:=`(source = factor(source, levels = sources), rate_per1000 = rate_per100000 / 100)]
cbh_atomic_csv(as.data.frame(long[order(source, year)]), file.path(out, "annual_rates_tensor.csv"))

# Figure (styling and labelling rules as in R_cbh/reporting/11_burden_comparison_figure.R).
paper <- theme_minimal(base_size = 14) + theme(panel.grid.minor = element_blank(), text = element_text(size = 16), axis.title = element_text(size = 16),
  axis.text = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 13), legend.position = "bottom",
  plot.tag = element_text(size = 22, face = "bold"), plot.margin = margin(8, 16, 8, 8))
cr <- nat[year == 2024][, `:=`(model_rate = 1000 * model_deaths / under5_person_years, ihme_rate = ihme_malaria_rate_per100000 / 100)]
rate_limits <- function(...) c(0, ceiling(1.05 * max(...)))
lim_a <- rate_limits(cr$model_rate, cr$ihme_rate); lim_b <- rate_limits(states$model_rate, states$ihme_rate)
rA <- ggplot(cr, aes(ihme_rate, model_rate)) + geom_abline(slope = 1, intercept = 0, colour = "grey65", linetype = 2) +
  geom_point(colour = "#215E91", size = 2.6, alpha = .85) +
  ggrepel::geom_text_repel(data = cr[iso3 %in% c("COD", "NGA", "AGO", "UGA", "TZA") | ihme_rate == max(ihme_rate) | abs(log(model_rate / ihme_rate)) >= .5 | abs(model_rate - ihme_rate) >= 1],
    aes(label = iso3), size = 3.8, seed = 20260915, max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey60", box.padding = .35) +
  coord_equal(xlim = lim_a, ylim = lim_a, expand = FALSE) +
  labs(x = "IHME malaria deaths per 1,000 child-years\nbefore age 5, 2024", y = paste0(model_label, "\ndeaths per 1,000 child-years before age 5, 2024")) + paper
zones <- list(`North Central` = c("Benue", "FCT (Abuja)", "Kogi", "Kwara", "Nasarawa", "Niger", "Plateau"), `North East` = c("Adamawa", "Bauchi", "Borno", "Gombe", "Taraba", "Yobe"),
  `North West` = c("Jigawa", "Kaduna", "Kano", "Katsina", "Kebbi", "Sokoto", "Zamfara"), `South East` = c("Abia", "Anambra", "Ebonyi", "Enugu", "Imo"),
  `South South` = c("Akwa Ibom", "Bayelsa", "Cross River", "Delta", "Edo", "Rivers"), `South West` = c("Ekiti", "Lagos", "Ogun", "Ondo", "Osun", "Oyo"))
zone_of <- setNames(rep(names(zones), lengths(zones)), unlist(zones)); stopifnot(setequal(names(zone_of), states$state))
states[, `:=`(zone = factor(zone_of[state], levels = names(zones)), label = sub(" [(]Abuja[)]$", "", state))]
rB <- ggplot(states, aes(ihme_rate, model_rate, colour = zone)) + geom_abline(slope = 1, intercept = 0, colour = "grey65", linetype = 2) + geom_point(size = 2.6) +
  ggrepel::geom_text_repel(aes(label = label), size = 3.3, seed = 20260917, max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey70", show.legend = FALSE, point.padding = .3, box.padding = .4) +
  coord_equal(xlim = lim_b, ylim = lim_b, expand = FALSE) +
  scale_colour_manual(values = c(`North Central` = "#E69F00", `North East` = "#D55E00", `North West` = "#CC79A7", `South East` = "#009E73", `South South` = "#0072B2", `South West` = "#56B4E9")) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3.5))) +
  labs(x = "IHME malaria deaths per 1,000 child-years\nbefore age 5, 2024", y = paste0(model_label, "\ndeaths per 1,000 child-years before age 5, 2024")) + paper
rC <- ggplot(long, aes(year, rate_per1000, colour = source, linetype = source, linewidth = source)) + geom_line() +
  scale_colour_manual(values = setNames(c("#16747C", "#253746", "#CC6A30", "#9AA5AE"), sources)) +
  scale_linetype_manual(values = setNames(c("solid", "longdash", "dotdash", "dotted"), sources)) +
  scale_linewidth_manual(values = setNames(c(1.15, 1.15, 1.15, .9), sources)) +
  scale_x_continuous(breaks = c(2000, 2005, 2010, 2015, 2020, 2024), limits = c(2000, 2024), expand = expansion(mult = c(.015, .025))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2), linewidth = guide_legend(nrow = 2)) +
  labs(x = "Year", y = "Under-five malaria deaths\nper 1,000 child-years") + paper + theme(legend.key.width = grid::unit(1.2, "cm"))
fig <- (rA | rB) / rC + plot_layout(heights = c(1, .7)) + plot_annotation(tag_levels = "A")
ggsave(file.path(out, "fig3_burden_comparison_tensor.png"), fig, width = 14, height = 13, dpi = 240, device = ragg::agg_png, bg = "white")

f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
t24 <- nat[year == 2024, .(m = sum(model_deaths), i = sum(ihme_malaria_deaths), u = sum(who_cacode_deaths))]
pc <- function(y0, y1) 100 * (tens[year == y1]$rate_per100000 / tens[year == y0]$rate_per100000 - 1)
caption <- paste(
  "Tensor-product sensitivity version of Figure 3: malaria mortality before age 5 from the PfPR-ACM model with PfPR[2–10] and calendar time modelled jointly (te(PfPR, year), survey and country random intercepts, no survey-region random intercept), compared with IHME and UN IGME cause-specific malaria estimates, as deaths per 1,000 under-five child-years.",
  sprintf("(A) National rates for 2024 in the %d countries (model %s deaths against IHME %s and UN IGME %s; %d countries above IHME).", nrow(cr), f0(t24$m), f0(t24$i), f0(t24$u), cr[, sum(model_deaths > ihme_malaria_deaths)]),
  sprintf("(B) The 36 Nigerian states and the Federal Capital Territory in 2024 (model %s deaths summed over states against %s IHME malaria deaths; %d of 37 states above IHME).", f0(sum(states$model_deaths)), f0(sum(states$ihme_malaria_deaths)), states[, sum(model_deaths > ihme_malaria_deaths)]),
  sprintf("(C) Annual under-five malaria mortality, 2000–2024, pooled across the same 42 countries on the common IHME-implied person-year denominator; the primary (time-constant) PfPR-ACM model is shown dotted for reference. With the tensor product the modelled rate falls %.1f%% from 2000 to 2015 and %.1f%% from 2015 to 2024.", -pc(2000, 2015), -pc(2015, 2024)),
  "Deaths equal IHME all-cause deaths in each of seven age bands multiplied by 1 − exp[f_g(0, t) − f_g(P, t)], where P is population-weighted MAP PfPR for the country or state and year and t is calendar time at mid-year, held at December 2023 (the last observed band entry) for 2024. The time-varying hazard ratio is an exploratory specification: its steepening at low PfPR at 24–59 months is not separable from between-survey differences (see the sensitivity report). Point estimates only.")
writeLines(c("# Figure 3, tensor-product sensitivity: caption", "", caption), file.path(out, "CAPTION.md"))
inputs <- c(unname(paths), "data/derived_cbh/models/tensor_pfpr_year_dhsmics_map_gamma2_v1/pfpr_year_components.rds", "R_cbh/sensitivity/tensor_pfpr_year/04_figure3.R", "R_cbh/sensitivity/tensor_pfpr_year/components.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")), file.path(out, "provenance.csv"))
print(t24); message("Tensor Figure 3 written: ", file.path(out, "fig3_burden_comparison_tensor.png"))

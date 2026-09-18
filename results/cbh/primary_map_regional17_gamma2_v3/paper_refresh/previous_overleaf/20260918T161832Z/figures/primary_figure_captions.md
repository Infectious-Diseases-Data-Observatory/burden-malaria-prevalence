# Main-paper figures and draft captions

Figures contain no title, subtitle or embedded caption. Axis labels, age/year facet labels and legends remain. File numbering here does not alter the order of figures in TeX.

## Figure 1 — survey map and timing

File: `survey_map_and_timing.png`.

Geographic coverage and timing of the 95 surveys in 34 countries contributing to the primary MAP analysis. Country shading indicates the number of included surveys. Timeline point area indicates the number of regions contributing analysis records in each survey; survey types are read from the registry. Region counts are the union across the seven fitted age-band samples. The figure contains 916 survey–region pairs; these are not counts of geographically distinct regions across survey years.

All surveys in the current selected sample are DHS surveys. Country outlines are dissolved from the most recent available DHS boundary files among the included surveys. Survey year describes fieldwork, not the calendar year assigned to each child's band entry. The included surveys and region counts follow the current primary complete-case sample after the declared covariate substitutions.

## Figure 2 — PfPR and mortality by age

File: `fig2_pfpr_mortality_by_age.png`.

PfPR-ACM model: age-specific associations between MAP PfPR[2–10] and all-cause mortality, estimated in seven separate age-band models at gamma=2. Curves show log mortality hazard ratios relative to PfPR=20% across the central 95% of each band's exposure distribution. Shading denotes conditional pointwise 95% intervals. Models use cr splines for PfPR (k=5) and calendar year (k=6), 17 survey-region/national annual covariates, including maternal age at first birth, urban percentage, DTP3/measles, delivery, birth interval, WASH, electricity, wasting and stunting, and survey/country/region random effects. Intervals condition on fitted smoothing parameters, exposure, the fixed posterior-median child HIV imputation and filled regional covariates. Sex, multiple births, birth order, Hib3, PCV, rotavirus and exclusive breastfeeding are excluded.

## Figure 3 — malaria-attributable fraction by age

File: `fig3_malaria_attributable_fraction_by_age.png`.

PfPR-ACM model: estimated share of all-cause deaths attributable to malaria within each completed-month age band, at MAP PfPR[2–10] levels of 10%, 20%, 30% and 40%. Points show 100 × {1 − exp[f_g(0) − f_g(P)]}, using the seven separate primary MAP gamma=2 models; lines connect discrete age-band estimates and do not imply a continuous-age fit. The denominator is all-cause deaths within each band, not all under-five deaths pooled across bands. The hazard-based attributable fraction follows the primary burden convention of holding annual person-time fixed; no country-specific IHME inputs are used here. Zero prevalence lies below observed exposure support, so the counterfactual requires extrapolation. Point estimates are plotted; conditional pointwise intervals and support flags are retained in attributable_fraction_by_age.csv.

## Figure 4 — comparison with IHME and UN IGME: countries (A), Nigerian states (B) and annual trend (C)

File: `fig4_burden_comparison.png`. Panel sources: `burden/country_totals.csv` (2024), `nigeria_states/state_totals_2024.csv`, `annual_comparison/figure5_data.csv`. The separate three-year country scatter (`burden/country_vs_ihme.png`), state figure (`nigeria_states/`) and annual figure (`annual_comparison/`) remain saved results with their own captions but are no longer paper figures.

Malaria-attributable deaths before age 5 from the PfPR-ACM model compared with IHME and UN IGME cause-specific malaria death estimates. (A) National estimates for 2024 in the 42 countries with both estimates available; each point is a country and the dashed line denotes equality. Both axes use a base-10 logarithmic scale from 100 to 250,000 deaths; 6 countries with fewer than 100 deaths on either axis lie outside the displayed range and are retained in the source table. (B) The same comparison for the 36 Nigerian states and the Federal Capital Territory in 2024, with axes from 200 to 30,000 deaths; colours group states into Nigeria's six geopolitical zones for orientation only. Summed over states the model gives 200,479 deaths against 132,138 IHME malaria deaths (ratio 1.52). (C) Annual under-five malaria mortality, 2004–2024, pooled across the same 42 countries covered by the national burden analysis. For each source the rate is the sum of national under-five malaria deaths divided by the same under-five person-years implied by the IHME all-cause death counts and rates, so the three series share one denominator. The UN IGME series is the CA-CODE 2026 release. In all panels, PfPR-ACM model deaths equal IHME all-cause deaths in each of seven age bands multiplied by 1 − exp[f_g(0) − f_g(P)], where P is population-weighted MAP PfPR[2–10] for the country or state and year and f_g are the seven fitted age-band effects; IHME early and late neonatal deaths share the <1-month effect and the IHME 2–4 year group is divided equally across the three oldest bands. All values are point estimates; joint uncertainty across countries and age bands is not available and marginal intervals are not summed. Model-attributable all-cause reductions and cause-specific malaria deaths are different estimands, and the model shares IHME all-cause inputs, so agreement is not independent validation.

[Source paths and hashes](manifest.csv). Copy these four images into the manuscript's `figures/` folder. The former fig4_country_estimates_vs_ihme.png, fig5_annual_malaria_mortality.png and fig6_nigeria_states_vs_ihme.png are superseded by the combined Figure 4; existing TeX references to them need updating by the author. No TeX edits are needed to prepare or copy the images; TeX references/captions are managed separately.

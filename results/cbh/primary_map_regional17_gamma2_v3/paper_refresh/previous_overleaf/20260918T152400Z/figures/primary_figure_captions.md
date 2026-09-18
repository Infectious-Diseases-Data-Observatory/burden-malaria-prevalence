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

## Figure 4 — country mortality estimates versus IHME

File: `fig4_country_estimates_vs_ihme.png`.

PfPR-ACM model malaria-attributable deaths before age 5 versus IHME cause-specific malaria deaths, by country in 2005, 2015 and 2024. Both axes use a base-10 logarithmic scale with identical limits starting at 1,000 deaths; each point represents a country and the dashed line denotes equality. The underlying estimates include all 42 estimable countries per year, with no pseudocount; points below 1,000 deaths on either axis are outside the displayed range. The three countries without model estimates remain in the source tables. Model-attributable deaths equal national IHME all-cause deaths multiplied by 1 − exp[f_g(0) − f_g(P_country)], summed over age bands. National MAP prevalence is population weighted using the saved 2020 population weights. Ages 2–4 share the IHME mortality rate and divide deaths/person-time equally. Country totals are point estimates; cross-age covariance and source/design/imputation uncertainty are not propagated. Model-attributable all-cause reductions and IHME cause-specific mortality are different estimands.

## Figure 5 — annual under-five malaria mortality

File: `fig5_annual_malaria_mortality.png`.

Annual malaria mortality before age five, 2004–2024, pooled across the same 42 countries covered by the national burden analysis. For each source, the plotted rate is 100,000 times the sum of national under-five malaria deaths divided by the sum of annual under-five person-years implied by the primary IHME all-cause death counts and rates. Rates therefore use a common population denominator and are not averages of country rates or probabilities per live birth. The PfPR-ACM model applies the seven separate primary MAP gamma=2 PfPR effects to annual age-specific IHME all-cause deaths, using the zero-PfPR counterfactual. National annual MAP PfPR[2–10] is weighted by population counts using a fixed GPW 2020 spatial distribution; these exposure weights are separate from the annual mortality denominator. IHME ages 2–4 share a mortality rate and divide deaths/person-time equally across the model's three annual bands. Signed age-specific contributions are retained. The UN IGME comparator is the CA-CODE 2026 series from the UN IGME portal, retrieved from UNICEF's CME_CAUSE_OF_DEATH dataflow on 17 September 2026 (under five, both sexes, malaria, deaths); it is not the World Malaria Report all-age series or a fixed proportion of it. All curves show point estimates. Joint uncertainty across countries and age bands is not available; marginal interval endpoints are not summed. Model-attributable all-cause mortality reductions and the comparator cause-specific estimates are different estimands. The PfPR-ACM model shares IHME all-cause inputs and this comparison is not independent validation. Zero exposure requires extrapolation. Cape Verde, Lesotho and São Tomé and Príncipe are excluded because usable MAP prevalence is unavailable; the other countries absent from the original national input set are not added.

## Figure 6 — Nigerian state estimates versus IHME

File: `fig6_nigeria_states_vs_ihme.png`.

PfPR-ACM model malaria-attributable deaths before age 5 versus IHME cause-specific malaria deaths, by Nigerian state (36 states and the Federal Capital Territory) in 2024. Both axes use a base-10 logarithmic scale with identical limits; the dashed line denotes equality. Colours group states into Nigeria's six geopolitical zones for orientation; all estimates are made state by state. Model-attributable deaths equal IHME state all-cause deaths in each of seven age bands multiplied by 1 − exp[f_g(0) − f_g(P_state)], summed over bands, where P_state is population-weighted MAP PfPR[2–10] for 2024 and f_g are the seven separate primary MAP gamma=2 age-band effects. IHME early and late neonatal deaths share the <1-month effect; the IHME 2–4 year group shares one rate and divides deaths and person-time equally across the 24–35, 36–47 and 48–59 month bands. IHME malaria deaths are the state malaria death rate applied to the under-5 person-years implied by the all-cause count and rate. Summed over states the model gives 200,479 deaths against 132,138 IHME malaria deaths (ratio 1.52); the state sum differs from the national estimate evaluated at national mean prevalence by -0.3%, because the age-band curves are nonlinear. State totals are point estimates; conditional age-band intervals are retained in the source table and are not summed. Model-attributable all-cause reductions and IHME cause-specific malaria deaths are different estimands, and the model shares IHME all-cause inputs, so agreement or disagreement is not independent validation. Zero prevalence lies below observed exposure support.

[Source paths and hashes](manifest.csv). Copy these six images into the manuscript's `figures/` folder. The former fig3_country_estimates_vs_ihme.png is superseded by the Figure 4 filename; retain a legacy copy if existing TeX references still need it. No TeX edits are needed to prepare or copy the images; TeX references/captions are managed separately.

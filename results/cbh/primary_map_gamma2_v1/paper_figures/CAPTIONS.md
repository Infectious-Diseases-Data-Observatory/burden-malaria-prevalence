# Main-paper figures and draft captions

Figures contain no title, subtitle or embedded caption. Axis labels, age/year facet labels and legends remain. File numbering here does not alter the order of figures in TeX.

## Figure 1 — survey map and timing

File: `survey_map_and_timing.png`.

Geographic coverage and timing of the 105 surveys in 34 countries contributing to the primary MAP analysis. Country shading indicates the number of included surveys. Timeline point area indicates the number of regions contributing analysis records in each survey; symbols distinguish DHS and MIS. Region counts are the union across the seven fitted age-band samples. The figure contains 1015 survey–region pairs; these are not counts of geographically distinct regions across survey years.

Country outlines are dissolved from the most recent available DHS boundary files among the included surveys. Survey year describes fieldwork, not the calendar year assigned to each child's band entry. The historical figure showed 120 surveys/36 countries; this version shows the 105-survey/34-country primary complete-case sample.

## Figure 2 — PfPR and mortality by age

File: `fig2_pfpr_mortality_by_age.png`.

Age-specific associations between MAP PfPR[2–10] and all-cause mortality, estimated in seven separate age-band models at gamma=2. Curves show log mortality hazard ratios relative to PfPR=20% across the central 95% of each band's exposure distribution. Shading denotes conditional pointwise 95% intervals. Models use cr splines for PfPR (k=5) and calendar year (k=6), the prespecified adjustment set, and survey/country/region random effects. Intervals condition on fitted smoothing parameters, exposure and the fixed posterior-median child HIV imputation.

## Figure 3 — country mortality estimates versus IHME

File: `fig3_country_estimates_vs_ihme.png`.

Primary-model malaria-attributable deaths before age 5 versus IHME cause-specific malaria deaths, by country in 2005, 2015 and 2024. Both axes use a base-10 logarithmic scale with identical limits; each point represents a country and the dashed line denotes equality. All 42 estimable countries per year are included, with no pseudocount. The three countries without model estimates remain in the source tables. Model-attributable deaths equal national IHME all-cause deaths multiplied by 1 − exp[f_g(0) − f_g(P_country)], summed over age bands. National MAP prevalence is population weighted using the saved 2020 population weights. Ages 2–4 share the IHME mortality rate and divide deaths/person-time equally. Country totals are point estimates; cross-age covariance and source/design/imputation uncertainty are not propagated. Model-attributable all-cause reductions and IHME cause-specific mortality are different estimands.

[Source paths and hashes](manifest.csv). Copy only these three images into the manuscript's `figures/` folder. No TeX edits are needed to prepare or copy the images; TeX references/captions are managed separately.

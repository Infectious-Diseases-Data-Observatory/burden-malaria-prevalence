# Updated primary results

Figures 1–4 (Figure 4 combines the 2024 country comparison, the Nigerian state
comparison and the 2004–2024 annual trend as panels A–C, all as deaths per 1,000
under-five child-years; its death-count version is a supplementary figure) and the
supplementary inclusion flow use the 17-variable
regional-adjustment MAP gamma=2 models in `primary_map_regional17_gamma2_v3`. The analysis contains
1,686,004 children, 5,465,305 child-band records and
75,726 deaths in 916 survey-regions,
95 surveys and 34 countries.

The revised adjustment removes sex, multiple births and birth order; replaces
maternal age at each birth with regional mean age at first birth; and adds
regional wasting and stunting prevalence. All retained DHS predictors are regional
summaries; the four existing national annual variables remain.

No TeX file has been edited. Update manuscript text and captions manually using:

- [Age-band table source](tables/age_band_results.latex.txt), with death
  percentages summing to 100.0% and both prevalence contrasts.
- [Country table source](tables/country_comparison_2024.latex.txt).
- [Main figure captions](figures/primary_figure_captions.md).
- [Inclusion-flow caption](<Supplementary Figures/sfig_study_flow_caption.md>).
- [Annual totals](tables/annual_totals_2004_2024.csv) and
  [Nigerian state totals](tables/nigeria_state_totals_2024.csv).

The 2024 totals are 565,616 PfPR-ACM deaths,
428,147 IHME deaths and
440,123 UN IGME deaths across the same 42 countries.
These are conditional model-based estimates; source and imputation uncertainty
are not fully propagated. Sensitivity fits remain historical and are not
relabelled as this revised primary specification.

Previous figure versions are backed up in the analysis project. `.latex.txt`
files contain source for manual insertion; they do not change the compiled
manuscript automatically. Figure placement and TeX references remain under the
author's control.

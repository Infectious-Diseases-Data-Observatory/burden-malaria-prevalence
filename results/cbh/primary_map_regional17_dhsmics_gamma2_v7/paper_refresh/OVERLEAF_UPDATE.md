# Updated primary results (DHS and MICS surveys)

Figures 1–4 (Figure 3 combines the 2024 country comparison, the Nigerian state
comparison and the 2000–2024 annual trend as panels A–C, all as deaths per 1,000
under-five child-years, with its death-count version as a supplementary figure;
Figure 4 is the probability of dying before age 5 by country, all causes and with
malaria transmission removed) and the
supplementary inclusion flow use the 17-variable
regional-adjustment MAP gamma=2 models in `primary_map_regional17_dhsmics_gamma2_v7`. The analysis contains
2,325,130 children, 7,607,122 child-band records and
103,987 deaths in 1,227 survey-regions,
135 surveys and 37 countries.

It combines 98 DHS and 37 UNICEF MICS surveys. The DHS part is the DHS-only complete-case sample (primary_map_regional17_gamma2_v3) plus Liberia's 2007, 2013 and 2019-20 DHS (5,573,968 records, 77,431 deaths, 932 survey-regions, 98 surveys, 35 countries); the MICS surveys add 2,033,154 records, 26,556 deaths, 295 survey-regions and 2 countries (CAF, GNB). None of the 5 MIS surveys passes complete-case selection.

The revised adjustment removes sex, multiple births and birth order; replaces
maternal age at each birth with regional mean age at first birth; and adds
regional wasting and stunting prevalence. All 13 survey-level predictors are survey-region summaries: for DHS surveys, published regional indicators or weighted recode summaries; for MICS surveys, values computed from the microdata with the DHS definitions and the same WUENIC and within-survey fallbacks. The four national annual
variables are unchanged.

No TeX file has been edited. Update manuscript text and captions manually using:

- [Manuscript text updates](MANUSCRIPT_TEXT_UPDATES.md) for the DHS and MICS
  sample (copied from `docs/MANUSCRIPT_UPDATE_DHS_MICS.md` in the analysis project).
- [Age-band table source](tables/age_band_results.latex.txt), with death
  percentages summing to 100.0% and both prevalence contrasts.
- [Country table source](tables/country_comparison_2024.latex.txt).
- [Source comparison table source](tables/source_comparison.latex.txt) and its
  [values](tables/source_comparison_2000_2024.csv).
- [Main figure captions](figures/primary_figure_captions.md) and the
  [Figure 3 caption](figures/fig3_burden_comparison_caption.md); Figure 3 is
  `figures/fig3_burden_comparison.png`.
- [Inclusion-flow caption](<Supplementary Figures/sfig_study_flow_caption.md>).
- [Annual totals](tables/annual_totals_2000_2024.csv) and
  [Nigerian state totals](tables/nigeria_state_totals_2024.csv).

The 2024 totals are 620,072 PfPR-ACM deaths,
428,147 IHME deaths and
440,123 UN IGME deaths, summed over the
42 countries in the national burden comparison (36 of the 37 survey countries; ZAF contributes survey data but is outside the burden set).
These are conditional model-based estimates; source and imputation uncertainty
are not fully propagated.

The supplementary sensitivity figures are refitted on the DHS and MICS surveys:
subgroups (`subgroups_dhsmics_map_gamma2_v3`) and no nutrition covariates (`nutrition_adjustment_dhsmics_map_gamma2_v3`)
on this sample, and imputed covariates on the larger MAP-eligible sample
(`primary_map_regional17_dhsmics_imputed_gamma2_v8`, which retains every MAP-eligible record: 8,797,963 records, 123,419 deaths, 166 surveys and 40 countries). Their captions are exported to
`Supplementary Figures/`. The DHS-only sensitivity fits are kept in the analysis
project as history.

Previous figure versions are backed up in the analysis project. `.latex.txt`
files contain source for manual insertion; they do not change the compiled
manuscript automatically. Figure placement and TeX references remain under the
author's control.

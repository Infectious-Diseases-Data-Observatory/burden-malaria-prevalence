# Primary study flow

Run from the repository root:

```sh
Rscript R_dhs/11_study_flow.R
```

The original entry point now calls `R_cbh/reporting/01_study_flow.R`. It reads the aggregate birth-history, eligibility and model-selection ledgers of both the DHS build (`data/derived_cbh/`, with the DHS exclusion attribution named by `settings$attribution_dir`; for v7 `results/cbh/planned17_covariate_missingness_lbr/`, written by `Rscript R_cbh/covariates/14_exclusion_attribution.R --liberia`) and the MICS build (`data/derived_mics/cbh_build/`, with the MICS registry, inventory and `results/mics_inventory/exclusion_attribution_mics_by_survey.csv`), and checks them against the seven full-sample MAP `gamma=2` fits. It does not refit models or read individual birth records. Every reporting script defaults to the current primary version (`CBH_PRIMARY_VERSION=regional_mics`, which is v7 since 24 September 2026; `regional_mics_v5` is the v5 history).

Outputs are in [`results/cbh/primary_map_regional17_dhsmics_gamma2_v7/study_flow/`](../../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/study_flow/): the PNG, separate caption, counts, survey exclusions and input hashes. The 300-dpi figure uses approximately 16–19-point text at its native 11-inch width. Eligibility definitions and modelling details are kept in the caption.

The checks also reconcile the current primary sample summary and survey coverage, and verify the prepared-data hash against the selected fit manifest. The registry holds 124 DHS/MIS and 93 MICS surveys; 46 MICS surveys have no complete birth history and 5 Lesotho surveys have no MAP surface, leaving 166 processed surveys in 40 countries. The final v7 sample is 2,325,130 distinct children contributing 7,607,122 records and 103,987 deaths across 135 surveys (98 DHS, 37 MICS) in 37 countries. The flows of v5 (2,292,089 children, 7,498,459 records, 102,282 deaths, 132 surveys in 36 countries) and of the DHS-only `primary_map_regional17_gamma2_v3` are preserved in those versions' `study_flow/`. The caption distinguishes unique children from child–age-band records and from all recorded births. The diagram shows inclusion only; no model or downstream-results boxes are added.

The supplementary manuscript image is `Supplementary Figures/sfig_study_flow.png` in the Overleaf project. Copying it is a separate action, never a TeX edit; the separate caption is saved alongside it as `sfig_study_flow_caption.md`.

`Rscript R_cbh/reporting/02_results_index.R` (also available through `R_dhs/33_key_results.R`) indexes authoritative MAP gamma=2 outputs and separates historical/supplementary figures. It does not copy older figures over current results.

Manuscript model naming: source `R_cbh/reporting/labels.R` and use `cbh_paper_model_label()` (currently “PfPR-ACM model”) for every legend, axis or caption that names the proposed method. This applies to main and supplementary figures; see the naming rule in `docs/ANALYSIS_PLAN.md`. Keep outcome-only labels and omit figure titles/subtitles.

The current runner is `Rscript R_cbh/primary/run_regional.R --report-only`. It regenerates only primary reporting defined in the analysis plan. Tables are written as CSV/Markdown and `.latex.txt` source, with no TeX files changed.

`Rscript R_cbh/reporting/13_under5_death_probability.R` (runner stage `u5_probability`, before `paper_manifest`; added 24 September 2026) draws Figure 4: the probability of dying before age 5 in 2024 by country, per 1,000 live births. All-cause values are synthetic-cohort probabilities from IHME age-band all-cause death rates, 5q0 = 1 − exp(−Σ H_g) over the seven bands (the neonatal band from the IHME early and late neonatal rates); the part caused by malaria is 5q0 minus 5q0 with each band hazard multiplied by the model's hazard ratio for PfPR 0 versus the country's current PfPR, with intervals from independent draws of each band's log hazard ratio. Outputs are in `<primary>/under5_probability/` (`under5_death_probability.csv`, `under5_death_probability_2024.png`, `CAPTION.md`, `provenance.csv`); `04_paper_figures.R` lists it as the fourth main figure and `export_paper.py` copies it to Overleaf as `figures/fig4_under5_death_probability.png`. Against UN IGME 2024, the IHME-based all-cause values have median ratio 1.07 and correlation 0.84.

Figure 1 (`03_survey_map.R`) labels timeline rows with full country names (DRC, CAF and RC abbreviated), grouped by UN M49 region, alphabetical within region, with the number of mothers (n) per country.

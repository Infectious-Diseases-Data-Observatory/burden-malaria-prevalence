# Primary study flow

Run from the repository root:

```sh
Rscript R_dhs/11_study_flow.R
```

The original entry point now calls `R_cbh/reporting/01_study_flow.R`. It reads the aggregate birth-history, eligibility and model-selection ledgers of both the DHS build (`data/derived_cbh/`) and the MICS build (`data/derived_mics/cbh_build/`, with the MICS registry, inventory and `results/mics_inventory/exclusion_attribution_mics_by_survey.csv`), and checks them against the seven full-sample MAP `gamma=2` fits. It does not refit models or read individual birth records. Every reporting script defaults to the current primary version (`CBH_PRIMARY_VERSION=regional_mics`).

Outputs are in [`results/cbh/primary_map_regional17_dhsmics_gamma2_v5/study_flow/`](../../results/cbh/primary_map_regional17_dhsmics_gamma2_v5/study_flow/): the PNG, separate caption, counts, survey exclusions and input hashes. The 300-dpi figure uses approximately 16–19-point text at its native 11-inch width. Eligibility definitions and modelling details are kept in the caption.

The checks also reconcile the current primary sample summary and survey coverage, and verify the prepared-data hash against the selected fit manifest. The registry holds 124 DHS/MIS and 93 MICS surveys; 46 MICS surveys have no complete birth history and 5 Lesotho surveys have no MAP surface, leaving 166 processed surveys in 40 countries. The final sample is 2,292,089 distinct children contributing 7,498,459 records and 102,282 deaths across 132 surveys (95 DHS, 37 MICS) in 36 countries. The DHS-only flow of `primary_map_regional17_gamma2_v3` is preserved in that version's `study_flow/`. The caption distinguishes unique children from child–age-band records and from all recorded births. The diagram shows inclusion only; no model or downstream-results boxes are added.

The supplementary manuscript image is `Supplementary Figures/sfig_study_flow.png` in the Overleaf project. Copying it is a separate action, never a TeX edit; the separate caption is saved alongside it as `sfig_study_flow_caption.md`.

`Rscript R_cbh/reporting/02_results_index.R` (also available through `R_dhs/33_key_results.R`) indexes authoritative MAP gamma=2 outputs and separates historical/supplementary figures. It does not copy older figures over current results.

Manuscript model naming: source `R_cbh/reporting/labels.R` and use `cbh_paper_model_label()` (currently “PfPR-ACM model”) for every legend, axis or caption that names the proposed method. This applies to main and supplementary figures; see the naming rule in `docs/ANALYSIS_PLAN.md`. Keep outcome-only labels and omit figure titles/subtitles.

The current runner is `Rscript R_cbh/primary/run_regional.R --report-only`. It regenerates only primary reporting defined in the analysis plan. Tables are written as CSV/Markdown and `.latex.txt` source, with no TeX files changed.

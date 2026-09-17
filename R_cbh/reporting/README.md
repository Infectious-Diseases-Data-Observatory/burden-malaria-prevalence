# Primary study flow

Run from the repository root:

```sh
Rscript R_dhs/11_study_flow.R
```

The original entry point now calls `R_cbh/reporting/01_study_flow.R`. It reads aggregate birth-history, eligibility and model-selection ledgers and checks them against the seven full-sample MAP `gamma=2` fits. It does not refit models or read individual birth records.

Outputs are in [`results/cbh/primary_map_gamma2_v1/study_flow/`](../../results/cbh/primary_map_gamma2_v1/study_flow/): the PNG, separate caption, counts, survey exclusions and input hashes. The 300-dpi figure uses approximately 16–19-point text at its native 11-inch width. Eligibility definitions and modelling details are kept in the caption.

The checks also reconcile the current primary sample summary and survey coverage, and verify the prepared-data hash against the selected fit manifest. The final sample is 1,817,912 distinct children contributing 5,885,022 records and 82,415 deaths across 105 surveys in 34 countries. The caption distinguishes unique children from child–age-band records and from all recorded births. The diagram shows inclusion only; no model or downstream-results boxes are added.

The supplementary manuscript image is `Supplementary Figures/sfig_study_flow.png` in the Overleaf project. Copying it is a separate action, never a TeX edit; the separate caption is saved alongside it as `sfig_study_flow_caption.md`.

`Rscript R_cbh/reporting/02_results_index.R` (also available through `R_dhs/33_key_results.R`) indexes authoritative MAP gamma=2 outputs and separates historical/supplementary figures. It does not copy older figures over current results.

Manuscript model naming: source `R_cbh/reporting/labels.R` and use `cbh_paper_model_label()` (currently “PfPR-ACM model”) for every legend, axis or caption that names the proposed method. This applies to main and supplementary figures; see the naming rule in `docs/ANALYSIS_PLAN.md`. Keep outcome-only labels and omit figure titles/subtitles.

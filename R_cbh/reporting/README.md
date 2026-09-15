# Primary study flow

Run from the repository root:

```sh
Rscript R_dhs/11_study_flow.R
```

The original entry point now calls `R_cbh/reporting/01_study_flow.R`. It reads aggregate birth-history, eligibility and model-selection ledgers and checks them against the seven full-sample MAP `gamma=2` fits. It does not refit models or read individual birth records.

Outputs are in [`results/cbh/map_snow_gamma2_v1/study_flow/`](../../results/cbh/map_snow_gamma2_v1/study_flow/): the PNG, separate caption, counts, survey exclusions and input hashes. The 300-dpi figure uses approximately 16–19-point text at its native 11-inch width. Eligibility definitions and modelling details are kept in the caption.

`Rscript R_cbh/reporting/02_results_index.R` (also available through `R_dhs/33_key_results.R`) indexes authoritative MAP gamma=2 outputs and separates historical/supplementary figures. It does not copy older figures over current results.

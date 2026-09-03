# Audit of the DHS household covariates (improved water, improved sanitation), 2026-09-03

`imp_water`, `imp_sanit` and `elec_dhs` in the analysis dataset are household-deduplicated,
design-weighted shares of births-recode value labels (v113, v116, v119) classified by a
regular expression in `R_dhs/03_build_analysis_dataset.R` (`household_fraction`). This
folder holds an audit of that rule against every survey's actual labels and against DHS
StatCompiler's own indicators (`WS_SRCE_H_IMP`, `WS_TLET_H_IMP`).

- `audit_household_labels.R` tabulates every survey's labels (`v113_labels_by_survey.csv`,
  `v116_v119_labels_by_survey.csv`, distinct-label tables) and recomputes region shares.
- `fetch_statcompiler_wash.R` pulls the StatCompiler indicators (`dhs_api_water_*.csv`);
  `compare_*.csv` set them against the pipeline; `wash_api_covariates.csv` is a
  StatCompiler-based version of both covariates for every survey-region (939 regional
  matches, 189 national fallbacks).
- `label_classes.json` is a classification of all 303 distinct text labels against JMP
  definitions (LLM judge panel with adversarial verification); `apply_label_classes.R`
  applies it survey by survey (`label_classes_*_by_survey.csv`).
- `wash_sensitivity_refit.R` refits the six age-band person-time models with the
  StatCompiler covariates (`wash_sensitivity.csv`): PfPR slopes move by <= 0.3 points.

Findings: two surveys (MW4EFL, CG51FL) carry bare numeric codes for water and four
(SN4AFL, SN5AFL, CD51FL, ST51FL) for sanitation, giving shares near zero; DHS-IV era
vocabulary ("forage", "covered well", "manual pumped water", location-style piped labels)
and negations ("non protected") are misread; "no slab" is not excluded for sanitation
(Uganda 86-92% against DHS 27-35%); "not de jure resident" sits in the denominator; and
delivered water is unimproved in the rule but improved in StatCompiler. Nothing in the
pipeline has been changed.

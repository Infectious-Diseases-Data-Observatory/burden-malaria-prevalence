# Key results

Curated figures from the DHS/MIS malaria prevalence and child mortality analysis. Every file here is a copy of a figure in `results/dhs_rebuild/`; the source script and the time the figure was generated are given for each. Regenerate with `Rscript R_dhs/33_key_results.R` after re-running the producing scripts (or `Rscript R_dhs/run_all.R`, which does both).

The primary analysis uses the 12-month mortality window and MAP prevalence in the survey year. Post-neonatal mortality is the primary outcome; neonatal mortality is the negative control.

| File | What it shows | Source | Generated |
|---|---|---|---|
| `subgroup_curves_postneonatal.png` | Dose-response curves for post-neonatal mortality: the malaria-attributable fraction against MAP PfPR2-10 by era (before / from 2013), region (West / Central & East) and the 2 x 2 cells. Bayesian (brms) fits with 95% credible intervals; each curve covers only the prevalence range its own subgroup observes. | `31_brms_subgroup_fits.R` (`figure15_brms_subgroup_curves.png`) | 2026-09-02 07:55 |
| `subgroup_af_10_30_50_postneonatal.png` | Post-neonatal attributable fraction at 10, 30 and 50% prevalence in each subgroup, with 95% credible intervals and the number of surveys and survey regions behind each estimate. | `31_brms_subgroup_fits.R` (`figure13_brms_subgroup_af.png`) | 2026-09-02 07:55 |
| `subgroup_curves_neonatal.png` | The same dose-response curves for neonatal mortality, the negative control. Not clipped at zero, so a null effect can sit on the zero line. | `31_brms_subgroup_fits.R with BRMS_OUTCOME=neonatal` (`figure15_brms_subgroup_curves_neonatal.png`) | 2026-09-02 07:55 |
| `subgroup_af_10_30_50_neonatal.png` | The same 10, 30 and 50% anchors for neonatal mortality, the negative control. | `31_brms_subgroup_fits.R with BRMS_OUTCOME=neonatal` (`figure13_brms_subgroup_af_neonatal.png`) | 2026-09-02 07:55 |
| `survey_map_and_timing.png` | Countries contributing to the panel and when each DHS/MIS survey was fielded; point size gives the number of survey regions. | `28_survey_map.R` (`figure10_survey_map.png`) | 2026-09-02 07:55 |
| `mortality_12_vs_60_month_window.png` | Region-level mortality estimated over the 12-month window (the primary analysis) against the 60-month window, for post-neonatal and neonatal mortality, with 95% delete-one-cluster jackknife intervals. | `27_horizon_lag_selection.R` (`figure8b_horizon_12_vs_60.png`) | 2026-09-02 07:55 |
| `map_vs_measured_prevalence.png` | MAP modelled PfPR2-10 against parasitaemia measured directly in the DHS/MIS surveys that tested children 6-59 months (RDT, or microscopy where available, converted to a microscopy-equivalent and age-standardised to 2-10 years), by era of fieldwork. | `21_map_vs_measured_prevalence.R` (`map_vs_measured_scatter.png`) | 2026-09-02 07:54 |

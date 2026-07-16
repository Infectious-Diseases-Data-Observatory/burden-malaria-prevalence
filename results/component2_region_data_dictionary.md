# Data dictionary — `component2_region_data.csv`

The **Component 2 analysis dataset**: one row per **DHS/MIS survey-region**
(600 rows; 44 surveys across 23 sub-Saharan African countries). This is the
model-ready table read by [`R/04_component2_dhs.R`](../R/04_component2_dhs.R)
(mixed models, GAM count model, exposure-response spline) and by
[`R/05_prediction_10_to_30.R`](../R/05_prediction_10_to_30.R).

**Grain:** a survey-region = one first-level administrative region (`v024`/`hv024`;
`sstate`/`shstate` for Nigeria) within one survey. Regions are matched between the
household/malaria recode (PR, prevalence + covariates) and the birth recode (BR,
mortality) on `svkey` + `regkey`.

**How it is built.** Prevalence and region covariates come from the PR recode via
`prev_cov_by_region()`; all-cause mortality from the BR recode via
`DHS.rates::chmort()` inside `mort_by_region()`; GDP and DTP3 are joined from World
Bank series at each survey's country and nearest year; the exposure is put on a common
footing (microscopy-equivalent → age-standardised PfPR₂₋₁₀). Rows are retained only
where `u5mr` is finite and `> 5`, `m1mo5y` is finite and `> 0`, and `pfpr2_10` is finite.

## Columns (in file order)

| # | Column | Type | Units | Description & derivation |
|---|--------|------|-------|--------------------------|
| 1 | `svkey` | chr | — | Survey key linking PR↔BR: 2-letter DHS country code + version (e.g. `AO62FL`). |
| 2 | `regkey` | chr | — | Normalised region key (lower-case alphanumerics of the region label) used for the PR↔BR join. |
| 3 | `region` | chr | — | Region label as recorded in the survey (`v024`/`hv024`, or `sstate`/`shstate` for Nigeria). |
| 4 | `rdt` | num | % | RDT parasite prevalence in tested children ~6–59 mo (`hml35` positive), design-weighted by `hv005` within the region. |
| 5 | `mic` | num | % | Microscopy parasite prevalence (`hml32` positive), same weighting. `NA` where microscopy was not measured in that survey. |
| 6 | `pct_urban` | num | % | Share of the region's PR sample in urban clusters (`hv025` = urban), design-weighted. |
| 7 | `stunting` | num | % | Share of under-5s stunted (`hc70` height-for-age Z < −2 SD), design-weighted. `NA` where `hc70` unavailable. |
| 8 | `survey` | chr | — | DHS PR recode file base for the survey (e.g. `AOPR62FL`). |
| 9 | `cc` | chr | — | DHS 2-letter country code (e.g. `AO`). |
| 10 | `country` | chr | — | Country name (from the DHS survey universe). |
| 11 | `iso3` | chr | — | ISO-3166 alpha-3 country code (from `country` via `countrycode`); the key for GDP/DTP3 joins. |
| 12 | `year` | int | calendar year | Survey year, from the median interview date (`hv008` CMC → `1900 + floor((cmc−1)/12)`). |
| 13 | `u5mr` | num | deaths per 1,000 live births | All-cause under-5 mortality (₅q₀) for the region, `DHS.rates::chmort` rate `R` for the `U5MR` row. **Primary outcome.** |
| 14 | `m1mo5y` | num | deaths per 1,000 live births | 1-month–5-year (neonatal-excluded) mortality = `u5mr − NNMR` from the same `chmort` fit. **Secondary outcome.** |
| 15 | `exposure` | num | weighted birth-exposure | `chmort`'s weighted N (`WN`) for the U5MR estimate — the birth-exposure denominator; used as the `log(exposure)` offset (and to reconstruct death counts) in the GAM count model. |
| 16 | `log_gdp` | num | log(US$) | Natural log of GDP per capita, current US$ (World Bank `NY.GDP.PCAP.CD`), matched to `iso3` at the year nearest `year`. |
| 17 | `dtp3` | num | % | DTP3 immunisation coverage (World Bank `SH.IMM.IDPT`), matched to `iso3` at the year nearest `year`. |
| 18 | `mic_eq` | num | % | Microscopy-equivalent prevalence: `mic` where measured, else `k × rdt` with `k` = through-origin slope (microscopy = k·RDT) from Component 3 (≈ 0.74). |
| 19 | `pfpr2_10` | num | % | Age-standardised prevalence: `mic_eq` (clamped to [0,100]) converted 0.5–5y → 2–10y via `malariaAtlas::convertPrevalence` (Smith et al. 2007). Comparable to MAP PfPR₂₋₁₀. |
| 20 | `pfpr10` | num | per 10 PfPR points | `pfpr2_10 / 10`. The model's prevalence term, so its coefficient reads as the effect **per +10 PfPR₂₋₁₀ points**. |
| 21 | `year_c` | num | years | Centred calendar year: `year − round(mean(year))`. Enters the LMM linearly and the GAM as `s(year_c)`. |

## Notes
- **Prevalence lineage:** `rdt`/`mic` (raw, region-weighted) → `mic_eq` (common test footing) → `pfpr2_10` (common age footing) → `pfpr10` (model unit). Models use `pfpr10`/`pfpr2_10`, not the raw `rdt`/`mic`.
- **Missingness:** `stunting` and `mic` can be `NA` (variable not collected in a given survey); `mic_eq` fills missing microscopy from RDT so the exposure chain is complete. Models using `stunting` drop rows where it is `NA` (complete-case).
- **Weighting:** prevalence/covariate shares use the DHS sample weight `hv005`; mortality rates come from `DHS.rates::chmort`, which applies the birth-history weights internally.
- **Upstream precursor:** the raw prevalence-only table (before the mortality join and age-standardisation) is `data/dhs_prevalence_by_region.csv` (git-ignored, under `data/`).

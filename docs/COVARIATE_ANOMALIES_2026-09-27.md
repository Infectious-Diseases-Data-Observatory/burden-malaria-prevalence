# Regional covariate anomalies (27 September 2026)

## Resolution (v9, 28 September 2026)

Each finding below, including the scan, was checked by an independent reviewer, who reproduced the v7 MICS values from the microdata (exactly or within 0.1 point, apart from one fallback region) and validated the DHS recode method against the published values of other survey rounds (Sierra Leone 2008 and 2019, within about 2.3 points; aggregates only). All were confirmed except the ML8AFL Kidal join, several more widely than reported (corrections below). The user then decided to fix every verified problem in one new primary version, **`primary_map_regional17_dhsmics_gamma2_v9`** (settings version `regional_mics` in `R_cbh/primary/settings.R`). v7 is history as `regional_mics_v7`; its covariate inputs are archived in `data/derived_cbh/history/v7_inputs_2026-09-28/`. The code changes are commits 7340e5c (MICS), 4fa4510 (DHS), 89cc4e8 (v9 wiring) and 8e0eb3e (South Sudan 2010); the rules are in `R_mics/README.md` (step 10), `R_cbh/covariates/README.md` ("v9 covariate corrections") and `docs/ANALYSIS_PLAN.md` Section 2.4. MICS national values below are the weighted national values in `results/mics_inventory/covariates/national_values_by_survey.csv` (v7: the file at commit 07bf0e9); the original report quotes means over regions, which differ slightly.

### DHS: fixed

- **UG7BFL wealth (Section 1).** `R_cbh/covariates/regional.R` decodes `v190` with both label sets and stops on any other label; `05_validate.R` requires wealth to be observed for at least 95% of eligible mothers in every survey with `v190`. Wealth is now observed for all 10,237 eligible mothers in 15 regions (v7: 1,616 of 8,572 in 13), with regional scores 1.32–4.95.
- **UG7BFL Buganda (Section 1, related).** `R_cbh/config/region_overrides.csv` matches south buganda to Central 1 and north buganda to Central 2, with boundary-versioned ids and their own recode means. UG7BFL's model-ready rows rise from 75,366 to 89,861 (887 deaths). This is the only change to the complete-case sample: 14,495 records, 4,258 children, 156 deaths and 2 survey-regions.
- **Published joins (scan).** Survey-specific rows in `R_cbh/covariates/published_region_overrides.csv` join the 4 confirmed v7 survey-regions, GA61FL Estuaire and Ogooué-Maritime, UG52FL North and ZW62FL Harare (14,680 v7 records), and 4 survey-regions outside the complete-case sample that enter the imputed-covariate sensitivity: ZW52FL Harare, AO62FL mesoendemic unstable, and UG5AFL East Central and Mid Western (the AO62FL and UG5AFL published rows have exactly the recode regions' household and non-first-birth denominators). 03 and 13 now stop on any partial published-region gap not listed in `published_region_gaps_reviewed.csv` (10 reviewed survey-regions, each with a reason).
- **SL61FL (scan).** `published_survey_exclusions.csv` drops the survey's mis-assigned published subnational rows (national rows kept), and `17_recode_replacements.R` computes all 9 published variables for its 4 regions from the recodes with DHS definitions: electricity, improved water, improved sanitation, stunting and wasting from the household-member recode SLPR61FL, downloaded for this purpose; DTP3 (all three DTP doses), measles, facility delivery and short birth interval from SLBR61FL. The national recode values reproduce the published national values within 0.1 point (stunting and wasting: the denominator-weighted mean of the published regional partition, since the regional API cache holds no national nutrition row): electricity 13.5, water 60.9, sanitation 49.3, DTP3 77.85 (published 77.9), measles 78.6, facility delivery 54.4, short birth interval 16.2 (16.1), stunting 37.9, wasting 9.3. Regional ranges: electricity 3.0–57.5, water 44.7–93.8, sanitation 40.3–78.2, DTP3 69.1–86.4.
- **Fallbacks.** In the DHS overlays, regional-mean or WUENIC fallback flags switch off in 11 survey-regions (the 8 joins, SL61FL Western and UG7BFL Central 1 and 2) and on in none; WUENIC substitution stops in 6 of them.
- **ML8AFL Kidal: unchanged.** DHS does not tabulate Kidal in the 2023 survey, so the within-survey regional-mean fallback stays (1,171 records); it is a reviewed gap.

### MICS: fixed

- **Vaccination (Section 2 and the first four scan rows).** DHS tabulation rules: every living child aged 12–23 months with a completed under-5 interview; only positive card or recall evidence counts; never received, blanks or 97–99 on a seen card, recall counts 0–2, don't know and no information count as not vaccinated. DTP3 uses every DTP/DTCoq and pentavalent card and recall column; measles uses the first-dose measles/VAR and MR/RR card columns and every measles-containing recall item, including ROR, VAS and MR1 labels. Survey-specific columns (`vacc_extra`): MRT2011, COG2014, BEN2021, GHA2017, GNB2018, MDG2018, TGO2017, MOZ2008 and ZWE2009. The seen-card code is taken by MICS round, not from labels. Excluded from the denominator (user decisions): COD2017 cards kept at the health centre (179 children) and SSD2010 children outside the module's age filter with a blank module (254).
- **Education (scan).** A reviewed per-survey level map (`R_mics/education_level_map.csv`); never-attended women from the item before the level question get 0 years; Koranic, Mahadra, non-formal and literacy levels 0; vocational levels a flat P + S/2; the MICS6 grade attended is converted to completed (minus 1 when WB7 says not completed; 16 surveys).
- **Facility delivery (scan).** Any public, private-medical, NGO/mission or sector-unknown facility, by code range (by label in MOZ2008 and SSD2010), so mis-encoded public-hospital labels and private maternity homes now count as facilities. COD2017 and TCD2019 code 96 ('AUTRE PRIVE MEDICAL', labels apparently shifted by one) stays other (status quo, documented in the code; not yet a user decision and open to revision); as a facility it would add about 1 point in COD2017.
- **Water and sanitation (scan).** MRT2011 traditional wells unimproved; COD2017 delivered water, MOZ2008 water code 32, flush to an unknown place and MWI2006 'w. slab' latrines improved; MDG2018 non-washable slab kept improved and GHA2017 toilets 24/61 left out (status quo, documented).
- **Mozambique 2008 (Section 2).** Delivery place routed from MN7_A or MN8; recall count IM16 and measles card IM6D used; households without a toilet, which skip the facility question, unimproved.

| Survey | Variable | v7 | v9 | Main correction |
|---|---|---:|---:|---|
| MC_MOZ2008 | DTP3 / measles | 100.0 / 84.7 | 73.3 / 73.7 | Seen-card blanks as 0; recall IM16; measles card IM6D; denominator |
| MC_MOZ2008 | Facility delivery | 14.6 | 59.4 | MN7_A/MN8 routing |
| MC_MOZ2008 | Improved water / sanitation | 41.5 / 31.3 | 45.9 / 17.4 | Water code 32; no-toilet households |
| MC_MRT2011 | DTP3 / measles | 32.3 / 49.5 | 55.6 / 70.5 | Pentavalent columns; mis-encoded measles recall label; denominator |
| MC_MRT2011 | Facility delivery / improved water | 51.8 / 84.2 | 66.5 / 71.8 | Code-range classification; traditional wells |
| MC_COG2014 | DTP3 | 41.9 | 65.3 | Penta3 card column; DTC and Penta recall |
| MC_BEN2021 | Measles | 38.4 | 53.6 | Measles-rubella card and recall |
| MC_NGA2016 | DTP3 / measles | 65.4 / 66.7 | 33.0 / 41.2 | Denominator |
| MC_CAF2018 | DTP3 | 57.5 | 32.3 | Denominator |
| MC_MDG2018 | DTP3 / measles | 79.7 / 74.7 | 57.6 / 54.7 | Denominator; ROR recall |
| MC_COD2017 | DTP3 | 69.2 | 47.5 | Denominator; cards kept at the health centre excluded |
| MC_CMR2014 | DTP3 / measles | 90.6 / 87.0 | 72.0 / 72.0 | Seen-card blanks as 0; denominator |
| MC_KEN2011NYA | DTP3 / measles | 93.4 / 96.2 | 84.2 / 79.1 | Seen-card blanks as 0; denominator |
| MC_GNB2018 | Maternal education (years) | 5.3 | 2.7 | Never-attended women; level map |
| MC_COD2017 | Maternal education (years) | 7.5 | 6.0 | Never-attended women; upper-secondary map |
| MC_ZWE2019 | Maternal education (years) | 8.0 | 9.2 | Lower-secondary grades |
| MC_TGO2017 | Facility delivery | 71.3 | 80.0 | Code-range classification |

On the complete-case sample the scaling moves accordingly (v7 → v9 `covariate_scaling.csv`): DTP3 mean 71.1 → 68.8 (SD 23.3 → 23.9), measles 71.6 → 70.0, facility delivery 61.2 → 62.0, maternal education 4.24 → 4.22 years, wealth score 2.786 / 0.714 → 2.785 / 0.719.

### v9 sample and results

**7,621,617 child-band records, 2,329,388 children, 104,143 deaths, 1,229 survey-regions, 135 surveys and 37 countries** (v7: 7,607,122 / 2,325,130 / 103,987 / 1,227 / 135 / 37). DHS part: 5,588,463 records, 1,723,303 children, 77,587 deaths, 934 survey-regions, 98 surveys, 35 countries; MICS part: 2,033,154, 606,085, 26,556, 295, 37, 24, the same records as v7. MAP-eligible records: 8,812,458 (6,372,297 DHS/MIS). All seven fits converged without restarts. The 40%→20% hazard ratios change by at most 0.011 (48–59 months, 0.919 → 0.908) and the 20%→0% ratios by at most 0.008 (24–35 months, 0.470 → 0.462); national attributable deaths are 1,045,865 (2005), 736,151 (2015) and 629,136 (2024), 1.3–1.5% above v7 (`results/cbh/primary_map_regional17_dhsmics_gamma2_v9/`: `comparison_contrasts.csv`, `burden/year_summary.csv`). The sensitivities are being refitted on v9 (subgroups v4, no-nutrition v4, imputed covariates v10, SMC v3; tensor, linear 1–30%, no-region-effect, bam-versus-gam and reference-floor v2).

### Corrections to the original report

- **Join failures (scan, Medium row).** The 21,860-record figure does not reproduce. The six fallback survey-regions behind the row were GA61FL Estuaire and Ogooué-Maritime, UG52FL North, ZW62FL Harare, ML8AFL Kidal and SL61FL Western; only the first four are join failures, with 14,680 v7 records (3,057 + 2,289 + 9,100 + 234), and in each of them all 9 published variables were fallbacks, not only the household ones.
- **ML8AFL Kidal** was not a join failure: the 2023 survey publishes no Kidal estimate (the '....Kidal' row belongs to the 2018 survey, ML7AFL), and its regional-mean fallback follows the plan's rule.
- **SL61FL (scan, High row).** The problem covered all 9 published variables in all 4 regions (68,076 v7 records, not 26,609): the SL2013 API subnational rows are assigned to the wrong areas (no Western row; 'Southern' matches Western Area plus Pujehun and 'North Western' Bo, Bonthe and Moyamba in the recode; the selected 'Northern (before 2017)' row does not match the old Northern Province, whereas the post-2017-labelled 'Northern' row does). A national-consistency check cannot detect this, because the mis-assignment preserves the national totals.
- **MC_COG2014 (scan, High row).** Besides the split card column, v7 used only the pentavalent recall item (IM15D) and ignored the DTC recall item (IM12), asked of the same children. v9 uses both card columns and both recall items (65.3; pentavalent recall alone would give 59.4).
- **Values quoted as current** in the scan are means over v7 regions (record-weighted unless stated), not national values: for example MRT2011 30.6 / 48.7 (national 32.3 / 49.5), COG2014 41.7 (41.9), BEN2021 38.2 (38.4), GNB2018 education 4.9 (5.3) and COD2017 education 7.0 (the unweighted region mean; record-weighted 6.9; national 7.5).
- **Wider than reported** (found during verification and fixed in v9): doses missing from a seen card stored as 99 (the whole date missing) in the four Kenya sub-national surveys of 2011–2013 and CMR2014; measles recall items labelled ROR, VAS or MR1 missed in MDG2018, TGO2017, GHA2017 and GNB2018; education grade-coding errors also in CAF2018, MLI2015, GHA2017, GHA2011, MWI2019 and SLE2017, and the MICS6 grade recorded as attended rather than completed; private maternity homes scored as home deliveries in 13 surveys; ZWE2009 (outside v7) missing its DPT-HepB3 card column; MOZ2008 households without a toilet dropped from the sanitation denominator. The denominator rule alone does not explain MRT2011 and COG2014 (split card and recall columns).

### Open items (not changed in v9)

- MC_SOM2006, MC_SOM2011NE and MC_SOM2011SL (outside the complete-case sample, in the imputed-covariate sensitivity): vaccination under the new rules is unreviewed (national DTP3 71.3 → 16.6, 33.3 → 9.3 and 29.4 → 12.8; SOM2011NE has about 5% item nonresponse counted as not vaccinated).
- MC_SOM2006 water codes 52–54 (roof top, berkad, balli; about 18% of households) and MC_MDG2012S toilet codes 24/25 stay unclassified, as in v7.
- Candidate sample changes not applied, each needing a boundary review: ZW72FL recode 'harare' as boundary Harare Chitungwiza (at most 4,062 records, 35 deaths) and UG5AFL recode 'south' as South Western (at most 1,819 records, 18 deaths).
- The legacy individual wealth decode in `R_cbh/R/child_bands.R` (lines 102–103) still reads only *poorest…richest*; v9 does not use it.
- Reporting: the counts of WUENIC-substituted and regional-mean-filled survey-regions quoted in reports and the paper need refreshing from the v9 outputs, and the v9 figures have not yet been exported.

The original report follows unchanged, as the record of what was found on 27 September.

---

Found through the supplementary covariate pairs figure (`R_cbh/reporting/17_covariate_pairs.R`) and investigated read-only at the user's request: two root-cause investigations, each checked by an independent reviewer who re-derived the key numbers, and a scan of all 135 v7 surveys. **Nothing has been changed**; the v7 primary, its sensitivities and the paper exports still use the values described here. Record counts are v7 child-band records (7,607,122 records, 103,987 deaths). Only national values and ranges are given; region-level values stay in `data/`.

## 1. Confirmed: Uganda DHS 2016 (UG7BFL) wealth-quintile score is 3.00 in every region

- **Root cause (label decoding).** The UG7B recode stores `v190` (the national combined wealth index) as text with the DHS-7 labels *lowest, second, middle, fourth, highest*. `R_cbh/covariates/regional.R:18` decodes only *poorest, poorer, middle, richer, richest*, and `cbh_decode_category()` falls back to `as.numeric()`, which gives NA for words. Only *middle* decodes (to 3); the other four quintiles are silently dropped as missing, so each region averages only its middle-quintile mothers. Wealth is observed for 1,616 of 8,572 eligible mothers (18.9%); every other DHS survey is at 98% or more. Label order checked against the factor score `v191` (means −1.01, −0.57, −0.22, 0.32, 1.65).
- **Only this survey.** Of the 120 DHS/MIS recodes, 110 use the *poorest…richest* labels and decode fully; 9 have no `v190` at all (the documented whole-survey gaps); UG7BFL alone uses *lowest…highest*. MICS is unaffected (numeric `windex5`).
- **Affected:** 13 survey-regions, 75,366 records, 731 deaths (about 1% of records).
- **Corrected values:** regional scores 1.32 (Karamoja) to 4.95 (Kampala), consistent with UG52FL (1.70–4.92) and UG61FL (1.47–4.87); national mean 2.99 with about 20% of mothers in each quintile. The v7 scaling of the wealth score would change from mean 2.786 / SD 0.714 to 2.783 / 0.718.
- **Fix:** extend the labels at `regional.R:18` to `c("poorest","poorer","middle","richer","richest","lowest","second","fourth","highest")` with codes `c(1:5,1,2,4,5)`, and add a guard that stops on any undecoded non-missing `v190` label; add a `05_validate.R` check that wealth is observed for at least 95% of eligible mothers when `v190` exists. The same decode appears in the legacy `R_cbh/R/child_bands.R:102-103` (superseded individual factor; changing it forces a rebuild of every child-band shard) and `R_dhs/03_build_analysis_dataset.R:171-174`.
- **Related, separate:** UG7BFL's *north buganda* and *south buganda* regions are unmatched in `data/derived_cbh/region_crosswalk.csv` (the boundary file calls them Central 1 and Central 2), so about 14,495 records are not model-ready. `R_cbh/covariates/published_region_aliases.csv:10-11` already maps southbuganda→central1 and northbuganda→central2; adding the same pairing to the recode-level region overrides would admit them. This changes the sample.

## 2. Confirmed: Mozambique MICS 2008 (MC_MOZ2008) DTP3 = 100% in every region, and facility delivery and measles are also wrong

- **DTP3 root cause (two defects in `R_mics/09_regional_covariates.R`).** (1) A dose missing from a seen card is left blank, with no 0 code; `card_dose()` (line 60) returns 1 or NA, never 0, so the 316 children with a seen card and no DPT3 (most with DPT1 or DPT2 recorded) drop out. (2) The recall-count regex (line 168) needs the vaccine name in the label, but the Portuguese label of IM16 is *Quantas vezes recebeu?*, so no recall variable is found and recall is NA for everyone. Only positive card records remain: 100% on n = 1,742.
- **Measles:** the card-day variable (IM6D, label *Sarampo*) is not matched by the name or label patterns (lines 173–174), so measles rests on recall for 270 of 2,399 children.
- **Facility delivery:** the survey splits place of delivery into intended place (MN7_A), whether she delivered there (MN7_B) and actual place only if different (MN8). `R_mics/01_inventory.R:99` selects MN8 and line 136–138 codes only MN8, so the indicator covers the 1,207 women whose delivery place changed, mostly to home: 14.6% nationally against 54% and 62% in the 2003 and 2011 DHS.
- **Affected:** 11 survey-regions, 66,350 records, 1,396 deaths.
- **Corrected values** (blank on a seen card = not given; IM16 as DPT recall with "never received" = 0; IM6D as measles card date; delivery place routed from MN7_A or MN8): DTP3 74.9% nationally (regions 57–91; WUENIC 75), measles 75.1% (62–94; WUENIC 77), facility delivery 59.4% (41–94), all in line with the Mozambique DHS rounds. The all-children denominator gives DTP3 73.3%.
- **Fix:** survey-specific changes in `09_regional_covariates.R` (the recall override must be placed after line 168, which otherwise overwrites it; the measles card-date fallback should be guarded to this survey) and the routed delivery variable, preferably recorded in `01_inventory.R`. A simpler alternative for vaccination only is adding the survey to `card_only` (flat WUENIC 75/77), which leaves facility delivery wrong.

## 3. Scan of all surveys (single reviewer; not yet independently verified)

The reviewer re-derived values from the microdata, reproducing the current pipeline exactly before applying corrections. Beyond the two issues above:

| Severity | Surveys | Variable | Finding | Records |
|---|---|---|---|---:|
| High | 34 MICS surveys not on the WUENIC list; 5+ point bias in NGA2016, MDG2018, CAF2018, MLI2015, COD2017, BEN2021, CIV2016, MRT2015, BEN2014, CMR2014, TGO2017 | DTP3 (and measles where recall is unmatched) | Recall is scored only from the dose-count item, which is skipped when "ever received" is No, so unvaccinated children without a card drop out of the denominator. MICS DTP3 sits a median 7.9 points above the nearest DHS (region-matched); re-derived national values fall, e.g. NGA2016 65→34 (WUENIC 53), CAF2018 58→34 (38), MDG2018 80→59 (74). Record-weighted mean shift −7.1 points (DTP3) and −4.8 (measles); 83 regions move by 10 points or more | 841,608 (5+ point surveys); 1,827,552 (all 34) |
| High | MC_MRT2011 | DTP3, measles | Card uses DTCoq3 though doses are in the Penta3 column; measles recall label mis-encoded. 30.6 / 48.7 against 69 / 79 in MICS 2015; corrected 57 / 72 (WUENIC 75 / 67) | 51,041 |
| High | MC_COG2014 | DTP3 | Penta3 dates in a different column; 305 children with Penta3 scored 0. 41.7 against WUENIC 80; corrected 66.7 | 51,604 |
| High | MC_BEN2021 | Measles | Measles-rubella doses ignored; 38.2 against 68–76 in neighbouring rounds; corrected 54.3 (WUENIC 53) | 78,127 |
| High | MC_GNB2018, MC_COD2017 | Maternal education | Never-schooled women dropped (ever-attended and level labels fail the patterns): GNB2018 4.9 years against 2.2 in 2014 (corrected about 2.5); COD2017 7.0 → 5.4, distorting the regional ranking | 161,664 |
| High | SL61FL | Electricity, water, sanitation | Cached DHS API rows for Southern look mislabelled (electricity 45% in a 16%-urban region; 6–8% in 2008 and 2019); Western Area has regional-mean fallbacks (electricity 18.6% in a 91%-urban capital region) | 26,609 |
| Medium | MC_MRT2011, MC_TGO2017, MC_MWI2006, MC_NGA2016 | Facility delivery | Label patterns miss mis-encoded or local facility names; "Private maternity home" matches "home" | 371,535 |
| Medium | MC_ZWE2019, MC_NGA2016, MC_CMR2014, MC_ZWE2014 | Maternal education | Grade coding: level×10+grade and a lower-secondary cap misconvert years (ZWE2019 about −1.5 years) | 286,627 |
| Medium | UG52FL, ZW62FL, GA61FL (6 regions) | Published household variables | Published values exist but did not join, so regional-mean fallbacks were used (e.g. Harare electricity 30.7 against published 80.7) | 21,860 |
| Medium | MC_MRT2011 | Improved water | Traditional wells unclassified and dropped from the denominator (84 → 72–74%) | 51,041 |

Judged genuine or documented (no change proposed): the WUENIC surveys MC_GIN2016, MC_COM2022 and MC_TCD2019; the within-sample wealth ranking of the sub-national Kenya and Dakar MICS; MICS6 off-grid electricity (ZWE2019); Congo's 0% urban departments; Liberia 2007's separate Monrovia region; published DHS values that differ between rounds (Rwanda 2005 water and sanitation, Mauritania 2019–21 wasting, Mali 2018 Kidal vaccination from a tiny weighted sample, Cameroon 2018 South-West urban-only fieldwork). Clean checks: every covariate is constant within every survey-region; prepared values equal the overlays exactly; no impossible values; DHS recode variables are at least 98% observed except UG7BFL wealth.

## What a correction would involve

The scaling changes, so the plan (`docs/ANALYSIS_PLAN.md`, scaling rule) requires a new analysis version rather than overwriting v7. In order: re-run the covariate extraction and overlays (`R_cbh/covariates/` 02–05 and 13; `R_mics/09_regional_covariates.R` and the MICS merge), re-prepare and refit the primary as a new version, re-run the v8 imputation (the flawed values were predictors and donors), every v7-based sensitivity (subgroups, nutrition, SMC, tensor, linear 1–30%, no region RE, bam vs gam, reference floor, Sahel), the reporting stages (covariate forest, pairs, Figures 1–4 and the tables) and the paper export. The complete-case sample is unchanged unless the Buganda regions are admitted. The effect on the PfPR hazard ratios has not been estimated; the affected surveys carry survey random intercepts, which absorb survey-wide level shifts but not within-survey distortions.

Scripts and aggregate outputs of the investigation are in the session scratchpad (`anomalies/`), not in the repository.

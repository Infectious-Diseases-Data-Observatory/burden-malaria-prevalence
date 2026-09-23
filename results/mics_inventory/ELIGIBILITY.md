# Can MICS surveys enter the primary analysis?

Assessment of the 93 MICS surveys in `data/MICS_Datasets/` against the requirements of the primary child age-band analysis. Code: `R_mics/01_inventory.R` (inventory), `R_mics/02_inventory_tables.R` (tables), `R_mics/03_eligibility_report.py` (this page). Details: [inventory](INVENTORY.md), [eligibility table](eligibility.csv).

Verdicts come from two independent routes: a scripted inventory with direct checks, and per-country audit agents each followed by an adversarial re-check. The two agreed on 38 of the 40 surveys both covered; the two corrections are applied below. Malawi, Mozambique and Zimbabwe were checked directly only.

## Outcome (23 September 2026)

All 47 surveys with a complete birth history were added to the analysis dataset. Under the rule then set for the primary, include surveys with anthropometry, complete-case selection kept **37 MICS surveys**, which enter `primary_map_regional17_dhsmics_gamma2_v5` alongside the 95 DHS surveys. That rule admits the six sub-national surveys with anthropometry (five Kenyan and Dakar city), whose wealth quintiles are ranked within their own sample. Guinea-Bissau and the Central African Republic enter with analysis regions that merge each capital with a neighbour. The verdicts below are the pre-integration assessment, corrected for Zimbabwe 2009.

## Summary

- **46 of 93 surveys have no complete birth history** (all 24 MICS2, 12 of 16 MICS3, 10 of 19 MICS4). They carry only summary questions (children ever born, children dead, date of last birth) and cannot support the age-band design.
- **47 surveys in 28 countries have a complete birth history**, structurally equivalent to the DHS Births Recode. Seven use non-standard names (Portuguese in Mozambique 2008, lowercase in Zimbabwe 2009, month-year dates or unit-coded ages elsewhere) but are complete.
- **25 surveys in 18 countries can enter with the existing method** (18,898 under-five deaths among births in the 5 years before interview), **6 in 5 countries after extra work** (4,369 deaths), and **16 cannot** (6,040 deaths). The eligible surveys would add roughly a third to the 75,726 deaths in the current primary sample, including two countries with no DHS survey (Central African Republic, Guinea-Bissau).

## What each requirement needs for MICS

- **Birth history.** Complete in all 47. Variable names differ across rounds and must be mapped (BH4C or CCDOB or month and year for birth date; BH9C, or units and number, for age at death; WDOI, WM6C or month and year for interview).
- **Regional covariates.** For DHS, 4 of the 13 came from the recodes and 9 from published DHS survey-region indicator tables. MICS has no such tables, so all 13 must be computed from the microdata. Every source exists in the eligible surveys: water, toilet type and electricity in `hh.sav`; DTP3, measles and anthropometry in `ch.sav` or `who_z.sav`; place of delivery and education in `wm.sav`; age at first birth and birth intervals from `bh.sav`.
- **Definition differences to resolve.** MICS records education as level and grade, not completed years, so a country-specific mapping is needed. It asks about place of delivery only for the last birth in the 2 years before interview, against the DHS 5-year window. MICS3 to MICS5 supply both old NCHS-referenced and WHO z-scores; only the WHO versions (`HAZ2`, `WHZ2`) match. MICS6 electricity has an off-grid category whose coding must match the DHS definition. Guinea 2016, Comoros 2022 and Chad 2019 lack usable recall doses, so DTP3 and measles are card-based there.
- **Regions and MAP prevalence.** Most MICS regions match an existing DHS boundary file in the same country name for name, sometimes after merging MICS regions into a DHS grouping (Chad, Mauritania, Togo) or aggregating districts to regions (Malawi 2013). Guinea-Bissau and the Central African Republic have no DHS boundary, so new admin-1 polygons are needed. Nigeria's DHS surveys in the study use the six zones, which the MICS zone variable matches.
- **National series.** Complete for every eligible country. They block Somalia (health expenditure absent before 2013), South Sudan (political stability absent before 2011, health expenditure before 2017) and Sao Tome (no child HIV incidence series), as they block the equivalent DHS records.

## Can enter with the existing method (25 surveys)

| Survey | Year | Regions | Under-5 deaths, last 5 y | Boundary | Note |
|---|---:|---:|---:|---|---|
| BEN_2014_MICS5 | 2014 | 12 | 1,028 | Matches BJ2017DHS and BJ2012DHS (LEVELNA 'Departments', 12 polygons) one-to-one: Atlantique = A | National survey. The birth history is complete and CMC-dated. All 13 covariates can be computed, reusing DHS department boundaries, and the national series and MAP surfac |
| BEN_2021_MICS6 | 2021 | 12 | 787 | Matches BJ2017DHS and BJ2012DHS (12 'Departments') one-to-one: ATLANTIQUE = Atlantic, OUEME = O | National survey with a complete CBH, all covariates computable and existing department boundaries. No DHS for Benin after 2017 is in the registry, so no duplication. Rout |
| CIV_2016_MICS5 | 2016 | 11 | 778 | Matches CI2012DHS (11 polygons, LEVELNA 'Groups of Regions') name for name: Centre, Centre-Est, | National survey with a complete CBH. All covariates can be computed (age at first birth from the birth history), and CI2012DHS boundaries can be reused as they are. Does  |
| CMR_2014_MICS5 | 2014 | 12 | 552 | Matches the 12 polygons of CM2011DHS and CM2018DHS (adamaoua, centre, douala, est, extreme-nord | Complete birth history. Entry years 2009-2014, between the CM2011 and CM2018 DHS windows; a separate survey, not a duplicate. All 13 covariates are present or derivable,  |
| COD_2017_MICS6 | 2018 | 26 | 1,186 | Matches the 26 provinces of CD2023DHS | All covariates present with WHO z-scores; national series complete. |
| COG_2014_MICS5 | 2015 | 12 | 467 | Matches CG2011DHS 12 departments name for name | All covariates present; national series complete. |
| COM_2022_MICS6 | 2022 | 3 | 132 | Matches KM2012DHS: DHSREGEN Mwali, Ndzuwani, Ngazidja (LEVELNA Autonomous Islands; REGNAME mohe | The birth history is complete, all 13 covariates are present or derivable (DTP3 slightly understated because recall dose counts are unusable), the 3 islands reuse the KM2 |
| GHA_2011_MICS4 | 2011 | 10 | 546 | Matches the 10 regions of GH2008DHS and GH2014DHS. One synonym is needed: MICS 'Asante' = DHS ' | Complete birth history. Entry years 2006-2011. All 13 covariates are present or derivable, with WHO z-scores that reproduce the published national figures. The 10 regions |
| GHA_2017_MICS6 | 2017 | 10 | 362 | The names (WESTERN ... UPPER WEST) match the 10 regions of GH2014DHS and GH2008DHS exactly afte | Complete birth history. Entry years 2012-2017. All 13 covariates are present or derivable, with WHO z-scores that reproduce the published figures. Regions match GH2014DHS |
| GMB_2018_MICS6 | 2018 | 8 | 489 | Exact match to the 8 LGA polygons in GM2019DHS and GM2013DHS (Banjul, Kanifing, Brikama, Mansak | The birth history is complete with near-complete date reporting, all 13 covariates are available or routinely derivable (WHO z-scores supplied), and the regions match exi |
| MDG_2018_MICS6 | 2018 | 22 | 641 | Matches MD2008DHS exactly: 22 regions, Analamanga to Sava (MICS AMORON'I MANIA vs DHSREGEN Amor | Large, nationally representative survey with a complete birth history, every covariate present or routinely derivable (WHO z-scores, card plus recall vaccination), and 22 |
| MLI_2015_MICS5 | 2015 | 8 | 1,184 | Matches 8 of the 9 region polygons in ML2018DHS (and ML2006DHS/ML2023DHS) by name. REGNAME spel | The birth history is complete, all 13 covariates are available or routinely derivable, WHO z-scores are supplied, and the 8 admin-1 regions match existing Mali DHS polygo |
| MOZ_2008_MICS3 | 2008 | 11 | 1,259 | Matches the 11 provinces of MZ2003/2011/2022 DHS | Portuguese variable names; WHO z-scores in who_z.sav. |
| MRT_2015_MICS5 | 2015 | 13 | 510 | MR2020DHS (14 polygons) can be reused. Dissolve Nouakchott Nord/Ouest/Sud into one Nouakchott ( | The birth history is complete, all 13 covariates are available or routinely derivable with WHO z-scores, and the regions map onto existing MR2020DHS polygons after a rout |
| MWI_2006_MICS3 | 2006 | 3 | 2,075 | Matches the 3 regions of every Malawi DHS | WHO z-scores in who_z.sav; national series complete. |
| MWI_2013_MICS5 | 2014 | 3 | 1,200 | 31 districts aggregate to the 3 regions used by Malawi DHS | All covariates present. |
| MWI_2019_MICS6 | 2020 | 3 | 657 | Matches the 3 regions of every Malawi DHS | All covariates present. |
| NGA_2016_MICS5 | 2016 | 6 | 2,671 | Zone labels match the 6 'Groups of States' in NG2013DHS, NG2018DHS and NG2024DHS exactly (North | Complete birth history. All 13 covariates are present or derivable, including WHO z-scores. Entry years 2011-2017. National series are complete and a MAP surface exists.  |
| SLE_2017_MICS6 | 2017 | 4 | 834 | HH7 matches the 4 province polygons in SL2013DHS and SL2008DHS by name (EAST=Eastern, etc.). It | The birth history is complete with near-complete dates, all 13 covariates are available or routinely derivable with WHO z-scores, and the 4 regions match the existing SL2 |
| SWZ_2010_MICS4 | 2010 | 4 | 216 | Matches SZ2006DHS exactly: DHSREGEN Hhohho, Manzini, Shiselweni, Lubombo (LEVELNA Regions). | The birth history is complete, all 13 covariates are available or routinely derivable (WHO z-scores present), the 4 admin-1 regions match the existing SZ2006DHS boundary, |
| SWZ_2014_MICS5 | 2014 | 4 | 128 | Matches SZ2006DHS exactly (Hhohho, Manzini, Shiselweni, Lubombo). | The birth history is complete with an imputed CMC and age at death in months, all 13 covariates are available or derivable from standard variables, the regions reuse the  |
| SWZ_2021_MICS6 | 2021 | 4 | 73 | Matches SZ2006DHS (names in upper case: HHOHHO, MANZINI, SHISELWENI, LUBOMBO). | The birth history is complete, all covariates are present (WHO z-scores, vaccination by card plus recall), there are 4 admin-1 regions matching SZ2006DHS, national series |
| TGO_2017_MICS6 | 2017 | 7 | 267 | TG2013DHS has 6 'Regions': Centrale, Kara, Plateaux, Savanes, 'Maritime (Sans Agglomération de  | National survey with a complete CBH and all covariates computable. Regions map onto the existing TG2013DHS polygons once the two Lomé strata are merged. Does not duplicat |
| ZWE_2014_MICS5 | 2014 | 10 | 571 | Matches ZW DHS 10 provinces | All covariates present. |
| ZWE_2019_MICS6 | 2019 | 10 | 285 | Matches ZW DHS 10 provinces | All covariates present. |

## Can enter after extra work (6 surveys)

| Survey | Year | Regions | Under-5 deaths, last 5 y | Boundary | Note |
|---|---:|---:|---:|---|---|
| CAF_2018_MICS6 | 2019 | 7 | 633 | No DHS boundary for CAF; dissolve prefecture polygons from an admin-1 source into the 7 regions | New country for the study; needs a new boundary source. |
| GIN_2016_MICS5 | 2016 | 8 | 485 | Matches GN2018DHS and GN2012DHS (8 'Regions') one-to-one: Boké, Conakry, Faranah, Kankan, Kindi | National survey with a complete CBH, existing DHS region boundaries and a MAP surface, and no duplication of DHS 2012 or 2018. One substantive step: decide how to handle  |
| GNB_2014_MICS5 | 2014 | 9 | 504 | No DHS boundary for Guinea-Bissau (no GNB DHS in the study). Needs an admin-1 source such as GA | National survey with a complete CBH, all 13 covariates computable, complete national series and a MAP surface. It adds a country with no DHS. The substantive step is sour |
| GNB_2018_MICS6 | 2019 | 9 | 303 | No DHS boundary for Guinea-Bissau. Same 9 regions as GNB 2014; needs GADM admin-1 or an equival | National survey with a complete CBH, all covariates computable, complete national series and a MAP surface. Needs new admin-1 polygons for Guinea-Bissau, shared with GNB  |
| MRT_2011_MICS4 | 2011 | 12 | 517 | MR2020DHS has 14 polygons: Nouakchott Nord/Ouest/Sud and a combined 'Tiris Zemour et Inchiri'.  | The birth history is complete and all 13 covariates are derivable, including WHO z-scores. However, Tiris Zemmour was surveyed without Inchiri, and the existing MR2020DHS |
| TCD_2019_MICS6 | 2019 | 22 | 1,927 | TD2014DHS has 21 regions; merge MICS regions into its combined polygons (e.g. Borkou/Tibesti) | Regions merge cleanly into TD2014DHS (Ennedi Est+Ouest to Ennedi; Borkou to Borkou/Tibesti). But recall vaccination answers were largely not entered (coded 'donnee non sa |

## Cannot enter the complete-case primary (16 surveys)

| Survey | Year | Regions | Under-5 deaths, last 5 y | Boundary | Note |
|---|---:|---:|---:|---|---|
| KEN(Bungoma County)_2013_MICS5 | 2013 | 1 | 22 | County polygons exist in KE2014DHS | Sub-national sample (one county, one former province or informal settlements), not nationally representative; excluded for consistency with the national-survey design. |
| KEN(Kakamega County)_2013_MICS5 | 2013 | 1 | 39 | County polygons exist in KE2014DHS | Sub-national sample (one county, one former province or informal settlements), not nationally representative; excluded for consistency with the national-survey design. |
| KEN(Mombasa Informal Settlements)_2009_MICS4 | 2009 | 1 | 39 | County polygons exist in KE2014DHS | Sub-national sample (one county, one former province or informal settlements), not nationally representative; excluded for consistency with the national-survey design. |
| KEN(Nyanza Province)_2011_MICS4 | 2011 | 6 | 377 | All six county names match KE2022DHS county polygons | Sub-national (6 of 47 counties), excluded under the national-survey rule. It is the strongest sub-national candidate if that rule changes (377 recent deaths, exact county |
| KEN(Turkana County)_2013_MICS5 | 2013 | 1 | 46 | County polygons exist in KE2014DHS | Sub-national sample (one county, one former province or informal settlements), not nationally representative; excluded for consistency with the national-survey design. |
| LSO_2018_MICS6 | 2018 | 10 | 172 | HH7A districts match the LS2004/2009/2014/2023 DHS boundaries (10 districts; spelling differenc | The data are otherwise complete and district-matchable, but Lesotho has no Malaria Atlas Project prevalence surface, so there is no exposure value for any region. |
| MDG(South)_2012_MICS4 | 2012 | 4 | 207 | The 4 regions (Anosy, Androy, Atsimo Atsinanana, Atsimo Andrefana) match MD2008DHS region polyg | Two of the 13 required covariates (wasting, stunting) are entirely absent because the survey had no anthropometry, so it cannot enter the complete-case primary analysis.  |
| NGA_2021_MICS6 | 2021 | 6 | 2,149 | zone matches the 6 zones in NG2018DHS and NG2024DHS exactly, and the state-to-zone mapping is s | The birth history, regions and the other 11 regional covariates are usable. Stunting and wasting cannot be computed because no heights or weights were collected, so the s |
| SEN(Dakar City)_2015_MICS5 | 2015 | 4 | 159 | No match. The Senegal DHS boundaries in the study (SN2014-SN2019DHS) are 4 groups of regions (C | This is a city-only, sub-national sample with admin-2 regions that cannot be matched to existing boundaries and no urban variable. It overlaps the Senegal continuous DHS  |
| SOM(Northeast Zone)_2011_MICS4 | 2011 | 3 | 214 | No DHS boundary | Sub-national zone sample; no anthropometry; health expenditure absent before 2013. |
| SOM(Somaliland)_2011_MICS4 | 2011 | 5 | 372 | No DHS boundary | Sub-national zone sample; no anthropometry; health expenditure absent before 2013. |
| SOM_2006_MICS3 | 2006 | 3 | 643 | No DHS boundary for Somalia | Health expenditure series absent before 2013, so every band entry fails complete-case selection; also needs new boundaries. |
| SSD_2010_MICS4 | 2010 | 10 | 997 | No DHS boundary for South Sudan | Political stability absent before 2011 and health expenditure before 2017, so band entries in 2005-2010 fail complete-case selection. |
| STP_2014_MICS5 | 2014 | 4 | 80 | ST2008DHS has the same 4 groups | No UNAIDS child HIV incidence series for Sao Tome, the same reason the Sao Tome DHS is excluded; usable only in the imputed-covariate sensitivity. |
| STP_2019_MICS6 | 2019 | 5 | 24 | ST2008DHS 4 groups after merging the two central districts | No UNAIDS child HIV incidence series for Sao Tome, the same reason the Sao Tome DHS is excluded. All other covariates are present (place of delivery is MN20, labelled in  |
| ZWE_2009_MICS3 | 2009 | 10 | 500 | Matches the 10 provinces of ZW2005/2010/2015 DHS | Zimbabwe's health expenditure series starts in 2010, so every band entry (2004-2009) fails complete-case selection. Missed in the first assessment and found when the surv |

## Caveats

- Adding MICS changes the primary sample definition from DHS and MIS to DHS, MIS and MICS, so it would be a new primary version rather than a sensitivity. The survey random effect absorbs level differences between programmes, but not differences in covariate definitions.
- This assessment excluded sub-national MICS (five Kenyan county, province or settlement samples, Dakar city, southern Madagascar, two Somali zones) for consistency with the national-survey design. Their wealth quintiles are also ranked within the sub-national sample, not nationally, so that covariate would have to be rebuilt against a national reference. Nyanza 2011 was the strongest candidate. The rule later set for the primary admits those with anthropometry, with the within-sample ranking retained (see Outcome).
- Several MICS surveys overlap DHS surveys in the same country and years (for example Gambia 2018 with DHS 2019, Sierra Leone 2017 with DHS 2019). They are independent samples, not duplicates, so both can enter.
- Nigeria 2021 (no anthropometry) and Sao Tome 2014 and 2019 (no child HIV series) fit the methods the imputed-covariate sensitivity already used, so they could enter that version. Somalia and South Sudan would need several years of national series extrapolated. The imputed-covariate version on the DHS and MICS sample (`primary_map_regional17_dhsmics_imputed_gamma2_v6`) now includes all of them, with the Somali and South Sudanese national series extrapolated by GAMs.
- The under-five death counts here are among births in the 5 years before interview. The primary analysis counts deaths in eligible complete age bands, which is a similar but not identical quantity.


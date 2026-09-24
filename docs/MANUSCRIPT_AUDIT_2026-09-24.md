# Manuscript audit: main.tex against the DHS + MICS pipeline (24 September 2026)

**Scope and ground rules.**
- Nothing in the Overleaf folder was edited: not `main.tex`, not `ref.bib`, not `SText1_model_specification.tex`. Every change below is a suggestion for the author to paste in.
- The source of truth is the current DHS + UNICEF MICS primary, `results/cbh/primary_map_regional17_dhsmics_gamma2_v5/`. Its sensitivity analyses are:
  - subgroups: `subgroups_dhsmics_map_gamma2_v2/`
  - no-nutrition refit: `nutrition_adjustment_dhsmics_map_gamma2_v2/`
  - imputed covariates: `primary_map_regional17_dhsmics_imputed_gamma2_v6/`
- The SMC before/after analysis (`smc_dhsmics_map_gamma2_v1/`) is **not** part of the paper and was not used.
- File paths below are relative to `results/cbh/primary_map_regional17_dhsmics_gamma2_v5/` unless given in full.
- **The figures are current but the text is not.** Every figure image in Overleaf is the current pipeline export: Figures 1–3 match `paper_figures/manifest.csv`, and Figures S1–S5 match the v5/v2/v6 PNGs by MD5. Most of the text, Tables 1, 2 and S1, and the captions still carry the DHS-only (v3) numbers. The manuscript therefore currently contradicts its own figures.
- Each section was checked by an auditor and then a verifier. The verifier's input was cut off partway through the Methods. Methods lines 210–229 and all of lines 234–333 (Methods continued, Supplement, Table S1, SText1) were therefore checked directly while writing this report, against the files named under each item.

## Summary

This audit covered the manuscript line by line, including the abstract, captions and table cells. About 180 statements and numbers were found correct. **90 findings** were confirmed or corrected; none of the verified findings was rejected.

Findings by category:

| Category | Count | Meaning |
|---|---:|---|
| Outdated | 32 | Correct for the DHS-only version but not for v5 |
| Wrong | 12 | Wrong statement or number, including unfilled placeholders that need a number |
| Internally inconsistent | 5 | Contradicts another part of the manuscript |
| Clarity, typo, citation or placeholder | 41 | Wording, formatting, missing references |

A further 6 items are external literature claims that the pipeline cannot check. They are listed separately at the end and are not counted as errors.

Findings by section:

| Section | Lines | Findings | Unverifiable |
|---|---|---:|---:|
| Front matter, abstract, introduction | 1–72 | 17 | 3 |
| Results: sample and curves | 73–101 | 17 | 1 |
| Results: comparison with IHME and UN IGME | 102–169 | 14 | 0 |
| Discussion | 170–192 | 11 | 2 |
| Methods: data and covariates | 193–233 | 17 | 0 |
| Methods continued, supplement | 234–333 | 14 | 0 |

**The most important corrections:**
1. **Counts, abstract (L51) and Results (L76).** Replace the DHS-only counts (120 surveys, 36 countries; 95 surveys, 34 countries; 1,686,004 children; 75,726 deaths) with the v5 counts:
   - processed: 166 surveys in 40 countries
   - analysed: 132 surveys (95 DHS + 37 MICS) in 36 countries
   - sample: 2,292,089 children and 102,282 deaths
2. **Attributable fractions at 20% prevalence (L51).** Change to 5% neonatal, 30% at 6–11 months and 42–53% at 1–4 years. The peak at 40% prevalence (L91) is 60%, not 56%.
3. **Burden rates and declines (L52).** The rate falls from 10.6 to 4.5 to 3.4 per 1,000 child-years: a 57% decline to 2015, then a further 24%. Fill the XX placeholders with 49% (IHME) and 60% (UN IGME). Change "about a third higher" to about 40% (43% vs IHME, 39% vs UN IGME).
4. **Table 1 (L120–131) and Table 2 (L154–165).** Replace the model column in both. Table 2's top ten changes: DR Congo in, Uganda out.
5. **Text about Nigeria and Chad (L105, L137).**
   - Nigeria: 215,000 deaths, 56% above IHME and 30% above UN IGME.
   - Chad: 4.0-fold IHME.
   - States: the model exceeds IHME in 36 of 37 states (all except Lagos); median state ratio 1.7 in the northern zones.
6. **MICS is missing everywhere.** The heading, Figure 1 caption and Methods describe a DHS/MIS-only sample, and L76 says "the analysed sample is DHS only". The Methods need a MICS paragraph (L197–201). MICS also needs to be added to the weights (L222), the vaccination fallback (L225), the data-sharing statement (L257) and `ref.bib`.
7. **Methods L247 sensitivity counts.**
   - Subgroups: Sahel 31 surveys, 8 countries, 18,919 deaths; Eastern Africa 51 surveys, 12 countries, 31,589 deaths; median survey year 2014 (70 and 62 surveys).
   - Imputed analysis: 8,797,963 records, 166 surveys, 40 countries, 123,419 deaths.
8. **Supplement Table S1 (L319–327).** It is entirely DHS-only; replace it from `tables/age_band_results.latex.txt`.
9. **SText1 describes the wrong model.** Methods L238 points to SText1 for "full model equations". SText1 actually describes the archived region-level negative-binomial model (600 survey regions, 23 countries), not the seven cloglog age-band models.
10. **Other wrong statements and inconsistencies:**
    - L76: "within the previous 5 years" misstates the eligibility rule.
    - L85: the curves no longer "decrease" at high prevalence.
    - L100: the Sahel curve is steeper from 6 months, not 1 year.
    - L105: "mostly agreed" is not supported; the model is above IHME in 35 of 42 countries.
    - L142: the Figure 3 caption does not match the figure (2000–2024, three series).
    - L209: cites `igme` for the IHME data.
    - L225: Liberia and São Tomé do have survey data; they were excluded because they lack an HIV incidence series.
    - L242: estimates run "2004–2024"; they actually run 2000–2024.
    - L276: Figure S1 caption says HIV prevalence; the model uses HIV incidence.

---

## 1. Front matter, abstract and introduction (lines 1–72)

About 17 statements were checked and found correct, including:
- 7 age bands to 5 years; surveys since 2000; MAP as the \textit{Pf}PR source
- weak neonatal association (HR 20%→0% 0.95, 0.90–1.00); strong association from 6 months
- 17 covariates; 42 countries; the 24% decline 2015–2024 (23.75%)
- IHME and UN IGME post-2015 declines of 8% and 6%
- more attributable deaths than either source (611,440 vs 428,147 and 440,123)

### L51: survey counts and sources (outdated)
- **Current:** "We pooled all available Demographic and Health Surveys with complete birth histories since 2000 (120 surveys, 36 countries)"
- **Problem:** These are the DHS-only processed counts, and MICS is missing. v5 processes 166 surveys in 40 countries: 115 DHS, 5 MIS and 46 MICS, after dropping 5 Lesotho surveys. The registry has 171 surveys with complete birth histories in 41 countries, so "all available" needs a qualifier. L76 has the same outdated counts.
- **Evidence:** `study_flow/flow_counts.csv`: registry_with_cbh 171 (41 countries), omitted 5, built 166, built_countries 40, surveys 132, countries 36. The DHS-only v3 file has built 120 and built_countries 36.
```latex
We pooled all available Demographic and Health Surveys (DHS) and UNICEF Multiple Indicator Cluster Surveys (MICS) with complete birth histories conducted since 2000 in malaria-endemic countries (166 surveys, 40 countries) and matched individual child trajectories, discretised into 7 age bands until 5 years, with population-weighted \textit{Pf}PR$_{2-10}$ estimates from the Malaria Atlas Project for each survey region (mostly first-level administrative units).
```

### L51: children and deaths (outdated)
- **Current:** "in data from 1,686,004 children under 5 and 75,726 deaths"
- **Problem:** These are the DHS-only sample counts (95 surveys, 34 countries).
- **Evidence:** `primary_sample.csv`: 2,292,089 children; 7,498,459 records; 102,282 deaths; 132 surveys; 36 countries. `comparison_sample.csv` gives the previous 1,686,004 children and 75,726 deaths.
```latex
After adjusting for a series of potential confounding variables, a robust age-dependent association between \textit{Pf}PR$_{2-10}$ and all-cause mortality was observed in data from 2,292,089 children under 5 and 102,282 deaths in 132 surveys (95 DHS and 37 MICS) from 36 countries.
```

### L51: attributable fractions at 20% (outdated)
- **Current:** "the model attributes 6\% of neonatal deaths, 25\% of deaths at 6-11 months and 42-47\% of deaths at 1-4 years"
- **Evidence:** `attributable_fraction_by_age.csv` at pfpr_pct = 20:

| Age band | v5 fraction | v3 fraction (current text) |
|---|---:|---:|
| <1 | 0.053 | 0.062 |
| 1–5 | 0.136 | |
| 6–11 | 0.304 | 0.245 |
| 12–23 | 0.456 | 0.422 |
| 24–35 | 0.527 | 0.471 |
| 36–47 | 0.457 | |
| 48–59 | 0.421 | |

  The v3 range for 1–4 years was 0.422–0.471, which is where the current 42–47% comes from.
```latex
At 20\% prevalence, the model attributes 5\% of neonatal deaths, 30\% of deaths at 6-11 months and 42-53\% of deaths at 1-4 years to malaria.
```

### L52: burden rates and declines (outdated)
- **Current:** "fell from 10.1 to 4.2 deaths per 1,000 child-years between 2000 and 2015 (a 58\% decline), and then to 3.2 deaths per 1000 in 2024 (24\% decline)"
- **Problem:** These are the DHS-only rates (10.06, 4.20, 3.19). The sentence also lacks "child-years" after the 2024 figure and does not name IHME as the all-cause baseline.
- **Evidence:** `burden_comparison/source_comparison_2000_2024.csv`, model rate: 10.614 (2000), 4.519 (2015), 3.446 (2024). The 2015–2024 decline is 23.75%, and 1 − 4.519/10.614 = 57.4%.
- **Rounding note:** Table 1 prints 3.45. Writing "3.45" in the abstract instead of "3.4" would avoid an apparent rounding mismatch.
```latex
Applied to IHME estimates of national all-cause mortality in 42 sub-Saharan African countries, we estimate that under-5 malaria mortality fell from 10.6 to 4.5 deaths per 1,000 child-years between 2000 and 2015 (a 57\% decline), and then to 3.4 deaths per 1,000 child-years in 2024 (a further 24\% decline).
```

### L52: XX placeholders; declines refer to rates (placeholder / clarity)
- **Current:** "In contrast, IHME and UN IGME estimate substantial declines up until 2015  (XX and XX\% decrease, respectively) but little change since (8 and 6\% decrease)."
- **Problem:**
  - The 2000–2015 rate declines are 48.9% (IHME) and 60.1% (UN IGME).
  - UN IGME's pre-2015 decline is *larger* than the model's (57%), so "in contrast" applies only after 2015.
  - The 8% and 6% are rate declines. Death counts rose over 2015–2024 (IHME +2.9%, UN IGME +5.6%) because person-time grew by 12%.
- **Evidence:** `source_comparison_2000_2024.csv`, rates:
  - IHME: 5.132 → 2.624 → 2.413 (decline_from_2015 8.07)
  - UN IGME: 6.595 → 2.629 → 2.480 (decline_from_2015 5.65)
  - Deaths decline_from_2015: IHME −2.89, UN IGME −5.59
```latex
IHME and UN IGME estimate similarly large declines in the mortality rate up until 2015 (49\% and 60\% decrease, respectively) but, in contrast to our estimates, little change in the rate since (8\% and 6\% decrease).
```

### L52: "about a third higher" (outdated)
- **Problem:** "About a third" fitted the DHS-only ratios (1.32 and 1.29). The v5 ratios are 1.43 and 1.39, for both deaths and rates.
- **Evidence:** `source_comparison_2000_2024.csv`, 2024: model 611,440 deaths; IHME 428,147; UN IGME 440,123.
```latex
We estimate malaria mortality in 2024 was about 40\% higher than either source in the same set of countries (43\% higher than IHME and 39\% higher than UN IGME).
```

### L70: "extremely robust based on 25 years of DHS surveys" (outdated)
- **Problem:**
  - The sample is DHS plus MICS: 132 surveys in 36 countries, fielded 2003–2024.
  - "Extremely robust" overstates the evidence. The sensitivity analyses agree in direction but differ in size: the Sahel 40%→20% HRs are 0.70–0.78 against 0.82–0.91 overall, and Eastern Africa shows a null neonatal effect.
  - The association is age-dependent, so "consistent across age bands" would also be wrong.
- **Evidence:**
  - `survey_map/survey_timeline_all.csv`: included surveys fielded 2003–2024; 95 DHS + 37 MICS
  - `subgroups_dhsmics_map_gamma2_v2/REPORT.md`: HR table
  - `nutrition_adjustment_dhsmics_map_gamma2_v2/CONTRASTS.md`: largest change 0.008
```latex
We show that the relationship between malaria transmission as measured by \textit{Pf}PR$_{2-10}$ and age-specific all-cause mortality is strong from six months of age and consistent in direction across sensitivity analyses, based on 132 DHS and MICS surveys conducted in 36 African countries between 2003 and 2024.
```

### L70: "considerably greater than previously reported" (clarity)
- **Problem:** The claim holds for the absolute decline but not fully for the relative one:
  - 2004–2024 absolute decline: model −5.36 per 1,000 child-years, against −2.35 (IHME) and −2.96 (UN IGME).
  - 2004–2024 relative decline: 60.9% vs 49.3% and 54.4%, a modest gap.
  - 2000–2015 relative decline: the model's 57.4% is *smaller* than UN IGME's 60.1%.
  - The large relative gap is after 2015: 23.8% vs 8.1% and 5.7%.
  - Grammar: "in last 2 decades".
- **Evidence:** `source_comparison_2000_2024.csv`, rate rows: decline_from_2004 60.87 / 49.29 / 54.43; decline_from_2015 23.75 / 8.07 / 5.65.
```latex
These data suggest that the absolute reduction in childhood malaria mortality over the last two decades, and the relative reduction since 2015, have been greater than currently estimated, and that the number of malaria-attributable deaths is higher.
```

### L50 and L68: "directly measurable quantities" / "most reliable inputs" (clarity)
- **Problem:** Neither input is used as a direct measurement:
  - Exposure is MAP model-based \textit{Pf}PR$_{2-10}$, from annual rasters weighted by GPW density.
  - The national burden applies the fitted effects to IHME *modelled* all-cause deaths.
  - Survey-measured mortality is used only to fit the model.
- **Evidence:** `docs/ANALYSIS_PLAN.md` §2.2; `R_mics/06_polygons_pfpr_registry.R`; `annual_comparison/CAPTION.md`.
```latex
% L50
... using two quantities that are routinely measured in nationally representative household surveys: the prevalence of infection in children (\textit{Pf}PR$_{2-10}$, a proxy for transmission intensity, here taken from Malaria Atlas Project model-based estimates) and all-cause childhood mortality.
% L68, last sentence
Because prevalence of infection and all-cause childhood mortality are measured frequently and at wide geographic scale in representative surveys, they are more reliable inputs than cause-of-death assignment by verbal autopsy.
```

### L51: "admin level-1 ... population weight prevalence" (clarity)
- **Problem:**
  - "Population weight" should be "population-weighted".
  - The units are survey regions, not strictly admin-1. Nigeria uses the six geopolitical zones; Bissau is merged with Biombo, and Bangui with Ombella-M'Poko.
- **Evidence:** `docs/ANALYSIS_PLAN.md` §2.2; `data/derived_mics/region_map.csv` (notes on merged capitals).
```latex
... with population-weighted \textit{Pf}PR$_{2-10}$ estimates from the Malaria Atlas Project for each survey region (mostly first-level administrative units; the six geopolitical zones in Nigeria).
```

### L65: UN IGME expansion (clarity)
- **Current:** "the openly available geospatial models of the Institute for Health Metrics and Evaluation (IHME) and the United Nations (UN IGME)"
- **Problem:** UN IGME is the Inter-agency Group for Child Mortality Estimation. The comparator is the WHO/UNICEF CA-CODE 2026 series, which L210 and L242 already call "UN IGME (CA-CODE)". Whether these models are "geospatial" is not assessed here.
- **Evidence:** `annual_comparison/who_source.json`, source_label "CA-CODE / UN IGME (WHO and UNICEF)".
```latex
The most influential and oft cited estimates are from the World Health Organization (WHO) yearly World Malaria Reports and the cause-of-death models of the Institute for Health Metrics and Evaluation (IHME; Global Burden of Disease study) and of the WHO and UNICEF Child and Adolescent Causes of Death Estimation (CA-CODE) group, disseminated by the United Nations Inter-agency Group for Child Mortality Estimation (UN IGME).
```

### L68: [REF] placeholder and "\textit{Pf}PR$_{2-10}$>20\%" (placeholder)
- **Problem:**
  - `[REF]` is unfilled.
  - The `>` in text mode renders correctly under T1 font encoding, but math mode is cleaner.
```latex
In addition to chronic anaemia, malaria increases the risk of bacterial sepsis \cite{scott} and impairs child growth \cite{<add reference>}. ... In areas of high transmission (e.g.\ \textit{Pf}PR$_{2-10}$ $>$20\%), ...
```

### L68: transmission proxy stated backwards (clarity)
- **Current:** "Transmission intensity, a proxy for inoculation rates, can be tracked by the proportion of children who are infected"
- **Problem:** The EIR is the direct measure of transmission, and \textit{Pf}PR$_{2-10}$ is the proxy. L50 has this the right way round.
```latex
Transmission intensity, usually quantified by the entomological inoculation rate, can be tracked by the proportion of children who are infected at any given point, ...
```

### L50, L68: hyphenation (typo)
- "all cause" is used at L50 and L68; "all-cause" is used at L51 and L70 and in all pipeline captions. Use "all-cause childhood mortality" throughout.

### L70: tautology (typo)
```latex
We propose an alternative approach to estimating the number of deaths caused by malaria, combining direct and indirect effects.
```

### L47: stale word count (clarity)
- **Current:** "% 210 words"
- **Problem:** Lines 50–52 contain 316 tokens, about 310 words once markup is discounted. That is above the 300-word limit for a Lancet-style abstract. Update the comment and trim the abstract.

### L28, L38: author list and affiliation (placeholder)
- **Problem:**
  - The author list contains a `[...]` placeholder.
  - "Robert W Snow" lacks the full stop used in the other names.
  - The Oxford unit name cannot be checked from the pipeline; please confirm it.
```latex
James A. Watson$^{1,2,3,*}$, Dhruv Darji$^{1,2}$, <complete author list>, Robert W. Snow$^{2,4}$
\small\textbf{2} Centre for Tropical Medicine and Global Health, Nuffield Department of Medicine, University of Oxford, Oxford, UK;\\
```

---

## 2. Results: sample and curves (lines 73–101)

About 24 statements were checked and found correct, including:
- 7 bands and the band definitions; 37% of deaths are neonatal (36.8%)
- 17 covariates; weak neonatal association
- monotone increase at 1–11 months
- at 10% PfPR: 2.5% in neonates, 27–35% at 12–59 months ("approximately 30%")
- no-nutrition refit essentially unchanged (largest change 0.008)
- all figure files match the v5 export

### L74: subsection heading (outdated)
- **Problem:** The analysed sample is 95 DHS + 37 MICS surveys. It contains no MIS surveys: all 5 were excluded.
- **Evidence:** `survey_map/survey_timeline_all.csv`.
```latex
\subsection*{Prevalence and all-cause mortality in DHS and MICS surveys}
```

### L76: registry and processed counts (outdated)
- **Current:** "microdata from Demographic and Health Surveys (DHS) and Malaria Indicator Surveys (MIS) ... 120 surveys across 36 countries, representing 1,114 survey regions"
- **Evidence:**
  - `study_flow/flow_counts.csv`: registry 171 surveys in 41 countries (124 DHS/MIS + 47 MICS); built 166 in 40 countries.
  - `primary_map_regional17_dhsmics_imputed_gamma2_v6/prepared_sample.csv` (all MAP-eligible records of those 166 surveys): 1,457 regions. The DHS-only analogue was 1,113.
```latex
We pooled all available microdata from Demographic and Health Surveys (DHS), Malaria Indicator Surveys (MIS) and UNICEF Multiple Indicator Cluster Surveys (MICS) conducted in sub-Saharan Africa since the year 2000 which recorded complete birth histories (124 DHS/MIS and 47 MICS surveys in 41 countries). Excluding Lesotho, which is malaria free and has no MAP prevalence estimates, this gave a total of 166 surveys (120 DHS/MIS and 46 MICS) across 40 countries, representing 1,457 survey regions, Figure \ref{fig:map}.
```

### L76: "admin level-1" units and "estimated admin level-1 PfPR estimates" (wrong)
- **Problem:** The units are each survey's own regions:
  - Nigeria: 6 zones in every survey
  - Senegal continuous DHS 2012–2016, 2018 and 2019: 4 macro-regions (2017 has 14)
  - Tanzania: 8–9 zones
  - Kenya 2003–2014: 8 provinces
  - some Kenya MICS surveys: a single region each

  The \textit{Pf}PR value is also not a MAP-published admin-1 estimate. The pipeline computes a GPW-weighted mean of the annual raster over each region polygon, for the year of band entry.
- **Evidence:** `survey_map/survey_coverage.csv`; `docs/ANALYSIS_PLAN.md` §2.2.
```latex
... representing 1,457 survey regions (the subnational regions of each survey, usually the first administrative level but the six geopolitical zones in Nigeria), Figure \ref{fig:map}. ... Each age band was then linked to the population-weighted mean \textit{Pf}PR$_{2-10}$ of its survey region from the annual Malaria Atlas Project (MAP) surfaces for the calendar year in which the child entered the band \cite{map2025}.
```

### L76: "events ... included if they occurred within the previous 5 years" (wrong)
- **Problem:** This is not the rule the pipeline uses:
  - A band is eligible if the child *entered* it within the 60 months before interview.
  - The whole band must end by interview, for deaths and survivors alike.
  - The entry year must be 2000–2024.

  Children born more than 5 years before interview still contribute their later bands.
- **Evidence:**
  - `R_cbh/R/build.R`: `band_end_cmc <= interview_cmc`
  - `study_flow/flow_counts.csv`: entered 10,657,524; incomplete 1,347,434; outside_year 418,661
```latex
A child contributed each age band that it entered alive within the 60 months before interview, provided the whole band had ended by the date of interview and the band began in 2000--2024.
```

### L76: "compete birth histories from the DHS/MIS" (typo)
- **Problem:**
  - "Compete" should be "complete".
  - MICS birth histories are also used.
  - The band split is not how age at death is recorded.
```latex
Each child trajectory from birth was divided into discrete age bands of $<$1 month (neonatal), 1--5 months, 6--11 months, and then yearly up to 59 months, reflecting the precision with which age at death is recorded in the complete birth histories of the DHS, MIS and MICS (days in the first month, months up to age 2 years and years thereafter).
```

### L76: "95 surveys in 34 countries" and complete-case wording (outdated; the wording issue was missed by the auditor)
- **Current:** "restricting the analysis to complete cases with identified confounders recorded, a total of 95 surveys in 34 countries were included"
- **Problem:** 95 surveys in 34 countries is the DHS-only sample; v5 has 132 surveys (95 DHS + 37 MICS) in 36 countries (`prepared_sample_dhs_mics.csv`). Complete-case selection is applied only *after* several substitutions:
  - national WUENIC values replace missing regional DTP3 and measles values
  - the survey's mean of other regions fills remaining regional gaps
  - HIV incidence uses a fixed imputation
- **Evidence:** `study_flow/CAPTION.md`: "1,299,504 lack at least one required covariate after the declared HIV, vaccination and available-region substitutions".
```latex
Following the merge with MAP prevalence estimates and restricting the analysis to records with all 17 adjustment covariates available (after substituting national WHO/UNICEF coverage for missing regional DTP3 and measles estimates, the mean of the survey's other regions for remaining regional gaps, and a model-based imputation of national child HIV incidence), a total of 132 surveys (95 DHS and 37 MICS) in 36 countries were included in the primary analysis (Figure \ref{fig:flow}).
```

### L76: "All five MIS surveys and 20 DHS surveys were excluded ... so the analysed sample is DHS only" (outdated)
- **Problem:** 34 surveys contribute no records, not 25: 9 MICS surveys are also excluded. The analysed sample is no longer DHS only.
- **Evidence:** `study_flow/CAPTION.md`: "25 DHS or MIS surveys (including all 5 MIS surveys ...) and 9 MICS surveys ... combines 95 DHS and 37 MICS surveys".
```latex
Thirty-four surveys were excluded because a required covariate was unavailable for every region of the survey: 25 DHS/MIS surveys, including all five MIS surveys (which do not publish the facility-delivery and anthropometry indicators), and nine MICS surveys, so the analysed sample combines 95 DHS and 37 MICS surveys.
```

### L76: children and deaths (outdated)
- **Evidence:** `prepared_sample_dhs_mics.csv`: 7,498,459 records; 2,292,089 children; 102,282 deaths; 1,211 regions.
```latex
The included surveys contributed 7,498,459 child age-band records from 2,292,089 children in 1,211 survey regions, with 102,282 deaths recorded.
```

### L76: "65\% in the first year of life, and 14\% between 1 and 2 years" (outdated)
- **Problem:**
  - The first-year share is 65.8%, which rounds to 66%.
  - The 12–23 month share is 12.8%, which rounds to 13%. Even in the DHS-only version the value was 13.49%, so 14% came from double rounding.
- **Evidence:** `tables/age_band_results.csv`: first-year deaths 37,623 + 15,586 + 14,126 = 67,335 of 102,282; 12–23 months 13,138.
```latex
37\% of these deaths were in neonates (first month of life); 66\% in the first year of life, and 13\% between 1 and 2 years of age.
```

### L81: Figure 1 caption (outdated)
- **Problem:** The figure itself is current (it matches the manifest MD5 616ccc01…). The caption is not:
  - It mentions neither MICS nor the included/excluded symbols.
  - It says countries are labelled with ISO3 codes, but the map has no labels; ISO3 codes label only the timeline rows.
- **Evidence:** `paper_figures/CAPTIONS.md` (Figure 1); `survey_map/survey_timeline_all.csv`.
```latex
\caption{Geographic coverage and timing of the 124 DHS and MIS surveys and 47 MICS surveys with complete birth histories conducted since 2000 in 41 countries of sub-Saharan Africa. Circles are DHS and MIS surveys and triangles MICS surveys; filled symbols are the 95 DHS and 37 MICS surveys that contribute to the primary analysis (36 countries) and open symbols the 39 excluded surveys: 5 in Lesotho, which is malaria free and has no MAP prevalence estimates, and 34 for which a required covariate was unavailable for every region (20 DHS, all 5 MIS and 9 MICS surveys). Map shading gives the number of included surveys in each country (grey: none); timeline rows are labelled with ISO3 country codes. The 46 MICS surveys without a complete birth history are not shown.}
```

### L85: "17 socio-economic and health-system covariates" (clarity)
- **Problem:** The count of 17 is right, but the set also includes nutritional measures (wasting, stunting), governance (political stability) and child HIV incidence. The sentence also omits the calendar-year spline and the random effects.
- **Evidence:** `model_formula.txt`.
```latex
adjusted for 17 survey-region and national covariates (socio-economic, health-system, nutritional, governance and child HIV incidence measures), a smooth function of calendar year, and survey, country and region random effects.
```

### L85: "plateaued above ... 30\% and even decreased for children over 3 years" (outdated)
- **Problem:** The v5 curves do not decline at high prevalence:

| Age band | Log HR vs 20% prevalence |
|---|---|
| 36–47 months | peaks at 0.179 at 48.5%; 0.170 at 63% |
| 48–59 months | 0.097 at 40%; peaks at 0.098 at 44%; 0.095 at 63% |
| 12–23 and 24–35 months | still rising above 30%, only more slowly |

  The decline described in the text belongs to the v3 fit (36–47 months: 0.178 → 0.150).
- **Evidence:** `pfpr_curves.csv` over the central support. The Figure S2 image in Overleaf is the v5 export and shows the same shapes.
```latex
However, above prevalences of around 30\% the increase in risk of death in older children slowed markedly, and for children over 3 years of age the curves were essentially flat above prevalences of about 40--50\%.
```

### L91: "peaking at 56\% at 2–3 years" (outdated)
- **Problem:**
  - The v5 attributable fraction at 40% prevalence peaks at 59.8% (54.0–64.8%) at 24–35 months.
  - At 48–59 months it is 47.5%, below half, so "1 and 4 years" needs narrowing to 12–47 months.
- **Evidence:** `attributable_fraction_by_age.csv` at 40%: 0.553, 0.598, 0.542, 0.475 for 12–23 through 48–59 months. v3 gave 0.564 at 24–35 months.
```latex
For an average \textit{Pf}PR$_{2-10}$ of 40\% (upper end of prevalence estimates in sub-Saharan Africa now), there is a sharper age-dependent increase in malaria deaths, with more than half of all deaths between 1 and 4 years of age (12--47 months) attributable to malaria, peaking at 60\% at 2--3 years.
```

### L91: bare \ref and \textit{PfPR} style (typo)
- The numbers in this sentence are correct for v5.
```latex
Under this model, for an average \textit{Pf}PR$_{2-10}$ of 10\%, approximately 30\% of all deaths in children 1 year and older are caused by malaria, compared with less than 5\% in neonates (Figure~\ref{fig:main_result2}).
```

### L95: Figure 2 caption (clarity)
- **Problem:** The caption says "mean" estimates and omits the interval bars, the four prevalence levels, and the fact that the zero counterfactual lies outside the observed data (extrapolation). The figure matches the manifest MD5 (89794f9b…).
```latex
\caption{Estimated malaria-attributable share of all-cause deaths within each age band, $1-\exp\{f_g(0)-f_g(P)\}$, at \textit{Pf}PR$_{2-10}$ ($P$) of 10\%, 20\%, 30\% and 40\%, from the seven separate age-band models. Points are point estimates and vertical bars conditional pointwise 95\% intervals; lines connect the discrete age bands. Zero prevalence lies below the observed exposure range, so the counterfactual involves extrapolation.}
```

### L100: subgroup results (wrong)
- **Problem:**
  - The Sahel curves are steeper from 6 months, not from 1 year. At 1–5 months the Sahel matches the full sample (0.91 vs 0.92).
  - Eastern Africa is flatter in neonates (1.00, 0.96–1.04) and at 36–47 months (0.95).
  - "Sahel versus East Africa" suggests the two are complementary subsets; they are not.
  - The period split is at median survey year 2014.
- **Evidence:** 40%→20% HRs from `subgroups_dhsmics_map_gamma2_v2/REPORT.md`, with the split year from `subgroup_definitions.csv`:

| Age | Full sample | Sahel | Eastern Africa |
|---|---:|---:|---:|
| <1 | 0.94 | 0.90 | 1.00 |
| 1–5 | 0.92 | 0.91 | 0.90 |
| 6–11 | 0.83 | 0.73 | 0.85 |
| 12–23 | 0.82 | 0.70 | 0.81 |
| 24–35 | 0.85 | 0.73 | 0.84 |
| 36–47 | 0.84 | 0.78 | 0.95 |
| 48–59 | 0.91 | 0.74 | 0.88 |
```latex
We conducted a series of sensitivity analyses looking at subgroups of the surveys (Sahel, Eastern Africa, and surveys conducted up to or after the median survey year, 2014). In general the subgroup fits gave very similar results, apart from survey regions in the Sahel (boundary centroid at or north of 12$^\circ$N and west of 36$^\circ$E, excluding the Horn of Africa), where transmission is highly seasonal and the estimated relationship between prevalence and all-cause mortality was steeper for children between 6 months and 5 years of age (hazard ratios for a reduction in \textit{Pf}PR$_{2-10}$ from 40\% to 20\% of 0.70--0.78, against 0.82--0.91 in the full sample), and Eastern Africa, where the curves were flatter in neonates and at 36--47 months (Figure \ref{fig:subgroups}).
```

### L100: multiple-imputation sentence (clarity; missed by the auditor)
- **Problem:** "Very similar" is fair: the largest change in the 40%→20% HR is 0.040, at 48–59 months. But the sentence undersells the analysis:
  - It retains 166 surveys in 40 countries.
  - The figure shows the point imputation.
  - The Rubin's-rules pooling is a check, and it changed the HRs by at most 0.002.
- **Evidence:** `primary_map_regional17_dhsmics_imputed_gamma2_v6/REPORT.md`.
```latex
A sensitivity analysis that imputed the missing covariates, retaining all 166 surveys in 40 countries (8,797,963 child age-band records and 123,419 deaths), gave very similar results to the complete-case primary analysis: hazard ratios for a reduction in \textit{Pf}PR$_{2-10}$ from 40\% to 20\% changed by at most 0.04, and pooling fits to 10 multiply imputed datasets with Rubin's rules changed them by at most 0.002 (Figure \ref{fig:MI}).
```

---

## 3. Results: comparison with IHME and UN IGME (lines 102–169)

About 58 statements were checked and found correct, including:
- the method description at L104
- the 42 countries
- all IHME and UN IGME cells in Tables 1 and 2, and both comparator totals
- the Nigeria comparators (138,111 and 165,329)
- southern-zone median state ratio 1.2

### Figure 3 is v5 but the text and tables are v3 (internal inconsistency)
- **Problem:** The Figure 3 image is the v5 export (MD5 25e87587… matches `paper_figures/manifest.csv`), so the text and tables contradict it:
  - Panel C runs from about 10.6 to about 3.4; Table 1 says 10.06 and 3.19.
  - Panel B shows only Lagos below the equality line; the text says 35 of 37 states.
- **Fix:** Apply the replacements below together with the abstract fix (L52) and the Methods fix (L242, 2004 → 2000).

### L104: "the 42 countries with MAP coverage" (clarity)
- **Problem:** The 42 are the 45 sub-Saharan locations in the IHME export, minus Cape Verde, Lesotho and São Tomé and Príncipe. South Africa is covered by MAP but is not in the IHME export. Panel A compares the model with IHME only.
- **Evidence:** `annual_comparison/CAPTION.md`; `annual_comparison/included_countries.csv` (42 rows).
```latex
... population-weighted national \textit{Pf}PR$_{2-10}$ in the 42 sub-Saharan African countries in the IHME export with MAP estimates, for the year 2024, and compared the results with the cause-specific under-5 malaria death estimates of IHME and of UN IGME, expressed per 1,000 under-5 child-years (Figure~\ref{fig:IHME_v_model}A; Table~\ref{tab:country}).
```

### L104: zero-prevalence counterfactual is an extrapolation (clarity; missed by the auditor)
- **Problem:** Zero prevalence lies below the observed range in every age band. The pipeline flags this everywhere; `main.tex` never mentions it.
- **Evidence:**
  - `REPORT.md` line 56: "Zero PfPR is just below observed support in every band"
  - `attributable_fraction_by_age.csv`: zero_below_observed_support = TRUE in every row
```latex
... can be estimated by applying the age-band hazard ratios for a reduction from the observed prevalence to zero to national all-cause death counts (zero prevalence lies just below the range of survey-region \textit{Pf}PR$_{2-10}$ observed in every age band, so this counterfactual involves some extrapolation).
```

### L105: "mostly agreed" (wrong)
- **Problem:**
  - The model rate is above IHME's in 35 of 42 countries (median ratio 1.56; ratio >1.5 in 22 countries, >2 in 13). Only 9 countries lie within 0.8–1.25.
  - These are malaria mortality rates, not under-5 mortality rates.
  - "PfPR ACM" should be "PfPR-ACM model".
- **Evidence:** `burden_comparison/country_rates_2024.csv`. The Spearman correlation of 0.91 was computed from this file and is not a saved output; drop it if every number must trace to a file.
```latex
Although national under-5 malaria mortality rates from IHME and our approach (PfPR-ACM model) were strongly correlated across the 42 countries (Spearman $\rho = 0.91$), the model estimate exceeded IHME's in 35 of the 42 countries (median ratio 1.56), with the largest absolute differences in several populous, high-burden countries (Table~\ref{tab:country}).
```

### L105: Nigeria (outdated)
- **Evidence:** `tables/country_comparison_2024.csv`, NGA: model 214,938; IHME 138,111; UN IGME 165,329. Ratios 1.556 and 1.300. The v3 model value was 200,988.
```latex
In Nigeria, IHME estimates approximately 138 thousand deaths and UN IGME approximately 165 thousand deaths, whereas our approach estimates 215 thousand deaths (56\% and 30\% higher, respectively).
```

### L105: Chad (outdated)
- **Evidence:** TCD: model 19,407; IHME 4,884; UN IGME 14,345. Ratios 3.97 and 1.35. The v3 ratio was 3.58.
```latex
In Chad we estimate a 4.0-fold higher total burden than IHME (1.4-fold higher than UN IGME).
```

### L105: bare \ref and "IGME" (typo)
- Write `(Table~\ref{tab:country})` and "UN IGME". Both are already covered if you use the replacement sentences above.

### L120–124: Table 1, deaths block, model column (outdated)
- **Evidence:** `burden_comparison/source_comparison_2000_2024.csv`; `burden_comparison/source_comparison.latex.txt`. The current values are the v3 values exactly.
```latex
\quad 2000 & 1,165,686 & 563,657 & 724,276  \\
\quad 2015 & 716,520 & 416,122 & 416,815  \\
\quad 2024 & 611,440 & 428,147 & 440,123  \\
\quad Change 2000--2024 & $-47.5\%$ & $-24.0\%$ & $-39.2\%$  \\
\quad Change 2015--2024 & $-14.7\%$ & $+2.9\%$ & $+5.6\%$ \\
```

### L127–131: Table 1, rates block, model column (outdated)
```latex
\quad 2000 & 10.61 & 5.13 & 6.59  \\
\quad 2015 & 4.52 & 2.62 & 2.63  \\
\quad 2024 & 3.45 & 2.41 & 2.48  \\
\quad Change 2000--2024 & $-67.5\%$ & $-53.0\%$ & $-62.4\%$  \\
\quad Change 2015--2024 & $-23.8\%$ & $-8.1\%$ & $-5.7\%$  \\
```

### L126 and after L133: rate unit and table note (clarity)
- **Problem:**
  - The rate unit leaves out the denominator.
  - There is no note on the shared denominator or on the different estimands.
  - For IHME and UN IGME, counts rose over 2015–2024 while rates fell, because person-time grew by 11.9% (158.6 million to 177.5 million child-years). The model's counts and rates both fell.
- **Evidence:** `annual_comparison/annual_totals_2000_2024.csv`, under5_person_years.
```latex
\multicolumn{4}{@{}l}{\textit{Malaria mortality rate (per 1,000 under-5 child-years)}} \\
% after \end{tabular*}:
\par\vspace{0.5em}\begin{minipage}{\linewidth}\footnotesize All three columns are deaths before age five in the same 42 countries. Rates use the same annual under-five person-years, implied by the IHME all-cause death count and rate. Model deaths are IHME all-cause deaths in each age band multiplied by the estimated malaria-attributable fraction; IHME and UN IGME (CA-CODE 2026) are cause-specific malaria death estimates, so these are different estimands. IHME and UN IGME death counts rose between 2015 and 2024 while their rates fell, because under-five person-time in these countries grew by 12\% over that period.\end{minipage}
```

### L137: Nigerian states (outdated)
- **Problem:**
  - The model exceeds IHME in 36 of 37 units; Lagos is the exception (ratio 0.70). v3 had 35 of 37.
  - The northern-zone median ratio is 1.73, which rounds to 1.7, not 1.6.
  - The Federal Capital Territory (FCT) is not a state.
- **Evidence:** `nigeria_states/README.md`: "exceeds IHME in 36 of 37 states ... median 1.73 ... 1.24"; state sum 214,127 vs 132,138 (ratio 1.62).
```latex
Within Nigeria, applying the same age-band effects to state-level IHME all-cause deaths by age and state \textit{Pf}PR$_{2-10}$ gives malaria mortality estimates greater than IHME's in every one of the 36 states and the Federal Capital Territory except Lagos (Figure~\ref{fig:IHME_v_model}B); summed over states, the model gives 214,127 deaths against 132,138 IHME malaria deaths (ratio 1.62). The differences are concentrated in the northern states: the median state ratio is 1.7 in the three northern zones and 1.2 in the three southern zones.
```

### L142: Figure 3 caption (wrong)
- **Problem:** The caption does not match the figure:
  - Panel C covers 2000–2024 (not "since 2004") and shows three series.
  - The rates are pooled across the 42 countries, not averaged over sub-Saharan Africa.
  - The caption says "prevalence-ACM model" instead of "PfPR-ACM model", and omits UN IGME.
- **Evidence:** `paper_figures/CAPTIONS.md` (Figure 3); image MD5 matches the manifest.
```latex
\caption{Malaria mortality before age 5 from the PfPR-ACM model compared with IHME and UN IGME (CA-CODE 2026) cause-specific estimates, expressed as deaths per 1,000 under-five child-years. A: model against IHME national rates in 2024 in the 42 countries; the dashed line denotes equality. B: model against IHME for the 36 Nigerian states and the Federal Capital Territory in 2024, coloured by geopolitical zone; summed over states the model gives 214,127 deaths against 132,138 IHME malaria deaths (ratio 1.62). C: annual under-five malaria mortality, 2000--2024, for the PfPR-ACM model, IHME and UN IGME, pooled across the same 42 countries. In every panel each source's deaths are divided by the under-five person-years implied by the IHME all-cause death counts and rates, so the sources share a denominator. Model deaths are IHME all-cause deaths in each age band multiplied by $1-\exp\{f_g(0)-f_g(P)\}$. Model-attributable reductions in all-cause mortality and cause-specific malaria deaths are different estimands, and the model shares the IHME all-cause inputs.}
```

### L154–165: Table 2 (outdated)
- **Problem:** The model column and the top-ten selection are v3.
  - DR Congo enters at rank 6.
  - Uganda falls out (to rank 21).
  - The order changes below rank 6.
  - Mali (difference 7,309) and DR Congo (7,305) are only 4 deaths apart, so their order is fragile.
- **Evidence:** `tables/country_comparison_2024.latex.txt`.
```latex
Nigeria & 214,938 & 138,111 & 165,329 \\
Niger & 36,990 & 20,343 & 23,877 \\
Chad & 19,407 & 4,884 & 14,345 \\
Mozambique & 18,975 & 10,285 & 10,315 \\
Mali & 20,254 & 12,945 & 6,984 \\
DR Congo & 69,080 & 61,774 & 105,372 \\
South Sudan & 11,086 & 4,779 & 3,980 \\
Cameroon & 22,345 & 16,295 & 6,286 \\
Angola & 20,520 & 14,652 & 4,734 \\
Central African Republic & 9,138 & 3,865 & 5,220 \\
\midrule
All 42 countries & 611,440 & 428,147 & 440,123 \\
```

### After L167: Table 2 note (clarity; missed by the auditor)
- **Evidence:** The generated `tables/country_comparison_2024.latex.txt` includes this note.
```latex
\par\smallskip\begin{minipage}{\linewidth}\footnotesize UN IGME denotes the under-five CA-CODE 2026 cause-specific malaria series. Totals include all 42 countries, not only the ten shown. Death counts are rounded after estimation.\end{minipage}
```

---

## 4. Discussion (lines 170–192)

About 9 statements were checked and found correct, including:
- the plateau at about 600,000 (WMR global 578,000 in 2015 and 610,000 in 2024)
- "predominantly young children"
- DRC prevalence above 30% (36.1%)
- 89% of recorded deaths occurring before age five (88.8%)
- a larger impact than malaria-specific deaths suggest (611,440 vs 428,147 and 440,123)

### L174: "approximately a XX percent decrease" (placeholder)
- **Problem:** The placeholder is unfilled, the source and measure are unstated, and "sub-Saharan African" should be "sub-Saharan Africa".
- **Evidence:** `source_comparison_2000_2024.csv`: WMR African Region 804,000 → 579,000 (−28.0%); IHME −24.0%; UN IGME −39.2%.
```latex
Compared with the year 2000, WHO estimates of malaria deaths (all ages) in the WHO African Region for 2024 represent a decrease of approximately 28\% (from 804,000 to 579,000) \cite{wmr}, and IHME and UN IGME estimates of under-5 malaria deaths in the 42 sub-Saharan African countries analysed here fell by 24\% and 39\%, respectively (Table~\ref{tab:source-comparison}).
```

### L175 (commented out): wrong numbers if restored
- **Problem:**
  - 580,000 and 28% are WHO African Region all-age figures, measured from 2000 (from 2004 the decline is 22.9%).
  - Under-5 person-time grew by 62%, not 85% (109.8 to 177.5 million child-years).
  - The rates are per child-year, not per 1,000 births.
  - `rowe2006` estimates the burden in 2000; it is not a source for these figures.
```latex
%UN IGME and IHME estimate that under-5 malaria deaths in the 42 countries fell from 724,276 and 563,657 in 2000 to 440,123 and 428,147 in 2024 (39\% and 24\% decreases) [REF: IHME GBD; UN IGME CA-CODE]. Over the same period under-5 person-time grew from about 110 to 177 million child-years (62\% increase), so the malaria mortality rate per 1,000 child-years fell by 62\% (UN IGME) and 53\% (IHME).
```

### L176: the "challenge" is phrased as the narrative itself (clarity)
```latex
This narrative of disappointing progress has recently been challenged: have we really failed to roll back malaria \cite{white:22}?
```

### L182 (commented out): DRC one third applies to all under-5 deaths
- **Evidence:** `annual_comparison/country_age_estimates_2000_2024.csv`, COD 2024:
  - all under-5 deaths: 69,080 of 207,594 attributable (33.3%)
  - post-neonatal deaths: 63,436 of 150,513 (42.1%)
  - prevalence 36.1%
```latex
%In the Democratic Republic of the Congo, where population-weighted \textit{Pf}PR$_{2-10}$ was 36\% in 2024, the model attributes one third of all under-5 deaths (69,080 of the 207,594 estimated by IHME) and 42\% of post-neonatal under-5 deaths (63,436 of 150,513) to malaria.
```

### L173: `gbd` citation cannot support a plateau to 2024 (citation)
- **Problem:** The `gbd` entry is GBD 2021 (1990–2021).
- **Evidence:** `ref.bib`.
```latex
... has plateaued at approximately 600,000 per year (610,000 in 2024) \cite{wmr}.
```

### L188: "Mechanistically we think that higher the malaria transmission, the greater the absolute number of deaths" (clarity)
- **Problem:** The grammar is off ("higher the"), and the claim is strictly monotone. The fitted curves level off at high prevalence in the older bands (downturns under 0.01 on the log-HR scale).
- **Evidence:** `pfpr_curves.csv`.
```latex
Mechanistically, we expect all-cause mortality to increase with malaria transmission intensity, although in older children the fitted relationship levels off at high prevalence (Figure~\ref{fig:main_result}).
```

### L188: "reduced all cause mortality by around 20\%" (imprecise; checked against the source)
- **Problem:** Pryce et al. report 17% for both outcomes, so "around 20%" is loose and "substantially" overstates a 17% prevalence reduction.
- **Evidence:** `data/Pryce_ITN_systematic_review.pdf`, ITNs vs no nets:
  - all-cause child mortality: RR 0.83 (0.77–0.89)
  - \textit{P. falciparum} prevalence: RR 0.83 (0.71–0.98)
```latex
This is supported by large randomised trials of insecticide-treated nets, which reduced \textit{P. falciparum} prevalence and all-cause child mortality each by about 17\% compared with no nets \cite{pryce}.
```

### L188: "neonatal deaths (one third of all deaths)" (internal inconsistency)
- **Problem:** The Results (L76, L85) say 37%, and the v5 value is 36.8%. The sentence also uses the \textit{PfPR} style.
- **Evidence:** `tables/age_band_results.csv`: HR 40%→20% 0.944 (0.916–0.973).
```latex
The weak association between \textit{Pf}PR$_{2-10}$ and neonatal death (37\% of deaths in the analysed sample; hazard ratio 0.94 for a reduction in \textit{Pf}PR$_{2-10}$ from 40\% to 20\%) suggests that residual confounding is likely to be small.
```

### L188: "because DHS data are the most reliable for children under 5" (outdated)
```latex
We limited our analysis to children under 5 because this is the most important subpopulation, but also because birth-history data from DHS and MICS surveys are the most reliable for children under 5.
```

### L188: "and the 89\% of all reported deaths" (typo)
- **Evidence:** 660,139 of 743,122 recorded deaths (88.8%). This combines v3 `study_flow/recorded_death_ages_by_survey.csv` (120 DHS/MIS surveys) with `results/mics_inventory/survey_inventory.csv` (46 built MICS surveys).
```latex
Age at death is recorded in days for neonates and in months up to age two, and 89\% of all deaths recorded in the DHS, MIS and MICS birth histories occurred before age five.
```

### L191: draft notes appear in the compiled PDF (clarity)
- **Problem:** L190 is commented out but L191 is not, so the PDF shows:
  - a lowercase fragment
  - the undefined abbreviation U5M
  - an uncited historical claim
```latex
Malaria infection is a risk factor for child death that extends beyond a single assigned cause. Both infection prevalence and all-cause under-5 mortality can be measured with little ambiguity in representative household surveys. Our estimates imply that reducing exposure to malaria parasites prevents more child deaths than cause-specific estimates suggest: in 2024 the model attributes 611,440 under-5 deaths in 42 countries to malaria, compared with 428,147 (IHME) and 440,123 (UN IGME) malaria deaths, consistent with historical malaria control and elimination experiments in Africa more than 50 years ago [REF].
```

---

## 5. Methods: data sources and covariates (lines 193–233)

About 42 statements were checked and found correct, including:
- the 7 band definitions
- the 4 recode and 9 API covariates for DHS
- the WUENIC fallback (L201)
- the 13 regional and 4 national covariates listed at L219
- log transforms; centring and scaling; entry-year matching of national covariates
- regional summaries computed before complete-case selection
- fixed posterior-median HIV imputation

The verifier's input was cut off after L209, so items from L210 onward were re-checked directly (evidence cited).

### L197–198: data sources (outdated)
- **Problem:** MICS appears nowhere in `main.tex` or `ref.bib`, and no MIS survey is in the analysed sample.
- **Evidence:** `survey_map/survey_timeline_all.csv`; `study_flow/CAPTION.md`.
```latex
Complete birth histories were obtained from the Births Recodes of the Demographic and Health Surveys and Malaria Indicator Surveys (DHS/MIS) \cite{dhs} and from UNICEF Multiple Indicator Cluster Surveys (MICS) that recorded a complete birth history \cite{mics}.
```
(`ref.bib` needs a `mics` entry; `docs/MANUSCRIPT_UPDATE_DHS_MICS.md` has a draft to check.)

### L201: DHS-only covariate sources; no MICS methods (outdated)
- **Problem:** The paragraph describes only DHS sources, and "comprised from the Births Recodes" is ungrammatical. There is no paragraph on MICS.
- **Evidence:**
  - `results/mics_inventory/ELIGIBILITY.md`: 93 MICS assessed, 47 with a complete birth history
  - `data/derived_mics/region_map.csv`: 366 labels, of which 291 match exactly and 75 use overrides
  - `R_mics/09_regional_covariates.R`: `card_only` lists MC_GIN2016, MC_COM2022 and MC_TCD2019
  - `study_flow/excluded_surveys.csv`: 4 Lesotho DHS surveys + MC_LSO2018
```latex
For DHS/MIS surveys, child, maternal and household characteristics comprised, from the Births Recodes, maternal age at first live birth, ... [rest unchanged].

Of the 93 MICS surveys from sub-Saharan Africa that we assessed, 47 (rounds 3 to 6, fieldwork 2006--2022, 28 countries) included a complete birth history; the other 46 recorded only summary birth histories and cannot support the age-band design. Each MICS birth history was converted to the DHS Births Recode layout and processed with the same age-band construction and eligibility rules as the DHS/MIS data. The 366 MICS region labels were mapped to analysis regions, defined as unions of DHS boundary polygons for the same country, by exact name (291) or a reviewed correspondence table (75); Guinea-Bissau and the Central African Republic have no DHS survey, so their regions were built from admin-1 polygons [REF]. Annual regional \textit{Pf}PR$_{2-10}$ was extracted for these regions by the same method as for DHS regions. MICS does not publish survey-region indicator tables, so all 13 regional covariates were computed from the MICS microdata with the DHS definitions, using the MICS women's, household and children's sampling weights. Education, recorded as level and grade, was converted to completed years using each country's primary and secondary school durations; facility delivery refers to the last birth in the two years before interview. Recall vaccination doses were unusable in Guinea 2016, Comoros 2022 and Chad 2019, so the national WUENIC DTP3 and measles estimates for the survey year were used for all regions of these surveys. Lesotho (four DHS surveys and one MICS survey) was excluded because it is malaria free and has no MAP prevalence surface.
```

### L204: regional weighting not described; en dash (clarity)
- **Problem:** The regional extracts use GPW density weights *without* cell area, while national values use density × area (L206). The text does not say how the regional values were weighted. "2000 - 2024" should use an en dash.
- **Evidence:** `R_mics/06_polygons_pfpr_registry.R` lines 4–6; `R_cbh/burden/README.md`.
```latex
\emph{P.\ falciparum} prevalence estimates were taken from the Malaria Atlas Project (MAP) annual modelled surfaces of the parasite rate standardised to ages 2 to 10 years (\textit{Pf}PR$_{2-10}$) for 2000--2024 \cite{map,map2025}. Survey-region values were means of each annual surface over the region's boundary polygons, weighted by the Gridded Population of the World version 4 (GPWv4) 2020 population density (without grid-cell area) \cite{gpwv4}.
```

### L204: MAP release not cited (citation)
- **Problem:** `map` and `map2025` cover 2000–17 and 2000–22. The rasters come from MAP dataset `Malaria__202508_Global_Pf_Parasite_Rate`.
- **Evidence:** `R_dhs/00_config.R` line 70.
```latex
... for 2000--2024 (MAP release \texttt{Malaria\_\_202508\_Global\_Pf\_Parasite\_Rate}, accessed via the malariaAtlas R package) \cite{map,map2025,map_release}.
```

### L205: ungrammatical; "estimated start" is vague (clarity)
```latex
Each age band in a child's trajectory (from birth until interview or death) was assigned the \textit{Pf}PR$_{2-10}$ estimate for the mother's survey region of residence at interview and the calendar year in which the child entered that band. To be eligible for the primary analysis, entry into the band had to fall within the 60 months before the interview.
```

### L206: Nigerian state exposure extracted differently (clarity)
- **Evidence:** `nigeria_states/README.md`: "MAP admin-1 boundaries, GPW 2020 density weights without cell area".
```latex
For national burden calculations, country-level prevalence was calculated as population-weighted means, with weights derived from the Gridded Population of the World version 4 population density surface for 2020 multiplied by grid-cell area \cite{gpwv4}. State-level prevalence for Nigeria was calculated over MAP admin-1 boundaries with the same 2020 density weights (without grid-cell area).
```

### L209: wrong citation (wrong)
- **Problem:** `igme` is the UN IGME 2024 child mortality report, but the data are IHME GBD exports.
- **Evidence:** `ref.bib`; `R_cbh/burden/README.md`: IHME all-cause export dated 2026-09-09.
```latex
National age-specific all-cause death counts and rates were obtained from the IHME Global Burden of Disease data exports and used to calculate malaria-attributable mortality \cite{ihme}.
```

### L210: "only country for which subnational estimates are available" (wrong)
- **Problem:** GBD publishes subnational estimates for several other countries, including Ethiopia, Kenya and South Africa. Nigeria is simply the only country analysed subnationally here. The IHME state malaria rates used as the comparator are also not mentioned.
- **Evidence:** `data/external/nigeria_states/source.json` lists the all-cause export and "GBD NG Malaria mortality incidence prevalence rates 2011-2024.csv", with the note "the 2024 death rate is the cause-specific comparator".
```latex
In addition, we used IHME estimates of all-cause deaths and death rates by age group, and of under-five malaria death rates, for the 36 states and the Federal Capital Territory of Nigeria, the only country we analysed subnationally \cite{ihme}.
```

### L212–213: sources uncited (citation)
- **Problem:** The UNAIDS/UNICEF HIV workbook, World Bank WDI and WGI series are uncited, and the political-stability indicator is not named.
- **Evidence:** `results/cbh/hiv_incidence/REPORT.md`; `R_cbh/00_config.R` (`wb_gdp_pc.csv`, `wb_hexp_pc.csv`, `wb_polstab.csv`).
```latex
... were obtained from the UNAIDS 2025 estimates distributed by UNICEF \cite{unaids2025}. ... National GDP per capita and current health expenditure per capita were obtained from the World Bank World Development Indicators \cite{wdi}, and political stability from the Worldwide Governance Indicators estimate of political stability and absence of violence \cite{wgi}.
```
(New bib entries needed.)

### L217: "LLIN coverage (measured in DHS/MIS surveys)" (outdated)
```latex
We did not adjust for LLIN coverage (measured in DHS, MIS and MICS surveys) because ...
```

### L217: "with and without measured of malnutrition" (typo)
- **Evidence:** `nutrition_adjustment_dhsmics_map_gamma2_v2/`, a refit on the same sample.
```latex
To assess the potential impact of this, we carried out a sensitivity analysis omitting both measures of malnutrition (wasting and stunting prevalence), fitted to the same analysis sample.
```

### L219: "regional summaries describe conditions at the time of the survey" (partly wrong)
- **Problem:** Facility delivery, birth interval and the maternal summaries refer to retrospective windows, not to the time of the survey.
- **Evidence:** `docs/ANALYSIS_PLAN.md` line 90: "describe survey-time conditions or the indicator's stated retrospective window"; lines 73–77 (five-year and 60-month windows).
```latex
We note that regional summaries describe conditions at the time of the survey, or over each indicator's own recall window (for example, births in the five years before the survey for facility delivery and birth intervals), rather than values at the time each child entered an age band.
```

### L222: weights are DHS-only (outdated)
- **Evidence:** `R_mics/09_regional_covariates.R` lines 103, 146–149, 161–163 (`wmweight`, `hhweight`/`hhweightMICS`, `chweight`/`chweightMICS`).
```latex
Recode-derived DHS regional summaries used DHS women's sampling weights. Published regional coverage indicators retained their stated denominators and reference periods. MICS regional summaries were computed from the microdata with the MICS women's, household and children's sampling weights.
```

### L225: vaccination fallback omits the MICS override (outdated)
- **Problem:** Three MICS surveys use national WUENIC values for every region because of data quality, not missingness. This sentence also repeats L201, where the source is called "WHO/UNICEF (WUENIC)" rather than "UNICEF".
- **Evidence:** `R_mics/09_regional_covariates.R` line 16; `docs/ANALYSIS_PLAN.md` line 92.
```latex
Missing regional DTP3 and measles coverage was replaced by the corresponding WHO/UNICEF national estimate (WUENIC) for the survey year (\url{https://immunizationdata.who.int/}); the same national values were used for all regions of three MICS surveys (Guinea 2016, Comoros 2022 and Chad 2019) whose recall vaccination doses were unusable.
```

### L225: "Liberia and Sao Tome had no data and were excluded" (wrong)
- **Problem:** Both countries have survey data: 4 DHS/MIS surveys in Liberia; 1 DHS and 2 MICS surveys in São Tomé. They were excluded because UNAIDS publishes no child or adolescent HIV incidence series for them.
- **Evidence:**
  - `study_flow/survey_selection.csv`: pfpr_available_rows LB51FL 31,692; LB5AFL 22,821; LB6AFL 44,084; LB7AFL 32,887; ST51FL 10,839; MC_STP2014 11,856; MC_STP2019 10,970 (total 165,149); complete_case_rows 0 for all
  - `R_cbh/hiv/README.md` line 49
  - `study_flow/flow_counts.csv`: missing_covariates 1,299,504 of with_map 8,797,963
```latex
Records with remaining missing required covariates were excluded: 1,299,504 of 8,797,963 MAP-eligible records (14.8\%), including all records from 34 surveys (25 DHS/MIS, among them all five MIS surveys, and 9 MICS). Liberia and S\~{a}o Tom\'{e} and Pr\'{i}ncipe have no UNAIDS child or adolescent HIV incidence series, so their surveys (four DHS/MIS surveys in Liberia; one DHS and two MICS surveys in S\~{a}o Tom\'{e} and Pr\'{i}ncipe; 165,149 child age-band records with MAP prevalence) were excluded from the primary analysis.
```

### L228–229: causal diagram shows "HIV prevalence" (internal inconsistency)
- **Problem:** The model adjusts for child HIV incidence (L212, L219). The diagram has no malnutrition node, although L217 discusses malnutrition.
- **Evidence:** Text extracted from `figures/Causal_diagram.pptx` includes "HIV prevalence" and "Malaria inoculation rate"; `model_formula.txt` uses `z_log_hiv_incidence`.
- **Fix:** Relabel the node "Child HIV incidence" and optionally add a "Malnutrition (wasting, stunting)" node, or add a sentence to the caption (see the next item).

### L229: "the malaria inoculation rate" (clarity)
- L68 uses "entomological inoculation rate". Relabel the figure node to match.
```latex
\caption{Assumed causal structure relating malaria transmission intensity (the entomological inoculation rate, proxied by \textit{Pf}PR$_{2-10}$) to all-cause child mortality, used to guide covariate selection. Solid arrows denote direct causes and dashed arrows effect modifiers. The HIV node is measured by national child (0--14 years) HIV incidence. LLIN: long-lasting insecticide-treated bednets; SES: socio-economic status.}
```

---

## 6. Methods continued and supplement (lines 234–333; checked directly)

About 30 statements were checked and found correct:
- **L236:** unit of observation; inclusion rules; completion rule for deaths and survivors; 7 records for a survivor, 3 for a death at 8 months; binomial equivalence (calendar year is fractional, by entry month: `R_cbh/R/child_bands.R`)
- **L238:** seven separate models; cloglog link; `offset(log(band_years))`; cr splines with k = 5 and k = 6; 17 linear covariates; survey, country and region random intercepts; `bam` fREML with `discrete=TRUE`; γ = 2; mgcv 1.9 (1.9-4 in `session_info.txt`); HR formula
- **L242:** AF formula; mapping of the IHME neonatal groups and the 2–4 year group; application to all 42 countries; point estimates only; shared denominator; Nigeria (36 states + FCT) and the reconciliation
- **L247:** Sahel and M49 definitions; knots at subset quantiles; PMM with 10 imputations; 2001 interpolation; Zimbabwe 2000–2009; Liberia and São Tomé via the latent adolescent series
- **L257:** the GitHub URL matches the repository remote
- **Figures S2–S5:** all match the current pipeline PNGs

### L238: SText1 is the wrong model (internal inconsistency)
- **Current:** "Full model equations are given in the Supplement (SText~1)."
- **Problem:** `SText1_model_specification.tex` describes the archived "Method 2" model:
  - a negative-binomial GAM of region-level post-neonatal death counts with a births offset
  - country random slopes; M1–M4 prevalence forms
  - 600 survey regions in 23 countries

  None of this is the seven-band cloglog model. SText1 needs to be rewritten; `docs/MANUSCRIPT_UPDATE_DHS_MICS.md` has a draft Overview. Its formula matches `model_formula.txt`, but check it before use.
- **Evidence:** SText1 lines 29–54 and 74–80; `model_formula.txt`.
```latex
% SText1 Overview (replace lines 29-54; delete or relabel the archived Method-2 sections)
The PfPR-ACM model relates annual survey-region \pfpr\ at entry into each age band to all-cause mortality in seven child age bands ($<$1, 1--5, 6--11, 12--23, 24--35, 36--47 and 48--59 completed months), using complete birth histories from Demographic and Health Surveys and UNICEF Multiple Indicator Cluster Surveys. The primary sample contains 7,498,459 child--age-band records from 2,292,089 children, with 102,282 deaths, in 1,211 survey-regions, 132 surveys (95 DHS and 37 MICS) and 36 countries. For child $i$ in band $g$ of width $\Delta_g$ years, $\mathrm{cloglog}\{\Pr(D_{ig}=1)\}=\log\Delta_g+\alpha_g+f_g(P_{r(i),y_{ig}})+h_g(t_{ig})+\mathbf{X}^{\top}\boldsymbol\beta_g+u_{s(i),g}+v_{c(i),g}+b_{r(i),g}$, where $f_g$ and $h_g$ are penalised cubic regression splines (basis dimensions 5 and 6), $\mathbf X$ holds the 17 standardised covariates, and $u$, $v$, $b$ are Gaussian random intercepts for survey, country and survey-region, fitted separately by band with \texttt{mgcv::bam} (fast REML, discretised covariates, $\gamma=2$).
```

### L238: what the intervals condition on (clarity)
- **Problem:** The pipeline's intervals also condition on exposure and on the filled regional covariates.
- **Evidence:** footnote to `tables/age_band_results.latex.txt`.
```latex
Confidence intervals used the coefficient covariance conditional on the estimated smoothing parameters, the exposure values, the single HIV imputation and the filled regional covariates.
```

### L242: "Estimates were produced annually for 2004-2024" (internal inconsistency)
- **Problem:** The abstract, Table 1 and Figure 3C all start in 2000.
- **Evidence:** `annual_comparison/annual_totals_2000_2024.csv` has rows from 2000 (42 countries; 1,165,686 model deaths).
```latex
Estimates were produced annually for 2000--2024, in the 42 countries with national MAP estimates.
```

### L245–246: heading and description of the sensitivity analyses (clarity)
- **Problem:**
  - `\subsection` is numbered, while every other heading uses `\subsection*`.
  - "Paramerisation" is a typo.
  - "Seasonal versus non seasonal, West Africa vs East Africa" does not describe what was fitted. The subsets are the Sahel and Eastern Africa; they overlap and are not complementary.
  - No parameterisation sensitivity is reported in the paper.
  - The no-nutrition refit is not listed.
- **Evidence:** `subgroups_dhsmics_map_gamma2_v2/subgroup_definitions.csv`; `nutrition_adjustment_dhsmics_map_gamma2_v2/`.
```latex
\subsection*{Sensitivity analyses}
We carried out sensitivity analyses to assess whether restricting the data to subsets of survey regions or surveys, imputing rather than excluding missing covariates, or removing the nutritional covariates changed the estimated relationship between prevalence and age-specific all-cause mortality.
```

### L247: subgroup sizes (outdated)
- **Current:** "23 surveys, 7 countries, 14,545 deaths ... 38 surveys, 12 countries, 23,583 deaths ... median survey year of 2013 (50 and 45 surveys)"
- **Evidence:** `subgroups_dhsmics_map_gamma2_v2/sample_summary.csv`:

| Subset | Surveys | Countries | Deaths |
|---|---:|---:|---:|
| Sahel | 31 | 8 | 18,919 |
| Eastern Africa | 51 | 12 | 31,589 |
| Up to median year | 70 | | 55,805 |
| After median year | 62 | | 46,477 |

  `subgroup_definitions.csv`: median_survey_year 2014.
```latex
We refitted the 7 age-band models in four subsets of the primary sample: Sahelian survey regions (boundary centroid at or north of 12°N and west of 36°E, excluding the Horn of Africa; 31 surveys, 8 countries, 18,919 deaths), Eastern Africa (UN M49 delimitation of Eastern Africa: 51 surveys, 12 countries, 31,589 deaths), and surveys conducted up to and including versus after the median survey year of 2014 (70 surveys and 55,805 deaths; 62 surveys and 46,477 deaths). Subset fits used the same specification and covariate scaling as the primary analysis, with knots at quantiles of the subset's own predictor values.
```

### L247: imputed-covariate sensitivity (outdated)
- **Current:** "This retained all 6,357,802 eligible records (120 surveys, 36 countries, 90,938 deaths)"
- **Problem:**
  - The sample sizes are DHS-only.
  - The list of imputed components leaves out Somalia and South Sudan.
  - The Rubin's-rules refits used fixed smoothing parameters, and the curves in Figure S4 come from the point imputation (the mean of the 10 imputations).
- **Evidence:** `primary_map_regional17_dhsmics_imputed_gamma2_v6/REPORT.md`: 8,797,963 records; 2,656,472 children; 123,419 deaths; 1,457 survey regions; 166 surveys; 40 countries. `multiple_imputation/REPORT.md`: smoothing parameters fixed at the point fit.
```latex
In a second sensitivity analysis we imputed every remaining covariate gap rather than excluding records: whole-survey gaps in regional indicators by chained-equation multiple imputation across the combined DHS and MICS survey regions (predictive mean matching, 10 imputations; the point analysis used their mean); the missing 2001 political-stability round by interpolation; health expenditure for Zimbabwe 2000--2009, Somalia 2000--2012, South Sudan before 2017 and all countries in 2024, and South Sudan's GDP before 2008 and political stability before independence, from generalised additive models on the observed national panel; and child HIV incidence for Liberia and S\~{a}o Tom\'{e} and Pr\'{i}ncipe from the incidence model with a latent adolescent series. This retained all 8,797,963 MAP-eligible records (2,656,472 children, 123,419 deaths, 1,457 survey regions, 166 surveys in 40 countries). As a check on imputation uncertainty, we refitted the models in each of the ten imputed datasets, with smoothing parameters fixed at those of the point fit, and pooled the results with Rubin's rules.
```

### L251: Contributors placeholder (placeholder)
- "XXXXX" is unfilled.

### L257: data sharing omits most sources (outdated)
```latex
All input data are available from their sources: DHS and MICS microdata on registration from the DHS Program (\url{https://dhsprogram.com}) and UNICEF (\url{https://mics.unicef.org}); IHME estimates from the GBD results tool; UN IGME (CA-CODE) cause-specific estimates from the UN IGME portal (\url{https://childmortality.org}); MAP prevalence surfaces from the Malaria Atlas Project; and UNAIDS, WHO/UNICEF (WUENIC) and World Bank series from their publishers. Analysis code and derived data tables are available at \url{https://github.com/Infectious-Diseases-Data-Observatory/burden-malaria-prevalence}.
```

### L276: Figure S1 caption (outdated)
- **Problem:** The caption describes child HIV *prevalence* (UNAIDS case counts ÷ World Bank population), but the model uses child HIV *incidence*. It also gives no counts and does not mention MICS. The image is current (MD5 77cb1109… = `study_flow/study_flow_diagram.png`).
- **Evidence:** `study_flow/CAPTION.md`; `study_flow/flow_counts.csv`; `model_formula.txt`.
```latex
\caption{Data sources and sample inclusion. The survey registry contains 124 DHS/MIS and 93 MICS surveys in 43 countries; 46 MICS surveys without a complete birth history and 5 surveys in Lesotho (malaria free; no MAP prevalence estimates) were excluded, leaving 166 surveys in 40 countries. Of 10,657,524 child age-band entries within 60 months of interview, 1,347,434 had a band not complete by interview and 418,661 an entry year outside 2000--2024, leaving 8,891,429 eligible records; 93,466 lacked a regional MAP prevalence value and 1,299,504 lacked at least one required covariate, leaving 7,498,459 records (2,292,089 children; 102,282 deaths) from 132 surveys (95 DHS and 37 MICS) in 36 countries. External (non-survey) data comprise MAP \textit{Pf}PR$_{2-10}$ (population-weighted with a gridded population surface), WHO/UNICEF national vaccine coverage, UNAIDS child HIV incidence (ages 0--14, per 1,000 uninfected) and World Bank GDP, health expenditure and political stability series.}
```

### L282: Figure S2 caption (typo / clarity)
- **Problem:**
  - "Log hazard ration relative for" is a typo.
  - The curves are point estimates, not "mean values".
  - The curves span the central 95% of exposure, and the caption does not say so.
- **Evidence:** image MD5 eea85640… = `pfpr_splines.png`.
```latex
\caption{Malaria prevalence and all-cause age-band-specific child mortality in the primary sample (7,498,459 child age-band records, 102,282 deaths, 132 DHS and MICS surveys). The y-axis shows the log hazard ratio of age-specific mortality relative to a \textit{Pf}PR$_{2-10}$ of 20\%. Point estimates (pointwise 95\% conditional intervals) are shown by the blue lines (shaded areas), over the central 95\% of each age band's exposure distribution.}
```

### L290: Figure S3 caption "Sensitivity analyses: subgroups" (clarity)
- **Evidence:** image MD5 = `subgroups_dhsmics_map_gamma2_v2/sfig_pfpr_splines_by_subgroup.png`; `CAPTION.md`.
```latex
\caption{Sensitivity of the \textit{Pf}PR$_{2-10}$-mortality relationship to the analysis subset. Log mortality hazard ratios relative to \textit{Pf}PR$_{2-10}$ = 20\% by age band from the full primary sample (black, with its pointwise 95\% interval shaded) and from separate refits in Sahelian survey regions (31 surveys, 8 countries, 18,919 deaths), Eastern Africa (51 surveys, 12 countries, 31,589 deaths), and surveys conducted up to (70 surveys, 55,805 deaths) and after (62 surveys, 46,477 deaths) the median survey year of 2014. Each curve is drawn over the central 95\% of its own exposure distribution. The subsets overlap the full sample and one another, so their differences are descriptive.}
```

### L297: Figure S4 caption "Results from multiply imputed datasets" (clarity)
- **Problem:** The figure compares the complete-case fit with the point-imputed fit. Multiple imputation with Rubin's rules was only a check.
- **Evidence:** image MD5 bc567c26… = v6 `sfig_pfpr_splines_imputed_covariates.png`; v6 `REPORT.md`.
```latex
\caption{Sensitivity of the \textit{Pf}PR$_{2-10}$-mortality relationship to the treatment of missing covariates: complete-case primary sample (7,498,459 records, 102,282 deaths, 132 surveys in 36 countries) versus all MAP-eligible records after imputing every remaining covariate gap (8,797,963 records, 123,419 deaths, 166 surveys in 40 countries). Log mortality hazard ratios relative to \textit{Pf}PR$_{2-10}$ = 20\%, with pointwise 95\% conditional intervals, over the central 95\% of each curve's exposure distribution. Hazard ratios for \textit{Pf}PR$_{2-10}$ 40\% to 20\% changed by at most 0.04; refitting in each of 10 imputed datasets and pooling with Rubin's rules changed them by at most 0.002.}
```

### L304: Figure S5 caption (clarity)
- **Evidence:** image MD5 = `nutrition_adjustment_dhsmics_map_gamma2_v2/sfig_pfpr_splines_without_nutrition.png`; `CONTRASTS.md` gives a largest change of 0.008, at 1–5 months (0.92 → 0.91).
```latex
\caption{Sensitivity of the \textit{Pf}PR$_{2-10}$-mortality relationship to removing regional wasting and stunting prevalence from the adjustment set. Log mortality hazard ratios relative to \textit{Pf}PR$_{2-10}$ = 20\% by age band from the primary model (17 covariates) and from a refit on the identical sample (7,498,459 records, 102,282 deaths) with 15 covariates. The largest change in the hazard ratio for a reduction in \textit{Pf}PR$_{2-10}$ from 40\% to 20\% is 0.008.}
```

### L310–327: Table S1 (outdated)
- **Problem:** Every data row and the total are the DHS-only values (966,553 records, 27,788 deaths and so on).
- **Evidence:** `tables/age_band_results.latex.txt`. You can also paste that file whole; it adds a footnote on rounding and the conditional intervals.
```latex
% caption (L310): ... based on the primary analysis dataset (132 DHS and MICS surveys, 36 countries). Each child can contribute to multiple age band records.
$<1$ & 1,322,149 & 37,623 (36.8\%) & 2.34 & 0.94 (0.92--0.97) & 0.95 (0.90--1.00) \\
1--5 & 1,189,209 & 15,586 (15.2\%) & 2.25 & 0.92 (0.88--0.95) & 0.86 (0.81--0.93) \\
6--11 & 1,143,915 & 14,126 (13.8\%) & 2.90 & 0.83 (0.79--0.87) & 0.70 (0.63--0.77) \\
12--23 & 989,211 & 13,138 (12.9\%) & 3.46 & 0.82 (0.77--0.88) & 0.54 (0.48--0.62) \\
24--35 & 980,344 & 11,338 (11.1\%) & 3.63 & 0.85 (0.80--0.91) & 0.47 (0.41--0.54) \\
36--47 & 949,032 & 6,650 (6.5\%) & 3.34 & 0.84 (0.78--0.91) & 0.54 (0.47--0.63) \\
48--59 & 924,599 & 3,821 (3.7\%) & 3.16 & 0.91 (0.83--0.99) & 0.58 (0.49--0.68) \\
\midrule
Total & 7,498,459 & 102,282 (100.0\%) & -- & -- & -- \\
```
Table S1 shows 12.9% at 12–23 months (largest-remainder rounding of 12.84%), while the Results text says 13%. Both are consistent with the file.

**Optional:** `Supplementary Figures/sfig_burden_comparison_counts.png` (the death-count version of Figure 3) is exported but not referenced in `main.tex`.

---

## Unverifiable from the pipeline (external literature; not counted as errors)

Check each of these against its cited source:
- **L65:** "around half a million preventable deaths in children" (`wmr`). The pipeline's comparators give about 430,000–440,000 under-5 deaths in 2024 (WMR African Region all-age 579,000, of which about 75% are under 5). A safer wording is "more than 400,000 deaths in children under 5".
- **L50, L65:** claims about verbal autopsy methods: single cause per death, poor sensitivity and specificity, few surveillance sites, the era of failing treatment (`todd1994`, `snowva1992`, `white_misdiagnosis`, `scott`, `gates2017ihme`).
- **L68:** EIR magnitudes, acquisition of immunity, severe anaemia admissions, and the `snowsevere` relationship.
- **L85:** passive immunity in neonates and causes of neonatal death (`WHO:14`, `who_newborn_mortality_2024`). There is also a minor redundancy ("Most ... are predominantly caused").
- **L174:** "the start of the Roll Back Malaria Partnership" in 2000. RBM is usually dated to its 1998 launch; 2000 is the year of the Abuja Declaration.
- **L188:** "In lower transmission areas, the mortality burden is less concentrated in early childhood" has no citation. The pipeline covers only ages 0–59 months.

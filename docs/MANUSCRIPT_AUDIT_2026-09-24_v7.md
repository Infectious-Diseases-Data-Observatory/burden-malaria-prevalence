# Manuscript audit against primary v7 (DHS + MICS + Liberia), 24 September 2026

- **Manuscript audited:** `~/Dropbox/Apps/Overleaf/malaria-burden-reassessment/main.tex`. This is the current file, 328 lines, last modified 24 Sep 12:17. The audit covers the abstract (L47–52), Results (L74–166) and Methods (L192–242). Supplementary items cited from those sections are covered briefly at the end.
- **No TeX was edited.** Nothing in the Overleaf folder or in `ref.bib` was changed. Every replacement below is text to paste by hand.
- **Source of truth:** the current primary, `results/cbh/primary_map_regional17_dhsmics_gamma2_v7/`. It combines 98 DHS and 37 MICS surveys, plus Liberia's three DHS surveys. Liberia's child HIV incidence is derived from UNAIDS counts (`results/cbh/hiv_incidence/liberia_aidsinfo/REPORT.md`).
- **This report supersedes `docs/MANUSCRIPT_AUDIT_2026-09-24.md`,** which audited against v5.
- **Verification:** each item was proposed by a section auditor and re-checked by a verifier. Rejected items have been dropped. The synthesis step then spot-checked the key numbers against the v7 CSVs and the figure checksums.
- **Three items were checked by the synthesis step alone,** because the verifier's output for the Methods-model section was cut off: L231 grammar, L237 age mapping and reconciliation, and L240 numbering. Two further items are new from this step: L233 on the interval conditioning, and Supplementary Table S1.
- **Dash style:** suggested replacements use `--` for numeric ranges, as the generated table fragments do. The body text currently uses hyphens for age ranges ("1-5 months"). Pick one style and apply it throughout.

## Summary

### Status of the v7 sensitivity refits (checked 24 Sep 13:23 BST)

| Analysis | Output | Status | What the text can rely on now |
|---|---|---|---|
| Subgroups (Sahel, Eastern Africa, early, late) | `results/cbh/subgroups_dhsmics_map_gamma2_v3/` | **Complete.** All 28 fits converged; REPORT.md was written at 13:17. | The v3 numbers quoted below are final. |
| No wasting/stunting | `results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v3/` | **Complete.** All 7 fits converged. | Largest change in HR 40%→20%: 0.008 (1–5 months: 0.916 → 0.908); in HR 20%→0%: 0.031 (1–5 months: 0.862 → 0.831). |
| Imputed covariates plus multiple imputation | `results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v8/` | **Complete.** All 7 point fits and 70 multiple-imputation refits done. | Sample unchanged from v6: 8,797,963 records, 123,419 deaths, 166 surveys, 40 countries. Largest change against v7 in HR 40%→20%: 0.028 (48–59 months: 0.919 → 0.948); in HR 20%→0%: 0.018 (12–23 months). Rubin pooling changed HRs by at most 0.002; between-imputation share of variance at most 2.6%. |

### Counts

About 300 statements and table cells were checked. 198 were confirmed correct. The flagged items below often cover several numbers each; a whole table block counts as one item.

| Section | Confirmed correct | Flagged items | Outdated | Wrong statement / number | Citation | Internal inconsistency | Clarity / typo / omission | Pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Abstract (L47–52) | 22 | 10 | 5 | 1 | 0 | 0 | 4 | 0 |
| Results: sample and curves (L74–97) | 39 | 16 | 6 | 4 | 0 | 1 | 5 | 0 |
| Results: IHME/UN IGME comparison, Tables 1–2 (L99–166) | 56 | 14 | 7 | 2 | 0 | 1 | 4 | 0 |
| Methods: data and covariates (L192–225) | 35 | 16 | 6 | 2 | 3 | 1 | 4 | 0 |
| Methods: model, burden, sensitivity (L229–242) | 46 | 16 | 5 | 3 | 0 | 1 | 7 | 0 |
| **Total in the assigned range** | **198** | **72** | **29** | **12** | **3** | **4** | **24** | **0** |
| Supplementary material cited from these sections | – | 5 | | | | | | |
| Unverifiable literature claims (listed at the end) | – | 4 | | | | | | |

### The 10 most important corrections

1. **Table 1 model column (L117–121, L124–128).** The death counts and rates come from the old DHS-only fit (v3). v7 gives 1,170,371 / 725,613 / 620,072 deaths and 10.66 / 4.58 / 3.49 per 1,000 child-years. The current values contradict both Figure 3C and the abstract.
2. **Table 2 (L151–162).** The model column and the top-ten list are v3. In v7, DR Congo and Côte d'Ivoire enter the list and Central African Republic and Uganda leave it. The all-country total is 620,072.
3. **L102, Nigeria and Chad.** Nigeria is 218 thousand deaths (58% above IHME, 32% above UN IGME), not 200 thousand (46% and 22%). Chad is 4.0-fold IHME, not 3.6-fold.
4. **L102, "mostly agreed".** The model rate exceeds IHME's in 36 of 42 countries, and the median ratio is 1.58.
5. **Sample counts (L51, L76).** 132 surveys (95 DHS) in 36 countries, 2,292,089 children and 102,282 deaths should read 135 surveys (98 DHS, 37 MICS) in 37 countries, 2,325,130 children and 103,987 deaths. DHS covariate exclusions are 17, not 20.
6. **Abstract headline numbers (L51–52).** Attributable fractions: 6%, 31% and 43–53%, not 5%, 30% and 42–53%. Rates: 10.7 → 4.6 → 3.5, not 10.6 → 4.5 → 3.4. The 2024 gap is 45% above IHME and 41% above UN IGME, not "about 40%".
7. **Liberia and São Tomé (L220, L207–208).** "Liberia and Sao Tome had no data and were excluded" is wrong. Liberia's three DHS surveys are now in the sample, using child incidence derived from UNAIDS counts, and the Methods must describe that derivation. São Tomé has surveys but no HIV incidence series.
8. **Curve shape (L88, L86).** At 40% prevalence the attributable fraction peaks at 59% (24–35 months), not 56%. The "plateau and even decreased" description of the curves in older children does not match the v7 curves.
9. **Sensitivity analyses (L242, L97, L241).** The subset and imputed-sample sizes are DHS-only values (23/7/14,545; 38/23,583; median 2013, 50/45; 6,357,802). The Results description of the subgroup findings is wrong: the Sahel is steeper from 6 months, not from 1 year, and Eastern Africa is flatter in neonates and at 36–47 months. The v3 subgroup numbers are now final.
10. **Stale figures and the Figure 3 caption (L79–93, L137–141, S1–S3).** Figures 1–3 in Overleaf are v5 exports (checksums match v5), and S1–S3 are v5 or v2 exports. The Figure 3 caption says panel C runs "since 2004" and names only IHME; the Methods (L237) also say 2004. Estimates run from 2000.

Other substantive items: the age-band inclusion rule at L76 is misstated; L204 cites `igme` for IHME data; SText1 (cited at L233) describes the archived negative-binomial model; Nigeria is 36 of 37 states, not 35 (L134); and Supplementary Table S1 (L303–325) still holds an old sample (5,465,305 records, 75,726 deaths).

---

## Abstract (L47–52)

A consolidated replacement for abstract paragraphs 2–3, with all v7 corrections applied, comes after the individual items.

### L47: word-count comment (clarity)
- **Current:** `% 210 words`
- **Problem:** The abstract body (L50–52) is about 333 words, and the replacements add about 10 more. That is over a 300-word Lancet-style limit. The file header (L4) promises a structured abstract, but there are no Background/Methods/Findings/Interpretation headings.
- **Evidence:** `sed -n '50,52p' main.tex | wc -w` = 333; main.tex L4 says "structured abstract, IMRaD".
```latex
% ~335 words -- trim to <=300 (Lancet) and consider Background/Methods/Findings/Interpretation headings
```

### L50, L52: style (clarity)
- **Current:** "all cause childhood mortality" (L50); "3.4 deaths per 1000" (L52)
- **Problem:** L50 writes "all cause" where L51 writes "all-cause". L52 mixes "per 1,000 child-years" and "per 1000". Hyphenated ranges ("6-11 months") match the body style at L76, so they can stay.
- **Evidence:** main.tex L50, L51, L52.
```latex
... and all-cause childhood mortality.
% and on L52 write "per 1,000 child-years" throughout
```

### L51: survey counts (outdated)
- **Current:** "Demographic and Health Surveys and Multiple Indicator Cluster Surveys (MICS) ... (132 surveys, 36 countries)"
- **Problem:** These are the v5 counts. v7 adds Liberia's 2007, 2013 and 2019–20 DHS surveys. The abstract also defines MICS but not DHS.
- **Evidence:** v7 `primary_sample.csv`: surveys 135, countries 37. `comparison_sample.csv`: 132→135 surveys, 36→37 countries. `survey_map/survey_timeline_all.csv`: 98 DHS and 37 MICS included; all 5 MIS excluded.
```latex
We pooled all available Demographic and Health Surveys (DHS) and Multiple Indicator Cluster Surveys (MICS) that recorded complete birth histories since 2000 and potential confounding variables (135 surveys: 98 DHS, 37 MICS; 37 countries)
```

### L51: "admin level-1" exposure (wrong statement)
- **Current:** "with the admin level-1 (subnational first level administrative boundaries) population weighted prevalence estimates from the Malaria Atlas Project"
- **Problem:** PfPR is assigned by survey region, and survey regions are not always admin-1 units. Nigeria uses its six geopolitical zones; Malawi 2013 is aggregated to 3 regions; several capitals are merged with a neighbouring unit. L76 already says "survey region (mostly first-level administrative units)".
- **Evidence:** `docs/ANALYSIS_PLAN.md` §2.2; v7 `prepared_sample_dhs_mics.csv`: 1,227 survey regions; main.tex L76.
```latex
with annual population-weighted \textit{Pf}PR$_{2-10}$ estimates from the Malaria Atlas Project for each survey region (mostly first-level administrative units).
```

### L51: children and deaths (outdated)
- **Current:** "in data from 2,292,089 children under 5 of whom 102,282 died"
- **Problem:** These are the v5 counts.
- **Evidence:** v7 `primary_sample.csv`: distinct_children 2,325,130, deaths 103,987. `study_flow/flow_counts.csv` gives the same.
```latex
a robust age-dependent association between \textit{Pf}PR$_{2-10}$ and all-cause mortality was observed in data from 2,325,130 children under 5 of whom 103,987 died.
```

### L51: attributable fractions at 20% (outdated)
- **Current:** "the model attributes 5\% of neonatal deaths, 30\% of deaths at 6-11 months and 42-53\% of deaths at 1-4 years"
- **Problem:** These are the v5 values (5.3%, 30.4%, 42.1–52.7%).
- **Evidence:** v7 `pfpr_20_to_zero_contrasts.csv` and `attributable_fraction_by_age.csv` at 20%: <1 month 0.0576; 6–11 months 0.3114; 12–23 0.4666; 24–35 0.5297; 36–47 0.4480; 48–59 0.4261. The 1–4 year range is therefore 42.6–53.0%.
```latex
At 20\% prevalence, the model attributes 6\% of neonatal deaths, 31\% of deaths at 6-11 months and 43-53\% of deaths at 1-4 years to malaria.
```

### L52: pooled rates (outdated)
- **Current:** "fell from 10.6 to 4.5 deaths per 1,000 child-years between 2000 and 2015 (a 57\% decline), and then to 3.4 deaths per 1000 in 2024 (24\% decline)"
- **Problem:** The rates are v5 (10.61, 4.52, 3.45); v7 gives 10.66, 4.58 and 3.49. The percentage declines are unchanged: 57.1% and 23.6%. The all-cause baseline (IHME) is not named.
- **Evidence:** v7 `burden_comparison/source_comparison_2000_2024.csv` (PfPR-ACM model, per 1,000 child-years): 10.656 (2000), 4.576 (2015), 3.494 (2024); decline_from_2015 23.64; 1 − 4.576/10.656 = 57.1%. The file covers 42 countries.
```latex
Applied to IHME estimates of national all-cause mortality in 42 sub-Saharan African countries, we estimate that under-5 malaria mortality fell from 10.7 to 4.6 deaths per 1,000 child-years between 2000 and 2015 (a 57\% decline), and then to 3.5 deaths per 1,000 child-years in 2024 (a further 24\% decline).
```

### L52: "In contrast" (clarity)
- **Current:** "In contrast, IHME and UN IGME estimate substantial declines up until 2015 (49 and 60\% decrease, respectively) but little change since (8 and 6\% decrease)."
- **Problem:** All four numbers are correct. But the model's 57% decline to 2015 lies between IHME's 49% and UN IGME's 60%, so the contrast holds only after 2015. All four figures are rate declines; the IHME and UN IGME death counts rose after 2015 (+2.9% and +5.6%).
- **Evidence:** v7 `source_comparison_2000_2024.csv`. IHME rates: 5.132 → 2.624 → 2.413 (declines 48.9% and 8.07%). UN IGME rates: 6.595 → 2.629 → 2.480 (declines 60.1% and 5.65%).
```latex
IHME and UN IGME estimate similarly large declines in the mortality rate up until 2015 (49\% and 60\% decrease, respectively) but, in contrast to our estimates, little change in the rate since (8\% and 6\% decrease).
```

### L52: "about 40\% higher" (outdated)
- **Current:** "We estimate malaria mortality in 2024 was about 40\% higher than either source in the same set of countries."
- **Problem:** v7 is 44.8% above IHME and 40.9% above UN IGME; v5 was 43% and 39%. The model estimand is a malaria-attributable reduction in all-cause deaths, not cause-specific deaths.
- **Evidence:** v7 `source_comparison_2000_2024.csv` 2024 deaths: model 620,072; IHME 428,147; UN IGME 440,123. Ratios 1.448 and 1.409. v7 REPORT.md: "predicted all-cause reductions under zero PfPR".
```latex
We estimate that malaria-attributable under-5 mortality in 2024 was 45\% higher than the IHME estimate and 41\% higher than the UN IGME estimate in the same set of countries.
```

### L52: closing sentence (clarity)
- **Current:** "Recent progress against childhood malaria has been substantially greater than currently reported, but there are more deaths caused by malaria than currently reported."
- **Problem:** "Than currently reported" appears twice. "Deaths caused by malaria" suggests cause-specific counting, whereas the pipeline captions call the model and cause-of-death estimates "different estimands". The substance is supported: the model's 2015–2024 rate decline is 23.6%, against 8.1% (IHME) and 5.7% (UN IGME).
- **Evidence:** `burden_comparison/source_comparison.latex.txt` note ("These are different estimands"); v7 REPORT.md L60.
```latex
Recent progress against childhood malaria has been substantially greater than currently reported, yet the malaria-attributable burden of child mortality remains larger than current cause-of-death estimates suggest.
```

### Consolidated replacement for abstract paragraphs 2–3 (L51–52)
```latex
    We pooled all available Demographic and Health Surveys (DHS) and Multiple Indicator Cluster Surveys (MICS) that recorded complete birth histories since 2000 and potential confounding variables (135 surveys: 98 DHS, 37 MICS; 37 countries) and matched individual child trajectories, discretised into 7 age bands from birth until 5 years, with annual population-weighted \textit{Pf}PR$_{2-10}$ estimates from the Malaria Atlas Project for each survey region (mostly first-level administrative units). After adjusting for a series of potential confounding variables, a robust age-dependent association between \textit{Pf}PR$_{2-10}$ and all-cause mortality was observed in data from 2,325,130 children under 5 of whom 103,987 died. \textit{Pf}PR$_{2-10}$ only weakly predicted neonatal death, but was a strong predictor of death from six months and older. At 20\% prevalence, the model attributes 6\% of neonatal deaths, 31\% of deaths at 6-11 months and 43-53\% of deaths at 1-4 years to malaria.
    Applied to IHME estimates of national all-cause mortality in 42 sub-Saharan African countries, we estimate that under-5 malaria mortality fell from 10.7 to 4.6 deaths per 1,000 child-years between 2000 and 2015 (a 57\% decline), and then to 3.5 deaths per 1,000 child-years in 2024 (a further 24\% decline). IHME and UN IGME estimate similarly large declines in the mortality rate up until 2015 (49\% and 60\% decrease, respectively) but, in contrast to our estimates, little change in the rate since (8\% and 6\% decrease). We estimate that malaria-attributable under-5 mortality in 2024 was 45\% higher than the IHME estimate and 41\% higher than the UN IGME estimate in the same set of countries. Recent progress against childhood malaria has been substantially greater than currently reported, yet the malaria-attributable burden of child mortality remains larger than current cause-of-death estimates suggest.
```

---

## Results: sample and curves (L74–97)

Confirmed correct: 171 surveys (119 DHS, 5 MIS, 47 MICS) in 41 countries; the five Lesotho surveys; death shares of 37% neonatal, 66% in the first year and 13% at 12–23 months (36.8%, 65.9%, 12.9%); the neonatal share at L86; the 10% attributable fractions at L88 (28–36% from 12 months; 2.8% in neonates); and monotone post-neonatal infant curves.

### L74: subsection heading (outdated)
- **Current:** `\subsection*{Prevalence and all cause mortality in DHS/MIS surveys}`
- **Problem:** The analysed sample is DHS plus MICS, and no MIS survey is included.
- **Evidence:** `survey_map/survey_timeline_all.csv`: 98 DHS and 37 MICS included; MIS 0 of 5.
```latex
\subsection*{Prevalence and all-cause mortality in DHS and MICS surveys}
```

### L76: age-band inclusion rule (wrong statement)
- **Current:** "the band was completed (by survival to its end or by death within it) before the interview"
- **Problem:** Under this wording, a death inside a band that would have ended after the interview counts as "completed". The pipeline excludes such deaths: the whole potential band must end by the interview date, for deaths and survivors alike. Band entry must also fall in 2000–2024. (Methods L231 states the rule correctly.)
- **Evidence:** `docs/ANALYSIS_PLAN.md` §2.3 rule 3; v7 `study_flow/flow_counts.csv`: incomplete 1,347,434 records, outside_year 418,661.
```latex
An age band was included if the child entered it alive within the 60 months before the interview, in 2000--2024, and the whole band would have ended by the date of interview, for deaths and survivors alike (deaths within bands that would otherwise have been incomplete were excluded).
```

### L76: included surveys (outdated)
- **Current:** "Restricting the analysis to complete cases with the identified confounders recorded, a total of 132 surveys (95 DHS, 37 MICS) in 36 countries were included"
- **Problem:** These are the v5 counts. The sentence also omits that selection happens after the declared substitutions described at L220.
- **Evidence:** v7 `study_flow/flow_counts.csv`: surveys 135, countries 37; `survey_timeline_all.csv`: 98 DHS and 37 MICS included.
```latex
Restricting the analysis to records with all 17 adjustment covariates available (after the substitutions described in the Methods), a total of 135 surveys (98 DHS, 37 MICS) in 37 countries were included in the primary analysis (Figure \ref{fig:flow});
```

### L76: covariate exclusions (outdated)
- **Current:** "all five MIS surveys, 20 DHS surveys and 9 MICS surveys were excluded because of missing data on identified confounders"
- **Problem:** Liberia's three DHS surveys now enter the sample, so 17 DHS surveys are excluded, not 20. Total covariate exclusions are 31.
- **Evidence:** `survey_timeline_all.csv`: excluded for covariates, DHS 17, MIS 5, MICS 9. `study_flow/CAPTION.md`: "31 processed surveys contribute no records".
```latex
all five MIS surveys, 17 DHS surveys and 9 MICS surveys were excluded because a required covariate was unavailable for every region of the survey.
```

### L76: children and deaths (outdated)
- **Current:** "The included surveys contributed data from a total of 2,292,089 children with 102,282 deaths recorded."
- **Problem:** These are the v5 counts.
- **Evidence:** v7 `primary_sample.csv`: 2,325,130 children, 7,607,122 records, 103,987 deaths; `prepared_sample_dhs_mics.csv`: 1,227 survey regions.
```latex
The included surveys contributed 7,607,122 child age-band records from 2,325,130 children in 1,227 survey regions, with 103,987 deaths recorded.
```

### L76, L86: wording (clarity)
- **Current:** "population-weighted \textit{Pf}PR$_{2-10}$ prevalence estimate" (L76); "and all cause mortality" (L86)
- **Problem:** "PfPR prevalence" is redundant. "All-cause" needs a hyphen.
- **Evidence:** Text only.
```latex
population-weighted \textit{Pf}PR$_{2-10}$ estimate ... and all-cause mortality
```

### L79–84: Figure 1 image and caption (outdated)
- **Current:** `figures/fig1_survey_map_and_timing.png`; the caption says twice that rows are grouped by region.
- **Problem:** The Overleaf image is the v5 export, in which Liberia's surveys are drawn open (excluded) and the mothers' counts are old. The caption repeats the grouping statement, gives no key for the map shading or the exclusion reasons, and has no total number of mothers.
- **Evidence:** The md5 of the Overleaf image is `efbbd5dc…`, identical to v5 `survey_map/survey_map_and_timing.png`; the v7 image is `4768438a…`. v7 `survey_timeline_all.csv` gives 171 surveys in 41 countries. `mothers_by_survey.csv` sums to 1,123,556. `country_summary.csv` gives LBR 16,892 mothers.
- **Action:** copy the v7 `survey_map/survey_map_and_timing.png` into `figures/`, then use:
```latex
\caption{DHS/MIS (circles) and MICS (triangles) surveys with complete birth histories conducted in sub-Saharan Africa since 2000 (124 DHS/MIS and 47 MICS surveys in 41 countries). Filled symbols are the 98 DHS and 37 MICS surveys included in the primary analysis (37 countries); open symbols are the 36 excluded surveys: five in Lesotho, which is not malaria endemic, and 31 (17 DHS, all 5 MIS and 9 MICS) for which a required covariate was unavailable for every survey region. Map shading gives the number of included surveys per country (grey: none). Timeline rows are labelled with ISO3 country codes and grouped by UN sub-region (Central Africa: UN Middle Africa); $n$ is the number of mothers (complete birth histories) contributing at least one child to the primary analysis, summed over the country's included surveys (1,123,556 in total).}
```

### L86: covariate description (clarity)
- **Current:** "adjusted for a set of 17 socio-economic and health-system covariates"
- **Problem:** The count of 17 is correct (13 regional, 4 national). The label is too narrow: the set also includes nutrition, governance and child HIV incidence. The calendar-year spline and the random effects are not mentioned.
- **Evidence:** v7 `covariate_scaling.csv` (17 rows); `model_formula.txt`.
```latex
adjusted for 17 survey-region and national covariates (socio-economic, health-system, nutritional, governance and child HIV incidence measures), a smooth function of calendar year, and survey, country and region random effects.
```

### L86: neonatal causes sentence (typo)
- **Current:** "Most neonatal deaths in these areas are predominantly caused by ... although malaria in primigravidae pregnancies increase risk of low birth \cite{Walker:14}"
- **Problem:** "Most … predominantly" is redundant; "increase" should be "increases"; "low birth" is missing "weight".
```latex
Most neonatal deaths in these areas are caused by birth-related complications (preterm birth, birth asphyxia or trauma), bacterial and viral infections and congenital anomalies \cite{who_newborn_mortality_2024}, although malaria in primigravid pregnancies increases the risk of low birth weight \cite{Walker:14}.
```

### L86: negative control (clarity; overstated)
- **Current:** "neonatal mortality provides an approximate negative control group and the weak association suggests no major residual confounding."
- **Problem:** The neonatal HR is small, but its interval excludes 1. The same paragraph cites a direct pathway (malaria in pregnancy leading to low birth weight), so neonates are not a clean negative control. In the Sahel subset the neonatal HR is 0.90 (0.85–0.96).
- **Evidence:** v7 `tables/age_band_results.csv` <1 month: HR 40→20 0.945 (0.917–0.974); HR 20→0 0.942 (0.895–0.992).
```latex
As such, neonatal mortality serves as an approximate negative control: the association was small (hazard ratio 0.94, 95\% interval 0.92--0.97, for a reduction in \textit{Pf}PR$_{2-10}$ from 40\% to 20\%), consistent with limited residual confounding, although malaria in pregnancy provides a plausible direct pathway.
```

### L86: plateau and "even decreased" (wrong statement)
- **Current:** "risk of death in older children plateaued above prevalences of around 30\% and even decreased for children over 3 years of age for prevalences above 30\% and 40\%"
- **Problem:** At 12–23 and 24–35 months the curves keep rising above 30%, only more slowly. The 36–47 month curve peaks at about 45% and the 48–59 month curve at about 35%, so the thresholds are in the reverse order to the text. The declines after the peak are 0.014 and 0.009 on the log-hazard scale, far smaller than their standard errors (0.04–0.07).
- **Evidence:** v7 `pfpr_curves.csv` (log HR relative to 20%). 12–23: 0.138 at 30%, 0.189 at 40%, 0.223 at 63%. 24–35: 0.122 rising to a maximum of 0.171 at 60.5%. 36–47: maximum 0.167 at 45.5% (se 0.040), 0.153 at 63% (se 0.053). 48–59: maximum 0.085 at 35.5% (se 0.040), 0.076 at 63% (se 0.067).
```latex
In older children, however, the increase in risk slowed markedly above prevalences of around 30\%, and for children aged 3 and 4 years the curves were essentially flat above prevalences of about 35--45\%, with at most a slight decline at higher prevalence that is well within the uncertainty.
```

### L88: 10% and "upper end" sentence (clarity)
- **Current:** "(\ref{fig:main_result2}). For an average \textit{PfPR}$_{2-10}$ of 40\% (upper end of prevalence estimates in sub-Saharan Africa now)"
- **Problem:** The bare `\ref` needs "Figure~". `\textit{PfPR}` does not match the `\textit{Pf}PR` style used elsewhere. 40% is above every 2024 national estimate (maximum 36.4%, Benin) and every Nigerian state estimate (maximum 34.7%). The 10% figures are correct.
- **Evidence:** v7 `annual_comparison/national_pfpr_2000_2024.csv` 2024: BEN 36.4, COD 36.1; `nigeria_states/state_totals_2024.csv` maximum 34.7; `attributable_fraction_by_age.csv` at 10%: 0.028 (<1 month) and 0.279–0.356 (12–59 months).
```latex
Under this model, for an average \textit{Pf}PR$_{2-10}$ of 10\%, approximately 30\% (28--36\%) of all deaths in children 1 year and older are caused by malaria, compared with less than 5\% in neonates (Figure~\ref{fig:main_result2}). For an average \textit{Pf}PR$_{2-10}$ of 40\% (close to the highest national estimates for 2024, 36\% in Benin and the Democratic Republic of the Congo), ...
```

### L88: peak attributable fraction (wrong number)
- **Current:** "more than half of all deaths between 1 and 4 years attributable to malaria, peaking at 56\% at 2–3 years"
- **Problem:** 56% is the 12–23 month value. The peak is 59% at 24–35 months. At 48–59 months the fraction is 47%, so "more than half" holds only for 12–47 months.
- **Evidence:** v7 `attributable_fraction_by_age.csv` at 40%: 12–23 0.558; 24–35 0.595 (0.537–0.646); 36–47 0.531; 48–59 0.472.
```latex
For an average \textit{Pf}PR$_{2-10}$ of 40\%, there is a sharper age-dependent increase in malaria deaths, with more than half of all deaths between 1 and 4 years of age (12--47 months) attributable to malaria, peaking at 59\% (95\% interval 54--65\%) at 2--3 years (24--35 months).
```

### L90–94: Figure 2 image and caption (outdated)
- **Current:** `figures/fig2_malaria_attributable_fraction_by_age.png`; caption begins "Mean model estimates of the malaria attributable share ..."
- **Problem:** The image is the v5 export. The caption says "mean" where these are point estimates, does not name the four prevalence levels or the interval bars, does not say that zero prevalence is an extrapolation, and uses `\textit{PfPR}`.
- **Evidence:** The Overleaf md5 is `89794f9b…`, identical to v5; the v7 `attributable_fraction_by_age.png` is `b100b37a…`. See v7 `paper_figures/CAPTIONS.md`, Figure 2.
- **Action:** copy the v7 `attributable_fraction_by_age.png`, then use:
```latex
\caption{Estimated malaria-attributable share of all-cause deaths within each age band at \textit{Pf}PR$_{2-10}$ of 10\%, 20\%, 30\% and 40\%, from the seven separate age-band models. Points are point estimates and vertical bars conditional pointwise 95\% intervals; lines connect the discrete age bands. Zero prevalence lies below the observed exposure range, so the counterfactual involves extrapolation.}
```

### L97: subgroup sensitivity description (wrong statement; v3 is final)
- **Current:** "(seasonal areas in the Sahel versus perennial areas in East Africa, earlier versus later surveys) ... steeper for children between 1 and 5 years of age"
- **Problem:** The Sahel curves are steeper from 6 months (6–11 months: 0.73 against 0.83), not from 1 year; at 1–5 months the Sahel matches the full sample. Eastern Africa is flatter in neonates and at 36–47 months, which the text omits. The subsets are defined geographically, not as "perennial areas", and the period split is at the median survey year, 2014.
- **Evidence:** `subgroups_dhsmics_map_gamma2_v3/pfpr_contrasts.csv` and REPORT.md (all 28 fits converged). HR 40→20 by age band:

  | Series | Hazard ratios, <1 month to 48–59 months |
  |---|---|
  | Sahel | 0.90, 0.91, 0.73, 0.70, 0.73, 0.78, 0.74 |
  | Eastern Africa | 1.00, 0.90, 0.85, 0.81, 0.84, 0.95, 0.88 |
  | Full sample | 0.94, 0.92, 0.83, 0.83, 0.86, 0.85, 0.92 |

  `sample_summary.csv`: Sahel 31 surveys in 8 countries; Eastern Africa 51 in 12; early 72; late 63. `subgroup_definitions.csv`: median 2014.
```latex
We conducted a series of sensitivity analyses in subsets of the primary sample: survey regions in the Sahel (boundary centroid at or north of 12$^{\circ}$N; 31 surveys in 8 countries), where transmission is highly seasonal; countries in Eastern Africa (51 surveys in 12 countries); and surveys conducted up to (72 surveys) or after (63 surveys) the median survey year, 2014. In general the subset fits gave very similar results, apart from survey regions in the Sahel, where the estimated relationship between prevalence and all-cause mortality was steeper for children between 6 months and 5 years of age (hazard ratios for a reduction in \textit{Pf}PR$_{2-10}$ from 40\% to 20\% of 0.70--0.78, against 0.83--0.92 in the full sample), and Eastern Africa, where it was flatter in neonates (1.00 against 0.94) and at 36--47 months (0.95 against 0.85) (Figure \ref{fig:subgroups}).
```

### L97: "very similar" and the 20%→0 contrasts (internal inconsistency)
- **Current:** "In general all the subgroup fits gave very similar results"
- **Problem:** This holds for the 40%→20% contrasts. The 20%→0 contrasts, which drive the burden estimates in the next subsection, vary substantially between subsets.
- **Evidence:** v3 `pfpr_contrasts.csv`, HR 20%→0. At 12–23 months: full sample 0.533, Eastern Africa 0.302, Sahel 0.700, early surveys 0.465. At 36–47 months: full sample 0.552, late surveys 0.364.
- **Suggested addition,** after the Sahel/Eastern Africa sentence:
```latex
Hazard ratios for a reduction from 20\% to zero prevalence, which involve extrapolation below the observed exposure range, varied more across subsets (at 12--23 months, 0.30 in Eastern Africa and 0.70 in the Sahel, against 0.53 in the full sample).
```

### L97: multiple imputation sentence (confirmed; v8 is final)
- **Current:** "A sensitivity analysis using multiply imputed datasets gave very similar results to the primary analysis restricted to complete cases"
- **Status:** **Confirmed.** v8 changes the 40%→20% hazard ratios by at most 0.028 against v7 (48–59 months: 0.919 → 0.948) and the 20%→0% hazard ratios by at most 0.018 (12–23 months). Rubin pooling over the 10 imputations changes them by at most 0.002, and the between-imputation share of variance is at most 2.6%. The current wording is supported; the suggested text adds the numbers.
- **Evidence:** `primary_map_regional17_dhsmics_imputed_gamma2_v8/REPORT.md`, `comparison_contrasts.csv` and `multiple_imputation/REPORT.md`.
- **Suggested text:**
```latex
A sensitivity analysis that imputed the missing covariates, retaining all 166 surveys in 40 countries (8,797,963 child age-band records and 123,419 deaths), gave very similar results to the complete-case primary analysis: hazard ratios for a reduction in \textit{Pf}PR$_{2-10}$ from 40\% to 20\% changed by at most 0.03, and pooling fits to 10 multiply imputed datasets with Rubin's rules changed them by at most 0.002 (Figure \ref{fig:MI}).
```

### L97: no-nutrition sentence (confirmed; v3 is final)
- **Current:** "we refitted the models on the same sample without adjusting for wasting and stunting; the estimated prevalence-mortality relationship was essentially unchanged"
- **Status:** **Confirmed.** All seven v3 refits converged. The largest absolute change in the 40%→20% hazard ratio is 0.008, at 1–5 months (0.916 → 0.908), and in the 20%→0% hazard ratio 0.031, also at 1–5 months (0.862 → 0.831). The current wording is supported.
- **Evidence:** `nutrition_adjustment_dhsmics_map_gamma2_v3/CONTRASTS.md` and `comparison_contrasts.csv`.
- **Optional addition:**
```latex
(largest change in the hazard ratio for a reduction in \textit{Pf}PR$_{2-10}$ from 40\% to 20\%: 0.008)
```

---

## Results: comparison with IHME and UN IGME, Tables 1–2 (L99–166)

Confirmed correct: every IHME and UN IGME cell in Tables 1 and 2; the 42-country scope; and the Nigeria comparators (138,111 IHME deaths; 165,329 UN IGME deaths). Every PfPR-ACM model value in this section's text and tables is stale: v3 (DHS-only) in the text and tables, v5 in the Figure 3 image.

### L101: extrapolation to zero prevalence (clarity)
- **Current:** "applying the age-band hazard ratios for a reduction from the observed prevalence to zero to national all-cause death counts"
- **Problem:** Zero PfPR lies below the observed exposure range in every band. The text does not say so.
- **Evidence:** v7 REPORT.md L56; `attributable_fraction_by_age.csv`: zero_below_observed_support is TRUE in all 28 rows.
```latex
... can be estimated by applying the age-band hazard ratios for a reduction from the observed prevalence to zero to national all-cause death counts (zero prevalence lies just below the range of survey-region \textit{Pf}PR$_{2-10}$ observed in every age band, so this counterfactual involves some extrapolation).
```

### L101: scope and figure citation (internal inconsistency)
- **Current:** "... for the year 2024, and compared the results with the cause-specific under-5 malaria death estimates of IHME and of UN IGME, expressed per 1,000 under-5 child-years (Figure~\ref{fig:IHME_v_model}A)."
- **Problem:** Panel A plots the model against IHME only. UN IGME appears only in Table 2 (counts) and panel C (pooled). The text limits the analysis to 2024, but Table 1 and panel C cover 2000–2024.
- **Evidence:** v7 `burden_comparison/fig4_burden_comparison.png` and `CAPTION.md`; `annual_comparison/annual_totals_2000_2024.csv`.
```latex
We applied the fitted age-band effects to IHME all-cause deaths by age band and population-weighted national \textit{Pf}PR$_{2-10}$ in the 42 countries with MAP coverage for each year from 2000 to 2024, and compared the results with the cause-specific under-5 malaria death estimates of IHME and of UN IGME, expressed per 1,000 under-5 child-years (national rates in 2024 against IHME in Figure~\ref{fig:IHME_v_model}A; national counts in Table~\ref{tab:country}; pooled trends in Table~\ref{tab:source-comparison} and Figure~\ref{fig:IHME_v_model}C).
```

### L102: "mostly agreed" (wrong statement)
- **Current:** "Although the country specific under 5 mortality rates mostly agreed between IHME and our new approach (\textit{Pf}PR ACM), there were substantial differences in several populous and high burden countries (\ref{tab:country})."
- **Problem:** In 2024 the model rate exceeds IHME's in 36 of 42 countries (all except BDI, GNQ, KEN, MRT, RWA and UGA). The median ratio is 1.58, and only 9 of 42 countries fall within 0.8–1.25. These are malaria rates, not under-5 all-cause rates. The bare `\ref` prints without "Table", and "PfPR ACM" should read "PfPR-ACM".
- **Evidence:** v7 `burden_comparison/country_rates_2024.csv`, recomputed in this audit: 36 of 42 above IHME, median ratio 1.582, Spearman ρ = 0.909. ρ is computed here, not a saved output; drop it if only saved outputs may be quoted.
```latex
Although national under-5 malaria mortality rates from IHME and our approach (PfPR-ACM model) were strongly correlated across the 42 countries (Spearman $\rho = 0.91$), the model estimate exceeded IHME's in 36 of the 42 countries (median ratio 1.58), with the largest absolute differences in several populous, high-burden countries (Table~\ref{tab:country}).
```

### L102: Nigeria (outdated)
- **Current:** "IHME estimate approximately 140 thousand deaths, and UN IGME estimate approximately 165 thousand deaths whereas our approach estimates 200 thousand deaths (46 and 22\% increase respectively)."
- **Problem:** 200,988 deaths, +46% and +22% are v3 values.
- **Evidence:** v7 `tables/country_comparison_2024.csv`, NGA: model 217,542; IHME 138,111; UN IGME 165,329. Ratios 1.575 and 1.316.
```latex
In Nigeria, IHME estimates approximately 138 thousand deaths and UN IGME approximately 165 thousand deaths, whereas our approach estimates 218 thousand deaths (58\% and 32\% higher, respectively).
```

### L102: Chad (outdated)
- **Current:** "In Chad we estimate a 3.6 fold higher total burden."
- **Problem:** 3.6 is the v3 ratio (17,472/4,884). The sentence also does not name the comparator.
- **Evidence:** v7 TCD: model 19,780; IHME 4,884; UN IGME 14,345. Ratios 4.05 and 1.38.
```latex
In Chad we estimate a 4.0-fold higher total burden than IHME (1.4-fold higher than UN IGME).
```

### L117–121: Table 1 death counts (outdated)
- **Current:** model column 1,104,990 / 665,559 / 565,616, with changes −48.8% and −15.0%
- **Problem:** These are v3 values. The IHME and UN IGME cells are correct.
- **Evidence:** v7 `burden_comparison/source_comparison_2000_2024.csv` and `source_comparison.latex.txt`: 1,170,371 / 725,613 / 620,072; declines 47.0% and 14.5%.
```latex
\quad 2000 & 1,170,371 & 563,657 & 724,276  \\
\quad 2015 & 725,613 & 416,122 & 416,815  \\
\quad 2024 & 620,072 & 428,147 & 440,123  \\
\quad Change 2000--2024 & $-47.0\%$ & $-24.0\%$ & $-39.2\%$  \\
\quad Change 2015--2024 & $-14.5\%$ & $+2.9\%$ & $+5.6\%$ \\
```

### L123, L131: rate unit and table note (clarity)
- **Current:** `\multicolumn{4}{@{}l}{\textit{Malaria mortality rate (per 1,000)}}`
- **Problem:** The unit leaves out the denominator. The table has no note on the shared denominator, on the different estimands, or on why IHME and UN IGME counts rose after 2015 while their rates fell.
- **Evidence:** v7 `annual_comparison/annual_totals_2000_2024.csv`: under-5 person-years 158.55 M (2015) and 177.45 M (2024), an increase of 11.9%.
```latex
\multicolumn{4}{@{}l}{\textit{Malaria mortality rate (per 1,000 under-5 child-years)}} \\
% after \end{tabular*}:
\par\vspace{0.5em}\begin{minipage}{\linewidth}\footnotesize All three columns are deaths before age five in the same 42 countries. Rates use the same annual under-five person-years, implied by the IHME all-cause death count and rate. Model deaths are IHME all-cause deaths in each age band multiplied by the estimated malaria-attributable fraction; IHME and UN IGME (CA-CODE 2026) are cause-specific malaria death estimates, so these are different estimands. IHME and UN IGME death counts rose between 2015 and 2024 while their rates fell, because under-five person-time in these countries grew by 12\% over that period.\end{minipage}
```

### L124–128: Table 1 rates (outdated)
- **Current:** model column 10.06 / 4.20 / 3.19, with changes −68.3% and −24.1%
- **Problem:** These are v3 values. They contradict Figure 3C and the abstract. The IHME and UN IGME cells are correct.
- **Evidence:** v7 `source_comparison_2000_2024.csv`: 10.656 / 4.576 / 3.494; declines 67.2% and 23.6%.
```latex
\quad 2000 & 10.66 & 5.13 & 6.59  \\
\quad 2015 & 4.58 & 2.62 & 2.63  \\
\quad 2024 & 3.49 & 2.41 & 2.48  \\
\quad Change 2000--2024 & $-67.2\%$ & $-53.0\%$ & $-62.4\%$  \\
\quad Change 2015--2024 & $-23.6\%$ & $-8.1\%$ & $-5.7\%$  \\
```
The generated fragment `burden_comparison/source_comparison.latex.txt` can replace the whole table.

### L105–131: Table 1 is never cited and the trend is missing from the Results (clarity)
- **Problem:** main.tex contains no `\ref{tab:source-comparison}`. The temporal comparison appears only in the abstract.
- **Evidence:** v7 `source_comparison_2000_2024.csv`. Rate declines 2000–2015: model 57.1%, IHME 48.9%, UN IGME 60.1%. Rate declines 2015–2024: 23.6%, 8.1%, 5.7%. 2024 death ratios: 1.448 and 1.409.
- **Suggested addition** after L102:
```latex
Pooled across the 42 countries, our estimated under-5 malaria mortality rate fell from 10.7 to 4.6 deaths per 1,000 child-years between 2000 and 2015 (a 57\% decline) and to 3.5 in 2024 (a further 24\% decline), whereas the IHME and UN IGME rates fell by 49\% and 60\% to 2015 but by only 8\% and 6\% thereafter (Table~\ref{tab:source-comparison}; Figure~\ref{fig:IHME_v_model}C). In 2024 the model gives 620,072 under-5 malaria deaths in these countries, 45\% more than IHME (428,147) and 41\% more than UN IGME (440,123).
```

### L134: Nigerian states (outdated)
- **Current:** "greater than IHME's estimates in 35 of the 37 states ... the median state ratio is 1.6 in the three northern zones and 1.2 in the three southern zones."
- **Problem:** These are v3 values (35 of 37; medians 1.63 and 1.16). The 37 units include the Federal Capital Territory. "PfPR ACM" should read "PfPR-ACM".
- **Evidence:** v7 `nigeria_states/README.md`: "exceeds IHME in 36 of 37 states"; medians 1.76 (north) and 1.27 (south); state sum 216,626 against IHME 132,138 (ratio 1.64). Lagos is the only state below IHME (ratio 0.72).
```latex
Within Nigeria, applying the same age-band effects to state-level IHME all-cause deaths by age and state \textit{Pf}PR$_{2-10}$ gives malaria mortality estimates greater than IHME's in every one of the 36 states and the Federal Capital Territory except Lagos (Figure~\ref{fig:IHME_v_model}B); summed over states, the model gives 216,626 deaths against 132,138 IHME malaria deaths (ratio 1.64). The differences are concentrated in the northern states: the median state ratio is 1.8 in the three northern zones and 1.3 in the three southern zones.
```

### L138: Figure 3 image (outdated)
- **Current:** `figures/fig3_burden_comparison.png`
- **Problem:** The image is the v5 export. Its sidecar `figures/fig3_burden_comparison_caption.md` still carries v5 state totals (214,127; ratio 1.62). There is also a stale duplicate, `figures/fig4_burden_comparison_caption.md`.
- **Evidence:** The Overleaf md5 is `25e87587…`, identical to v5 `burden_comparison/fig4_burden_comparison.png`; the v7 image is `21bb035c…`.
- **Action:** copy v7 `burden_comparison/fig4_burden_comparison.png` to `figures/fig3_burden_comparison.png` and refresh the sidecar. No TeX change is needed.

### L139: Figure 3 caption (wrong statement)
- **Current:** "... between the prevalence-ACM model and IHME. ... C: average predicted mortality in children under 5 in sub-Saharan Africa since 2004."
- **Problem:** Panel C covers 2000–2024, not "since 2004". It shows three series (model, IHME, UN IGME), each a pooled 42-country rate rather than an "average" for sub-Saharan Africa. The caption omits UN IGME and the shared denominator, and says "prevalence-ACM".
- **Evidence:** the v7 image; `burden_comparison/CAPTION.md`; `R_cbh/reporting/11_burden_comparison_figure.R` asserts years 2000:2024.
```latex
\caption{Malaria mortality before age 5 from the PfPR-ACM model compared with IHME and UN IGME (CA-CODE 2026) cause-specific estimates, expressed as deaths per 1,000 under-five child-years. A: model against IHME national rates in 2024 in the 42 countries; the dashed line denotes equality. B: model against IHME for the 36 Nigerian states and the Federal Capital Territory in 2024, coloured by geopolitical zone; summed over states the model gives 216,626 deaths against 132,138 IHME malaria deaths (ratio 1.64). C: annual under-five malaria mortality, 2000--2024, for the PfPR-ACM model, IHME and UN IGME, pooled across the same 42 countries. In every panel each source's deaths are divided by the under-five person-years implied by the IHME all-cause death counts and rates, so the sources share a denominator. Model-attributable reductions in all-cause mortality and cause-specific malaria deaths are different estimands, and the model shares the IHME all-cause inputs, so agreement is not independent validation.}
```

### L151–162: Table 2 body (outdated)
- **Current:** Nigeria 200,988 … Central African Republic 8,544 … Uganda 21,899 … All 42 countries 565,616
- **Problem:** The model column and the top-ten selection are v3. In v7, DR Congo (rank 5) and Côte d'Ivoire (rank 10) enter the list; Central African Republic (rank 11) and Uganda (rank 22) leave it. The IHME and UN IGME cells shown are correct. Two orderings are fragile: DR Congo versus Mali differ by 85 deaths, and Côte d'Ivoire versus CAR by 40.
- **Evidence:** v7 `tables/country_comparison_2024.csv` (rank_absolute_difference_from_ihme) and `country_comparison_2024.latex.txt`.
```latex
Nigeria & 217,542 & 138,111 & 165,329 \\
Niger & 37,567 & 20,343 & 23,877 \\
Chad & 19,780 & 4,884 & 14,345 \\
Mozambique & 19,296 & 10,285 & 10,315 \\
DR Congo & 69,518 & 61,774 & 105,372 \\
Mali & 20,604 & 12,945 & 6,984 \\
South Sudan & 11,251 & 4,779 & 3,980 \\
Cameroon & 22,637 & 16,295 & 6,286 \\
Angola & 20,839 & 14,652 & 4,734 \\
C\^{o}te d'Ivoire & 18,233 & 12,848 & 5,753 \\
\midrule
All 42 countries & 620,072 & 428,147 & 440,123 \\
```

### L164–165: Table 2 note (clarity)
- **Problem:** Table 2 has no note on the UN IGME series, on the total covering all 42 countries, or on rounding.
- **Evidence:** v7 `tables/country_comparison_2024.latex.txt`.
```latex
\par\smallskip\begin{minipage}{\linewidth}\footnotesize UN IGME denotes the under-five CA-CODE 2026 cause-specific malaria series. Totals include all 42 countries, not only the ten shown. Death counts are rounded after estimation.\end{minipage}
```

---

## Methods: data and covariates (L192–225)

Confirmed correct: the seven age bands; the list of 13 regional and 4 national covariates; the log transforms; standardisation; the within-survey regional-mean fill; fixed posterior-median HIV imputation; GPWv4 2020 density times cell area for national PfPR; and the IHME all-cause source.

### L194–195: birth-history sources (outdated)
- **Current:** "Complete birth histories were obtained from Demographic and Health Surveys and Malaria Indicator Surveys (DHS/MIS) Births Recodes \cite{dhs}."
- **Problem:** 37 of the 135 analysed surveys are MICS, which the Methods never mention. `ref.bib` has no MICS entry.
- **Evidence:** v7 `study_flow/CAPTION.md`: "combines 98 DHS and 37 MICS surveys".
```latex
Complete birth histories were obtained from the Births Recodes of the Demographic and Health Surveys and Malaria Indicator Surveys (DHS/MIS) \cite{dhs} and from UNICEF Multiple Indicator Cluster Surveys (MICS) that recorded a complete birth history \cite{mics}.
```

### L196: DHS-only covariate paragraph (outdated; grammar)
- **Current:** "Child, maternal and household characteristics obtained from the DHS/MIS surveys comprised from the Births Recodes, ..."
- **Problem:** "Comprised from" is ungrammatical. Nothing describes how MICS surveys were converted, how their regions were mapped, or how their covariates were derived.
- **Evidence:** `docs/ANALYSIS_PLAN.md` §2.1–2.4; `data/derived_mics/region_map.csv`: 366 labels (291 exact, 75 overrides); `R_mics/README.md`: 93 MICS assessed, 47 with a complete birth history.
```latex
For DHS/MIS surveys, child, maternal and household characteristics comprised, from the Births Recodes, maternal age at first live birth, ... [rest unchanged]. \par Of the 93 MICS surveys from sub-Saharan Africa that we assessed, 47 (28 countries) recorded a complete birth history; the other 46 recorded only summary birth histories and cannot support the age-band design. Each MICS birth history was converted to the DHS Births Recode layout and processed with the same age-band construction and eligibility rules as the DHS/MIS data. The 366 MICS region labels were mapped to analysis regions by exact name (291) or a reviewed correspondence table (75); analysis regions were unions of DHS boundary polygons for the same country where these exist, and admin-1 polygons otherwise (for example, Guinea-Bissau and the Central African Republic). MICS does not publish survey-region indicator tables, so all 13 regional covariates were computed from the MICS microdata using the DHS definitions; education, recorded as level and grade, was converted to completed years, and facility delivery refers to the last birth in the two years before interview.
```

### L199–201: MAP extraction (clarity)
- **Current:** "population-weighted prevalence estimates were taken from the Malaria Atlas Project (MAP) modelled surfaces ... for 2000 - 2024 \cite{map,map2025}."
- **Problem:** The survey-region values were computed by the authors from the annual MAP rasters, using GPWv4 2020 density weights without cell area. The text implies MAP supplied population-weighted values. The cited papers cover 2000–17 and 2000–22. "2000 - 2024" needs an en dash.
- **Evidence:** `R_dhs/00_config.R`: MAP release `Malaria__202508_Global_Pf_Parasite_Rate`; `ANALYSIS_PLAN.md` §2.2; `R_mics/06_polygons_pfpr_registry.R`.
```latex
\emph{P.\ falciparum} prevalence estimates were taken from the Malaria Atlas Project (MAP) annual modelled surfaces of the parasite rate standardised to ages 2 to 10 years (\textit{Pf}PR$_{2-10}$) for 2000--2024 (MAP release \texttt{Malaria\_\_202508\_Global\_Pf\_Parasite\_Rate}) \cite{map,map2025}. Survey-region values were means of each annual surface over the region's boundary polygons, weighted by the Gridded Population of the World version 4 (GPWv4) 2020 population density \cite{gpwv4}.
```

### L200: exposure assignment sentence (grammar)
- **Current:** "Each child trajectory, each of their recorded age-bands ... the estimated start needed to be within 60 months of the interview."
- **Problem:** The sentence is ungrammatical and "estimated start" is vague. The content is correct.
```latex
Each age band in a child's trajectory (from birth until interview or death) was assigned the \textit{Pf}PR$_{2-10}$ estimate for the mother's survey region of residence at interview and the calendar year in which the child entered that band. To be eligible for the primary analysis, entry into the band had to fall within the 60 months before the interview.
```

### L201: Nigerian state PfPR extraction (omission)
- **Current:** "For national burden calculations, country-level prevalence was calculated using population-weighted means ... and grid-cell areas \cite{gpwv4}."
- **Problem:** State prevalence for the Nigeria comparison uses MAP admin-1 boundaries and density weights without cell area. This extraction difference contributes to the state sum being 0.42% below the national estimate.
- **Evidence:** `nigeria_states/README.md`; `nigeria_states/national_reconciliation_2024.csv` note.
```latex
... and grid-cell areas \cite{gpwv4}. Nigerian state prevalence for 2024 used MAP admin-1 boundaries and the same density surface without grid-cell areas.
```

### L204: wrong citation for IHME all-cause data (citation)
- **Current:** "... obtained from the IHME Global Burden of Disease data exports and used to calculate malaria-attributable mortality \cite{igme}."
- **Problem:** The key `igme` is the UN IGME 2024 child mortality report. The data are IHME GBD exports, for which the key is `ihme`.
- **Evidence:** `ref.bib`, entries `igme` and `ihme`.
```latex
National age-specific all-cause death counts and rates were obtained from the IHME Global Burden of Disease data exports and used to calculate malaria-attributable mortality \cite{ihme}.
```

### L205: CA-CODE citation (citation)
- **Current:** "... the CA-CODE 2026 series disseminated through the UN Inter-agency Group for Child Mortality Estimation (UN IGME) portal \cite{ihme,Villavicencio:24}."
- **Problem:** Villavicencio:24 describes the CA-CODE portal with estimates for 2000–21. The release used here covers 2000–2024 and is documented as the BMJ study (doi 10.1136/bmj-2025-088686), which has no entry in `ref.bib`.
- **Evidence:** `ref.bib` Villavicencio:24; `R_cbh/reporting/07_annual_mortality_comparison.R`; `burden_comparison/CAPTION.md`.
```latex
Comparator estimates of under-five malaria deaths were obtained from IHME \cite{ihme} and from the CA-CODE 2026 release (estimates for 2000--2024) disseminated through the UN Inter-agency Group for Child Mortality Estimation (UN IGME) portal \cite{cacode2025,Villavicencio:24}.
```

### L205: "only country for which subnational estimates are available" (wrong statement)
- **Current:** "In addition we used IHME estimates for state specific all-cause mortality by age band in Nigeria (only country for which subnational estimates are available)."
- **Problem:** GBD also publishes subnational estimates for other sub-Saharan countries, including Ethiopia, Kenya and South Africa. This is a literature fact, not checkable in the project files. Nigeria is simply the one country analysed subnationally here. IHME supplies GBD age groups, not the model's bands, and the state malaria comparator is not mentioned.
- **Evidence:** `nigeria_states/ihme_source.json`; `state_totals_2024.csv` (37 rows).
```latex
In addition, we used IHME estimates of all-cause deaths and death rates by age group, and of under-five malaria death rates, for the 36 states and the Federal Capital Territory of Nigeria, the only country we analysed subnationally \cite{ihme}.
```

### L196, L208: missing data-source citations (citation)
- **Problem:** WUENIC (L196; only a URL at L220), the World Bank series (L208) and the UNAIDS/UNICEF and MICS sources have no citations or `ref.bib` entries.
- **Evidence:** a grep of `ref.bib` for unaids, aidsinfo, mics, world bank and wuenic finds nothing relevant.
```latex
... political stability indicators were obtained from the World Bank World Development Indicators and Worldwide Governance Indicators \cite{wdi,wgi}.
```

### L207–208: child HIV incidence source, Liberia derivation (outdated)
- **Current:** "National child HIV incidence estimates for ages 0-14 years, expressed as new infections per 1,000 uninfected population, were obtained from the UNAIDS 2025 estimates distributed by UNICEF."
- **Problem:** The UNICEF workbook has no Liberia series, and the Methods do not describe v7's derived Liberia rate. That rate is calculated from UNAIDS AIDSinfo counts in the 2026 estimates round.
- **Evidence:** `results/cbh/hiv_incidence/liberia_aidsinfo/REPORT.md`: the formula; derived/published median 0.988 (IQR 0.944–1.055) across 761 country-years in 34 countries with at least 200 new infections. (The v7 REPORT.md said 741 because of a hard-coded count; since 24 September it reads the count from the check file.)
```latex
National child HIV incidence (new infections at ages 0--14 years per 1,000 uninfected population, both sexes) was obtained from the UNAIDS 2025 estimates distributed by UNICEF \cite{unaids2025}; the corresponding adolescent incidence at ages 15--19 years informed the imputation of missing child incidence. The UNICEF workbook has no incidence series for Liberia, so its child incidence was calculated from the UNAIDS 2026 estimates of annual new HIV infections at ages 0--14 years published on AIDSinfo \cite{aidsinfo}, as $1{,}000\times$ new infections divided by the under-five person-years implied by the IHME all-cause death count and rate, less the children aged 0--4 years living with HIV. Applied to the other countries, this formula reproduces the published child rates (median ratio of derived to published rate 0.99 across 761 country-years with at least 200 new infections, in 34 countries). National GDP per capita, health expenditure per capita and political stability indicators were obtained from World Bank data series.
```

### L212: LLIN and malnutrition sentence (typo/clarity)
- **Current:** "We did not adjust for LLIN coverage (measured in DHS/MIS surveys) ... with and without measured of malnutrition (wasting and stunting)."
- **Problem:** "Measured of" is a typo. LLIN coverage is also measured in MICS. The sensitivity analysis omits both nutrition measures on the same sample.
- **Evidence:** `ANALYSIS_PLAN.md` §2.4.
```latex
We did not adjust for LLIN coverage (measured in DHS, MIS and MICS surveys) because ... To assess the potential impact of this, we carried out a sensitivity analysis omitting both measures of malnutrition (wasting and stunting prevalence), fitted to the same analysis sample.
```

### L214: reference period of the regional summaries (wrong statement)
- **Current:** "We note that regional summaries describe conditions at the time of the survey, not values in the past 60 months."
- **Problem:** Facility delivery and birth interval refer to births in the five years before the survey. Maternal summaries are computed among mothers with a birth in the preceding 60 months. What the summaries are not is values at the time each child entered a band.
- **Evidence:** `ANALYSIS_PLAN.md` §2.4 (table, and the sentence "Do not relabel them as reconstructed historical values at band entry").
```latex
We note that regional summaries describe conditions at the time of the survey, or over each indicator's own recall window (for example, births in the five years before the survey for facility delivery and birth intervals), rather than values at the time each child entered an age band.
```

### L217: sampling weights (outdated)
- **Current:** "Recode-derived regional summaries used DHS women's sampling weights."
- **Problem:** This covers DHS only. MICS summaries use the MICS women's, household and children's weights.
- **Evidence:** `ANALYSIS_PLAN.md` §2.4 (wmweight, hhweight, chweight).
```latex
DHS recode-derived regional summaries used DHS women's sampling weights, and published DHS regional indicators retained their stated denominators and reference periods. MICS regional summaries were computed from the microdata with the MICS women's, household and children's sampling weights. Regional summaries were calculated before complete-case selection, excluding missing individual responses from each variable's denominator.
```

### L220: vaccination substitution (outdated)
- **Current:** "Missing regional DTP3 and measles coverage was replaced by the corresponding national UNICEF estimate for the survey year"
- **Problem:** This repeats L196 and calls the source "UNICEF" instead of WHO/UNICEF (WUENIC). It omits three included MICS surveys that use national values for every region because their recall doses were unusable, not because values were missing.
- **Evidence:** `R_mics/09_regional_covariates.R`: `card_only <- c('MC_GIN2016','MC_COM2022','MC_TCD2019')`. All three are included (`survey_selection.csv`).
```latex
Missing regional DTP3 and measles coverage was replaced by the corresponding WHO/UNICEF national estimate (WUENIC) for the survey year (\url{https://immunizationdata.who.int/}); the same national values were used for all regions of three MICS surveys (Guinea 2016, Comoros 2022 and Chad 2019) whose recall vaccination doses were unusable.
```

### L220: "Liberia and Sao Tome had no data and were excluded" (outdated; wrong for Liberia)
- **Problem:** Liberia's 2007, 2013 and 2019–20 DHS surveys are in v7. Only its 2009 MIS is excluded, because MIS surveys do not report anthropometry or facility delivery. São Tomé and Príncipe does have surveys (one DHS, two MICS); they are excluded only because the UNAIDS/UNICEF workbook has no child or adolescent incidence series for it. The size of the complete-case exclusion is not stated.
- **Evidence:** v7 `study_flow/survey_selection.csv` complete-case rows: LB51FL 31,692; LB6AFL 44,084; LB7AFL 32,887; LB5AFL 0. The three STP surveys have 0 rows between them; their records sum to 33,665, matching the "child HIV incidence" count in `study_flow/CAPTION.md`. `flow_counts.csv`: with_map 8,797,963; missing_covariates 1,190,841 (13.5%). `CAPTION.md`: 31 surveys (22 DHS/MIS, 9 MICS) contribute no records.
```latex
Records with remaining missing required covariates were excluded: 1,190,841 of 8,797,963 MAP-eligible records (13.5\%), including all records from 31 surveys (22 DHS/MIS, among them all five MIS surveys, and 9 MICS). Imputation uncertainty was not propagated in the primary analysis. The UNAIDS/UNICEF workbook has no child or adolescent HIV incidence series for S\~{a}o Tom\'{e} and Pr\'{i}ncipe, so its three surveys (one DHS and two MICS; 33,665 child age-band records) were excluded. Liberia's 2007, 2013 and 2019--20 DHS surveys were included using the derived child incidence described above; its 2009 MIS was excluded because, like all MIS surveys, it does not report the anthropometry and facility-delivery indicators.
```
This replacement depends on the L207–208 text being adopted.

### L224 and `figures/Causal_diagram.pdf`: causal diagram (internal inconsistency)
- **Current:** "(the malaria inoculation rate, proxied by \textit{Pf}PR$_{2-10}$)"
- **Problem:** L68 uses the standard term "entomological inoculation rate". The diagram labels a node "HIV prevalence", but the model uses child HIV incidence. Malnutrition, discussed at L212, has no node.
- **Evidence:** text extracted from `Causal_diagram.pdf`; `model_formula.txt`: `z_log_hiv_incidence`.
```latex
\caption{Assumed causal structure relating malaria transmission intensity (the entomological inoculation rate, proxied by \textit{Pf}PR$_{2-10}$) to all-cause child mortality, used to guide covariate selection. Solid arrows denote direct causes and dashed arrows effect modifiers. The HIV node is measured by national child (0--14 years) HIV incidence. LLIN: long-lasting insecticide-treated bednets; SES: socio-economic status.}
```
Also relabel the diagram node as "Child HIV incidence".

---

## Methods: model, attributable mortality, sensitivity (L229–242)

Confirmed correct: the unit of observation; the band-inclusion rule at L231; the record-count examples; the binomial cloglog model with a band-width offset; basis dimensions 5 and 6; the three random intercepts; `bam`, fREML, discrete fitting, γ = 2 and mgcv 1.9 (v7 `session_info.txt`: mgcv_1.9-4); the HR and attributable-fraction formulas; point-estimate national totals; the shared IHME person-time denominator; the 36 states plus FCT; the 42 countries; and knots at quantiles of each subset's own values.

### L229: subsection heading (clarity)
- **Current:** `\subsection*{Statistical model and model selection}`
- **Problem:** No model-selection step is described; the specification is fixed. The AIC selection belongs to the archived model.
- **Evidence:** `R_cbh/primary/specification.R`. In `fit_manifest.csv`, selected_version only distinguishes original from strict_restart.
```latex
\subsection*{Statistical model}
```

### L231: inclusion-rule sentence (grammar; checked by the synthesis step)
- **Current:** "this inclusion criteria applied equally to children who died ... (e.g. if a child was 4 and half and alive at the time of the interview"
- **Problem:** "Criteria" is plural, and "4 and half" is missing "a". The rule itself is correct.
- **Evidence:** `R_cbh/R/child_bands.R` (window and complete-band conditions).
```latex
these inclusion criteria applied equally to children who died, so a death in a band cut short by the interview was not counted (e.g.\ a child aged four and a half years and alive at interview contributed the six younger age bands but not the 48--59 month band).
```

### L233: SText1 cross-reference (internal inconsistency)
- **Current:** "Full model equations are given in the Supplement (SText~1)."
- **Problem:** `SText1_model_specification.tex` still describes the archived region-level negative-binomial GAM: births offset, post-neonatal deaths, AIC over M1–M4, n = 600 survey-regions in 23 countries.
- **Evidence:** the SText1 file; v7 `model_formula.txt`; `R_cbh/primary/01_fit.R` (`binomial(link='cloglog')`, fREML, discrete, γ = 2).
- **Suggested body for SText1:**
```latex
For child $i$ entering age band $g$ (width $\Delta_g$ years) at calendar time $t_i$ in survey-region $r$ of survey $s$ and country $c$, $\Pr(y_i=1)=1-\exp\{-\exp(\eta_i)\}$, with $\eta_i=\log\Delta_g+\alpha_g+f_g(\mathrm{PfPR}_{r,t_i})+h_g(t_i)+\mathbf{x}_{i}^{\top}\boldsymbol{\beta}_g+u_{g,s}+v_{g,c}+w_{g,r}$, where $f_g$ and $h_g$ are penalised cubic regression splines (basis dimensions 5 and 6), $\mathbf{x}_i$ holds the 17 standardised covariates, and $u_{g,s}\sim N(0,\sigma^2_{u,g})$, $v_{g,c}\sim N(0,\sigma^2_{v,g})$, $w_{g,r}\sim N(0,\sigma^2_{w,g})$ are independent.
```
Until SText1 is replaced, delete the L233 sentence.

### L233: interval conditioning (clarity; new from the synthesis spot-check)
- **Current:** "Confidence intervals used the coefficient covariance conditional on the estimated smoothing parameters and the single HIV imputation."
- **Problem:** The intervals also condition on the exposure values and the filled regional covariates.
- **Evidence:** v7 `tables/age_band_results.latex.txt` note: "Intervals condition on fitted smoothing parameters, exposure, the fixed HIV imputation and filled regional covariates".
```latex
Confidence intervals used the coefficient covariance conditional on the estimated smoothing parameters, the exposure values, the single child HIV incidence imputation and the filled regional covariates.
```

### L237: IHME age mapping (wrong statement; checked by the synthesis step)
- **Current:** "using the 1-5 and 6-11 month groups directly,  and assuming the 2-4 year group's deaths were spread equally across the three years"
- **Problem:** The mapping leaves out the 12–23 month band, which is also taken directly. The 2–4 year group divides both deaths and person-time equally, so the three oldest bands share one baseline rate.
- **Evidence:** v7 `burden/country_age_estimates.csv` age_mapping. 12–23: "direct age match". 24–35, 36–47 and 48–59: "IHME 2-4 rate shared; deaths and person-time divided equally".
```latex
IHME age groups were mapped to the model bands by combining early and late neonatal deaths (0--27 days) for the $<$1 month band, using the 1--5 month, 6--11 month and 12--23 month groups directly, and dividing the deaths and person-time of the 2--4 year group equally across the 24--35, 36--47 and 48--59 month bands, which therefore share one baseline rate.
```

### L237: years covered (wrong number)
- **Current:** "Estimates were produced annually for 2004-2024, in the 42 countries with national MAP estimates."
- **Problem:** The estimates start in 2000; Table 1 and Figure 3C use 2000. The count of 42 is correct: 45 countries less CPV, LSO and STP, which have no MAP estimate.
- **Evidence:** v7 `annual_comparison/annual_totals_2000_2024.csv`: 2000–2024, with 42 countries in every year.
```latex
Estimates were produced annually for 2000--2024, in the 42 countries with national MAP estimates.
```

### L237: "reconciled" (clarity; checked by the synthesis step)
- **Current:** "with the state sum reconciled against the national estimate"
- **Problem:** "Reconciled" suggests the state estimates were adjusted to the national total; they were only compared. The state IHME inputs sum exactly to the national inputs, and the attributable sum is 0.42% below the national estimate.
- **Evidence:** `nigeria_states/national_reconciliation_2024.csv`: IHME under-5 deaths 733,038 (national) and 733,038 (state sum); attributable deaths 217,542 (national) and 216,626 (state sum); relative difference −0.0042.
```latex
... using state-level IHME all-cause deaths by age and state MAP prevalence; the state IHME inputs sum to the national inputs, and the summed state estimate (216,626 deaths) was within 0.5\% of the national estimate (217,542).
```

### L240: heading numbering (typo; checked by the synthesis step)
- **Current:** `\subsection{Sensitivity analyses}`
- **Problem:** This heading is numbered; the other Methods subsections (L192, L210, L229, L235) use `\subsection*`.
```latex
\subsection*{Sensitivity analyses}
```

### L241: description of the sensitivity analyses (wrong statement)
- **Current:** "how changes in the paramerisation of the models, or if subsetting the data (seasonal versus non seasonal, West Africa vs East Africa, earlier vs later surveys)"
- **Problem:** No West Africa subset and no non-seasonal complement were fitted. "Paramerisation" is a typo. The no-nutrition and imputed analyses are not named.
- **Evidence:** v3 `sample_summary.csv` lists the subsets map_full, sahel, east_africa, early and late.
```latex
We carried out a series of sensitivity analyses to assess whether the estimated relationship between prevalence and age-specific all-cause mortality changed when the models were refitted in subsets of the data (the highly seasonal Sahel, Eastern Africa, and earlier versus later surveys), when wasting and stunting were omitted from the adjustment set, or when missing covariates were imputed rather than the affected records excluded.
```

### L240–242: the no-nutrition analysis is not described (omission)
- **Problem:** L212 and L97 promise an analysis without wasting and stunting, but this subsection never describes it. The v3 refit found a largest change of 0.008 in the 40%→20% hazard ratio, at 1–5 months.
- **Evidence:** `nutrition_adjustment_dhsmics_map_gamma2_v3/dropped_covariates.csv`: wasting_pct and stunting_pct. The v2 CAPTION.md gives 15 covariates.
```latex
In a second sensitivity analysis we refitted the seven models on the same sample without regional wasting and stunting prevalence (15 covariates), because malnutrition may mediate as well as confound the relationship between transmission and death. In a third sensitivity analysis we imputed every remaining covariate gap rather than excluding records: ...
```

### L242: Sahel and Eastern Africa subset sizes (outdated; v3 is final)
- **Current:** "23 surveys, 7 countries, 14,545 deaths ... 38 surveys, 12 countries, 23,583 deaths"
- **Problem:** These are the DHS-only v1 subsets.
- **Evidence:** `subgroups_dhsmics_map_gamma2_v3/sample_summary.csv`. Sahel: 31 surveys, 8 countries, 18,919 deaths. Eastern Africa: 51 surveys, 12 countries, 31,589 deaths.
```latex
Sahelian survey regions (boundary centroid at or north of 12$^\circ$N and west of 36$^\circ$E, excluding the Horn of Africa; 31 surveys, 8 countries, 18,919 deaths), Eastern Africa (UN M49 delimitation of Eastern Africa: 51 surveys, 12 countries, 31,589 deaths),
```

### L242: period split (outdated)
- **Current:** "surveys before versus after the median survey year of 2013 (50 and 45 surveys, respectively)"
- **Problem:** On v7 the median is 2014, and surveys in the median year go in the early group, so "before versus after" misdescribes the split. The counts are 72 and 63.
- **Evidence:** v3 `subgroup_definitions.csv` (survey_year <= 2014) and `sample_summary.csv` (early 72, late 63).
```latex
and surveys conducted up to and including the median survey year of 2014 versus later surveys (72 and 63 surveys, respectively).
```

### L242: national imputations listed (outdated)
- **Current:** "health expenditure for Zimbabwe 2000–2009 and all countries in 2024 from a generalised additive model on the observed panel"
- **Problem:** The 40-country DHS + MICS panel also imputes values for Somalia and South Sudan.
- **Evidence:** `primary_map_regional17_dhsmics_imputed_gamma2_v6/REPORT.md` caption.
```latex
the missing 2001 political-stability round by linear interpolation; health expenditure for Zimbabwe 2000--2009, Somalia 2000--2012, South Sudan before 2017 and all countries in 2024, and South Sudan's GDP before 2008 and political stability before independence, from generalised additive models fitted to the observed national panel;
```

### L242: HIV imputation for Liberia (outdated)
- **Current:** "and child HIV incidence for Liberia and São Tomé from the incidence model with a latent adolescent series"
- **Problem:** In v8 (and in v7), Liberia uses the rate derived from UNAIDS counts.
- **Evidence:** `R_cbh/sensitivity/imputation_dhsmics/settings.R` v8 header; `05_propagate.R`.
```latex
and child HIV incidence for S\~{a}o Tom\'{e} and Pr\'{i}ncipe from the incidence model with a latent adolescent series (Liberia's rate was derived from UNAIDS counts, as in the primary analysis).
```

### L242: imputed sample size (outdated)
- **Current:** "This retained all 6,357,802 eligible records (120 surveys, 36 countries, 90,938 deaths)"
- **Problem:** These are the DHS-only v4 counts.
- **Evidence:** v6 `prepared_sample.csv` and the expected counts in the v8 `settings.R`: 8,797,963 records, 166 surveys, 40 countries, 123,419 deaths.
```latex
This retained all 8,797,963 MAP-eligible records (166 surveys, 40 countries, 123,419 deaths).
```

### L242: Rubin's rules (clarity)
- **Current:** "imputation uncertainty was propagated by refitting with each of the ten imputed datasets and pooling with Rubin's rules."
- **Problem:** The text omits three details: smoothing parameters were fixed at the point fit, the point fit used the mean of the ten imputations, and each refit used one posterior draw of child HIV incidence. In v8, pooling changed HRs by at most 0.002, with a between-imputation variance share of at most 2.6%. Liberia's HIV rate, derived from UNAIDS counts, is held fixed across imputations.
- **Evidence:** `05_propagate.R`; v6 `multiple_imputation/REPORT.md`.
```latex
Imputation uncertainty was propagated by refitting the seven models in each of the ten imputed datasets (each also using one posterior draw of imputed child HIV incidence), with smoothing parameters fixed at the point fit (which used the mean of the ten imputations), and pooling the log hazard ratios with Rubin's rules.
```

### Consolidated replacement for L240–242
```latex
\subsection*{Sensitivity analyses}
We carried out a series of sensitivity analyses to assess whether the estimated relationship between prevalence and age-specific all-cause mortality changed when the models were refitted in subsets of the data (the highly seasonal Sahel, Eastern Africa, and earlier versus later surveys), when wasting and stunting were omitted from the adjustment set, or when missing covariates were imputed rather than the affected records excluded.
We refitted the 7 age-band models in four subsets of the primary sample: Sahelian survey regions (boundary centroid at or north of 12$^\circ$N and west of 36$^\circ$E, excluding the Horn of Africa; 31 surveys, 8 countries, 18,919 deaths), Eastern Africa (UN M49 delimitation of Eastern Africa: 51 surveys, 12 countries, 31,589 deaths), and surveys conducted up to and including the median survey year of 2014 versus later surveys (72 and 63 surveys, respectively). Subset fits used the same specification, with covariates scaled as in the full sample and knots at quantiles of the subset's own predictor values.
In a second sensitivity analysis we refitted the seven models on the same sample without regional wasting and stunting prevalence (15 covariates), because malnutrition may mediate as well as confound the relationship between transmission and death.
In a third sensitivity analysis we imputed every remaining covariate gap rather than excluding records: whole-survey gaps in regional indicators by chained-equation multiple imputation (predictive mean matching, 10 imputations); the missing 2001 political-stability round by linear interpolation; health expenditure for Zimbabwe 2000--2009, Somalia 2000--2012, South Sudan before 2017 and all countries in 2024, and South Sudan's GDP before 2008 and political stability before independence, from generalised additive models fitted to the observed national panel; and child HIV incidence for S\~{a}o Tom\'{e} and Pr\'{i}ncipe from the incidence model with a latent adolescent series (Liberia's rate was derived from UNAIDS counts, as in the primary analysis). This retained all 8,797,963 MAP-eligible records (166 surveys, 40 countries, 123,419 deaths). Imputation uncertainty was propagated by refitting the seven models in each of the ten imputed datasets (each also using one posterior draw of imputed child HIV incidence), with smoothing parameters fixed at the point fit (which used the mean of the ten imputations), and pooling the log hazard ratios with Rubin's rules.
```

---

## Supplementary material cited from these sections (outside the assigned range)

### S1 (L270–271): study-flow figure and caption (outdated; internal inconsistency)
- **Problem:** `Supplementary Figures/sfig_study_flow.png` is the v5 diagram, showing 132 surveys, 36 countries and 7,498,459 records. The caption describes a "data sources" figure rather than sample inclusion. It says "child HIV prevalence (UNAIDS 0-14 case counts divided by the World Bank population estimates)"; the model uses child HIV incidence, and only Liberia's rate is count-derived, with an IHME-based denominator. MICS is not mentioned.
- **Evidence:** Overleaf md5 `77cb1109…` = v5; the v7 `study_flow/study_flow_diagram.png` is `e3b6d421…`. See v7 `study_flow/CAPTION.md` and `flow_counts.csv`.
- **Action:** copy the v7 diagram into place, then use:
```latex
\caption{Sample inclusion for the primary analysis, from 171 surveys with complete birth histories in 41 countries to 7,607,122 child age-band records with 103,987 deaths from 135 surveys (98 DHS, 37 MICS) in 37 countries. Counts below the survey level refer to child age-band records. External data: MAP \textit{Pf}PR$_{2-10}$ (population-weighted with a gridded population surface), WHO/UNICEF vaccine coverage, national child (0--14 years) HIV incidence (UNAIDS/UNICEF; derived from UNAIDS new-infection counts for Liberia) and World Bank series for GDP, health expenditure and political stability.}
```

### S2 (L276–278): PfPR spline figure (outdated; typo)
- **Problem:** `sfig_pfpr_mortality_by_age.png` is the v5 export. The caption contains "log hazard ration relative for" and calls the point estimates "Mean values".
- **Evidence:** Overleaf md5 `eea85640…` = v5; the v7 `pfpr_splines.png` is `9157a282…`.
```latex
\caption{Malaria prevalence and all-cause age-band-specific child mortality. The y-axis shows the log mortality hazard ratio relative to a \textit{Pf}PR$_{2-10}$ of 20\%. Point estimates (lines) and conditional pointwise 95\% intervals (shaded areas) are shown for each of the seven separate age-band models.}
```

### S3 (L284–286): subgroup figure (outdated; the v3 image is ready)
- **Problem:** The image is the v2 export (v5 base), and the caption is only "Sensitivity analyses: subgroups".
- **Evidence:** Overleaf md5 `889021f8…` = v2; the v3 `sfig_pfpr_splines_by_subgroup.png` is `d275d5db…`.
- **Action:** copy `results/cbh/subgroups_dhsmics_map_gamma2_v3/sfig_pfpr_splines_by_subgroup.png`, then use:
```latex
\caption{Sensitivity of the fitted \textit{Pf}PR$_{2-10}$--mortality relationship to the analysis subset: log mortality hazard ratios relative to 20\% prevalence from the full primary sample (black, with pointwise 95\% interval) and from refits in the Sahel (31 surveys), Eastern Africa (51 surveys), and surveys up to (72) and after (63) the median survey year, 2014. Subsets overlap, so differences are descriptive.}
```

### fig:MI and fig:without_malnutrition (L291–300): both ready
The v3 no-nutrition figure (`nutrition_adjustment_dhsmics_map_gamma2_v3/sfig_pfpr_splines_without_nutrition.png`) is exported to `Supplementary Figures/` by `export_paper.py`. The v8 imputed-covariate figure (`primary_map_regional17_dhsmics_imputed_gamma2_v8/sfig_pfpr_splines_imputed_covariates.png`) is exported the same way; both replace the existing images under the same file names.

### Table S1 (L303–325): primary age-band results (outdated; new from the synthesis spot-check)
- **Problem:** The table holds an old sample. For example, its total row reads 5,465,305 records and 75,726 deaths, and the <1 month row reads 966,553 records and 27,788 deaths. v7 has 7,607,122 records and 103,987 deaths, and the HRs differ: for example, 6–11 months 20→0 is 0.69 in v7, against 0.75 in the table. The table is never cited in the text.
- **Evidence:** v7 `tables/age_band_results.latex.txt`.
- **Action:** replace L303–325 with the v7 fragment. Its rows:
```latex
$<1$ & 1,341,316 & 38,245 (36.8\%) & 2.33 & 0.94 (0.92--0.97) & 0.94 (0.90--0.99) \\
1--5 & 1,206,366 & 15,873 (15.3\%) & 2.28 & 0.92 (0.88--0.95) & 0.86 (0.80--0.93) \\
6--11 & 1,160,005 & 14,384 (13.8\%) & 2.97 & 0.83 (0.79--0.88) & 0.69 (0.63--0.76) \\
12--23 & 1,003,279 & 13,389 (12.9\%) & 3.51 & 0.83 (0.78--0.88) & 0.53 (0.47--0.60) \\
24--35 & 994,653 & 11,488 (11.0\%) & 3.66 & 0.86 (0.81--0.92) & 0.47 (0.41--0.54) \\
36--47 & 963,051 & 6,732 (6.5\%) & 3.33 & 0.85 (0.79--0.91) & 0.55 (0.48--0.64) \\
48--59 & 938,452 & 3,876 (3.7\%) & 3.20 & 0.92 (0.84--1.01) & 0.57 (0.49--0.67) \\
\midrule
Total & 7,607,122 & 103,987 (100.0\%) & -- & -- & -- \\
```
Also add the fragment's footnote minipage.

---

## Unverifiable literature claims (not errors; check the citations)

1. **L50:** how IHME and WHO/MCEE use verbal autopsy in geospatial cause-of-death models, and VA's sensitivity, specificity, single-cause assumption and inability to capture indirect effects. These need citations in the Introduction.
2. **L86:** neonatal passive immunity and foetal haemoglobin (`WHO:14`); causes of neonatal death (`who_newborn_mortality_2024`); malaria in primigravidae and low birth weight (`Walker:14`); disease concentrated in the first years of life (`Kamau:22`); and the interpretation that surviving 3–4-year-olds have disease-controlling immunity.
3. **L205:** that Nigeria is the only country with GBD subnational estimates. This is probably wrong, and a correction is given above.
4. **L231:** "Following previous analyses of DHS/MIS birth histories \cite{Wakefield:19,Burstein:18}, each child's trajectory was divided into seven age bands". Wakefield et al. (2019) appear to use six bands (0, 1–11, 12–23, 24–35, 36–47, 48–59 months), and Burstein et al. (2018) is a summary-birth-history method, so neither used this seven-band scheme. Check this against the papers. The sample also includes MICS. A possible rewording:
```latex
Adapting the discrete-time age-band approach of previous analyses of complete birth histories \cite{Wakefield:19,Burstein:18}, which used six bands (0, 1--11, 12--23, 24--35, 36--47 and 48--59 months), we divided each child's trajectory into the seven age bands defined above, splitting the 1--11 month band at six months.
```

## New bibliography keys needed (do not add to `ref.bib` automatically)

These keys are cited in the suggestions above: `mics`, `unaids2025`, `aidsinfo`, `cacode2025` (doi 10.1136/bmj-2025-088686; take the authors and title from the DOI), `wdi`, `wgi` and, optionally, `wuenic`. The skeletons below use the source URLs; check them before pasting.
```bibtex
@misc{mics, author = {{UNICEF}}, title = {Multiple Indicator Cluster Surveys (MICS)}, howpublished = {\url{https://mics.unicef.org}}}
@misc{unaids2025, author = {{UNAIDS}}, title = {UNAIDS 2025 HIV estimates, distributed by UNICEF}, howpublished = {\url{https://data.unicef.org/topic/hivaids/}}, year = {2025}}
@misc{aidsinfo, author = {{UNAIDS}}, title = {AIDSinfo: UNAIDS 2026 HIV estimates}, howpublished = {\url{https://aidsinfo.unaids.org}}, year = {2026}}
@article{cacode2025, doi = {10.1136/bmj-2025-088686}, note = {Complete author, title and volume from the DOI}}
@misc{wdi, author = {{World Bank}}, title = {World Development Indicators}, howpublished = {\url{https://databank.worldbank.org/source/world-development-indicators}}}
@misc{wgi, author = {{World Bank}}, title = {Worldwide Governance Indicators}, howpublished = {\url{https://www.worldbank.org/en/publication/worldwide-governance-indicators}}}
@misc{wuenic, author = {{WHO and UNICEF}}, title = {WHO/UNICEF estimates of national immunization coverage (WUENIC)}, howpublished = {\url{https://immunizationdata.who.int/}}}
```

## Pipeline notes found during the audit (not manuscript text)
- The v7 `REPORT.md` said the Liberia HIV validation covers 741 country-years; `liberia_aidsinfo/REPORT.md` and `definition_check_by_country.csv` say 761. Fixed: `05_compare_regional.R` now reads the count from the check file.
- `figures/fig3_burden_comparison_caption.md` still holds v5 state totals, and `figures/fig4_burden_comparison_caption.md` is a stale duplicate.

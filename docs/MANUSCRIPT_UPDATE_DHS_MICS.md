# Manuscript changes for the DHS and MICS primary analysis

Written 23 September 2026 for `main.tex` and `SText1_model_specification.tex` in the Overleaf project
`malaria-burden-reassessment`. No `.tex` or `.bib` file has been edited: each item below gives the `main.tex` line number (as
of the 22 September version, 333 lines), a short quote of the current text and replacement text to paste. Every number
was checked against the results files listed at the end. The export copies this file to Overleaf as
`MANUSCRIPT_TEXT_UPDATES.md`.

> **Superseded for numbers (24 September 2026).** The primary is now `primary_map_regional17_dhsmics_gamma2_v7`: this
> file's v5 sample plus Liberia's 2007, 2013 and 2019–20 DHS, whose child HIV incidence is now derived from UNAIDS counts
> (135 surveys: 98 DHS, 37 MICS; 37 countries; 7,607,122 child age-band records; 103,987 deaths). The sensitivities are
> refitted on v7 (subgroups and no-nutrition v3, imputed covariates v8). Every number below is v5's and is kept as the
> v5 record; for v7 values and replacement text use `docs/MANUSCRIPT_AUDIT_2026-09-24_v7.md` in the analysis repository
> ([link](MANUSCRIPT_AUDIT_2026-09-24_v7.md)), which audits the current `main.tex` against v7.
> Items superseded outright are marked below (Line 225, Liberia; Line 247, sensitivity counts), and the paper now has a
> fourth main figure (see "New Figure 4").

The imputed-covariate sensitivity (`primary_map_regional17_dhsmics_imputed_gamma2_v6`) and its multiple-imputation check
are complete. The wording below that depends on them is taken from its `REPORT.md`, `CAPTION_sfig_imputed_covariates.md`
and `multiple_imputation/REPORT.md`.

## What changed in the sample

The primary analysis (`primary_map_regional17_dhsmics_gamma2_v5`) now pools DHS surveys and UNICEF Multiple Indicator
Cluster Surveys (MICS) with complete birth histories. The model specification is unchanged (seven age-band cloglog
models, 17 covariates, gamma = 2, fixed reference knots); covariates are re-scaled on the larger sample.

| Sample | DHS only (v3) | DHS and MICS (v5) | Added by MICS |
|---|---:|---:|---:|
| Child age-band records | 5,465,305 | 7,498,459 | 2,033,154 (27.1%) |
| Children | 1,686,004 | 2,292,089 | 606,085 |
| Deaths | 75,726 | 102,282 | 26,556 (26.0%) |
| Survey-regions | 916 | 1,211 | 295 |
| Surveys | 95 | 132 (95 DHS, 37 MICS) | 37 |
| Countries | 34 | 36 | 2 (Central African Republic, Guinea-Bissau) |

- **Registry.** 124 DHS/MIS surveys and 93 MICS surveys in 43 countries. 46 MICS surveys have no complete birth
  history. Of the 47 that do (rounds 3 to 6, fieldwork 2006-2022, 28 countries), Lesotho 2018 is excluded with the four
  Lesotho DHS surveys (malaria free, no MAP surface). That leaves 166 processed surveys (120 DHS/MIS and 46 MICS) in 40
  countries.
- **Complete-case exclusions.** 34 processed surveys contribute no records: 25 DHS/MIS surveys (20 DHS and all 5 MIS,
  as before) and 9 MICS surveys. No MIS survey is in the analysed sample.
- **DHS part unchanged.** The 95 DHS surveys reproduce the DHS-only sample exactly.
- **Effects.** The 40%→20% hazard ratios are 0.944, 0.915, 0.830, 0.823, 0.850, 0.843 and 0.907 across the seven
  bands. They are within 0.026 of the DHS-only fit; the largest change is at 24-35 months.
- **Burden.** Attributable deaths across the 42 countries are 1,165,686 (2000), 716,520 (2015) and 611,440 (2024),
  7.7% and 8.1% above the DHS-only values for 2015 and 2024. The rate falls 67.5% from 2000 to 2024 (DHS only: 68.3%).
- **Sensitivity analyses.** The subgroup and no-nutrition analyses are refitted on the combined sample (v2). The
  imputed-covariate version (v6) covers 8,797,963 records from 166 surveys. Its 40%→20% hazard ratios are within 0.040
  of v5, and pooling across the ten imputed datasets changes them by at most 0.002.

## Before pasting

Run `python3 R_cbh/reporting/export_paper.py --destination <Overleaf folder> --copy`; since 24 September 2026 it exports
from v7 by default and stops if the v8 imputed-covariate figure or caption, or any of the four main figures, is missing
(the v5 export ran on 24 September 2026). It also copies Figure 4 as `figures/fig4_under5_death_probability.png`. It refreshes every results figure `main.tex` includes (all except `figures/Causal_diagram.pdf`),
including `figures/fig3_burden_comparison.png` for Figure 3, and writes the table sources `tables/age_band_results.latex.txt`,
`tables/country_comparison_2024.latex.txt` and `tables/source_comparison.latex.txt`. Until the export runs, the images
in Overleaf are the DHS-only versions.

---

## Abstract

### Line 51 (sample and attributable fractions)

Current: "We pooled all available Demographic and Health Surveys with complete birth histories since 2000 (120 surveys, 36
countries) ... in data from 1,686,004 children under 5 and 75,726 deaths. ... At 20\% prevalence, the model attributes 6\%
of neonatal deaths, 25\% of deaths at 6-11 months and 42-47\% of deaths at 1-4 years to malaria."

Replace the whole line with:

```latex
    We pooled all available Demographic and Health Surveys (DHS) and UNICEF Multiple Indicator Cluster Surveys (MICS) with complete birth histories since 2000 (166 surveys, 40 countries) and matched individual child trajectories, discretised into 7 age bands until 5 years, with the admin level-1 (subnational first level administrative boundaries) population weight prevalence estimates from the Malaria Atlas Project. After adjusting for a series of potential confounding variables, a robust age-dependent association between \textit{Pf}PR$_{2-10}$ and all-cause mortality was observed in data from 2,292,089 children under 5 and 102,282 deaths (132 surveys, 36 countries). \textit{Pf}PR$_{2-10}$ only weakly predicted neonatal death, but was a strong predictor of death from six months and older. At 20\% prevalence, the model attributes 5\% of neonatal deaths, 30\% of deaths at 6-11 months and 42-53\% of deaths at 1-4 years to malaria.
```

Values: attributable fraction at 20% = 5.3% (<1 month), 30.4% (6-11 months), 45.6%, 52.7%, 45.7% and 42.1% (12-59
months).

### Line 52 (burden rates, declines, comparison)

Current: "fell from 10.1 to 4.2 deaths per 1,000 child-years between 2000 and 2015 (a 58\% decline), and then to 3.2
deaths per 1000 in 2024 (24\% decline). ... (XX and XX\% decrease, respectively) ... about a third higher than either
source".

Replace the first three sentences of the line with the following, keeping the final sentence ("Recent progress ...") as it is:

```latex
    Applied to national all-cause mortality in 42 sub-Saharan African countries, we estimate that under-5 malaria mortality fell from 10.6 to 4.5 deaths per 1,000 child-years between 2000 and 2015 (a 57\% decline), and then to 3.4 deaths per 1,000 in 2024 (a further 24\% decline). In contrast, IHME and UN IGME estimate substantial declines up until 2015 (49 and 60\% decrease, respectively) but little change since (8 and 6\% decrease). We estimate malaria mortality in 2024 was about 40\% higher than either source in the same set of countries (43\% above IHME and 39\% above UN IGME).
```

Values (deaths per 1,000 child-years, common IHME-implied denominator): model 10.61 (2000), 4.52 (2015), 3.45 (2024);
IHME 5.13, 2.62, 2.41; UN IGME 6.59, 2.63, 2.48. The IHME and UN IGME values do not depend on MICS; they fill the XX
placeholders. The 2024 ratios are 1.43 (model/IHME) and 1.39 (model/UN IGME), for both deaths and rates.

## Introduction

### Line 70

Current: "...is extremely robust based on 25 years of DHS surveys in Africa."

```latex
... is extremely robust based on 25 years of DHS and MICS surveys in Africa.
```

## Results

### Line 74 (subsection heading)

Current: `\subsection*{Prevalence and all cause mortality in DHS/MIS surveys}`

```latex
\subsection*{Prevalence and all-cause mortality in DHS and MICS surveys}
```

### Line 76 (survey counts, sample, deaths by age)

Current: "...Demographic and Health Surveys (DHS) and Malaria Indicator Surveys (MIS)..."; "120 surveys across 36
countries, representing 1,114 survey regions"; "95 surveys in 34 countries"; "All five MIS surveys and 20 DHS surveys were
excluded ..., so the analysed sample is DHS only"; "1,686,004 children with 75,726 deaths"; "65\% in the first year of
life, and 14\% between 1 and 2 years". The paragraph also has the typo "compete birth histories".

Replace the whole paragraph with:

```latex
We pooled all available microdata from Demographic and Health Surveys (DHS), Malaria Indicator Surveys (MIS) and UNICEF Multiple Indicator Cluster Surveys (MICS) conducted in malaria endemic countries since the year 2000 which recorded complete birth histories. Excluding Lesotho, which is malaria free and has no MAP prevalence estimates, this gave a total of 166 surveys (120 DHS/MIS and 46 MICS) across 40 countries, representing 1,457 survey regions at the admin level-1 (subnational first level administrative boundaries), Figure \ref{fig:map}. For each birth recorded, survival and death events were included if they occurred within the previous 5 years. Each child trajectory from birth was divided into discrete age bands of $<$1 month (neonatal), 1-5 months, 6-11 months, and then yearly, corresponding to how data are recorded in the complete birth histories of the DHS, MIS and MICS. Each recorded survival and death event was then merged with the estimated admin level-1 \textit{Pf}PR$_{2-10}$ prevalence estimates from the Malaria Atlas Project (MAP) at the time of start of the age band \cite{map2025}. Following the merge with MAP prevalence estimates and restricting the analysis to complete cases with identified confounders recorded, a total of 132 surveys (95 DHS and 37 MICS) in 36 countries were included in the primary analysis (Figure \ref{fig:flow}). Thirty-four surveys were excluded because a required covariate was unavailable for every region of the survey: 25 DHS/MIS surveys, including all five MIS surveys, and nine MICS surveys. These surveys contributed data from a total of 2,292,089 children with 102,282 deaths recorded in 1,211 survey regions. 37\% of these deaths were in neonates (first month of life); 66\% in the first year of life, and 13\% between 1 and 2 years of age.
```

Values: 36.8% neonatal; 67,335 of 102,282 deaths (65.8%) in the first year; 13,138 (12.8%, shown as 12.9% after
largest-remainder rounding in Table S1) at 12-23 months. 1,457 is the number of survey-regions with MAP-eligible records
in the 166 processed surveys.

### Line 81 (Figure 1 caption)

Current: "DHS and MIS surveys in malaria endemic (or historically endemic) countries of sub-Saharan Africa conducted since
2000 with complete birth histories recorded. Countries are labelled using their ISO3 codes."

The exported map now shows MICS surveys as triangles, excluded surveys as crosses and Lesotho as diamonds.

```latex
    \caption{Geographic coverage and timing of the 124 DHS and MIS surveys and 47 MICS surveys with complete birth histories in malaria endemic (or historically endemic) countries of sub-Saharan Africa conducted since 2000 (41 countries). Filled circles mark the 95 DHS surveys and filled triangles the 37 MICS surveys that contribute to the primary analysis (36 countries); point area indicates the number of contributing survey regions, and country shading the number of included surveys. Crosses mark the 34 surveys excluded because a required covariate was unavailable for every region (20 DHS, all 5 MIS and 9 MICS surveys) and open diamonds the 5 surveys in Lesotho, which is malaria free and has no MAP prevalence estimates. The 46 MICS surveys without a complete birth history are not shown. Countries are labelled using their ISO3 codes.}
```

### Line 85 (shape of the curves)

Current: "However, risk of death in older children plateaued above prevalences of around 30\% and even decreased for
children over 3 years of age for prevalences above 30\% and 40\%."

The combined-sample curves no longer decline materially. At 36-47 months the log hazard ratio (vs 20%) peaks at 0.179
at 48.5% prevalence and is 0.170 at 63%; at 48-59 months it is 0.097 at 40% and flat thereafter (0.095 at 63%). At
12-35 months the curves keep rising, more slowly. Replace the sentence with:

```latex
However, the increase in risk of death in older children slowed above prevalences of around 30\%, and for children over 3 years of age the curves were flat above prevalences of about 40-50\%.
```

The following sentence ("This pattern fits our understanding ...") still reads correctly. "(37\% of included deaths)"
earlier in the paragraph is unchanged (36.8%).

### Line 91 (attributable fractions)

Current: "...compared with less than 5\% in neonates (\ref{fig:main_result2}). ... peaking at 56\% at 2–3 years."

Only the peak changes (and the missing "Figure~"):

```latex
Under this model, for an average \textit{Pf}PR$_{2-10}$ of 10\%, approximately 30\% of all deaths in children 1 year and older are caused by malaria, compared with less than 5\% in neonates (Figure~\ref{fig:main_result2}). For an average \textit{Pf}PR$_{2-10}$ of 40\% (upper end of prevalence estimates in sub-Saharan Africa now), there is a sharper age-dependent increase in malaria deaths with more than half of all deaths between 1 and 4 years attributable to malaria, peaking at 60\% at 2-3 years.
```

Values: at 10%, 2.5% (<1 month) and 29.1%, 35.2%, 29.4%, 27.4% (12-59 months); at 40%, 55.3%, 59.8%, 54.2% and 47.5%
(12-59 months).

### Line 95 (Figure 2 caption)

Current: "Mean model estimates of the malaria attributable share of all-cause mortality as a function of age-band at
different prevalences of infection (\textit{PfPR}$_{2-10}$)."

The exported figure now draws 95% intervals.

```latex
\caption{Estimated malaria-attributable share of all-cause mortality within each age band at \textit{Pf}PR$_{2-10}$ of 10\%, 20\%, 30\% and 40\%, from the seven age-band models fitted to the DHS and MICS sample. Points are point estimates of $1-\exp\{f_g(0)-f_g(P)\}$ and vertical bars conditional pointwise 95\% intervals; lines connect the discrete age bands. Zero prevalence lies below the observed exposure range, so the counterfactual involves extrapolation.}
```

### Line 100 (sensitivity analyses in Results)

Current: "...steeper for children between 1 and 5 years of age (Figure \ref{fig:subgroups}). A sensitivity analysis using
multiply imputed datasets gave very similar results to the primary analysis restricted to complete cases (Figure
\ref{fig:MI}). ... essentially unchanged (Figure \ref{fig:without_malnutrition})."

On the combined sample the Sahel curves are steeper from 6 months onward (40%→20% hazard ratios at 6-59 months 0.70-0.78
in the Sahel, 0.82-0.91 in the full sample); at 1-5 months they match (0.91 vs 0.92). Eastern Africa departs from the
full sample in the other direction at two ages: 1.00 (0.96-1.04) against 0.94 (0.92-0.97) at <1 month and 0.95
(0.86-1.06) against 0.84 (0.78-0.91) at 36-47 months.

```latex
We conducted a series of sensitivity analyses looking at subgroups of the surveys (Sahel versus Eastern Africa, earlier versus later surveys). In general all the subgroup fits gave very similar results, apart from the subgroup of survey-regions in seasonal areas (Sahel), where the estimated relationship between prevalence and all-cause mortality was steeper for children between 6 months and 5 years of age (hazard ratios for a reduction in \textit{Pf}PR$_{2-10}$ from 40\% to 20\% of 0.70-0.78, against 0.82-0.91 in the full sample), and Eastern Africa, where the curves were flatter in neonates and at 36-47 months (hazard ratios 1.00 and 0.95, against 0.94 and 0.84; Figure \ref{fig:subgroups}). A sensitivity analysis that imputed the remaining covariate gaps, retaining all 8,797,963 MAP-eligible records from 166 surveys, gave similar results to the primary analysis restricted to complete cases (40\% to 20\% hazard ratios within 0.04, the largest difference at 48-59 months; pooling across the ten imputed datasets changed them by at most 0.002; Figure \ref{fig:MI}). Because malnutrition may mediate the relationship between malaria transmission and death as well as confound it, we refitted the models on the same sample without adjusting for wasting and stunting; the estimated prevalence-mortality relationship was essentially unchanged (largest change in the 40\% to 20\% hazard ratio 0.008; Figure \ref{fig:without_malnutrition}).
```

Values: the v6 point fits give 40%→20% hazard ratios of 0.940, 0.913, 0.825, 0.843, 0.864, 0.858 and 0.947, against
0.944, 0.915, 0.830, 0.823, 0.850, 0.843 and 0.907 in v5 (largest difference 0.040, at 48-59 months). Pooling the ten
imputed-dataset refits with Rubin's rules changes the 40%→20% hazard ratios by at most 0.0016 (0.002 across both
contrasts); the between-imputation share of variance is at most 2.0% for this contrast and 2.3% across both (v6
`REPORT.md` and `multiple_imputation/pooled_contrasts.csv`). The original wording "very similar" is kept for the
subgroups; "similar" is used for v6 because of the 0.04 difference at 48-59 months.

### Line 105 (Nigeria and Chad)

Current: "(\ref{tab:country}). In Nigeria, IHME estimate approximately 140 thousand deaths, and IGME estimate approximately
165 thousand deaths whereas our approach estimates 200 thousand deaths (46 and 22\% increase respectively). In Chad we
estimate a 3.6 fold higher total burden."

```latex
Although the country specific under 5 mortality rates mostly agreed between IHME and our new approach (PfPR-ACM), there were substantial differences in several populous and high burden countries (Table~\ref{tab:country}). In Nigeria, IHME estimate approximately 138 thousand deaths, and UN IGME estimate approximately 165 thousand deaths, whereas our approach estimates 215 thousand deaths (56 and 30\% higher, respectively). In Chad we estimate a 4.0-fold higher total burden than IHME.
```

Values: Nigeria 214,938 (model), 138,111 (IHME), 165,329 (UN IGME); Chad 19,407 against IHME 4,884 (ratio 3.97).

### Lines 120-131 (Table 1, `tab:source-comparison`)

Only the PfPR-ACM column changes; IHME and UN IGME are identical to the current table. Replace lines 120-124 and
127-131 with:

```latex
\quad 2000 & 1,165,686 & 563,657 & 724,276  \\
\quad 2015 & 716,520 & 416,122 & 416,815  \\
\quad 2024 & 611,440 & 428,147 & 440,123  \\
\quad Change 2000--2024 & $-47.5\%$ & $-24.0\%$ & $-39.2\%$  \\
\quad Change 2015--2024 & $-14.7\%$ & $+2.9\%$ & $+5.6\%$ \\
```

```latex
\quad 2000 & 10.61 & 5.13 & 6.59  \\
\quad 2015 & 4.52 & 2.62 & 2.63  \\
\quad 2024 & 3.45 & 2.41 & 2.48  \\
\quad Change 2000--2024 & $-67.5\%$ & $-53.0\%$ & $-62.4\%$  \\
\quad Change 2015--2024 & $-23.8\%$ & $-8.1\%$ & $-5.7\%$  \\
```

Old PfPR-ACM values: 1,104,990 / 665,559 / 565,616; −48.8% / −15.0%; 10.06 / 4.20 / 3.19; −68.3% / −24.1%.
Alternatively paste the whole table from `tables/source_comparison.latex.txt`. Its footnote adds that the sources share
one denominator and that under-five person-time in these countries grew 12% between 2015 and 2024 (158,554,858 to
177,450,738 person-years), which is why counts and rates move differently after 2015.

### Line 137 (Nigerian states)

Current: "...greater than IHME's estimates in 35 of the 37 states (Figure~\ref{fig:IHME_v_model}B). ... the median state
ratio is 1.6 in the three northern zones and 1.2 in the three southern zones."

```latex
Within Nigeria, applying the same age-band effects to state-level IHME all-cause deaths by age and state \textit{Pf}PR$_{2-10}$ gives malaria mortality estimates greater than IHME's estimates in 36 of the 37 states, all except Lagos (Figure~\ref{fig:IHME_v_model}B). Summed over states, the model gives 214,127 deaths against 132,138 IHME malaria deaths (ratio 1.62). The differences are concentrated in the northern states: the median state ratio is 1.7 in the three northern zones and 1.2 in the three southern zones.
```

Values: medians 1.73 (north) and 1.24 (south); Lagos ratio 0.70.

### Line 141 (Figure 3 image)

Include `figures/fig3_burden_comparison.png` (updated in Overleaf on 24 September 2026); the export no longer writes
`figures/fig4_burden_comparison.png`, and the old file and `figures/fig4_burden_comparison_caption.md` can be deleted.

### New Figure 4 (probability of dying before age 5), added 24 September 2026

The export writes `figures/fig4_under5_death_probability.png` (source
`results/cbh/primary_map_regional17_dhsmics_gamma2_v7/under5_probability/under5_death_probability_2024.png`): the
probability of dying before age 5 in 2024 by country, per 1,000 live births, from all causes (synthetic cohort from IHME
age-band all-cause death rates) and the part caused by malaria under the PfPR-ACM model. If `main.tex` has no figure
environment for it, add one; the draft caption is `under5_probability/CAPTION.md`, also in the paper-figure `CAPTIONS.md`.

### Line 142 (Figure 3 caption)

Current: "Comparison of national and subnational (Nigeria only) estimates of malaria mortality rates in children under 5
between the prevalence-ACM model and IHME. ... C: average predicted mortality in children under 5 in sub-Saharan Africa
since 2004."

Panel C also shows UN IGME and covers 2000-2024; panels A and B compare the model with IHME only.

```latex
\caption{Malaria mortality before age 5 from the PfPR-ACM model compared with IHME and UN IGME (CA-CODE 2026) cause-specific estimates, expressed as deaths per 1,000 under-five child-years. A: model against IHME national rates in 2024 in the 42 countries with both estimates; the dashed line denotes equality. B: model against IHME for the 36 Nigerian states and the Federal Capital Territory in 2024, coloured by geopolitical zone; summed over states the model gives 214,127 deaths against 132,138 IHME malaria deaths (ratio 1.62). C: annual under-five malaria mortality, 2000-2024, pooled across the same 42 countries. In every panel each source's deaths are divided by the under-five person-years implied by the IHME all-cause death counts and rates, so the sources share a denominator. Model deaths are IHME all-cause deaths in each age band multiplied by $1-\exp\{f_g(0)-f_g(P)\}$. Model-attributable reductions in all-cause mortality and cause-specific malaria deaths are different estimands, and the model shares the IHME all-cause inputs.}
```

### Lines 148-165 (Table 2, `tab:country`)

The top ten changes: DR Congo enters and Uganda drops out. Replace lines 154-165 with:

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

Or paste `tables/country_comparison_2024.latex.txt`, which adds a footnote (UN IGME is the CA-CODE 2026 series; totals
cover all 42 countries). Rows remain ordered by the absolute model-IHME difference (Nigeria 76,827 to Central African
Republic 5,273).

## Discussion

### Line 188 (limitations)

Current: "...because DHS data are the most reliable for children under 5. ... (one third of all deaths) ... the 89\% of
all reported deaths occurred before age five."

```latex
... because birth-history data from DHS and MICS surveys are the most reliable for children under 5.
```

"(one third of all deaths)" can become "(37\% of deaths in the analysed sample)". The 89% figure comes from
`R_cbh/audit/17_recorded_death_ages.R` and covers the DHS/MIS birth histories (513,952 of 578,905 recorded deaths,
88.8%, in the 120 processed DHS/MIS surveys). `R_mics/01_inventory.R` tallies MICS with the same definition (deaths with
age at death under 60 months over all recorded deaths): the 46 built MICS surveys give 146,187 of 164,217 (89.0%;
`results/mics_inventory/survey_inventory.csv`, `deaths_under5`/`deaths`, excluding Lesotho 2018). The combined share is
88.8% (660,139 of 743,122), so "89\%" can stay, e.g. "in the DHS, MIS and MICS birth histories, 89\% of all reported
deaths occurred before age five".

## Methods

### Lines 197-199 (data sources, first sentence)

Current: "Complete birth histories were obtained from Demographic and Health Surveys and Malaria Indicator Surveys (DHS/MIS)
Births Recodes \cite{dhs}."

```latex
Complete birth histories were obtained from Demographic and Health Surveys and Malaria Indicator Surveys (DHS/MIS) Births Recodes \cite{dhs} and from UNICEF Multiple Indicator Cluster Surveys (MICS) \cite{mics}.
```

The rest of line 199-200 ("Birth dates, survival status ... 4-5 years.") is unchanged.

### Line 201 (DHS covariates)

Current opening: "Child, maternal and household characteristics obtained from the DHS/MIS surveys comprised from the Births
Recodes, ..."

```latex
For DHS/MIS surveys, child, maternal and household characteristics comprised, from the Births Recodes, ...
```

The remainder of the line is unchanged.

### New paragraph after line 201 (MICS)

Checked against `R_mics/README.md`, `results/mics_inventory/ELIGIBILITY.md`, `docs/ANALYSIS_PLAN.md` Sections 2.1–2.4,
`R_mics/04_boundaries_gnb_caf.R`, `R_mics/09_regional_covariates.R`, `data/derived_mics/region_map.csv` and the v5 study
flow. Insert:

```latex
Of the 93 MICS surveys from sub-Saharan Africa that we assessed, 47 (rounds 3 to 6, fieldwork 2006-2022, 28 countries) included a complete birth history equivalent to the DHS Births Recode; the other 46 recorded only summary birth histories (children ever born and children who died) and cannot support the age-band design. Each MICS birth history was converted to the DHS Births Recode layout (survival status, multiple births, age at death and sampling weights) and processed with the same child age-band construction and eligibility rules as the DHS/MIS data. The 366 MICS region labels were mapped to analysis regions defined as unions of boundary polygons, reusing the DHS boundary file for the same country where the regions matched, by exact name (291 labels) or a reviewed correspondence table (75 labels). Guinea-Bissau and the Central African Republic have no DHS survey; their analysis regions were built from an admin-1 boundary file [REF] in which each capital is merged with a neighbouring region (Bissau with Biombo; Bangui with Ombella-M'Poko), giving 8 and 6 analysis regions respectively, and the unnamed MICS regions of the Central African Republic were assigned prefectures following the country's seven health regions. The Dakar city survey was assigned the combined Dakar and Thiès polygon of the same file. Somalia and South Sudan, which also have no DHS boundary file, were mapped to the regions of the same file (the three zones of the 2006 Somalia survey as approximate unions of regions); none of their four surveys entered the primary analysis. Annual regional MAP \textit{Pf}PR$_{2-10}$ was extracted for these analysis regions by the same method as for DHS regions.

MICS does not publish survey-region indicator tables, so all 13 regional covariates were computed from the MICS microdata with the DHS definitions: maternal age at first birth, education, wealth quintile and urban residence (among women with a birth in the 60 months before interview), facility delivery and short birth intervals using the women's sampling weights; improved water, improved sanitation and electricity using household weights; and DTP3 and measles coverage among children aged 12-23 months (card or mother's report) and wasting and stunting among children aged 0-59 months with valid WHO Child Growth Standards z-scores, using children's weights. Education, recorded as level and grade, was converted to completed years using each country's primary and secondary school durations, with women who never attended school counted as zero years; facility delivery refers to the last birth in the two years before interview, not all births in the previous five years; water and sanitation follow the 2017 WHO/UNICEF Joint Monitoring Programme classification; and the Dakar city survey, which has no urban variable, was coded as fully urban. The same national WUENIC and within-survey fallbacks as for DHS were applied (below); in addition, recall vaccination doses were unusable in Guinea 2016, Comoros 2022 and Chad 2019, so the national WUENIC DTP3 and measles estimates for the survey year were used for all regions of these surveys.

The complete-case rules were the same as for DHS/MIS, so a MICS survey entered the primary analysis only if anthropometry and the national covariate series were available. Thirty-seven MICS surveys (fieldwork 2006-2022, 24 countries) entered, contributing 2,033,154 child age-band records and 26,556 deaths in 295 survey regions. They include the only surveys from the Central African Republic and Guinea-Bissau, and six sub-national samples (five in Kenya and the Dakar city survey) whose wealth quintiles are ranked within their own sample rather than nationally. Nine MICS surveys contributed no records: Nigeria 2021 and southern Madagascar 2012 (no anthropometry); São Tomé and Príncipe 2014 and 2019 (no child HIV incidence series); Somalia 2006, the two Somali zone surveys of 2011 (which also lack anthropometry) and Zimbabwe 2009, whose band-entry years precede the national health-expenditure series; and South Sudan 2010, whose band-entry years precede its political-stability and health-expenditure series. Lesotho 2018 was excluded, as were the Lesotho DHS surveys, because Lesotho is malaria free and has no MAP prevalence estimates.
```

Notes for the author:
- `\cite{mics}` needs a new `ref.bib` entry (see References). `[REF]` is the source of the admin-1 file
  (`Africa_New_Admin.shp` in the Snow prevalence project), which has no citation in the repository. The same file
  supplies the Guinea-Bissau, Central African Republic, Dakar, Somalia and South Sudan regions
  (`data/derived_mics/region_map.csv`: 18, 7, 1, 11 and 10 labels). The four Somalia and South Sudan surveys are in the
  imputed-covariate sensitivity sample (line 247) but not the primary.
- The national series start in 2013 (Somalia health expenditure), 2010 (Zimbabwe health expenditure), and 2011 and 2017
  (South Sudan political stability and health expenditure).
- If the Methods need to be shorter, the second paragraph can end after the first sentence plus the WUENIC sentence,
  with the definition details moved to the Supplement.

### Line 217 (LLIN coverage)

Current: "(measured in DHS/MIS surveys)"

```latex
(measured in DHS, MIS and MICS surveys)
```

### Line 222 (weights)

Current: "Recode-derived regional summaries used DHS women's sampling weights. Published regional coverage indicators
retained their stated denominators and reference periods."

Append:

```latex
MICS regional summaries were computed from the microdata with the MICS women's, household and children's sampling weights, as described above.
```

### Line 225 (fallbacks; Liberia and São Tomé)

> **Superseded for Liberia (24 September 2026).** From v7 Liberia's 2007, 2013 and 2019–20 DHS are in the sample, with
> child HIV incidence derived from UNAIDS AIDSinfo counts; only its 2009 MIS is excluded (no anthropometry). The
> replacement below is wrong for Liberia; use the v7 audit's L207–208 and L220 items. São Tomé and Príncipe remains
> excluded for want of an incidence series.

Current: "Missing regional DTP3 and measles coverage was replaced by the corresponding national UNICEF estimate for the
survey year (\url{https://immunizationdata.who.int/}). ... Liberia and Sao Tome had no data and were excluded."

After the first sentence insert:

```latex
This national substitution was also applied to all regions of three MICS surveys (Guinea 2016, Comoros 2022 and Chad 2019) whose recall vaccination doses were unusable.
```

Replace the last sentence with:

```latex
Liberia and São Tomé and Príncipe have no UNAIDS child or adolescent HIV incidence series, so their surveys (four DHS/MIS surveys in Liberia; one DHS and two MICS surveys in São Tomé and Príncipe; 165,149 child age-band records with MAP prevalence) were excluded from the primary analysis.
```

The 165,149 records with MAP prevalence are in LB51FL, LB5AFL (MIS), LB6AFL, LB7AFL, ST51FL, MC_STP2014 and MC_STP2019
(`pfpr_available_rows` in `study_flow/survey_selection.csv`; 165,253 eligible records before the MAP merge). The same
total is the "HIV incidence" exclusion in the study-flow figure, which counts after the MAP merge.

### Line 247 (sensitivity analyses)

> **Superseded counts (24 September 2026).** The subgroup refits on v7 (`subgroups_dhsmics_map_gamma2_v3`) have 72 early
> and 63 late surveys (Sahel and Eastern Africa unchanged), and the imputed-covariate version v8 takes Liberia's HIV rate
> from the UNAIDS counts; see the v7 audit.

Current: "...Sahelian survey regions (...; 23 surveys, 7 countries, 14,545 deaths), Eastern Africa (...: 38 surveys, 12
countries, 23,583 deaths), and surveys before versus after the median survey year of 2013 (50 and 45 surveys,
respectively). ... This retained all 6,357,802 eligible records (120 surveys, 36 countries, 90,938 deaths); imputation
uncertainty was propagated by refitting with each of the ten imputed datasets and pooling with Rubin's rules."

Replace the whole line with:

```latex
We refitted the 7 age-band models in four subsets of the primary sample: Sahelian survey regions (boundary centroid at or north of 12°N and west of 36°E, excluding the Horn of Africa; 31 surveys, 8 countries, 18,919 deaths), Eastern Africa (UN M49 delimitation of Eastern Africa: 51 surveys, 12 countries, 31,589 deaths), and surveys conducted up to and including versus after the median survey year of 2014 (70 surveys and 55,805 deaths; 62 surveys and 46,477 deaths). Subset fits used the same specification and covariate scaling as the primary analysis, with knots at quantiles of the subset's own predictor values. In a second sensitivity analysis we imputed every remaining covariate gap rather than excluding records: whole-survey gaps in regional indicators by chained-equation multiple imputation across the combined DHS and MICS survey regions (predictive mean matching, 10 imputations); the missing 2001 political-stability round by interpolation; South Sudan's political stability before independence (2000-2010) and its GDP outside the observed 2008-2015 series from generalised additive models; health expenditure for Zimbabwe 2000-2009, Somalia 2000-2012, South Sudan outside 2017-2023 and all countries in 2024 from a generalised additive model on the observed panel; and child HIV incidence for Liberia and São Tomé and Príncipe from the incidence model with a latent adolescent series. This retained all 8,797,963 MAP-eligible records (2,656,472 children, 123,419 deaths, 1,457 survey regions, 166 surveys: 120 DHS/MIS and 46 MICS; 40 countries). Imputation uncertainty was propagated by refitting with each of the ten imputed datasets, with smoothing parameters fixed at those of the main imputed-covariate fit, and pooling with Rubin's rules.
```

Notes:
- Subgroup counts are from `results/cbh/subgroups_dhsmics_map_gamma2_v2` (Sahel 1,538,691 records and 239 survey-regions;
  Eastern Africa 2,680,402 and 447; early 3,399,915 and 549; late 4,098,544 and 662). The subset refits use the
  full-sample covariate scaling.
- South Sudan's pre-independence political stability and GDP are GAM fills, a modelling choice. The alternatives are
  whole-Sudan values, carrying the 2011 values back, or omitting MC_SSD2010 (57,509 records, 0.65% of the imputed
  sample).
- The Rubin's-rules sentence matches v6 `multiple_imputation/REPORT.md`: each of the ten datasets varies the regional
  imputations, the GAM draws of the national fills and one posterior draw of child HIV incidence for every imputed or
  censored country-year (1,283,638 records), with smoothing parameters fixed at the point fit.

## Data sharing

### Line 257

Current: "All input data are publicly available via the GBD client portal and the DHS website."

```latex
All input data are available from their sources: DHS and MICS microdata on registration from the DHS Program (\url{https://dhsprogram.com}) and UNICEF (\url{https://mics.unicef.org}); IHME estimates from the GBD results tool; UN IGME (CA-CODE) cause-specific estimates from the UN IGME portal (\url{https://childmortality.org}); MAP prevalence surfaces from the Malaria Atlas Project; and UNAIDS, WHO/UNICEF (WUENIC) and World Bank series from their publishers.
```

The code sentence ("Analysis code and derived data tables are available at ...") is unchanged. The CA-CODE series was
retrieved through UNICEF's `CME_CAUSE_OF_DEATH` dataflow (`R_cbh/reporting/07_annual_mortality_comparison.R`). The
admin-1 boundary file used for the Guinea-Bissau, Central African Republic, Dakar, Somalia and South Sudan regions should
be added to this list once its source is confirmed (see `[REF]` in the Methods).

## Supplementary figures and table

### Line 276 (Figure S1, study flow)

Current: "Data sources and analysis flow. External non-DHS data include MAP ..., child HIV prevalence (UNAIDS 0-14 case
counts divided by the World Bank population estimates) and World Bank economic series for GDP and health expenditure."

This has no MICS counts and describes HIV prevalence, not the incidence covariate used (Methods line 212).

```latex
\caption{Data sources and sample inclusion. The survey registry contains 124 DHS/MIS and 93 MICS surveys in 43 countries; 46 MICS surveys without a complete birth history and 5 surveys in Lesotho (malaria free; no MAP prevalence estimates) were excluded, leaving 166 surveys in 40 countries. Of 10,657,524 child age-band entries within 60 months of interview, 1,347,434 had a band not complete by interview and 418,661 an entry year outside 2000-2024, leaving 8,891,429 eligible records; 93,466 lacked a regional MAP prevalence value and 1,299,504 lacked at least one required covariate, leaving 7,498,459 records (2,292,089 children; 102,282 deaths) from 132 surveys (95 DHS and 37 MICS) in 36 countries. External (non-survey) data comprise MAP \textit{Pf}PR$_{2-10}$ (population-weighted with a gridded population surface), WHO/UNICEF national vaccine coverage, UNAIDS child HIV incidence (ages 0-14, per 1,000 uninfected) and World Bank GDP, health expenditure and political stability series.}
```

Covariate exclusions, attributed to the first missing covariate: wasting 578,080; health expenditure 260,385; political
stability 181,520; HIV incidence 165,149; electricity 83,796; wealth score 30,574. The full caption is exported as
`Supplementary Figures/sfig_study_flow_caption.md`.

### Line 282 (Figure S2, PfPR and mortality by age)

No number in this caption depends on the sample, but it has a typo ("log hazard ration relative for") and the exported
figure is now the combined-sample fit. Suggested:

```latex
\caption{Malaria prevalence and all-cause age band specific child mortality in the primary sample (7,498,459 child age-band records, 102,282 deaths, 132 DHS and MICS surveys). The y-axis shows the log hazard ratio of age-specific mortality relative to a prevalence of 20\%. Mean values (pointwise 95\% confidence intervals) are shown by the blue lines (shaded areas), over the central 95\% of each age band's exposure distribution.}
```

### Line 290 (Figure S3, subgroups)

Current: "Sensitivity analyses: subgroups"

```latex
    \caption{Sensitivity of the \textit{Pf}PR$_{2-10}$-mortality relationship to the analysis subset. Log mortality hazard ratios relative to \textit{Pf}PR$_{2-10}$ = 20\% by age band from the full primary sample (black, with its pointwise 95\% interval shaded) and from separate refits in Sahelian survey regions (31 surveys, 8 countries, 18,919 deaths), Eastern Africa (51 surveys, 12 countries, 31,589 deaths), and surveys conducted up to (70 surveys, 55,805 deaths) and after (62 surveys, 46,477 deaths) the median survey year of 2014. Each curve is drawn over the central 95\% of its own exposure distribution. The subsets overlap the full sample and one another, so their differences are descriptive.}
```

### Line 297 (Figure S4, imputed covariates)

Current: "Results from multiply imputed datasets."

From `CAPTION_sfig_imputed_covariates.md` in the v6 report (exported as
`Supplementary Figures/sfig_pfpr_splines_imputed_covariates_caption.md`), in LaTeX. Until the export runs, the image in
Overleaf is the DHS-only v4 figure.

```latex
    \caption{Sensitivity of the \textit{Pf}PR$_{2-10}$-mortality relationship to the treatment of missing covariates. Log mortality hazard ratios relative to \textit{Pf}PR$_{2-10}$ = 20\% by completed-month age band, with pointwise 95\% conditional intervals shaded, from the PfPR-ACM model fitted to the complete-case primary sample (7,498,459 child age-band records, 102,282 deaths, 132 surveys in 36 countries) and to every MAP-eligible record after imputing all remaining covariate gaps (8,797,963 records, 123,419 deaths, 166 surveys in 40 countries, including all five Malaria Indicator Surveys and all 46 MICS surveys with complete birth histories outside Lesotho). Whole-survey gaps in the regional indicators (wasting, stunting, facility delivery, electricity and the wealth score) were imputed by chained equations at the survey-region level (predictive mean matching, 10 imputations, point value their mean); the missing 2001 World Governance Indicators round by linear interpolation; health expenditure for Zimbabwe 2000-2009, Somalia 2000-2012, South Sudan before 2017 and all countries in 2024, and South Sudan's GDP before 2008 and political stability before independence, from generalised additive models on the observed national panel; and child HIV incidence for Liberia and São Tomé and Príncipe from the incidence model with a latent adolescent series. Observed values were never replaced. The specification, reference knots, basis dimensions and $\gamma=2$ are those of the primary analysis; covariates were re-standardised on the enlarged sample. Each curve is drawn over the central 95\% of its own exposure distribution. Hazard ratios for \textit{Pf}PR$_{2-10}$ 40\% to 20\% are 0.94, 0.92, 0.83, 0.82, 0.85, 0.84, 0.91 (complete case) and 0.94, 0.91, 0.82, 0.84, 0.86, 0.86, 0.95 (imputed) for the seven bands from youngest to oldest. Refitting the imputed version in each of the 10 imputed datasets and pooling with Rubin's rules changed these hazard ratios by at most 0.002; the between-imputation share of interval variance was at most 2.3\%. The two samples are nested, so the comparison is descriptive.}
```

### Line 304 (Figure S5, without wasting and stunting)

Current: "Results when removing malnutrition indicators from the list of covariates included for adjustment."

```latex
    \caption{Sensitivity of the \textit{Pf}PR$_{2-10}$-mortality relationship to removing regional wasting and stunting prevalence from the adjustment set. Log mortality hazard ratios relative to \textit{Pf}PR$_{2-10}$ = 20\% by age band from the primary model (17 covariates) and from a refit on the identical sample (7,498,459 records, 102,282 deaths) with 15 covariates. The largest change in the hazard ratio for a reduction in \textit{Pf}PR$_{2-10}$ from 40\% to 20\% is 0.008, at 1-5 months.}
```

### New Figure S6 (death-count version of Figure 3), optional

The export writes `Supplementary Figures/sfig_burden_comparison_counts.png`, but `main.tex` has no figure environment
for it. If it belongs in the supplement, add after the Figure S5 environment (line 306):

```latex
\begin{figure}
    \centering
    \includegraphics[width=0.95\linewidth]{Supplementary Figures/sfig_burden_comparison_counts.png}
    \caption{Death-count version of Figure~\ref{fig:IHME_v_model}: malaria-attributable deaths before age 5 from the PfPR-ACM model compared with IHME and UN IGME cause-specific malaria deaths. A: model against IHME national deaths in 2024 in the 42 countries with both estimates, on base-10 logarithmic axes from 100 to 250,000 deaths (6 countries with fewer than 100 deaths on either axis lie outside the displayed range). B: model against IHME for the 36 Nigerian states and the Federal Capital Territory in 2024, axes from 200 to 30,000 deaths. C: annual under-five malaria mortality, 2000-2024, per 100,000 child-years, pooled across the same 42 countries on the common IHME-implied person-year denominator. All values are point estimates.}
    \label{fig:burden_counts}
\end{figure}
```

The generated caption file (`sfig_burden_comparison_counts_caption.md`) still calls the main figure "Figure 4".

### Lines 319-327 (Table S1, `tab:primary-age-band-effects`)

Replace lines 319-327 with:

```latex
$<1$ & 1,322,149 & 37,623 (36.8\%) & 2.34 & 0.94 (0.92-0.97) & 0.95 (0.90-1.00) \\
1-5 & 1,189,209 & 15,586 (15.2\%) & 2.25 & 0.92 (0.88-0.95) & 0.86 (0.81-0.93) \\
6-11 & 1,143,915 & 14,126 (13.8\%) & 2.90 & 0.83 (0.79-0.87) & 0.70 (0.63-0.77) \\
12-23 & 989,211 & 13,138 (12.9\%) & 3.46 & 0.82 (0.77-0.88) & 0.54 (0.48-0.62) \\
24-35 & 980,344 & 11,338 (11.1\%) & 3.63 & 0.85 (0.80-0.91) & 0.47 (0.41-0.54) \\
36-47 & 949,032 & 6,650 (6.5\%) & 3.34 & 0.84 (0.78-0.91) & 0.54 (0.47-0.63) \\
48-59 & 924,599 & 3,821 (3.7\%) & 3.16 & 0.91 (0.83-0.99) & 0.58 (0.49-0.68) \\
\midrule
Total & 7,498,459 & 102,282 (100.0\%) & - & - & - \\
```

Or paste `tables/age_band_results.latex.txt` (en-dash ranges and a footnote on rounding and conditional intervals). The
caption (line 310) can add "(132 DHS and MICS surveys)" after "primary analysis dataset".

## SText1 (`SText1_model_specification.tex`)

The only MICS-specific wording is "across DHS/MIS survey-regions" (line 32), but the whole document describes the
archived Method-2 negative-binomial model (600 survey-regions, +8.5% per 10 points, 1% counterfactual). That is not the
model `main.tex` line 238 refers to ("Full model equations are given in the Supplement (SText~1)"). None of it can be
updated with v5 numbers. Suggested replacement for the title (lines 21-22) and the Overview (lines 29-54), using the
file's own `\pfpr` macro:

```latex
\title{\textbf{Supplementary Text 1.\ Specification of the
prevalence--mortality (PfPR-ACM) model}}
```

```latex
\section*{Overview}
The PfPR-ACM model relates annual survey-region \pfpr\ at entry into each age band to all-cause mortality in seven child age bands ($<$1, 1--5, 6--11, 12--23, 24--35, 36--47 and 48--59 completed months), using complete birth histories from Demographic and Health Surveys and UNICEF Multiple Indicator Cluster Surveys. The primary sample contains 7,498,459 child--age-band records from 2,292,089 children, with 102,282 deaths, in 1,211 survey-regions, 132 surveys (95 DHS and 37 MICS) and 36 countries.

For child $i$ in band $g$ of width $\Delta_g$ years, $D_{ig}\sim\mathrm{Bernoulli}(q_{ig})$ with $q_{ig}=1-\exp(-\Delta_g\lambda_{ig})$ and
\[
\log\lambda_{ig}=\alpha_g+f_g(P_{r(i),y_{ig}})+h_g(t_{ig})+\mathbf{X}_{r(i),s(i),y_{ig}}^{\top}\boldsymbol\beta_g+u_{s(i),g}+v_{c(i),g}+b_{r(i),g},
\]
equivalently $\mathrm{cloglog}(q_{ig})=\log\Delta_g+\log\lambda_{ig}$. Here $P$ is regional MAP \pfpr\ in the calendar year of band entry $y$; $t$ is fractional calendar year at entry; $f_g$ and $h_g$ are penalised cubic regression splines (basis dimensions 5 and 6, knots fixed at those of the reference fit); $\mathbf X$ holds the 13 survey-region and 4 national annual covariates, centred and scaled on the primary sample; and $u$, $v$ and $b$ are independent Gaussian random intercepts for survey, country and survey-region. Each band is fitted separately with an unweighted likelihood using \texttt{mgcv::bam} (fast REML, discretised covariates, $\gamma=2$), so every band has its own intercept, curves, coefficients and variance parameters. Hazard ratios between prevalences $X$ and $Y$ are $\exp\{f_g(Y)-f_g(X)\}$, and the malaria-attributable fraction at prevalence $P$ is $1-\exp\{f_g(0)-f_g(P)\}$.
```

Then delete the Model comparison (lines 56-100), Choice of prevalence floor (lines 102-111) and Conclusion (lines
113-136) sections, including `figures/sfig1_spec_comparison.png`, or move them to an appendix clearly labelled as the
archived earlier model. The formula matches `results/cbh/primary_map_regional17_dhsmics_gamma2_v5/model_formula.txt` and
`docs/ANALYSIS_PLAN.md` sections 3.1-3.2.

## References (`ref.bib`)

`ref.bib` has a DHS entry (`dhs`, the Guide to DHS statistics) but nothing for MICS. These entries are suggestions
written from memory: the author should check authors, volume, pages and DOI before pasting them.

```bibtex
@article{mics,
  author  = {Khan, Shane and Hancioglu, Attila},
  title   = {Multiple Indicator Cluster Surveys: delivering robust data on children and women across the globe},
  journal = {Studies in Family Planning},
  year    = {2019},
  volume  = {50},
  number  = {3},
  pages   = {279--286},
  doi     = {10.1111/sifp.12103}
}

@misc{mics_data,
  author       = {{UNICEF}},
  title        = {Multiple Indicator Cluster Surveys ({MICS})},
  howpublished = {\url{https://mics.unicef.org/surveys}},
  note         = {Microdata available on registration}
}

@techreport{jmp2017,
  author      = {{WHO and UNICEF}},
  title       = {Progress on drinking water, sanitation and hygiene: 2017 update and {SDG} baselines},
  institution = {World Health Organization and United Nations Children's Fund},
  address     = {Geneva},
  year        = {2017}
}
```

- `mics` is the key used in the Methods text above. `mics_data` is optional, for the Data sharing statement.
- `jmp2017` is optional, if the JMP classification is cited in the MICS covariate paragraph.
- A citation is still needed for the admin-1 boundary file used for Guinea-Bissau, the Central African Republic,
  Dakar, Somalia and South Sudan (`[REF]` in the Methods). Its original publisher has not been confirmed.

## Not caused by MICS, but in passages touched above

- **Line 174:** "approximately a XX percent decrease" is still a placeholder. For reference, WHO World Malaria Report
  2025 all-age deaths in the African Region fell 28% from 2000 to 2024 (804,000 to 579,000). Under-five malaria deaths
  in the 42 countries fell 24% (IHME) and 39% (UN IGME).
- **Line 242:** "Estimates were produced annually for 2004-2024" should read 2000-2024. Table 1 and Figure 3C start in
  2000, and all 42 countries have inputs for all 25 years. Optionally add the Nigeria reconciliation: state sum 214,127
  against the national estimate 214,938 (−0.4%).
- **Line 246:** typo "paramerisation".
- **Line 91, 95, 188:** `\textit{PfPR}$_{2-10}$` differs from the manuscript's usual `\textit{Pf}PR$_{2-10}$`; the
  replacements above use the latter.

## Results files used

All paths are in the analysis repository. These are the v5 sources; the v7 equivalents are in
`results/cbh/primary_map_regional17_dhsmics_gamma2_v7/`, `subgroups_dhsmics_map_gamma2_v3/`,
`nutrition_adjustment_dhsmics_map_gamma2_v3/` and `primary_map_regional17_dhsmics_imputed_gamma2_v8/`.

- `results/cbh/primary_map_regional17_dhsmics_gamma2_v5/`:
  - `prepared_sample.csv`, `REPORT.md`, `RESULTS.md`
  - `study_flow/` (`CAPTION.md`, `flow_counts.csv`, `figure_labels.csv`, `survey_selection.csv`)
  - `survey_map/` (`CAPTION.md`, `survey_timeline_all.csv`)
  - `attributable_fraction_by_age.csv`, `pfpr_curves.csv`
  - `tables/` (`age_band_results.latex.txt`, `country_comparison_2024.csv`, `country_comparison_2024.latex.txt`)
  - `burden_comparison/` (`CAPTION.md`, `CAPTION_counts.md`, `source_comparison.latex.txt`,
    `source_comparison_2000_2024.csv`)
  - `annual_comparison/annual_totals_2000_2024.csv`
  - `nigeria_states/` (`README.md`, `state_totals_2024.csv`, `national_reconciliation_2024.csv`)
- `results/cbh/subgroups_dhsmics_map_gamma2_v2/` (`sample_summary.csv`, `subgroup_definitions.csv`, `REPORT.md`,
  `CAPTION.md`)
- `results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v2/` (`CAPTION.md`, `CONTRASTS.md`)
- `results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/` (`prepared_sample.csv`, `imputation_record_counts.csv`,
  `imputation_record_counts_by_survey.csv`, `prepared_selection_by_survey.csv`, `pfpr_curves.csv`, `REPORT.md`,
  `CAPTION_sfig_imputed_covariates.md`, `comparison_contrasts.csv`, `multiple_imputation/REPORT.md`,
  `multiple_imputation/pooled_contrasts.csv`)
- `results/cbh/covariate_imputation_dhsmics_v6/` (`NATIONAL_IMPUTATION.md`, `REGIONAL_IMPUTATION.md`,
  `somalia_south_sudan_zimbabwe_series.csv`)
- `results/mics_inventory/` (`ELIGIBILITY.md`, `survey_inventory.csv` (also the MICS recorded-death ages),
  `exclusion_attribution_mics_by_survey.csv`)
- `results/cbh/primary_map_regional17_gamma2_v3/study_flow/recorded_death_ages_by_survey.csv` (the DHS/MIS 88.8% figure)
- `R_mics/README.md`, `R_mics/04_boundaries_gnb_caf.R`, `R_mics/09_regional_covariates.R`,
  `data/derived_mics/region_map.csv` (read only), `docs/ANALYSIS_PLAN.md` Sections 2.1–2.4 and 3.1

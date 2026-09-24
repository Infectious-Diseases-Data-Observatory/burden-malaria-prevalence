# Adversarial review of main.tex (Methods and Results) against the v7 pipeline

24 September 2026. Scope: Results (main.tex L74–182), Methods (L204–261), the supplementary items they cite (L279–345) and `SText1_model_specification.tex`, checked against the current primary `results/cbh/primary_map_regional17_dhsmics_gamma2_v7/` and the sensitivity folders refitted on it. The review is written as a hostile Lancet/NEJM referee report.

Eight reviewers worked independently on different aspects: causal identification, exposure, birth-history outcome data, statistics, burden estimation, text-versus-code consistency, missing analyses, and reporting. Each reviewer's findings then went to a separate skeptic who tried to refute them. An editor-level critic read everything that survived, looked for gaps and wrote the recommendation. No finding was refuted outright. Several were downgraded, and those are listed in their own section below so the case is not overstated.

The numbers marked ✔ were recomputed independently for this report from saved v7 outputs. The others were reproduced by the skeptic pass. main.tex was not edited.

## Verdict

**Reject in its current form.** The assembled data are valuable and the question matters. As written, however, the paper's two headline claims do not come from what the data identify:
- malaria deaths 45% above IHME and 41% above UN IGME;
- progress since 2015 about three times what is reported.

The level depends on extrapolating below the observed exposure and on specification choices made after seeing results. The trend is mostly IHME's all-cause decline. None of the headline numbers has an uncertainty interval, and the supplement describes a different model. A serious revision could turn this into a publishable paper; the requirements are in the last section.

## The strongest argument against the paper

The method defines malaria deaths as IHME all-cause deaths × [1 − HR(PfPR = 0 vs current)], where the hazard ratio comes from a cross-sectional, between-area association.

**Both headline results are driven by parts of this formula that the data do not inform.**

**Trend.**
- With PfPR frozen at each country's 2015 value, the modelled malaria mortality rate still falls 21.2 of the 23.6 percentage points it falls between 2015 and 2024 ✔. The pooled attributable fraction barely moves, from 0.248 to 0.235 ✔, while person-year-weighted national PfPR goes from 20.0% to 18.6% ✔.
- The "greater recent progress" is therefore IHME's falling all-cause envelope passed through a near-constant fraction. The model's own malaria-specific signal after 2015 is small and consistent with the stagnation it is contrasted with.

**Level.**
- The 45% excess over IHME comes mainly from the 0–5% PfPR segment of a single cubic spline piece. That piece runs from the lowest knot (0.016%) to the first interior knot (about 19%), and its bottom lies below the 2.5th percentile of exposure (about 1.7%).
- Moving the reference prevalence from 0% to 1%, 2% or 5% gives 589k, 558k and 466k deaths in 2024 instead of 620k ✔. At 5% the total is 1.09 × IHME.
- The alternative exposure surface the authors fitted and then dropped (Snow) gave 439k against 592k with MAP on the same earlier sample ✔.

**Causality.**
- The association comes mainly from between-country and between-region contrasts, with random (not fixed) intercepts and regional covariates measured at the survey date, not at band entry.
- The one falsification test offered, neonates as a "negative control", is non-null: HR for PfPR 0 vs 20% is 0.94 (0.90–0.99) ✔. It also adds 49,853 attributed neonatal deaths to the 2024 total ✔.
- An unreported Sahel seasonality analysis shows no malaria-season excess in high-PfPR regions (August–October rate ratio 1.02 at 12–23 months where PfPR ≥ 25%) ✔. There, the model attributes about half of deaths to malaria.

A reviewer can therefore argue that the paper cannot distinguish "malaria kills more children than estimated" from "low-transmission areas have lower child mortality for other reasons, and IHME's envelope is falling".

## Major flaws, ranked

### 1. The burden level rests on extrapolation below the data
- **Where:** L91, L104–107, L140, L254; Abstract L52.
- **Evidence:**
  - All seven bands flag `zero_below_observed_support`.
  - About 24–31% of each band's 0→20% log hazard ratio accrues between 0 and 5% PfPR.
  - The 0–5% segment accounts for about 154,500 deaths, roughly 80% of the 191,925-death excess over IHME.
  - About 18% of records lie below 5% PfPR, about 73% of them from five countries (Senegal, Kenya, Ethiopia, Mauritania, Rwanda).
  - Reference-prevalence sensitivity for 2024: 0% → 620k (1.45 × IHME), 1% → 589k (1.38), 2% → 558k (1.30), 5% → 466k (1.09) ✔.
  - The subgroup curves disagree most in exactly this segment. The attributable fraction at 12–23 months and 5% PfPR is 0.60 in Eastern Africa and 0.09 in the Sahel.
  - SText1, cited as the model specification, argues for a 1% floor and a curve flat below 10%.
- **Authors' best defence:** true zero transmission exists in the data only in Lesotho, which is excluded. The alternative is to change the estimand to "deaths attributable to transmission above X%".
- **Request:** report the burden at 0/1/2/5% references; a shape-constrained or monotone sensitivity; leave-one-country-out refits of the low segment; the Snow comparison.

### 2. The post-2015 trend claim is almost entirely arithmetic from the IHME envelope
- **Where:** Abstract L52, L106, Table `tab:source-comparison`, Discussion L185–189.
- **Evidence:**
  - Freezing PfPR at its 2015 level reproduces a −21.2% fall against the reported −23.6% ✔.
  - Holding all-cause deaths fixed, the PfPR change alone lowers deaths by only about 3%.
  - For 2000–2015, 49.0 of the 57.1 percentage points are reproduced with the fractions frozen.
  - The assumption that the PfPR curve is constant over 25 years is contradicted by the saved period subgroups. The 20%→0% HR at 36–47 months is 0.50 in early surveys and 0.36 in late ones.
- **Endpoint:**
  - 2024 is an IHME forecast year: the median relative width of the under-5 death uncertainty interval jumps from 0.06 to 0.33, and IHME malaria deaths rise 7% in one year.
  - With 2023 as the endpoint the contrast survives (model −24.9% vs IHME −13.5%, 2015–2023 ✔), but it is smaller.
  - For 2015–2019 the model and IHME fall by similar amounts (−13.2% vs −10.1% ✔). The divergence opens only after 2019, when IHME shows a COVID-era rise.
- **Request:**
  - Publish the decomposition (PfPR held fixed vs all-cause held fixed).
  - Test the multiplicative assumption with PfPR × period and PfPR × baseline-mortality interactions.
  - Report 2015–2019 and 2019–2023 separately, and use 2023 as the primary endpoint.
  - Otherwise, withdraw the "substantially greater progress" sentence.

### 3. The causal interpretation is not identified, and the paper's own falsification check fails
- **Where:** L88, L104, L199, L228–230, L250.
- **Evidence:**
  - **Where the variation is.** Of the variance in PfPR across the 9,688 survey-region × entry-year cells of the included surveys, 59.7% lies between countries, 30.0% between regions within a country and only 10.3% within a region over time.
  - **Random effects.** Country and region enter only as random intercepts, and gamma = 2 shrinks the region random effect almost to nothing (EDF 0.08 of 1,213 at 1–5 months and 0.0014 at 48–59 months).
  - **Adjustment.** Crude adjustment already removes 42–64% of the log association. For example, the crude RR at 12–23 months is 3.30 against an adjusted 1.94. The 13 regional covariates are joined by survey and region, so they describe conditions at the interview date, not at band entry.
  - **The authors' own documents say this is unresolved.** ANALYSIS_PLAN §3.3 says "the existing coding does not itself identify a within-region causal effect". OPEN_ANALYSIS_ISSUES.md says the 17 covariates are "not a validated sufficient adjustment set".
  - **Neonatal "negative control".** It is non-null, measured worse than later outcomes, and has a plausible causal route through malaria in pregnancy. So it cannot establish "no major residual confounding" (L88). It also contributes 8% of the 2024 total ✔.
  - **Covariate behaviour.** Several confounder coefficients behave non-causally. Wasting is protective in all seven bands (HR 0.92–0.95), and DTP3 coverage is associated with higher mortality at 6–11 months (HR 1.10).
- **Authors' best defence:**
  - The country random effects are only lightly shrunk (EDF 14–27 of 36), so the estimate is partly a within-country one.
  - The model *under*-predicts the ITN and SMC trial effects, which argues against strong positive confounding (see "weakened" below).
  - Measured survey prevalence gives a steeper slope than MAP.
- **Request:**
  - A within/between (Mundlak) decomposition.
  - Country and survey fixed effects.
  - Harmonised-region fixed effects with country-specific trends, i.e. a difference-in-differences on the 2000–2015 regional PfPR declines.
  - E-values or a quantitative bias analysis.
  - A future-PfPR placebo.

### 4. No uncertainty on any headline number
- **Where:** L52, L106, L109–112, Tables `tab:country` and `tab:source-comparison`, L250, L254.
- **Evidence:**
  - National totals, the ratios to IHME and UN IGME, "36 of 42", ρ = 0.91 and the trend contrasts are all point estimates (L254 says so).
  - The likelihood is unweighted and ignores clustering. Sampling weights, PSU, stratum and mother IDs are built in the child-band data but never used. There are 2.07 children per mother. The 1,227 survey-regions collapse to about 503 distinct country-region names.
  - MAP exposure uncertainty is not propagated. Using the 2024 MAP lower and upper rasters, which are in the repository but unused, gives a 237k–806k outer bound.
  - IHME uncertainty, HIV-imputation uncertainty and smoothing-parameter uncertainty are also not propagated. L250 uses the conditional covariance although the smoothing-corrected one is available.
  - Figure 4's intervals assume independent bands, which L254 explicitly rejects.
- **Request:** a survey- or PSU-cluster bootstrap with common replicates across all seven bands, with MAP draws and IHME uncertainty. It should give 95% intervals for national and country totals, the ratios to IHME/UN IGME and the 2000–2015 and 2015–2024 changes. Also state that the likelihood is unweighted.

### 5. The specification was chosen after seeing results, and the alternatives are unreported
- **Where:** L246–250 ("Statistical model and model selection"), L257–259.
- **Evidence — the 2024 total under successive specifications, against IHME's 428,147:**

  | Specification | 2024 deaths | × IHME |
  |---|---:|---:|
  | Joint seven-band model | 433,510 | 1.01 |
  | Separate models, gamma = 1 | 545,150 | 1.27 |
  | Separate models, gamma = 2 | 592,353 | 1.38 |
  | Snow exposure instead of MAP | 439,132 ✔ | |
  | 18-variable adjustment | 501,065 | |
  | 17-variable adjustment | 565,616 | |
  | Adding MICS (v5) | 611,440 | |
  | Adding Liberia (v7) | 620,072 | 1.45 |

- **Further evidence:**
  - The plan records these as decisions taken after reviewing results.
  - The "model selection" subsection describes no selection procedure.
  - L250 says gamma = 2 was used "to force smoother curves". In a like-for-like comparison, gamma = 2 left the PfPR effective degrees of freedom unchanged, collapsed the region random effect and raised the burden by 8.7%.
  - The spline knots are copied from an archived Snow-comparison fit and not disclosed; the subgroup fits use quantile knots instead.
- **Authors' best defence:** the path is not monotone (18 variables lowered the estimate), and MICS and Liberia are data additions, not modelling choices.
- **Request:** a pre-specified specification curve on the v7 sample, reported in full:
  - MAP vs Snow exposure;
  - gamma 1 vs 2;
  - the 11-, 17- and 18-covariate sets;
  - joint vs separate models;
  - quantile knots and k = 8;
  - a monotone spline.

  Disclose the decision history.

### 6. Heterogeneity is described as "very similar", and transport to other settings is untested
- **Where:** L99, L104 ("generalisable across geographies"), L254.
- **Evidence:**
  - The paper judges similarity on the 40%→20% contrast, but the burden uses 20%→0%.
  - On 20%→0% at 12–23 months: Eastern Africa 0.30, early surveys 0.47, pooled 0.53, late surveys 0.59, Sahel 0.70. On that contrast the Sahel is *flatter*, not steeper as L99 says.
  - Applying the subgroup curves to the 2024 inputs gives 602k–719k. That total is reassuring, but it isn't reported.
  - The curves are applied to all 42 countries, including those with no surveys in the sample.
- **All-cause envelope:**
  - Country results depend on the IHME all-cause envelope, which differs from UN IGME by a factor of 0.72–1.81 across countries (Nigeria 5q0 87.5 vs 115.6).
  - On the UN IGME envelope, DRC rises from about 70k to about 128k, while the 42-country total changes by only about 5%.
  - The UN IGME comparator is built on its own envelope, so the country-level comparison mixes two envelopes.
- **Request:** report 20%→0% contrasts and subgroup-curve burdens; test PfPR × region and PfPR × period interactions; leave-one-country-out refits (Nigeria supplies 16% of fitting deaths and 35% of the 2024 burden); country-held-out prediction; a burden version on the UN IGME envelope.

### 7. The estimand and the language overreach
- **Where:** L50, L71, L104–112, L178.
- **Evidence:**
  - **The estimand.** It is all-cause mortality attributable to transmission above zero. It includes indirect effects and residual confounding, and it is compared with single-cause estimates while being called "malaria deaths".
  - **The correlation.** ρ = 0.91 with IHME is largely built in. National MAP PfPR alone correlates 0.88 with the IHME malaria rate, and IHME's own African estimates use MAP.
  - **"Directly measurable".** The abstract calls PfPR₂₋₁₀ a "directly measurable quantity". It is actually a MAP model output fitted partly to the same surveys, and the 2023–24 layers go beyond the cited methods paper, which covers 2000–22. Agreement between MAP and measured prevalence falls from r = 0.89 overall to 0.71 in 2019–2024, with MAP 6.8 points higher.
  - **Age distribution.** 66% of the 2024 attributable deaths fall at 12–59 months, where IHME's age distribution of deaths is most model-dependent.
- **Request:**
  - Relabel the estimand, for example "all-cause under-5 deaths attributable to *P. falciparum* transmission".
  - Drop ρ as evidence of agreement, or say it follows by construction.
  - Describe the exposure as modelled, and document the MAP release used for each year.

### 8. The Methods and supplement do not describe the analysis that produced the Results
- **SText1:** cited at L250 for the "full model equations", but it describes an archived negative-binomial region-level model. That model had 600 regions, 23 countries, AIC selection among four models and a 1% counterfactual.
- **L259:** still gives the DHS-only numbers:
  - Sahel 23 surveys / 14,545 deaths; Eastern Africa 38 / 23,583;
  - median year 2013 (50/45 surveys);
  - imputed sample 6,357,802 records / 120 surveys;
  - Liberia's HIV rate "from the latent adolescent model".

  L99 has the correct v7 values, and v8 takes Liberia's rate from UNAIDS counts.
- **L258** names analyses that were never run ("West vs East Africa", "parameterisation"). The no-nutrition refit, cited at L100, is missing from the Methods.
- **Table S1 (L320–343):** is the DHS-only v3 table (5,465,305 records), is never cited, and its HRs differ from v7. For example, 6–11 months 20%→0% is 0.75 in the table and 0.69 in v7.
- **L234:** says "Liberia and São Tomé had no data and were excluded", but Liberia's three DHS surveys are in the sample. Its HIV derivation exists only in the commented-out L236.
- **Figure S1 caption:** describes the archived HIV *prevalence* covariate. Nigeria's and Comoros's child HIV series are entirely model-imputed, which the paper does not say.
- **L254:** says estimates cover "2004–2024"; the tables and figures run 2000–2024.
- **L112:** the Nigeria state results say 35 of 37 states with medians 1.6/1.2; the v7 state README says 36 of 37 (Lagos is the only state below IHME) with medians 1.76/1.27 ✔. The README also notes that the model-to-IHME ratio has Spearman correlation 0.10 with state PfPR, so the north–south gap reflects IHME's state age distributions more than transmission.
- **Figure 4:** `fig:cumulative_risk` is never cited, and neither the Methods nor the Results describe the figure.
- **Figure 1 caption:** still says the rows are "labelled using their ISO3 codes", and repeats the grouping sentence.
- **Unweighted likelihood:** never stated.
- **Bibliography:** `\cite{mics}`, `\cite{wdi,wgi}` have no ref.bib entries.
- **Typos:** L209 repeats "maternal age at first live birth"; L99 has "multipl".
- **Front and back matter:** the contributors section is a placeholder and the author list contains "[...]". There is no funding, role-of-funder or ethics statement. The first author's GiveWell affiliation sits alongside "no competing interests". The data statement omits MICS, MAP, GPW, UNAIDS/AIDSinfo, WUENIC, WDI/WGI and CA-CODE, and implies DHS data are freely redistributable.

## Findings that weakened under scrutiny (do not overstate these)

- **ITN trial calibration** (the Discussion cites the ITN trials as support):
  - For an ITN-sized 17% relative fall in prevalence, the fitted curves predict only a 1–6% fall in all-cause mortality. The trials observed 17–20%.
  - This cuts against the paper's use of trials as support. However, it points to under- rather than over-estimation, and PfPR saturates as a short-term marker. Use it as "the curve is not validated against trials", not as proof of confounding.
- **MAP circularity:** measured survey prevalence gave a *steeper* mortality slope than MAP, so shared data probably attenuate the association rather than inflate it. The real issues are the "directly measurable" wording and exposure uncertainty.
- **2024 forecast endpoint:** the qualitative conclusions survive with a 2023 endpoint (1.51 × IHME; −24.9% vs −13.5% for 2015–2023 ✔), but the "8% vs 24%" framing is endpoint-sensitive.
- **Birth-history quality:** age heaping at 12 months, recall and omission for births 3–9 years earlier, and DHS–MICS comparability are plausible biases. They are unquantified, not demonstrated.
- **High-HIV southern African anchor and survivor bias:** a real selection concern at the PfPR = 0 end, but minor next to flaw 1.
- **UN IGME envelope:** the 42-country total is robust (about +5%); individual countries are not.
- **GiveWell affiliation:** a disclosure issue, not evidence of bias.

## Missing analyses, ranked

| Priority | Analysis | What it tests | Feasibility with data in the repo |
|---|---|---|---|
| Essential | Reference-prevalence sensitivity (0/1/2/5%), with records and deaths per PfPR bin below 10% | Dependence of the excess on extrapolation | Immediate from saved curves (values above) |
| Essential | Trend decomposition, plus PfPR × period and PfPR × baseline-mortality models | Whether "greater progress" reflects malaria | Decomposition immediate; interactions need 7 refits each |
| Essential | Within-area identification: Mundlak, country FE, harmonised-region FE with country trends / difference-in-differences | Between-area confounding | Feasible; needs stable region identifiers across survey rounds |
| Essential | Joint uncertainty: survey/PSU bootstrap across all bands, plus MAP draws and IHME uncertainty | Whether the excess and trend gap exceed uncertainty | About 17 min per 7-band fit, so ~57 CPU-hours for 200 replicates |
| Essential | Falsification battery: seasonality × PfPR model; future-PfPR placebo; predicted vs observed ITN/SMC/RTS,S effects; report the SMC before/after (null); Lesotho as a true zero anchor | Whether the association behaves causally | Seasonality cells, SMC results and MAP rasters exist; placebo and Lesotho need refits |
| Essential | Specification curve on v7 (exposure, gamma, covariate set, joint/separate, knots/k, monotone) | Post hoc choice dependence | Code exists for every variant; ~20 min each |
| Important | Leave-one-country-out refits and burden (Nigeria, DRC and the low-PfPR anchors first) | Influence | 37 × 7 fits, ~10–12 CPU-hours |
| Important | Period- and region-specific curves carried to the burden and trend; country random slopes; country-held-out prediction | Transportability | Subgroup curves saved; random slopes need refits |
| Important | Survey-weighted fit and design-based variance | Unweighted, cluster-ignoring likelihood | Weights and strata in the child-band data |
| Important | Measured RDT/microscopy prevalence as exposure; MAP vs measured prevalence by era | Exposure error and circularity | Regional measured prevalence already extracted |
| Important | Burden on the UN IGME all-cause envelope | Envelope dependence of country results | UN IGME age blocks are coarse; approximate allocation needed |
| Important | Birth-history sensitivities: merged 6–23 month band; heaping redistribution; exclude imputed dates; DHS-only vs MICS-only; drop sub-national MICS | Measurement artefacts in the age pattern | Refits on existing data |
| Desirable | Life-table (g-formula) counterfactual for death counts | Double counting and survivor selection when band fractions are summed | Extend the Figure 4 script |
| Desirable | Sub-national burden aggregation with year-specific under-5 weights | National-mean PfPR in a concave curve (Jensen's inequality) | Rasters exist |
| Desirable | 2023 endpoint; 2015–2019 vs 2019–2023 | Forecast and COVID-period sensitivity | Immediate |
| Desirable | Individual-level covariates (sex, twin, birth order), migration restriction, lagged/cumulative exposure | Compositional confounding, exposure timing | Individual variables in the child-band data |
| Desirable | CHAMPS/MITS triangulation of age-specific attributable fractions | External validity | Needs external aggregates |

## What would change the verdict

1. **Within-area estimates.** Country fixed-effect, Mundlak and region fixed-effect or difference-in-differences estimates of the 40%→20% and 20%→0% contrasts are close to the pooled ones.
2. **A robust excess.** The excess over IHME and UN IGME has a 95% interval excluding 1:
   - across reference prevalences of 1–5%;
   - with an independent exposure (Snow, or measured prevalence);
   - with a design-based bootstrap and MAP/IHME uncertainty propagated.
3. **A full specification curve.** Gamma, covariate set, knots/basis, joint vs separate models, period- and region-specific curves, leave-one-country-out refits and the UN IGME envelope do not overturn the conclusions.
4. **The trend claim.** It is either decomposed and shown to be PfPR-driven under a tested model, or withdrawn.
5. **Pre-specified falsification tests that pass:**
   - a null future-PfPR placebo;
   - seasonal amplitude that increases with PfPR in the Sahel;
   - a neonatal association that is explained or bounded.
6. **Accurate reporting.** The Methods and supplement describe the fitted model. All fitted variants and the decision history are disclosed. The estimand is relabelled, and the funding, competing-interest, ethics and data statements are completed.

Short of this, a narrower paper could suit a specialist journal. It would present an adjusted, age-specific association with explicit uncertainty and support limits, and the burden calculation as a clearly labelled scenario rather than a revised global estimate.

## Reproducing the checked numbers

The ✔ values were recomputed from these v7 outputs:
- `annual_comparison/country_age_estimates_2000_2024.csv` and `pfpr_curves.csv`:
  - trend decomposition, attributable fractions and PfPR;
  - reference-prevalence totals, which reproduce the saved 2024 total to within 0.003%;
  - neonatal attributable deaths.
- `pfpr_curves.csv`, <1 month row at PfPR 0: the neonatal HR, 0.942 (0.895–0.992).
- `annual_comparison/annual_totals_2000_2024.csv`: the 2015–2019 and 2015–2023 changes.
- `results/cbh/map_snow_gamma2_v1/burden/year_summary.csv`: Snow vs MAP.
- `results/cbh/seasonality_sahel_v1/REPORT.md`: the Sahel seasonality rate ratio.

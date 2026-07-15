# Malaria prevalence and child mortality — Methods & Results summary

Three analyses relating *P. falciparum* prevalence to all-cause under-5 mortality in
sub-Saharan Africa. All figures/tables below are reproduced by `run_all.R` into `results/`.

**Data sources:** MAP PfPR₂₋₁₀ raster (2024); IHME/GBD under-5 malaria deaths (2025);
UN IGME all-cause child mortality + World Bank live births, **GDP per capita**
(`NY.GDP.PCAP.CD`) and **DTP3** (`SH.IMM.IDPT`); DHS/MIS microdata (PR + birth-history
recodes) across 23 SSA countries.

---

## Component 1 — Country level

**Methods.** Malaria's share of all-cause under-5 deaths = IHME under-5 malaria deaths ÷
IGME all-cause under-5 deaths (U5MR/1000 × live births), computed for **all under-5** and
**neonatal-excluded** (denominator = 1mo–5y = U5MR − neonatal). Exposure = population-weighted
national **PfPR₂₋₁₀** (MAP 2024). Linear model `share ~ PfPR + log(GDP p.c.) + DTP3`; outlier
countries flagged by studentized residual from the prevalence-only fit.

**Results** (42 SSA countries). Malaria's share of child deaths rises strongly with
prevalence: **≈ +0.9 percentage-points of the all-U5 share (≈ +1.3 pp of the 1mo–5y share)
per +1 PfPR₂₋₁₀ point** (adjusted for GDP and DTP3; adjusted R² ≈ 0.5–0.6). **Burundi** is a
striking outlier (studentized residual ≈ 8) — malaria a far higher share of child deaths than
prevalence predicts — followed by Uganda, Rwanda and Ghana. GDP per capita is not significant
once DTP3 is included.

![Component 1: malaria share of child deaths vs national PfPR](results/component1_share_vs_pfpr.png)

Tables: `results/component1_model_coefficients.csv`, `component1_outliers.csv`.

---

## Component 2 — DHS subnational + national

**Methods.** One observation per **survey-region** across all usable DHS/MIS surveys
(malaria biomarker + birth history). All-cause **U5MR** and **1mo–5y** mortality from the
birth histories via `DHS.rates::chmort`. Exposure = child malaria prevalence put on a common
footing: **microscopy-equivalent** (RDT→microscopy via the Component-3 factor where microscopy
is missing), then **age-standardised 0.5–5y → PfPR₂₋₁₀ (2–10y)** with
`malariaAtlas::convertPrevalence` (Smith et al. 2007 age–prevalence model). Two specifications:

1. **Primary — log-rate linear mixed model:**
   `log(mortality) ~ PfPR2-10 + DTP3 + log(GDP) + %urban + year + stunting + (1 + PfPR2-10 || country)`.
2. **Count — negative-binomial GAM (`mgcv`):** region under-5 deaths with a `log(exposure)`
   offset, a smooth `s(year_c)`, and country random intercept + slope.

**Results** (600 survey-regions, 44 surveys, 23 countries). Higher child malaria prevalence
predicts higher child mortality, adjusted for all covariates:

| | LMM (log-rate) | GAM (nb count) |
|---|---|---|
| **U5MR**, per +10 PfPR₂₋₁₀ pts | **+6.9%** (95% CI +2.6, +11.4) | +6.0% (+2.2, +9.9) |
| **1mo–5y**, per +10 pts | **+11.0%** (+5.0, +17.5) | +9.9% (+4.3, +15.8) |

The effect is **positive in every country** (country random slopes; figure below) and covariates
behave as expected — DTP3 and % urban protective, stunting harmful, a ~3%/yr secular decline,
GDP ~null after the others. **Sensitivity** restricting to mid-transmission regions (PfPR₂₋₁₀
5–50%) attenuates the effect (LMM: U5MR +5.3%, 1mo–5y +7.4%; GAM: +3.9% / +6.2%) — positive but
weaker/less precise over the narrower range.

![Component 2: adjusted country-specific prevalence–mortality slopes](results/component2_country_slopes.png)

**Exposure-response shape.** Replacing linear prevalence with a smooth `s(PfPR₂₋₁₀)`:
**U5MR is linear** across the range (edf = 1.0 — a constant proportional effect), but
**1mo–5y is nonlinear** (edf ≈ 5.4 — steep at low prevalence, a plateau ~10–25%, then rising).
Removing neonatal deaths from the denominator reveals curvature the all-U5 relationship hides.

![Component 2: smooth exposure-response s(PfPR2-10)](results/component2_exposure_response_spline.png)

Tables: `results/component2_model_coefficients.csv`, `component2_count_model_coefficients.csv`
(`sample` = full / pfpr_5_50), `component2_country_slopes.csv`,
`component2_exposure_response_edf.csv`. Sensitivity figure: `component2_sensitivity_5to50.png`.

---

## Component 3 — RDT vs microscopy (supporting)

**Methods.** Across all DHS/MIS survey-regions with both tests, regress microscopy prevalence on
RDT prevalence → a consensus conversion factor, used by Component 2 to impute microscopy where
only RDT exists.

**Results.** RDTs over-detect (persistent HRP2 antigen), so microscopy runs lower:
**microscopy ≈ 0.74 × RDT** (Pearson r = 0.90, 431 survey-regions).

![Component 3: RDT vs microscopy](results/component3_rdt_vs_microscopy.png)

Table: `results/rdt_microscopy_conversion.csv`.

---

## Caveats
- **Component 1 is ecological** (country-level); **Component 2** adjusts for the main confounders
  but is observational, and region-level U5MR from a single survey is noisy.
- MAP PfPR₂₋₁₀ is 2024 vs the 2025 IHME export (1-year offset).
- IHME (GBD) and WHO differ substantially on malaria mortality; Component 1's numerator is IHME.
- DHS prevalence is age-standardised to PfPR₂₋₁₀ via a model (Smith 2007), and RDT is converted to
  a microscopy footing; both are documented approximations.

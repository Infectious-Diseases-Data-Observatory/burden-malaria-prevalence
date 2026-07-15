# Malaria prevalence and child mortality

Relationship between *Plasmodium falciparum* prevalence and all-cause under-5
mortality in sub-Saharan Africa, at country and subnational level, using MAP
(prevalence), IHME/GBD (malaria-attributable deaths), UN IGME (all-cause child
mortality) and DHS/MIS surveys (subnational prevalence + mortality).

The analysis has **three components**:

1. **Country level** — malaria's share of all-cause under-5 deaths (with and
   without neonatal deaths) as a function of national PfPR₂₋₁₀, with a
   multivariable model adding GDP per capita and DTP3 coverage, and
   identification of outlier countries where malaria is a higher share of child
   deaths than prevalence predicts.
2. **DHS subnational + national** — association between child malaria prevalence
   and all-cause under-5 mortality across all usable DHS/MIS survey-regions,
   adjusted for DTP3, GDP per capita, urban/rural, calendar year and stunting, in
   a mixed-effects model with country random slopes; reported as the **% change
   in under-5 mortality per +10 PfPR₂₋₁₀ points**.
3. **RDT vs microscopy** — a consensus conversion factor between the two malaria
   tests across DHS/MIS surveys, used in Component 2 to put all surveys on a
   common (microscopy) footing.

## Repository layout
```
R/00_utils.R                 shared config + helper functions
R/01_fetch_data.R            fetch/prepare all inputs into data/ (idempotent)
R/02_component1_country.R    Component 1
R/03_component3_rdt_microscopy.R  Component 3 (runs before Component 2)
R/04_component2_dhs.R        Component 2
run_all.R                    runs the whole pipeline
results/                     figures (.png) + summary tables (.csv)  [tracked]
data/                        raw inputs                              [git-ignored]
archive/                     superseded exploratory scripts          [git-ignored]
```

## Running it
From the repository root:
```r
Rscript run_all.R            # or source the R/ files in order in an R session
```
R packages: `terra, sf, jsonlite, httr, countrycode, malariaAtlas, geodata,
DHS.rates, lme4, ggplot2, patchwork`.

### Data inputs
`R/01_fetch_data.R` fetches automatically (internet required): the MAP PfPR₂₋₁₀
raster, GPW population raster, World Bank series (IGME U5MR / neonatal / infant,
live births, **GDP per capita `NY.GDP.PCAP.CD`**, **DTP3 `SH.IMM.IDPT`**), and the
DHS usable-survey catalogue.

Two inputs are **access-controlled and git-ignored** — a collaborator supplies
their own copy:

| Input | Path | How to obtain |
|---|---|---|
| IHME/GBD under-5 malaria deaths | `data/ihme_malaria_u5_deaths_by_country.csv` | GBD Results / Data Explorer (free login): Deaths · Malaria · Under 5 · Both · Number + Rate · all countries. |
| DHS/MIS microdata (PR + BR recodes) | `data/dhs/*.rds` | DHS account + a registered project; `01_fetch_data.R` downloads via `rdhs` using the cached login in `~/.rdhs.json`. |

## Methodology

### Component 1 — country level
Malaria's share of under-5 deaths = **IHME** under-5 malaria deaths ÷ **IGME**
all-cause under-5 deaths, where all-cause U5 deaths = (U5MR ÷ 1000) × live births.
Computed two ways: all under-5, and **neonatal-excluded** (denominator uses
U5MR − neonatal mortality = deaths ages 1 month–5 years; malaria kills virtually
no neonates, so the numerator is unchanged). Exposure is population-weighted
national **PfPR₂₋₁₀** from the MAP raster. A linear model
`share ~ PfPR + log(GDP p.c.) + DTP3` is fit; **outliers** are countries with a
large positive studentized residual from the prevalence-only fit (malaria a
higher share of child deaths than prevalence predicts), reported before and after
covariate adjustment.

### Component 2 — DHS subnational + national
Each survey-region contributes one observation. All-cause **U5MR** and **1mo–5y**
mortality are computed directly from the DHS birth histories with
`DHS.rates::chmort` (validated: Nigeria 2018 national U5MR = 132, matching the
published value). Child malaria prevalence (6–59 months) is taken from the PR
recode. The exposure is built in two steps:

1. **Microscopy-equivalent prevalence.** PfPR₂₋₁₀ (MAP) is microscopy-calibrated,
   but DHS RDTs over-detect (persistent HRP2 antigenaemia). Where only RDT is
   available we convert to microscopy using the Component-3 factor
   (microscopy ≈ 0.74 × RDT).
2. **Age-standardisation to PfPR₂₋₁₀ (see below).**

A mixed-effects model is then fit on the **log** of mortality (so coefficients are
% changes):
```
log(mortality) ~ PfPR2-10 + DTP3 + log(GDP p.c.) + % urban + calendar year
                 + stunting + (1 + PfPR2-10 | country) + (1 | survey)
```
for both U5MR and 1mo–5y mortality. The headline is the adjusted **% change in
under-5 mortality per +10 PfPR₂₋₁₀ points**, with country-specific random slopes.

#### PfPR age-standardisation (DHS 6–59 mo → PfPR₂₋₁₀)
DHS measures infection prevalence in children **~6–59 months (0.5–5 y)**, whereas
the modelled reference standard (MAP PfPR) is prevalence in **ages 2–10**. Infection
prevalence follows a characteristic **age–prevalence curve**: it rises steeply
through infancy and early childhood and plateaus in older children, so at a fixed
transmission intensity the prevalence measured over one age window implies a
predictable value over another. We use `malariaAtlas::convertPrevalence(p, 0.5, 5,
2, 10)`, which applies the **Smith et al. (2007, *Malaria Journal*)** model of that
age–prevalence relationship to map the 0.5–5 y estimate onto the 2–10 y standard.
In practice the adjustment is modest (a 30% prevalence in 0.5–5 y maps to ≈31% in
2–10 y) but it makes the DHS exposure directly comparable to MAP PfPR₂₋₁₀ and to
Component 1.

### Component 3 — RDT vs microscopy
Across all DHS/MIS survey-regions reporting both tests, microscopy prevalence is
regressed on RDT prevalence (through-origin and with intercept) to give a
**conversion factor** (`results/rdt_microscopy_conversion.csv`), used by
Component 2 to impute microscopy where only RDT exists.

## Key results
_(reproduced by the pipeline into `results/`)_
- **Component 3:** microscopy ≈ **0.74 × RDT** (r ≈ 0.90, 431 survey-regions).
- **Component 1:** `results/component1_share_vs_pfpr.png`, `component1_model_coefficients.csv`, `component1_outliers.csv`.
- **Component 2** (600 survey-regions, 44 surveys, 23 countries; adjusted for DTP3,
  GDP p.c., % urban, year, stunting): **+6.7% under-5 mortality per +10 PfPR₂₋₁₀
  points** (95% CI +2.4 to +11.2), and **+11.7% for 1mo–5y mortality**
  (95% CI +5.3 to +18.4). **Sensitivity** restricting to mid-transmission regions
  (PfPR₂₋₁₀ 5–50%, 354 regions): attenuated to +5.1% (CI −0.7 to +11.2) for U5MR
  and +7.7% (CI +0.1 to +15.9) for 1mo–5y — positive but weaker/less precise over
  the narrower range. See `results/component2_country_slopes.png`,
  `component2_model_coefficients.csv` (`sample` column = `full` / `pfpr_5_50`),
  `component2_country_slopes.csv`.

## Caveats
- MAP PfPR₂₋₁₀ is available through 2024; the IHME export is 2025, so Component 1
  pairs 2024 prevalence with 2025 mortality (documented 1-year offset).
- Component 1 is ecological (country-level). Component 2 adjusts for the main
  confounders but remains observational; U5MR at the region level from a single
  survey is noisy (wide sampling error).
- IHME (GBD) and WHO differ substantially on malaria mortality; this project uses
  IHME for Component 1's numerator.

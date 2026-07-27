# Vaccine, paediatric HIV and SMC covariate sources

This note separates variables already implemented in the DHS/MIS rebuild from
external series that still require a data-access and modelling decision.

## Vaccination

### Direct DHS/MIS measures implemented

DHS does not ask a single question about the date a vaccine was introduced.
Instead, the vaccination module records whether an eligible child received
each dose, using the vaccination card where available and the mother's report
otherwise. The standard recodes are:

| Construct | Standard DHS recodes | Implemented variable |
|---|---|---|
| DPT-HepB-Hib (pentavalent) | `h51`, `h52`, `h53` | `pentavalent3_reg` |
| Pneumococcal conjugate vaccine | `h54`, `h55`, `h56` | `pcv3_reg` |
| Rotavirus vaccine | `h57`, `h58`, `h59` | `rotavirus_complete_reg` |

The implementation estimates weighted coverage among living children aged
12–23 months in each survey region. Pentavalent and PCV coverage require three
recorded doses. Rotavirus completion requires two doses unless `h59` is
observed for eligible children in that survey, in which case it requires three.
This follows the DHS approach of using the number of doses known to have been
received rather than requiring an uninterrupted sequence of dose fields.

These variables are present mainly in newer recodes. In the 936-row validation
panel, their missingness is 55.1%, 57.1% and 59.0%, respectively. They remain in
the analysis dataset and missingness report but are not imputed and do not enter
the main ridge block under the prespecified 5% rule. A missing field is not
interpreted as zero coverage or as evidence that the vaccine had not been
introduced.

Source: [DHS Guide to Statistics: Vaccination](https://dhsprogram.com/data/Guide-to-DHS-Statistics/Vaccination.htm).

### National WUENIC series implemented

The pipeline now extracts the WHO/UNICEF Estimates of National Immunization
Coverage (WUENIC) country-year series from
`data/fusion_GLOBAL_DATAFLOW_UNICEF_1.0_all.csv`:

| UNICEF indicator | Analysis variable | Interpretation |
|---|---|---|
| `IM_HIB3` | `hib3_wuenic` | Third dose of a Hib-containing vaccine |
| `IM_PCVC` | `pcv3_wuenic` | Completion of the national PCV schedule |
| `IM_ROTAC` | `rotac_wuenic` | Final recommended rotavirus dose |

These series are better suited to tracking national rollout over time than the
sparse direct DHS fields. `R_dhs/02b_extract_unicef_immunisation.R` retains
total-sex, age-12–23-month, percentage WUENIC estimates and writes a compact
country-year panel. `R_dhs/03_build_analysis_dataset.R` joins the panel by ISO3
and exact survey year.

Reported zeros are preserved. Years before a country's first WUENIC estimate,
and countries with no series for an antigen, are assigned zero as
not-yet-introduced; gaps after a series begins remain missing. Separate status
columns distinguish `wuenic_estimate`,
`pre_series_assumed_not_introduced`,
`no_series_assumed_not_introduced`, and
`missing_after_series_start`. In the validation analysis, all 936
survey-region-years have values for all three national series, so the variables
pass the missingness rule and enter the ridge block.

Source: [UNICEF immunization resources and WUENIC downloads](https://data.unicef.org/resources/immunization/).

## Paediatric HIV

### Derived child HIV prevalence implemented

Child (0-14) HIV prevalence is now a national ridge covariate. The source is the
UNAIDS 2025 estimates workbook distributed by UNICEF,
`data/HIV_Epidemiology_Children_Adolescents_2025.xlsx`. That workbook does not
publish a prevalence rate; its child series are counts and rates (people living
with HIV, annual AIDS deaths and new infections, AIDS-death and incidence
rates, and mother-to-child transmission). Prevalence is therefore derived:

| Component | Source | Age / sex |
|---|---|---|
| Numerator | Estimated number of people living with HIV | 0-14, both sexes |
| Denominator | World Bank `SP.POP.0014.TO` population | 0-14, total |

`R_dhs/03_build_analysis_dataset.R` (`attach_hiv_prevalence`) computes
`hiv_prev = 100 * PLHIV(0-14) / population(0-14)`, joined by ISO3 and exact
survey year, and writes the compact panel
`data/derived_dhs/hiv_prevalence_country_year.csv`. The values in the workbook
that read `"<500"` are parsed at their midpoint (250); every derived prevalence
is strictly positive.

Derived child prevalence spans roughly three orders of magnitude across
sub-Saharan Africa (about 0.02% in Madagascar to about 3.8% in Eswatini), so it
enters the model on the log scale as `log_hiv_prev`, alongside the other
log-scaled national continuous covariates. It is standardised inside the same
ridge block as the other covariates.

Nigeria and Comoros publish only a 15-19 UNAIDS series in this workbook, so they
have no derived child prevalence and stay missing. That is about 4% of the
primary sample, below the 5% threshold, so `log_hiv_prev` is retained and singly
imputed (country median, then overall median) with the `hiv_prev_status` flag
(`unaids_2025_estimate` versus `no_under15_series`) preserved for audit.

Source: [UNICEF HIV/AIDS data
page](https://data.unicef.org/resources/dataset/hiv-aids-statistical-tables/)
and the underlying [UNAIDS AIDSinfo estimates
dataset](https://aidsinfo.unaids.org/dataset).

### Alternatives and caveats

Child prevalence alone does not fully track the mortality effect of HIV:
effective ART can increase survival and therefore leave prevalence stable even
as AIDS mortality falls. The same workbook also carries an AIDS-related death
rate per 100,000 (0-14) for the same 33 countries, which is the closest single
measure of the changing paediatric mortality burden and is the recommended
alternative or companion covariate; a further option is child prevalence plus
paediatric ART/PMTCT coverage. UNAIDS uncertainty intervals (`Lower`/`Upper` in
the workbook) are available for sensitivity work.

## Seasonal malaria chemoprevention

The [WHO World malaria report
2025](https://www.who.int/publications/i/item/9789240117822) is the most
complete public summary currently identified:

- Figure 5.6 maps the subnational areas receiving SMC and the number of cycles
  delivered per district in 2024.
- Table 5.1 reports the average number of children treated by country and year,
  2012–2024.
- The report states that national malaria programmes supplied the data to
  LSHTM and Medicines for Malaria Venture, and that the SMC Alliance assembled
  the subnational data.

The report annex download page does not expose a machine-readable historical
district-year SMC panel. The best route to a complete analysis variable is
therefore a request to the SMC Alliance/MMV/LSHTM or WHO Global Malaria
Programme. The older [ACCESS-SMC collection at LSHTM Data
Compass](https://datacompass.lshtm.ac.uk/id/eprint/2136/) contains coverage and
implementation data for seven West and Central African countries in 2015–2016,
but the data files require a request and data-sharing agreement. The [Malaria
Atlas Project SMC coverage
project](https://malariaatlas.org/project-resources/modelling-coverage-of-seasonal-malaria-chemoprevention-smc/)
is collating country reports, treatment cycles and targeted subnational units,
but describes this work as ongoing.

The preferred region-year covariate is the fraction of eligible children
treated, with the number of cycles delivered retained as a second dimension.
A binary "introduced" flag is a weaker fallback. Districts should be mapped to
DHS survey regions using spatial overlap.

SMC is also an intervention responding to malaria risk and may be on the causal
path from malaria prevalence to mortality. Including it changes the estimand
from the total association with malaria prevalence to an association
conditional on SMC delivery. It should therefore be analysed as a
prespecified sensitivity or effect-modification analysis, not added
automatically to the primary confounder block.

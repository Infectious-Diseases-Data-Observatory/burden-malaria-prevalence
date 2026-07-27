# Vaccine, paediatric HIV and SMC covariate sources

This note separates variables already implemented in the DHS/MIS rebuild from
external series that would require a subsequent data-access and modelling
decision.

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

### Recommended longitudinal source

For adjustment across the full 2000–2025 study period, the preferable source
is the WHO/UNICEF Estimates of National Immunization Coverage (WUENIC)
country-year series:

- `HIB3`: third dose of a Hib-containing vaccine;
- `PCV3`: third dose of pneumococcal conjugate vaccine;
- `ROTAC`: final recommended dose of rotavirus vaccine.

These series are better suited to tracking national rollout over time than the
sparse direct DHS fields. Before merging, the current release's metadata should
be used to distinguish a true zero before introduction from a missing estimate.
The current analysis does not yet include WUENIC.

Source: [UNICEF immunization resources and WUENIC downloads](https://data.unicef.org/resources/immunization/).

## Paediatric HIV

The preferred primary source is the [UNAIDS AIDSinfo estimates
dataset](https://aidsinfo.unaids.org/dataset), which provides annual
country-level Spectrum estimates. UNAIDS recommends the unrounded estimates
for calculations and derived indicators. The [UNICEF HIV/AIDS data
page](https://data.unicef.org/resources/dataset/hiv-aids-statistical-tables/)
provides a convenient child-focused download and country dashboard.

Candidate country-year variables are:

- HIV prevalence among children aged 0–14;
- AIDS-related mortality among children aged 0–14;
- paediatric ART coverage and PMTCT coverage.

Child prevalence alone does not fully track the mortality effect of HIV:
effective ART can increase survival and therefore leave prevalence stable even
as AIDS mortality falls. For the child-mortality model, AIDS-related mortality
is the closest single measure of the changing mortality burden; a prespecified
alternative is child prevalence plus paediatric ART coverage. Annual values
should be joined to survey year, with interpolation only across short internal
gaps and with the UNAIDS uncertainty intervals retained for sensitivity work.

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

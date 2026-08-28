# DHS/MIS manuscript alignment notes

These notes compare the rebuilt DHS/MIS pipeline with the current statements in
`main.tex`. They cover only the DHS/MIS survey analysis; burden time series,
IHME/WHO comparisons and trial triangulation are outside this rebuild.

## Statements that remain supported

- The primary sample contains 921 survey-region-years from 105 surveys and 34
  countries, covering survey years 2000–2024.
- Higher MAP PfPR2-10 is associated with higher post-neonatal and all-under-5
  mortality.
- The ridge-linear post-neonatal estimate is an 8.65% increase per 10
  percentage-point increase in PfPR2-10 (95% CI 5.69% to 11.68%).
- The neonatal negative-control estimate remains null in the ridge-linear model:
  0.99% (95% CI -1.18% to 3.21%; p=0.372).
- An independent five-model neonatal comparison selects the linear
  no-interaction form, narrowly ahead of the linear time-interaction model
  (delta AIC 1.59). In the prespecified REML refit (post-neonatal-selected
  spline), the neonatal prevalence smooth is linear (edf 1.00) and null
  (p=0.37).
- All primary models use a negative-binomial likelihood, log-exposure offset,
  smooth calendar year and country random intercepts. The country-specific
  random PfPR slope was removed from the primary structure (see
  INCLUDE_COUNTRY_PFPR_SLOPE in 00_config.R) and is carried only as a
  structural sensitivity, in script 07 (arm `with_country_pfpr_slope`) and
  script 15 (arm "With country PfPR random slope").

## Statements requiring revision

- **Functional form:** `main.tex` says that the best fit was linear. The formal
  comparison now selects the spline without a time interaction. Its AIC is 4.91
  lower than the spline interaction, 10.19 lower than the linear model without
  interaction, 9.23 lower than the linear interaction and 11.21 lower than the
  full tensor surface.
- **Time interaction:** the spline interaction has p=0.16 and is 4.91 AIC
  units behind the selected additive spline. The linear interaction has
  p=0.0075 but fits substantially worse. The manuscript should treat temporal
  variation as a sensitivity rather than the primary specification.
- **Single linear effect:** the selected spline model does not have one
  prevalence coefficient valid across all prevalence values. The 8.65%
  estimate is from the prespecified ridge-linear, no-interaction
  summary model and should be labelled as such.
- **Attributable fractions:** the selected time-stable spline estimates 14.1%,
  34.4% and 40.2% at 10%, 30% and 50% PfPR, respectively,
  relative to a 1% PfPR counterfactual. These replace the current statements of
  approximately 7% at 10% and one third at 50%.
- **Covariate description:** the current text contains duplicated and placeholder
  covariate wording. The rebuilt candidate set is screened before imputation.
  Seventeen variables meet the at-most-5% missingness rule and enter one ridge
  block: regional urban residence, DTP3, measles vaccination, facility delivery,
  maternal education, exclusive breastfeeding, short birth interval, maternal
  age at first birth, improved water, improved sanitation and electricity; and
  exact-year national WUENIC Hib3, PCV-completion and final-dose rotavirus
  coverage and derived child (0–14) HIV prevalence (log scale); and national
  log GDP per capita and log health expenditure per capita. National DTP3,
  electricity access and health-expenditure-share-of-GDP were dropped as
  redundant with their regional or per-capita counterparts.
- **Direct vaccine measures:** the candidate set now also includes regional
  pentavalent-dose-3 (`h51`–`h53`), PCV-dose-3 (`h54`–`h56`) and completion of
  the survey-specific rotavirus schedule (`h57`–`h59`). These fields are
  concentrated in recent DHS recodes. Their missingness in the validation
  panel is 55.1%, 57.1% and 59.0%, respectively, so all three are flagged and
  excluded rather than imputed.
- **National vaccine measures:** the WUENIC series have no missing values after
  exact country-year matching and explicit pre-introduction zero coding. Adding
  these three ridge covariates changes the AIC-selected model from the spline
  interaction to the spline without interaction, while changing the
  cross-sectional effect and attributable-fraction estimates only modestly.
- **Paediatric HIV:** child (0–14) HIV prevalence is now included as a national
  ridge covariate. It is derived from the UNAIDS 2025 estimates workbook
  (`HIV_Epidemiology_Children_Adolescents_2025.xlsx`, estimated number of
  children 0–14 living with HIV) divided by the World Bank 0–14 population
  (`SP.POP.0014.TO`), joined by ISO3 and exact survey year and entered on the
  log scale as `log_hiv_prev`. Nigeria and Comoros publish only a 15–19 UNAIDS
  series, so they are missing (~4% of the sample) and singly imputed under the
  ≤5% rule with a status flag. The manuscript's HIV-adjustment paragraph should
  cite this derived measure and its limitation.
- **Other excluded covariates:** wealth (5.56% missing), stunting (25.2%),
  underweight (27.5%), wasting (27.5%) and political stability (100%) are also
  flagged and excluded rather than imputed.
- **Date range:** the scripts request/support 2000–2025 and explicitly flag a
  missing 2025 MAP surface rather than substitute 2024. The reproduced analysis
  currently ends in 2024 because that is the last paired survey-year in the
  validated aggregate panel.
- **Likelihood sensitivity:** the ridge log-Gaussian rate sensitivity gives a
  positive but less precise linear post-neonatal estimate of 5.73% per 10
  points (95% CI -1.76% to 13.78%; p=0.14). Its selected additive-spline
  prevalence term also has p=0.14. This sensitivity should not be described as
  giving identical precision to the negative-binomial analysis.

The legacy and rebuilt headline results pass all prespecified checks in
`results/dhs_rebuild/reproduction_check.csv`; this establishes substantive
reproduction while allowing the rebuilt specification and covariate policy to
differ transparently.

# DHS/MIS manuscript alignment notes

These notes compare the rebuilt DHS/MIS pipeline with the current statements in
`main.tex`. They cover only the DHS/MIS survey analysis; burden time series,
IHME/WHO comparisons and trial triangulation are outside this rebuild.

## Statements that remain supported

- The primary sample contains 921 survey-region-years from 105 surveys and 34
  countries, covering survey years 2000–2024.
- Higher MAP PfPR2-10 is associated with higher post-neonatal and all-under-5
  mortality.
- The ridge-linear post-neonatal estimate is an 8.75% increase per 10
  percentage-point increase in PfPR2-10 (95% CI 5.77% to 11.83%).
- The neonatal negative-control estimate remains null in the ridge-linear model:
  0.98% (95% CI -1.16% to 3.17%; p=0.373).
- An independent five-model neonatal comparison also selects the additive
  spline without a time interaction, narrowly ahead of the linear
  no-interaction model (delta AIC 0.97). Its REML prevalence smooth has edf
  1.00 and p=0.374.
- All primary models use a negative-binomial likelihood, log-exposure offset,
  smooth calendar year, country random intercepts and country random linear
  PfPR slopes.

## Statements requiring revision

- **Functional form:** `main.tex` says that the best fit was linear. The formal
  comparison now selects the spline without a time interaction. Its AIC is 1.81
  lower than the spline interaction, 8.16 lower than the linear model without
  interaction and 16.20 lower than the linear interaction.
- **Time interaction:** the spline interaction has p=0.078 and is 1.81 AIC
  units behind the selected additive spline. The linear interaction has
  p=0.0028 but fits substantially worse. The manuscript should treat temporal
  variation as a sensitivity rather than the primary specification.
- **Single linear effect:** the selected spline model does not have one
  prevalence coefficient valid across all prevalence values. The 8.75%
  estimate is from the prespecified ridge-linear, no-interaction
  summary model and should be labelled as such.
- **Attributable fractions:** the selected time-stable spline estimates 14.0%,
  34.0% and 40.1% at 10%, 30% and 50% PfPR, respectively,
  relative to a 1% PfPR counterfactual. These replace the current statements of
  approximately 7% at 10% and one third at 50%.
- **Covariate description:** the current text contains duplicated and placeholder
  covariate wording. The rebuilt candidate set is screened before imputation.
  Nineteen variables meet the at-most-5% missingness rule and enter one ridge
  block: regional urban residence, DTP3, measles vaccination, facility delivery,
  maternal education, exclusive breastfeeding, short birth interval, maternal
  age at first birth, improved water, improved sanitation and electricity; and
  exact-year national WUENIC Hib3, PCV-completion and final-dose rotavirus
  coverage; national DTP3, log GDP per capita, health expenditure as a share of
  GDP, log health expenditure per capita and electricity access.
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
- **Paediatric HIV:** the local UNICEF global dataflow does not contain child
  HIV prevalence. `HVA_EPI_LHIV` is a count of 10–19-year-olds living with HIV
  for 2010–2024, not prevalence or an under-five measure, so it has not been
  merged under a different label. A separate UNAIDS prevalence series is still
  required.
- **Other excluded covariates:** wealth (5.56% missing), stunting (25.2%),
  underweight (27.5%), wasting (27.5%) and political stability (100%) are also
  flagged and excluded rather than imputed.
- **Date range:** the scripts request/support 2000–2025 and explicitly flag a
  missing 2025 MAP surface rather than substitute 2024. The reproduced analysis
  currently ends in 2024 because that is the last paired survey-year in the
  validated aggregate panel.
- **Likelihood sensitivity:** the ridge log-Gaussian rate sensitivity gives a
  positive but less precise linear post-neonatal estimate of 5.75% per 10
  points (95% CI -1.73% to 13.80%; p=0.136). Its selected additive-spline
  prevalence term also has p=0.136. This sensitivity should not be described as
  giving identical precision to the negative-binomial analysis.

The legacy and rebuilt headline results pass all prespecified checks in
`results/dhs_rebuild/reproduction_check.csv`; this establishes substantive
reproduction while allowing the rebuilt specification and covariate policy to
differ transparently.

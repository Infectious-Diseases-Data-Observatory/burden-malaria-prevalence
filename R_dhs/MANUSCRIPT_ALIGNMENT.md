# DHS/MIS manuscript alignment notes

These notes compare the rebuilt DHS/MIS pipeline with the current statements in
`main.tex`. They cover only the DHS/MIS survey analysis; burden time series,
IHME/WHO comparisons and trial triangulation are outside this rebuild.

## Statements that remain supported

- The primary sample contains 921 survey-region-years from 105 surveys and 34
  countries, covering survey years 2000–2024.
- Higher MAP PfPR2-10 is associated with higher post-neonatal and all-under-5
  mortality.
- The ridge-linear post-neonatal estimate is a 9.04% increase per 10
  percentage-point increase in PfPR2-10 (95% CI 5.98% to 12.19%).
- The neonatal negative-control estimate remains null in the ridge-linear model:
  0.95% (95% CI -1.19% to 3.14%; p=0.388).
- All primary models use a negative-binomial likelihood, log-exposure offset,
  smooth calendar year, country random intercepts and country random linear
  PfPR slopes.

## Statements requiring revision

- **Functional form:** `main.tex` says that the best fit was linear. The formal
  four-model comparison selects the spline-plus-time-interaction specification.
  Its AIC is 2.17 lower than the linear interaction model, 3.30 lower than the
  spline without interaction, and 11.99 lower than the linear model without
  interaction.
- **Time interaction:** the selected spline interaction is borderline
  (p=0.057). The linear interaction model has a clearer interaction
  (p=0.0013), but is not the AIC minimum. The manuscript should report both the
  selection rule and this uncertainty rather than simply assert an interaction.
- **Single linear effect:** the selected spline-plus-interaction model does not
  have one prevalence coefficient valid across all prevalence values and years.
  The 9.04% estimate is from the prespecified ridge-linear, no-interaction
  summary model and should be labelled as such.
- **Attributable fractions:** at the 2014 reference year, the selected model
  estimates 14.6%, 35.6% and 43.0% at 10%, 30% and 50% PfPR, respectively,
  relative to a 1% PfPR counterfactual. These replace the current statements of
  approximately 7% at 10% and one third at 50%.
- **Covariate description:** the current text contains duplicated and placeholder
  covariate wording. The rebuilt candidate set is screened before imputation.
  Sixteen variables meet the at-most-5% missingness rule and enter one ridge
  block: regional urban residence, DTP3, measles vaccination, facility delivery,
  maternal education, exclusive breastfeeding, short birth interval, maternal
  age at first birth, improved water, improved sanitation and electricity; and
  national DTP3, log GDP per capita, health expenditure as a share of GDP, log
  health expenditure per capita and electricity access.
- **Direct vaccine measures:** the candidate set now also includes regional
  pentavalent-dose-3 (`h51`–`h53`), PCV-dose-3 (`h54`–`h56`) and completion of
  the survey-specific rotavirus schedule (`h57`–`h59`). These fields are
  concentrated in recent DHS recodes. Their missingness in the validation
  panel is 55.1%, 57.1% and 59.0%, respectively, so all three are flagged and
  excluded rather than imputed. Re-running the model suite therefore leaves
  the model comparison and headline estimates unchanged.
- **Other excluded covariates:** wealth (5.56% missing), stunting (25.2%),
  underweight (27.5%), wasting (27.5%) and political stability (100%) are also
  flagged and excluded rather than imputed.
- **Date range:** the scripts request/support 2000–2025 and explicitly flag a
  missing 2025 MAP surface rather than substitute 2024. The reproduced analysis
  currently ends in 2024 because that is the last paired survey-year in the
  validated aggregate panel.
- **Likelihood sensitivity:** the ridge log-Gaussian rate sensitivity gives a
  positive but less precise linear post-neonatal estimate of 7.13% per 10
  points (95% CI -0.60% to 15.47%; p=0.072). The selected nonlinear
  post-neonatal prevalence term has p=0.049, while its time interaction has
  p=0.599. This sensitivity should not be described as giving identical
  precision to the negative-binomial analysis.

The legacy and rebuilt headline results pass all prespecified checks in
`results/dhs_rebuild/reproduction_check.csv`; this establishes substantive
reproduction while allowing the rebuilt specification and covariate policy to
differ transparently.

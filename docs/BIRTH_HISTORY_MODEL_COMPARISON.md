# Complete birth histories: Adebayo (2005) versus Burstein (2018)

Reviewed 7 September 2026. This comparison concerns the **complete-birth-history (CBH) training model** in Burstein et al., not the later reconstruction of histories from summary data. It is a methods comparison, not a change to the agreed primary analysis or a report of new empirical fits.

Sources: Adebayo and Fahrmeir (2005), *Statistics in Medicine* 24:709–728, especially pp. 714–718; Burstein, Wang, Reiner and Hay (2018), *PLOS Medicine* 15:e1002687, especially PDF pp. 5–8, supplied as `Hay_PLOS_Medicine.pdf`. See the [Burstein paper](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002687), the [Adebayo specification](ADEBAYO_2005_MGCV.md) and the authors' code links below.

## Main comparison

Both methods fit logistic discrete-time survival models to observed birth histories. Their main differences are **the length of one observation interval**, how covariate effects vary with age, the pooling of information, and the purpose of the adjustment variables. Burstein's implementation already uses `mgcv::bam`.

| Feature | Adebayo and Fahrmeir: M7 | Burstein et al.: CBH model |
|---|---|---|
| Observation | One child-month at risk | One child-age-band with adequate potential follow-up |
| Response | Death in that month, conditional on reaching it alive | Death in the whole band, conditional on entering it alive |
| Likelihood/link | Bernoulli, logit | Bernoulli, logit |
| Age baseline | Smooth function of attained age | Separate intercept for each of seven bands |
| Age-dependent covariates | Smooth changes with attained age in breastfeeding and maternal-age-category associations | Separate two-dimensional birth-year/SDI smooth for each band |
| Other covariates | Selected individual and household factors; optional socioeconomic extension | One joint smooth of maternal proportion of children dead, birth order and maternal age at birth, shared across bands |
| Geographic structure | Adjacent-district MRF plus independent district effect | Independent survey and country-by-age-band intercepts; no adjacency-based spatial effect |
| Fitting | BayesX MCMC with variance hyperpriors | `mgcv::bam`, restricted maximum likelihood; Gaussian coefficient draws for uncertainty |
| Pooling | One model across monthly ages in the Nigerian study | One model across age bands within each broad geographic region |
| Covariate purpose | Explore mortality associations and geographic variation | Predict mortality using variables available for subsequent summary-history applications |

Neither paper estimates a causal effect of PfPR. Adebayo's covariate list also needs causal review rather than automatic adoption.

## Burstein's model, restricted to observed complete histories

For child i of mother m in country c, survey s and age band a:

\[
Y_{mia}\sim\operatorname{Bernoulli}(q_{mia}),
\]

\[
\operatorname{logit}(q_{mia})=
\beta_a+g_{1a}(\mathrm{birth\ year}_i,\mathrm{SDI})
+g_2\left(\frac{CD_m}{CEB_m},\mathrm{birth\ order}_i,
\mathrm{maternal\ age\ at\ birth}_i\right)
+u_s+v_{ca}.
\]

Here CD and CEB are the mother's total children dead and children ever born **reported at interview**. The paper describes the second argument of g₂ as children born by the index child's birth; the published code implements it as `birthorder`. The terms g₁ and g₂ are joint thin-plate regression splines, not separate additive smooths of each covariate. Only g₁ differs by age band; g₂ and the survey effect are shared across ages. Country-by-age effects are independent Gaussian intercepts with a shared variance component. There is no maternal random intercept in this specification.

Using descriptive variable names, the full-model formula can be expressed as:

```r
library(mgcv)

# d: correctly prepared CBH child-band data for one broad geographic region.
# Continuous predictors must use the preprocessing chosen for the model.
d$age_band <- factor(d$age_band)
d$survey <- factor(d$survey)
d$country_age <- interaction(d$country, d$age_band, drop = TRUE)

fit <- bam(
  death ~ 0 + age_band +
    s(birth_year, sdi, by = age_band, bs = "tp", k = 9) +
    s(maternal_dead_fraction, birth_order, maternal_age_birth,
      bs = "tp", k = 9) +
    s(survey, bs = "re") +
    s(country_age, bs = "re"),
  data = d,
  family = binomial(link = "logit"),
  weights = weight_normalized,
  method = "fREML",
  na.action = na.fail
)
```

This is an explanatory transcription, not a ready-to-run DHS pipeline or a numerical replication. The authors use `cdceb100` rather than the descriptive name above. Their full configuration sets both smooth basis dimensions to 9. The fitting script normalizes survey weights to sum to the number of retained child-band observations within each survey, then supplies them to `bam`. It leaves the fitting method at the `bam` default; the example makes fast REML explicit. These likelihood weights do not themselves provide design-based uncertainty for DHS clustering and stratification. [Full-model configuration, lines 47–49](https://github.com/royburst/sbh_agespecific_indirect_paper_code/blob/master/EXTVAL_launch_models_fitting.R#L47-L49); [fitting and weights, lines 191–210](https://github.com/royburst/sbh_agespecific_indirect_paper_code/blob/master/EXTVAL_run_models_fitting.R#L191-L210).

The seven published bands are birth through 28 days, 29 days through 5 months, 6–11, 12–23, 24–35, 36–47 and 48–59 months. These **do not match our six groups**: <1, 1–3, 4–11, 12–23, 24–35 and 36–59 completed months. Changing boundaries requires rebuilding the child-band outcomes. It cannot be done by relabelling coefficients.

There are small but material paper/code differences to settle before replication:

- The paper describes SDI at birth year. The fitting scripts join SDI to the rounded year of **entry into each age band**, while retaining birth year as the other smooth argument. [SDI join](https://github.com/royburst/sbh_agespecific_indirect_paper_code/blob/master/EXTVAL_run_models_fitting.R#L160-L163).
- The external-fitting configuration uses a one-month neonatal boundary; the cross-validation configuration uses 29/30 months. [External configuration](https://github.com/royburst/sbh_agespecific_indirect_paper_code/blob/master/EXTVAL_launch_models_fitting.R#L64-L71); [cross-validation configuration](https://github.com/royburst/sbh_agespecific_indirect_paper_code/blob/master/CROSSVAL_launch_models.R#L92-L99).
- The paper says all covariates are standardized. The inspected fitting script explicitly skips SDI and birth order at its scaling step. [Preprocessing](https://github.com/royburst/sbh_agespecific_indirect_paper_code/blob/master/EXTVAL_run_models_fitting.R#L138-L149).
- The current cross-validation launcher has the reduced `INDIV` model active; its `FULL` formula is present but commented out. The external-fitting launcher contains the full formula used above. [Cross-validation model choices](https://github.com/royburst/sbh_agespecific_indirect_paper_code/blob/master/CROSSVAL_launch_models.R#L33-L53).

A notation correction: the text following the paper's baseline logistic equation describes exp(βₐ) as a probability. Under the stated logit model, the probability is **plogis(βₐ)** and exp(βₐ) is the odds. The prediction code uses `plogis`. [Prediction function](https://github.com/royburst/sbh_agespecific_indirect_paper_code/blob/master/utils.R#L246-L286).

## The observation interval and censoring matter

Burstein's qₐ is the probability of dying anywhere within the band; it is not a monthly hazard and not deaths divided by person-time. If monthly hazards are hₜ, the corresponding band risk for a given covariate trajectory is

\[
q_a=1-\prod_{t\in a}(1-h_t).
\]

Consequently, replacing monthly rows with wide bands changes both the response scale and the covariate model. A logistic model for qₐ is not generally equivalent to a logistic model for hₜ. Neither requires an offset for its own Bernoulli trials. Converting a band risk to a constant rate using −log(1−qₐ)/band width requires an additional constant-hazard assumption.

The authors' CBH code requires the child's **potential age-band end** to precede or coincide with interview. It excludes the band when this is not true, including observations from children who died before the band's end. This is more specific than simply dropping the unfinished band of surviving children. [Eligibility calculation and accompanying comment, lines 167–174](https://github.com/royburst/sbh_agespecific_indirect_paper_code/blob/master/EXTVAL_run_models_fitting.R#L167-L174).

For example, suppose a child would have been 18 months old at interview. The 12–23-month band is excluded whether that child is alive at interview or died at 14 months. A child who could have reached 24 months by interview can contribute to the band. Earlier eligible bands remain available.

This gives a clear denominator for complete-band risks, but discards observed follow-up in the incomplete band. With our combined 36–59-month group, a direct copy would require the child to have been able to reach 60 months by interview before that group contributes. Monthly records can retain substantially more observed history and align exposure more closely to calendar time. Their administrative-censoring convention must still be explicit: keeping early deaths while excluding comparable survivors solely because their interval is unfinished would create differential inclusion. Use equal potential-follow-up eligibility within the chosen interval, or a likelihood that explicitly accommodates partial follow-up.

## Implications for the six age-specific PfPR effects

**My recommendation is to borrow the discrete-survival/GAM framework and scalable fitting approach, while defining the data intervals and adjustment set for our causal question.**

1. **Retain six separate fits.** Burstein estimates age-specific components inside a joint model. That is not equivalent to the user's confirmed choice that each age group has its own prevalence curve, covariate effects and other parameters. Monthly or finer split person-time observations can sit inside each separately fitted group; reporting groups do not need to be the observation intervals.
2. **Exclude the interview-time CD/CEB predictor from the causal PfPR adjustment set.** It contains the index child's death, siblings' deaths and potentially deaths after the exposure period. This is useful outcome information for the paper's prediction task. My causal assessment is that conditioning on it when estimating the effect of malaria exposure would condition on a quantity partly determined by the outcome and potentially by malaria itself. It could remove or distort the effect we seek. A prior sibling-mortality variable would be a different variable requiring its own justification; the paper does not make that substitution.
3. **Assign PfPR to the time actually at risk.** An effect of exposure during ages 36–59 months should not automatically use prevalence at birth. Splitting at relevant age and calendar boundaries permits annual MAP and external covariates to track that follow-up. For wide-band risk models, the exposure history within the band needs an explicit definition.
4. **Choose risk versus rate explicitly.** A monthly logistic model can produce age-band risks by multiplying survival probabilities and then standardizing individual predictions. The existing primary plan targets deaths per person-time and uses a negative-binomial log-link model with an offset. Burstein's likelihood does not resolve or reproduce that rate estimand automatically.
5. **Use predictive validation as a model check.** Holding out surveys or countries can reveal poor transportability. The paper's reported validation primarily assesses mortality reconstruction from held-out summaries; it does not demonstrate identification of a causal exposure effect.

For a risk analysis, use the same target population under PfPR = X and PfPR = Y, calculate individual age-band risks, then average them with the chosen population weights. Report risk ratios and differences from those standardized predictions. For a rate analysis, define the corresponding common person-time standardization. In either case, exponentiating a logit contrast produces an **odds ratio**, not automatically the desired mortality risk or rate ratio.

No hypothetical children, probability-of-birth weights or summary-history imputation are needed when the observed complete birth histories themselves are the analysis data.

## Review scope

Read and visually checked the model pages of the supplied PDF and inspected five public R scripts from the authors' repository. Downloaded scripts were read as text only. No external code was executed, no DHS data were accessed and no empirical models were fitted. The comparison uses the repository's `master` contents retrieved on 7 September 2026; the launchers contain multiple model configurations, so they are not proof of which archived run produced every published result.

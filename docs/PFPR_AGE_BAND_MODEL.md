# Proposed joint PfPR model using Burstein's age bands

Draft for the next model trial, 7 September 2026. The user now proposes Burstein's seven age bands and a full model with confounders X and an age-dependent PfPR effect. This develops that proposal as a joint model. The earlier six-band, separately fitted models remain historical specifications; the existing data-building and fitting scripts have not yet been changed.

The model retains the observed complete-birth-history survival structure. It excludes maternal death fraction, hypothetical children, probability-of-birth reconstruction and all other machinery used to apply the Burstein model to summary histories. See the [method comparison](BIRTH_HISTORY_MODEL_COMPARISON.md).

## 1. Age bands and observation unit

For an initial trial using the authors' month-based implementation:

| Band | Completed months | Interval width in years |
|---|---|---|
| Neonatal | <1 | 1/12 |
| Early post-neonatal | 1–5 | 5/12 |
| Late post-neonatal | 6–11 | 6/12 |
| Age 1 | 12–23 | 1 |
| Age 2 | 24–35 | 1 |
| Age 3 | 36–47 | 1 |
| Age 4 | 48–59 | 1 |

Burstein's text instead describes a day-based neonatal boundary (birth through day 28, followed by day 29 onward), and their two launchers differ. The table follows their external-fitting launcher's [month boundaries](https://github.com/royburst/sbh_agespecific_indirect_paper_code/blob/master/EXTVAL_launch_models_fitting.R#L64-L71). Preserve recorded neonatal days during preprocessing; adopting the paper's exact day boundary requires correspondingly defined first and second interval widths.

Each row represents an actual child in an age band reached alive. `death` is 1 if the child died in that band, otherwise 0 if it survived the band. For this initial complete-band implementation, require that the child's potential band end is on or before interview, applying the same rule to deaths and survivors. Do not include later bands after a death. This eligibility rule loses partial follow-up, as discussed in the comparison.

`band_years` is the full, predetermined width of the age band. **It is not the child's observed time to death.** The current deaths/person-time CSV does not contain these child-band trials or their binomial denominators and cannot be used directly for this specification. Aggregation is possible over children sharing all model predictors and the same observation interval, using integer deaths and survivors.

## 2. Likelihood and linear predictor

Let qᵢₐ be the probability of dying in band a conditional on entering it alive, and λᵢₐ its constant mortality hazard per person-year under the specified covariate profile. Then

\[
D_{ia}\sim\mathrm{Bernoulli}(q_{ia}),\qquad
q_{ia}=1-\exp\{-\Delta_a\lambda_{ia}\},
\]

\[
\log\{-\log(1-q_{ia})\}
=\log(\Delta_a)+\eta_{ia},
\]

\[
\eta_{ia}=\log\lambda_{ia}
=\alpha_a+f_a(P_{ia})+h_a(t_{ia})
+\mathbf X_{ia}^{\mathsf T}\boldsymbol\beta_a
+u_{s(i)}+v_{c(i),a}+b_{r(i)}.
\]

Here:

- **αₐ:** one baseline log rate per age band.
- **fₐ(P):** a separate nonlinear PfPR₂–₁₀ curve per band, with PfPR expressed in percentage points. No smoothness across age-band boundaries is imposed.
- **hₐ(t):** a separate smooth calendar-time trend per band. The time index must describe the band's follow-up period; it is not automatically birth year.
- **X′βₐ:** the chosen confounders, with coefficients allowed to differ by age band. A simpler X′β version shares these coefficients; it is a restriction to assess, not required by a joint fit.
- **uₛ:** a survey random intercept, shared across age bands.
- **v꜀ₐ:** a country-by-age-band random intercept.
- **bᵣ:** a regional random intercept, with region identifiers nested within country and harmonized across surveys where possible. This is an extension to Burstein's hierarchy because our exposure varies regionally.

The random effects are independent mean-zero Gaussians with separate variances for survey, country-by-age, and region effects. Region-specific residual effects could also vary by age in a later sensitivity. A district adjacency-based spatial effect is another possible extension, rather than a requirement of the first trial. Random intercepts represent residual clustering; they do not by themselves guarantee control of spatial confounding. A region fixed-effect or within-region analysis would probe reliance on between-region comparisons.

The complementary log–log link is the deliberate change from Burstein's logit. It connects band probabilities to a log mortality rate under a constant-hazard interpretation. Without a constant within-band hazard, −log(1−q)/Δ is an interval-average cumulative hazard per unit time; it should not silently be labelled the empirical deaths/person-time rate. Both the hazard assumption and the target population need to be stated when comparing with the existing rate analysis.

Because the interval width is fixed within each band, its offset can be absorbed into the band intercepts without changing fitted probabilities. Including it makes the rate units explicit; the link function supplies the hazard interpretation.

## 3. mgcv specification

```r
library(mgcv)

# d contains prepared, eligible child-band records.
d$age_band <- factor(d$age_band, ordered = FALSE)
d$survey <- factor(d$survey)
d$region <- factor(d$region)  # globally unique country/region key
d$country_age <- interaction(d$country, d$age_band, drop = TRUE)

fit <- bam(
  death ~ 0 + age_band +
    s(pfpr_pct, by = age_band, bs = "cr", k = 5) +
    s(calendar_year, by = age_band, bs = "cr", k = 6) +
    age_band:(X1 + X2 + X3) +
    s(survey, bs = "re") +
    s(country_age, bs = "re") +
    s(region, bs = "re") +
    offset(log(band_years)),
  family = binomial(link = "cloglog"),
  data = d,
  method = "fREML",
  discrete = TRUE,
  na.action = na.fail
)
```

`X1`, `X2`, `X3` illustrate named numeric or factor confounders; replace them with the reviewed adjustment set. Continuous confounders need not be linear: an age-specific `s(X1, by = age_band, ...)` can replace its linear interaction. Use a set of specified effects rather than one high-dimensional joint spline over all confounders.

An unordered factor `by` term creates one PfPR smooth per age band. The `age_band` main effect supplies the band intercepts required alongside these centred smooths. Omitting `id` allows separate smoothing parameters; sharing an `id` shares the smoothing parameter, not the fitted curve. The basis dimensions shown are initial choices to check against available variation and diagnostics. [mgcv factor smooths](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/factor.smooth.html).

The example shows the unweighted conditional likelihood. A DHS analysis also needs an explicit weighting and uncertainty strategy. If using normalized survey weights in `bam`, pass them through `weights`; they are likelihood weights, not a replacement for handling survey clustering and stratification. Retain child, mother, PSU and stratum identifiers in the internal data build so these choices remain possible. [bam weights and fitting](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/bam.html). The binomial family supports the complementary log–log link directly. [R family documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/family.html).

## 4. PfPR contrasts

Use p₀ and p₁ for the two PfPR values, to distinguish them from the confounder vector X. For band a, holding the other predictors fixed,

\[
RR_a(p_0\to p_1)=\frac{\lambda_a(p_1)}{\lambda_a(p_0)}
=\exp\{f_a(p_1)-f_a(p_0)\}.
\]

This is a **rate/hazard ratio under the model**, with percentage change 100(RRₐ−1). It is generally different from the ratio of band death probabilities qₐ(p₁)/qₐ(p₀). The constant-rate assumption permits either quantity to be calculated; report them separately.

This contrast concerns an exposure change during the band in a defined population alive at its start. An intervention changing malaria exposure from birth can also change who survives into later bands; its total effect requires a corresponding life-course prediction rather than interpreting each conditional contrast as a separate total effect.

For example, comparing 40% with 20% PfPR uses exp{fₐ(20)−fₐ(40)}. No global PfPR coefficient is needed alongside the seven smooths. Adding PfPR-by-confounder interactions would make the relative effect depend on X as well; none are assumed in this initial specification.

For absolute effects, predict both scenarios for the same target population and average the resulting rates with the chosen common person-time weights, or average the band probabilities with the chosen common entrant-population weights. Report, for example, deaths per 1,000 person-years and deaths per 1,000 children entering the band as distinct quantities. Rate-ratio contrasts should be restricted to supported exposure ranges within each age group. Calculate uncertainty from contrasts of the full prediction matrix so the covariance of fₐ(p₀) and fₐ(p₁) is retained.

## 5. Exposure timing and partial follow-up

The displayed simple model assumes a defined band-level exposure and covariate profile. It does not specify which annual MAP values to use. Fix that definition before fitting. An exposure summary must not depend on how early the child died; for example, averaging only the child's observed pre-death exposure can make the summary itself depend on the outcome.

If the scientific intervention changes annual PfPR during a band, the more explicit model is

\[
q_{ia}=1-\exp\left\{-\sum_y
\Delta_{iay}\exp\left[\alpha_a+f_a(P_{ry})+h_a(y)
+\mathbf X_{iay}^{\mathsf T}\boldsymbol\beta_a
+u_s+v_{ca}+b_r\right]\right\}.
\]

Here Δᵢₐᵧ is the portion of the potential band in calendar year y. The nonlinear annual hazards must be integrated; evaluating fₐ at an average PfPR is generally not the same model. When observed event timing permits splitting into smaller intervals, the same predictor can be used on those intervals. When only the age band containing death is known, retain that interval censoring rather than assigning a fictitious exact death month. Integrating over those unresolved intervals requires an extension beyond the simple single-row `bam` formula above.

If actual time at risk until death or censoring is available at the required resolution, a piecewise-exponential implementation instead uses a Poisson likelihood representation with an offset for **actual person-time** and the same log-rate predictor. This permits partial follow-up. It must not be confused with putting actual time-to-death into the complementary log–log binomial offset: the binomial offset above is the predetermined observation-interval width.

## 6. Relationship to previous specifications

This is one joint model with seven PfPR curves. It shares survey/region effects and random-effect variance parameters, so it is not identical to seven independent fits. It also differs from the currently saved six-band negative-binomial count models and from the Adebayo monthly-logit template. The new age bands and exclusion of maternal death fraction follow the user's latest direction; the complementary log–log link, age-specific confounder effects and additional regional term are recommendations for this trial.

## 7. Verification

The displayed R block was fitted using mgcv 1.9-4 to simulated birth histories with 7,039 eligible child-band observations and 513 deaths. The fit converged, retained all input rows and produced seven PfPR smooths. Rate ratios calculated from linear-predictor contrasts matched ratios of rates reconstructed from fitted probabilities. This checks implementation and interpretation on synthetic data; no DHS records or empirical estimates were used.

# Adebayo and Fahrmeir (2005): model specification in mgcv

The accompanying [R specification](../R_dhs/specifications/adebayo_2005_mgcv.R) implements the structure of models **M7 and M6**, with the optional socioeconomic extension in Section 4.3. M7 was the authors' preferred specification. This is a proposed model reference; it does not change the primary analysis in [ANALYSIS_PLAN.md](ANALYSIS_PLAN.md).

See also the [comparison with Burstein et al. (2018)](BIRTH_HISTORY_MODEL_COMPARISON.md), including their complete-birth-history likelihood, age-band eligibility rules and implications for the PfPR analysis.

Source: Adebayo SB, Fahrmeir L. *Analysing child mortality in Nigeria with geoadditive discrete-time survival models.* Statistics in Medicine 2005;24:709–728. DOI: [10.1002/sim.1842](https://doi.org/10.1002/sim.1842). Equations (1)–(10), pp. 714–716, define the model family; p. 718 specifies M6/M7; Tables IV and VII identify the fixed effects.

## Likelihood and data layout

For child i and month of life t, define

\[
h_{it}=P(T_i=t\mid T_i\geq t,\mathbf{x}_{it}),\qquad
y_{it}\mid\text{at risk}\sim\operatorname{Bernoulli}(h_{it}),\qquad
\operatorname{logit}(h_{it})=\eta_{it}.
\]

Use one row for each observed month at risk. Set `event = 1` in the month of death and `0` in preceding survived months. A right-censored child contributes only observed, survived months. Include no rows after death or censoring. The resulting product of Bernoulli likelihood contributions is the discrete survival likelihood; repeated rows are not separate children.

In the template, `month = 1` represents age [0,1) months, `month = 2` represents [1,2), and so on. Thus `completed_month = month - 1`. The paper follows children for up to 36 months; extending to 60 months is a project-specific change. **Time here means attained age, not calendar year.**

A death in an interval contributes one event trial even if the child dies before its end. For a survivor censored within an unfinished interval, do not silently count that interval as a complete survived month. Establish a consistent interval convention from the available DHS dates; the template expects that preprocessing to have been done. No person-time offset is needed for these equal one-month Bernoulli trials.

The existing `results/dhs_rebuild/person_time_model_data.csv` aggregates deaths and person-time across children, age segments and retrospective windows. It cannot supply this individual-history likelihood or these individual covariates directly. Exact aggregation is possible only over trials sharing **all** model predictors, using integer deaths and numbers of trials in `cbind(deaths, trials - deaths)`; aggregated person-months or weighted effective counts are not automatically binomial trial counts.

## M7 and its mgcv translation

\[
\eta_{it}=\alpha+f_0(t)+u_{s_i}+v_{s_i}
  +bft_{it}f_1(t)+ma1_i f_2(t)+ma2_i f_3(t)
  +\mathbf{z}_i^\top\boldsymbol\gamma.
\]

The paper absorbs the intercept into its baseline function. Here `mgcv` represents it with an intercept and a centred baseline smooth. The two district terms are a spatially structured effect u and an independent district effect v.

```r
library(mgcv)

# person_months must already contain the variables described below.
d <- person_months
d$ma1 <- as.numeric(d$maternal_age_birth < 22) -
  as.numeric(d$maternal_age_birth > 35)
d$ma2 <- as.numeric(d$maternal_age_birth >= 22 & d$maternal_age_birth <= 35) -
  as.numeric(d$maternal_age_birth > 35)
d$district <- factor(d$district, levels = names(nb))
d$district_iid <- d$district

fit_m7 <- gam(
  event ~
    s(month, bs = "ps", k = 10, m = c(2, 2)) +
    s(month, by = bft, bs = "ps", k = 10, m = c(2, 2)) +
    s(month, by = ma1, bs = "ps", k = 10, m = c(2, 2)) +
    s(month, by = ma2, bs = "ps", k = 10, m = c(2, 2)) +
    s(district, bs = "mrf", xt = list(nb = nb)) +
    s(district_iid, bs = "re") +
    urban_ec + male_ec + primary_or_less_ec + assisted_ec +
    hospital_ec + long_interval_ec + assisted_ec:antenatal_ec,
  data = d,
  family = binomial(link = "logit"),
  method = "REML",
  na.action = na.fail,
  drop.unused.levels = FALSE,
  knots = list(district = names(nb))
)
```

| Paper term | Meaning | Implementation |
|---|---|---|
| f₀(t) | Baseline log odds of death as a function of attained age | Intercept plus `s(month, bs = "ps")` |
| bft × f₁(t) | Breastfeeding association changes with the child's age | Numeric `by = bft`, with bft = 1 if breastfed in that month, otherwise 0 |
| ma1 × f₂(t), ma2 × f₃(t) | Maternal-age-category associations change with the child's age | Numeric effect-coded `by` terms defined above |
| uₛ | Adjacent districts have similar effects | Full district MRF using the adjacency graph |
| vₛ | Independent residual district heterogeneity | Gaussian district random intercept |
| z′γ | Remaining covariate effects are constant over the child's age | Parametric terms |

For mothers under 22, aged 22–35, and over 35, the maternal contributions are respectively f₂(t), f₃(t), and −f₂(t)−f₃(t). These are the paper's effect codes, rather than reference-category dummy variables. Keep `ma1`, `ma2` and `bft` **numeric**: numeric `by` smooths multiply a function by the supplied covariate and ordinarily retain its constant component. Additional main effects for these variables would duplicate components already represented by these smooths. Factor `by` smooths follow different centring rules. [mgcv model specification](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.models.html).

The spline choice `bs = "ps", m = c(2,2)` uses cubic B-splines with a second-difference penalty, corresponding to the second-order random-walk option in equation (9). `m = c(2,1)` gives the first-difference alternative. `k = 10` is a starting basis dimension for this implementation, not a recovered setting from the paper; assess it against the age support and fitted diagnostics. [mgcv P-splines](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.ps.smooth.spec.html).

Supply `nb` as a named symmetric adjacency list, with names matching district factor levels and entries identifying neighbouring districts. Its MRF penalty is the graph Laplacian D − W, corresponding to equation (10). Do not supply `k` for this term if retaining one coefficient per district, subject to identifiability constraints. The separate `district_iid` column is a copy of the same district factor. The template assumes a connected graph: islands and disconnected country graphs require an explicit choice of component constraints and intercepts before use. Interpret the combined district effect cautiously when the structured and independent components are weakly distinguished. [mgcv MRF smooths](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.mrf.smooth.spec.html).

## Covariates, M6 and the socioeconomic extension

The template uses +1 for the named category and −1 for its comparison category:

| Column | Category coded +1 |
|---|---|
| `urban_ec` | Urban residence |
| `male_ec` | Male child |
| `primary_or_less_ec` | Mother's education at most primary |
| `assisted_ec` | Assistance at delivery |
| `hospital_ec` | Hospital delivery |
| `long_interval_ec` | Long preceding birth interval, using the paper's category definition |
| `antenatal_ec` | Antenatal visit |

The assistance × antenatal term is implemented as the product of these two numeric effect codes. Table IV retains that interaction but omits the antenatal main effect; the formula follows that term structure. Confirm the original category and interaction construction before attempting numerical replication. These names specify model inputs, not a validated crosswalk to DHS recode fields.

**M6** replaces the two maternal-category varying-coefficient terms with a smooth of the mother's age **at the child's birth**:

```r
s(maternal_age_birth, bs = "ps", k = 10, m = c(2, 2))
```

It retains the age-varying breastfeeding term. M6 does not include both that maternal-age smooth and the two maternal-category smooths.

The optional Section 4.3 extension adds six effect-coded main effects: `working_ec`, `household_le5_ec`, `toilet_ec`, `electricity_ec`, `quality_house_ec`, and `water_residence_ec`. Their +1 categories are mother currently working, household size ≤5, has toilet, has electricity, high-quality wall and floor, and water in the residence.

```r
source("R_dhs/specifications/adebayo_2005_mgcv.R")
fit_m7 <- fit_adebayo_2005(person_months, nb)
fit_m6 <- fit_adebayo_2005(person_months, nb, model = "M6")
fit_m7_extended <- fit_adebayo_2005(person_months, nb,
                                  include_socioeconomic = TRUE)
```

## Estimation and relation to the PfPR question

The implementation matches the likelihood and additive structure, with the stated spline choices. **REML is a penalized-likelihood analogue, not a reproduction of the paper's fully Bayesian inference.** The authors use BayesX MCMC, inverse-gamma hyperpriors on variance components, posterior intervals and DIC. `gam(method = "REML")` estimates smoothing parameters and does not reproduce those hyperpriors, intervals or DIC. Neither the basic formula nor the wrapper provides DHS survey-design uncertainty automatically.

For our question, a possible extension is

\[
\operatorname{logit}(h_{it,g})=
\alpha_g+a_g(t)+f_g(PfPR_{it})+q_g(\text{calendar year}_{it})
+\mathbf{c}_{it}^{\top}\boldsymbol\beta_g+u_{s_i,g}+v_{s_i,g},
\]

fitted **separately** in the six confirmed groups: <1, 1–3, 4–11, 12–23, 24–35 and 36–59 completed months. For example, `s(pfpr_pct, bs = "ps", k = 5, m = c(2,2))` represents a nonlinear prevalence effect within a group. PfPR and calendar year refer to the interval at risk. This extension is not in the 2005 paper, and its adjustment set needs to follow the PfPR causal question; the paper's selected covariates, including breastfeeding, are not automatically suitable confounders.

Within the neonatal group there is only one monthly age value, so use an intercept without an age smooth. For 1–3 completed months, a three-level age factor is a straightforward baseline. Wider groups can support an age smooth with a suitable basis dimension. The paper's age-varying maternal and breastfeeding terms also cannot simply be copied with `k = 10` into these narrow groups.

To compare X% with Y%, predict the same target population under both exposure scenarios. For an individual profile and a band spanning monthly intervals a through b, its probability of dying in that band, **conditional on entering the band alive**, is

\[
Q_g(x)=1-\prod_{t=a}^{b}\{1-h_{it,g}(x)\}.
\]

Compute these products per profile before standardizing over the chosen population. A contrast `exp(eta(Y) - eta(X))` is a monthly **odds ratio**; it is not exactly a mortality risk ratio or the deaths-per-person-time rate ratio in the current primary plan. To retain the planned rate estimand, explicitly define the person-time standardization or retain a count/rate likelihood with a person-time offset. Adopting this paper's likelihood therefore requires an explicit choice of outcome scale as well as a finer dataset.

## Verification

Checked with R 4.6.0 and mgcv 1.9-4 on 7 September 2026. M7, M6 and the socioeconomic extension converged on simulated censored birth histories containing 7,535 child-months and 127 deaths. Checks confirmed the logit link, retention of all input rows, finite probabilities between zero and one, and rejection of factor-valued breastfeeding indicators and non-Bernoulli event counts. No DHS records were read and no empirical model was fitted. These checks establish that the specifications execute; they do not validate a DHS data build or reproduce the paper's estimates.

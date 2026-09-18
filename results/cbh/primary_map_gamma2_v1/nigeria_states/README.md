# Figure 6: Nigerian states versus IHME, 2024

Secondary subnational extrapolation of the primary MAP gamma=2 model to the 37 Nigerian states, using the same convention as the national burden: IHME all-cause deaths by age band multiplied by the age-specific attributable fraction at the state's MAP prevalence. Calculation code is [R_cbh/burden/05_nigeria_state_burden.R](../../../../R_cbh/burden/05_nigeria_state_burden.R); this figure and table are produced by [R_cbh/reporting/08_nigeria_state_comparison.R](../../../../R_cbh/reporting/08_nigeria_state_comparison.R).

## Input audit

- Primary model: seven verified saved MAP gamma=2 fits; no refitting.
- IHME all-cause by state and GBD age group, 2024: the seven constructed bands reproduce every state's under-5 total, and their sum over states reproduces the national 2024 age-band inputs used by `primary/02_effects.R` to within 1e-6, confirming the same export vintage.
- IHME malaria by state, 2024: under-5 death rate per 100,000, applied to all-cause implied person-years.
- MAP PfPR[2–10] by state, 2024: `data/pfpr_admin1_ng_cd_2024.csv`, an inherited extraction (MAP admin-1 boundaries, GPW 2020 density weights without cell area). All 37 states lie inside the central 95% of fitted prevalence; none requires extrapolation beyond observed support.
- Source metadata and hashes: [ihme_source.json](ihme_source.json), [input_provenance.csv](input_provenance.csv).

## Reconciliation with the national estimate

State sum 209,594 attributable deaths versus the national estimate 210,622 evaluated at the national mean PfPR of 24.70% (state person-year-weighted mean 25.21%): difference -0.49%. The national and state exposure extractions use different boundary sources and population weights (see the note in [national_reconciliation_2024.csv](national_reconciliation_2024.csv)). IHME state malaria deaths sum to 132,138, so the model attributes 1.59 times as many under-5 deaths to malaria as IHME's cause-specific estimate.

## Pattern of the discrepancy

The model exceeds IHME in 36 of 37 states. The ratio of model to IHME deaths has median 1.67 in the three northern zones and 1.18 in the three southern zones; its Spearman correlation with state PfPR is 0.10. IHME's malaria share of all-cause under-5 deaths ranges from 13.5% to 29.9% and is highest in the southeast, whereas the model's attributable fraction ranges from 12.9% to 32.9% and rises with prevalence. The two methods therefore disagree on geography more than on the national total.

## State table

Sorted by the ratio of model-attributable to IHME malaria deaths. Deaths are 2024 annual counts before age 5.

| State | Zone | PfPR (%) | IHME all-cause | PfPR-ACM model | Attributable fraction | IHME malaria | IHME malaria share | Model / IHME |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Sokoto | North West | 24.3 | 50,749 | 16,190 | 31.9% | 6,831 | 13.5% | 2.37 |
| Kebbi | North West | 29.9 | 27,830 | 9,089 | 32.7% | 4,097 | 14.7% | 2.22 |
| Zamfara | North West | 34.7 | 20,179 | 6,613 | 32.8% | 3,012 | 14.9% | 2.20 |
| Jigawa | North West | 23.3 | 37,811 | 11,856 | 31.4% | 5,606 | 14.8% | 2.11 |
| Bauchi | North East | 24.5 | 33,384 | 9,787 | 29.3% | 4,781 | 14.3% | 2.05 |
| Taraba | North East | 28.5 | 36,844 | 10,548 | 28.6% | 5,328 | 14.5% | 1.98 |
| Gombe | North East | 23.3 | 15,302 | 4,542 | 29.7% | 2,427 | 15.9% | 1.87 |
| Kaduna | North West | 25.5 | 41,097 | 12,429 | 30.2% | 6,656 | 16.2% | 1.87 |
| Katsina | North West | 25.7 | 44,458 | 14,363 | 32.3% | 7,840 | 17.6% | 1.83 |
| Ebonyi | South East | 26.0 | 9,228 | 2,744 | 29.7% | 1,502 | 16.3% | 1.83 |
| Yobe | North East | 26.5 | 16,062 | 4,451 | 27.7% | 2,610 | 16.2% | 1.71 |
| Nasarawa | North Central | 28.8 | 19,142 | 5,230 | 27.3% | 3,203 | 16.7% | 1.63 |
| Bayelsa | South South | 26.1 | 1,598 | 427 | 26.7% | 269 | 16.8% | 1.59 |
| Niger | North Central | 30.1 | 30,818 | 8,528 | 27.7% | 5,475 | 17.8% | 1.56 |
| Kano | North West | 19.8 | 91,761 | 27,423 | 29.9% | 18,251 | 19.9% | 1.50 |
| Osun | South West | 26.7 | 6,792 | 1,765 | 26.0% | 1,241 | 18.3% | 1.42 |
| Adamawa | North East | 24.3 | 10,734 | 2,739 | 25.5% | 1,965 | 18.3% | 1.39 |
| Kogi | North Central | 25.3 | 5,277 | 1,386 | 26.3% | 1,014 | 19.2% | 1.37 |
| FCT (Abuja) | North Central | 25.3 | 17,161 | 4,577 | 26.7% | 3,399 | 19.8% | 1.35 |
| Plateau | North Central | 19.8 | 34,127 | 7,722 | 22.6% | 5,833 | 17.1% | 1.32 |
| Benue | North Central | 29.6 | 23,741 | 5,681 | 23.9% | 4,362 | 18.4% | 1.30 |
| Borno | North East | 18.9 | 13,264 | 3,136 | 23.6% | 2,437 | 18.4% | 1.29 |
| Edo | South South | 28.1 | 7,317 | 1,965 | 26.9% | 1,553 | 21.2% | 1.27 |
| Oyo | South West | 26.4 | 8,565 | 2,352 | 27.5% | 1,895 | 22.1% | 1.24 |
| Ondo | South West | 32.0 | 6,212 | 1,971 | 31.7% | 1,611 | 25.9% | 1.22 |
| Abia | South East | 29.3 | 13,063 | 4,299 | 32.9% | 3,558 | 27.2% | 1.21 |
| Cross River | South South | 31.1 | 3,634 | 1,014 | 27.9% | 841 | 23.1% | 1.21 |
| Kwara | North Central | 24.0 | 4,644 | 1,158 | 24.9% | 967 | 20.8% | 1.20 |
| Ekiti | South West | 29.2 | 3,691 | 1,081 | 29.3% | 914 | 24.8% | 1.18 |
| Delta | South South | 20.0 | 9,654 | 2,207 | 22.9% | 1,868 | 19.3% | 1.18 |
| Ogun | South West | 23.1 | 29,614 | 6,292 | 21.2% | 5,359 | 18.1% | 1.17 |
| Rivers | South South | 26.6 | 8,576 | 2,408 | 28.1% | 2,099 | 24.5% | 1.15 |
| Imo | South East | 30.1 | 22,305 | 7,263 | 32.6% | 6,664 | 29.9% | 1.09 |
| Akwa Ibom | South South | 31.3 | 3,523 | 908 | 25.8% | 838 | 23.8% | 1.08 |
| Enugu | South East | 20.6 | 6,832 | 1,843 | 27.0% | 1,707 | 25.0% | 1.08 |
| Anambra | South East | 21.7 | 9,549 | 2,506 | 26.2% | 2,494 | 26.1% | 1.00 |
| Lagos | South West | 10.6 | 8,501 | 1,100 | 12.9% | 1,635 | 19.2% | 0.67 |

[State totals](state_totals_2024.csv) · [State-by-age estimates](state_age_estimates_2024.csv) · [Reconciliation](national_reconciliation_2024.csv) · [Model specification](model_specification.txt)

![Figure 6](fig6_nigeria_states_vs_ihme.png)

PfPR-ACM model malaria-attributable deaths before age 5 versus IHME cause-specific malaria deaths, by Nigerian state (36 states and the Federal Capital Territory) in 2024. Both axes use a base-10 logarithmic scale with identical limits; the dashed line denotes equality. Colours group states into Nigeria's six geopolitical zones for orientation; all estimates are made state by state. Model-attributable deaths equal IHME state all-cause deaths in each of seven age bands multiplied by 1 − exp[f_g(0) − f_g(P_state)], summed over bands, where P_state is population-weighted MAP PfPR[2–10] for 2024 and f_g are the seven separate primary MAP gamma=2 age-band effects. IHME early and late neonatal deaths share the <1-month effect; the IHME 2–4 year group shares one rate and divides deaths and person-time equally across the 24–35, 36–47 and 48–59 month bands. IHME malaria deaths are the state malaria death rate applied to the under-5 person-years implied by the all-cause count and rate. Summed over states the model gives 209,594 deaths against 132,138 IHME malaria deaths (ratio 1.59); the state sum differs from the national estimate evaluated at national mean prevalence by -0.5%, because the age-band curves are nonlinear. State totals are point estimates; conditional age-band intervals are retained in the source table and are not summed. Model-attributable all-cause reductions and IHME cause-specific malaria deaths are different estimands, and the model shares IHME all-cause inputs, so agreement or disagreement is not independent validation. Zero prevalence lies below observed exposure support.

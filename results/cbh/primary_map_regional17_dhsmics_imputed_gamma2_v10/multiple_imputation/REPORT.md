# Multiple-imputation check: imputed-covariate sensitivity, DHS and MICS

The point fit in `primary_map_regional17_dhsmics_imputed_gamma2_v10` uses one imputed dataset (mean of the mice imputations, GAM fitted means and posterior-median HIV incidence). This check refits all seven age-band models in each of 10 imputed datasets with the smoothing parameters fixed at the point fit, and pools the PfPR log hazard ratios with Rubin's rules (within-imputation variance from the conditional covariance, between-imputation variance across the 10 fits; Barnard–Rubin degrees of freedom).

Varied components per imputation: the mice imputations of the regional covariates; the GAM draws of health expenditure (Zimbabwe, Somalia, South Sudan, 2024), GDP and political stability (South Sudan before independence); and one posterior draw of the child HIV incidence series for every country-year whose value is imputed or censored (São Tomé and Príncipe among them; Liberia's rate is derived from UNAIDS counts and held fixed). Observed values never change. Records affected:

- regional covariates (any):   795,848 records
- health expenditure:   260,385 records
- GDP per capita:    33,840 records
- political stability (model fill):    57,509 records
- child HIV incidence (imputed or censored series): 1,152,154 records

![Pooled versus point curves](pooled_vs_point_curves.png)

### PfPR 40% to 20%

| Age (months) | Point imputation HR (95% interval) | Pooled MI HR (95% interval) | Between-imputation share of variance |
|---|---:|---:|---:|
| <1 | 0.94 (0.91–0.97) | 0.94 (0.91–0.96) | 0.9% |
| 1-5 | 0.91 (0.88–0.95) | 0.91 (0.88–0.95) | 4.5% |
| 6-11 | 0.82 (0.78–0.86) | 0.82 (0.78–0.86) | 3.6% |
| 12-23 | 0.84 (0.79–0.89) | 0.84 (0.79–0.89) | 1.4% |
| 24-35 | 0.86 (0.81–0.92) | 0.86 (0.81–0.92) | 1.3% |
| 36-47 | 0.85 (0.80–0.92) | 0.85 (0.79–0.91) | 0.9% |
| 48-59 | 0.94 (0.86–1.02) | 0.94 (0.86–1.02) | 0.5% |

### PfPR 20% to 0%

| Age (months) | Point imputation HR (95% interval) | Pooled MI HR (95% interval) | Between-imputation share of variance |
|---|---:|---:|---:|
| <1 | 0.93 (0.88–0.98) | 0.93 (0.88–0.98) | 0.5% |
| 1-5 | 0.87 (0.81–0.93) | 0.87 (0.81–0.93) | 2.7% |
| 6-11 | 0.69 (0.63–0.75) | 0.69 (0.63–0.75) | 1.8% |
| 12-23 | 0.55 (0.49–0.62) | 0.55 (0.49–0.62) | 1.1% |
| 24-35 | 0.46 (0.40–0.52) | 0.46 (0.40–0.52) | 1.9% |
| 36-47 | 0.55 (0.48–0.63) | 0.55 (0.48–0.63) | 0.9% |
| 48-59 | 0.57 (0.49–0.67) | 0.57 (0.49–0.66) | 0.9% |

The between-imputation share is (1 + 1/M) B / T. A small share means the covariate imputation adds little to the sampling uncertainty already in the point fit. Smoothing parameters are fixed, so smoothing uncertainty is not part of either interval. Survey design and exposure uncertainty remain unpropagated.

Reproduce: `Rscript R_cbh/sensitivity/imputation_dhsmics/05_propagate.R` after scripts 01–04 in the same folder.

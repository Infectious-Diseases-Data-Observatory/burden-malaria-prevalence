# Disjoint attribution of covariate exclusions

Of 6,357,802 MAP-eligible child-band records, 892,497 were excluded because at least one of the 17 covariates was unavailable after the HIV, UNICEF vaccination and available-region substitutions. Each excluded record is attributed here to the first missing covariate in the order: national annual series (child HIV incidence, health expenditure, political stability, GDP), then regional summaries (wasting, stunting, facility delivery, electricity, wealth, vaccination, birth interval, water, sanitation, urban, education, age at first birth). Counts are therefore disjoint and sum to the total; a record missing several covariates is counted once, under the first.

| Attributed reason | Records | Deaths | Surveys affected | Share of excluded |
|---|---:|---:|---:|---:|
| Child HIV incidence (national series unavailable) | 142,323 | 2,225 | 5 | 15.9% |
| Health expenditure per capita | 63,255 | 850 | 7 | 7.1% |
| Political stability | 173,291 | 3,852 | 25 | 19.4% |
| Wasting prevalence | 399,258 | 5,740 | 9 | 44.7% |
| Household electricity | 83,796 | 1,591 | 1 | 9.4% |
| Wealth-quintile score | 30,574 | 954 | 9 | 3.4% |

Overlapping per-variable counts are in covariate_missingness.csv; per-survey attribution in exclusion_attribution_by_survey.csv. Reproduce: `Rscript R_cbh/covariates/14_exclusion_attribution.R`.

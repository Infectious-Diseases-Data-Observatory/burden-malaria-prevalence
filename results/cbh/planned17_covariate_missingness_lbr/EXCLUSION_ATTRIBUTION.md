# Disjoint attribution of covariate exclusions

Of 6,357,802 MAP-eligible child-band records, 783,834 were excluded because at least one of the 17 covariates was unavailable after the HIV, UNICEF vaccination and available-region substitutions. Each excluded record is attributed here to the first missing covariate in the order: national annual series (child HIV incidence, health expenditure, political stability, GDP), then regional summaries (wasting, stunting, facility delivery, electricity, wealth, vaccination, birth interval, water, sanitation, urban, education, age at first birth). Counts are therefore disjoint and sum to the total; a record missing several covariates is counted once, under the first.

| Attributed reason | Records | Deaths | Surveys affected | Share of excluded |
|---|---:|---:|---:|---:|
| Child HIV incidence (national series unavailable) | 10,839 | 92 | 1 | 1.4% |
| Health expenditure per capita | 63,255 | 850 | 7 | 8.1% |
| Political stability | 173,291 | 3,852 | 25 | 22.1% |
| Wasting prevalence | 422,079 | 6,168 | 10 | 53.8% |
| Household electricity | 83,796 | 1,591 | 1 | 10.7% |
| Wealth-quintile score | 30,574 | 954 | 9 | 3.9% |

HIV panel with Liberia's child incidence derived from UNAIDS counts (DHS part of the v7 primary). Overlapping per-variable counts are in ../planned17_covariate_missingness/covariate_missingness.csv; per-survey attribution in exclusion_attribution_by_survey.csv. Reproduce: `Rscript R_cbh/covariates/14_exclusion_attribution.R --liberia`.

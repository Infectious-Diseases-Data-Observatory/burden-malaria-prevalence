# Missingness with pre-nationwide zeros and available-region means

Additional scenario: each missing regional covariate receives the arithmetic mean of finite available regional estimates within the same survey, calculated separately for each variable. Every donor region counts once; child-band row counts do not weight this mean. Surveys with no available donor region remain missing. National country-year covariates do not use this fallback. Existing values are preserved and donor counts/imputation flags are retained in the private wide overlay. This scenario does not change the adopted primary data policy.

Scenario dated 17 September 2026. Apply zero only to missing Hib3, PCV and rotavirus coverage strictly before the first WHO-reported nationwide introduction year. Preserve all observed coverage. Where introduction is not recorded and the latest status is No, apply zero through that last reported year (2025). This includes missing values during partial rollout by user specification; it is not a claim that true coverage was zero. Introduction-year gaps remain missing. Existing UNICEF DTP3/measles fallback is retained. No base datasets, model inputs or fits are overwritten.

Records mean child–age-band rows. The previous-primary denominator matches the earlier missingness tables; the full MAP-eligible denominator is the starting pool for the revised regional adjustment model. These count surveys with at least one missing record, not necessarily surveys completely excluded.

| Sample | Records | Remaining vaccine gaps | % of records | % of previously vaccine-missing records | Surveys with vaccine gaps | Countries |
|---|---:|---:|---:|---:|---:|---:|
| previous_primary | 5,885,022 | 198,852 | 3.38% | 4.93% | 17 | 10 |
| eligible_MAP | 6,357,802 | 206,318 | 3.25% | 4.61% | 18 | 11 |

## Surveys with remaining vaccine gaps: full MAP-eligible pool

Survey years are registry start years; surveys can span two calendar years. Kenya 2003 is additional to the previous-primary sample; the other countries/surveys are present in both samples.

| Country | Number of surveys | Survey start years | Missing vaccine(s) | Records with vaccine gaps |
|---|---:|---|---|---:|
| Angola | 1 | 2011 | Hib3 | 7,950 |
| Burundi | 1 | 2016 | Rotavirus | 17,143 |
| Gambia | 1 | 2013 | PCV | 9,725 |
| Guinea | 1 | 2012 | Hib3 | 9,127 |
| Kenya | 1 | 2003 | Hib3 | 7,404 |
| Malawi | 1 | 2015 | PCV | 23,225 |
| Namibia | 1 | 2013 | Hib3 | 5,242 |
| Niger | 1 | 2012 | Hib3 | 16,827 |
| Rwanda | 1 | 2010 | PCV | 11,878 |
| Senegal | 7 | 2012; 2014; 2015; 2016; 2017; 2018; 2019 | PCV, Rotavirus | 88,829 |
| Zambia | 2 | 2013; 2018 | PCV, Rotavirus | 8,968 |

## All 22 required covariates

| Sample | Records with any remaining gap | % of records | Surveys with any gap | Countries |
|---|---:|---:|---:|---:|
| previous_primary | 1,127,585 | 19.16% | 41 | 25 |
| eligible_MAP | 1,538,399 | 24.20% | 65 | 30 |

Countries/surveys with any remaining adjustment-variable gap in the full MAP-eligible pool:

| Country | Surveys with any gap | Survey start years |
|---|---:|---|
| Angola | 1 | 2011 |
| Benin | 2 | 2001; 2006 |
| Burkina Faso | 1 | 2003 |
| Burundi | 1 | 2016 |
| Cameroon | 1 | 2004 |
| Chad | 1 | 2004 |
| Congo | 1 | 2005 |
| Democratic Republic of the Congo | 1 | 2007 |
| Eswatini | 1 | 2006 |
| Ethiopia | 3 | 2000; 2005; 2024 |
| Gabon | 1 | 2000 |
| Gambia | 1 | 2013 |
| Ghana | 1 | 2003 |
| Guinea | 2 | 2005; 2012 |
| Kenya | 1 | 2003 |
| Liberia | 4 | 2007; 2009; 2013; 2019 |
| Madagascar | 1 | 2004 |
| Malawi | 4 | 2000; 2004; 2015; 2024 |
| Mali | 3 | 2001; 2006; 2023 |
| Mozambique | 2 | 2003; 2011 |
| Namibia | 3 | 2000; 2006; 2013 |
| Niger | 2 | 2006; 2012 |
| Nigeria | 3 | 2003; 2010; 2024 |
| Rwanda | 4 | 2000; 2005; 2008; 2010 |
| Sao Tome and Principe | 1 | 2008 |
| Senegal | 9 | 2005; 2008; 2012; 2014; 2015; 2016; 2017; 2018; 2019 |
| Uganda | 3 | 2000; 2006; 2009 |
| United Republic of Tanzania | 1 | 2004 |
| Zambia | 4 | 2002; 2013; 2018; 2024 |
| Zimbabwe | 2 | 2005; 2010 |

[Variable totals](summary.csv), [individual surveys and covariates](survey.csv), [country totals](country.csv), [complete-case selection](selection.csv), [annual missingness](annual.csv).

Starting missingness was reconciled with the established UNICEF-fallback audit for both samples and all covariates. Boundary checks enforce strictly before nationwide rollout and preservation of observed values. Source hashes are in input_manifest.csv.

Reproduce: `Rscript R_cbh/covariates/10_prenationwide_zero_scenario.R --region-mean`.

## Where the regional gaps occur

Full MAP-eligible pool, before regional-mean substitution. Counts across covariates overlap.

| Covariate | Regions fillable from other regions | Surveys with fillable regions | Surveys without any donor region |
|---|---:|---:|---:|
| mean_wealth_quintile | 0 | 0 | 9 |
| facility_delivery_pct | 7 | 6 | 5 |
| exclusive_breastfeeding_pct | 32 | 1 | 33 |
| short_birth_interval_pct | 10 | 8 | 0 |
| improved_water_pct | 10 | 8 | 0 |
| improved_sanitation_pct | 10 | 8 | 0 |
| electricity_pct | 10 | 8 | 1 |

Within-region means already exclude missing individual responses and invalid weights. A regional mean is missing only when it cannot be estimated from usable source values under the existing extraction/quality rules. Available-region imputation does not solve whole-survey absence, questionnaire-review exclusions affecting all donor regions, or missing national country-year inputs.

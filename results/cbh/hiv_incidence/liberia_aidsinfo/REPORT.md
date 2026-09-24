# Child HIV incidence for Liberia from UNAIDS counts

The UNAIDS/UNICEF workbook (`data/HIV_Epidemiology_Children_Adolescents_2025.xlsx`) has no incidence series for Liberia, so Liberia was excluded from the complete-case primary. AIDSinfo publishes Liberia's annual new HIV infections among children aged 0–14 (UNAIDS epidemiological estimates 2026; retrieved 24 September 2026 and saved in `data/unaids_aidsinfo/`).

**Definition check.** The workbook's child rate (new infections per 1,000 uninfected population, age 0–14) is reproduced by 1,000 × new infections aged 0–14 / (under-5 person-years − children aged 0–4 living with HIV): across 761 country-years in 34 countries with at least 200 new infections, derived/published has median 0.988 (IQR 0.944–1.055; 5–95% 0.864–1.237). The under-5 person-years are the IHME-implied denominators of the burden step. Dividing by the whole population aged 0–14 instead gives about 0.38 of the published rate, so the published child rate is effectively per uninfected child under 5.

**Liberia.** The same formula with Liberia's AIDSinfo counts, IHME-implied under-5 person-years and the workbook's children aged 0–4 living with HIV:

| Year | New infections 0–14 | Rate per 1,000 (range from the count bounds) |
|---|---:|---:|
| 2000 | 708 | 1.32 (1.07–1.59) |
| 2005 | 867 | 1.53 (1.31–1.73) |
| 2010 | 859 | 1.26 (1.07–1.47) |
| 2015 | 727 | 1.01 (0.80–1.22) |
| 2020 | 624 | 0.84 (0.56–1.12) |
| 2024 | 359 | 0.46 (0.22–0.77) |

Published rates in neighbouring countries for comparison: CIV 2005: 3.46; CIV 2010: 1.66; CIV 2015: 0.61; CIV 2020: 0.40; CIV 2024: 0.23; GIN 2005: 1.54; GIN 2010: 1.30; GIN 2015: 1.35; GIN 2020: 0.76; GIN 2024: 0.48; SLE 2005: 2.10; SLE 2010: 1.46; SLE 2015: 1.14; SLE 2020: 1.01; SLE 2024: 0.55.

**Panels.** `child_incidence_country_year_lbr.csv` is the frozen base panel plus Liberia 2000–2024 (status `derived_aidsinfo_counts`), used by the primary from v7; `child_incidence_country_year_extended_lbr.csv` is the extended panel with Liberia's latent-adolescent imputation replaced by these values, used by the imputed-covariate sensitivity. The frozen panels are unchanged. Caveats: the counts come from the 2026 estimates round while the other countries' rates come from the 2025 workbook; Liberia's rate is fixed (its range is shown here but not propagated).

[Liberia series](liberia_child_incidence.csv) · [definition check by country](definition_check_by_country.csv) · [comparison with the extended model's imputation](liberia_derived_vs_extended_model.csv)

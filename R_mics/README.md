# MICS eligibility assessment

Assesses whether UNICEF MICS surveys can join the primary child age-band analysis. Microdata stay in the
git-ignored `data/` tree; only aggregate tables are written to `results/mics_inventory/`.

1. Extract the SPSS files from `data/MICS_Datasets/*.zip` into `data/MICS_extracted/<survey>/` (flattened, nested zips opened).
2. `Rscript R_mics/00_dump_metadata.R` catalogues every variable name and label into `data/derived_mics/variable_metadata.csv`.
3. `Rscript R_mics/01_inventory.R` reads each survey and writes `results/mics_inventory/survey_inventory.csv`: birth-history
   completeness and counts, region variable and labels, and the source variable for each covariate. Variable names are
   matched case-insensitively, with English, French and Portuguese label patterns and round-specific names
   (for example DTP3 is `IM4C` in MICS3, `IM3D3D`/`IM3PENTA3D` in MICS4-5, `IM6PENTA3D` in MICS6).
4. `Rscript R_mics/02_inventory_tables.R` writes `INVENTORY.md`.
5. `results/mics_inventory/eligibility.csv` holds the per-survey verdicts, from direct checks reconciled with per-country
   audit agents and adversarial re-checks; `python3 R_mics/03_eligibility_report.py` renders `ELIGIBILITY.md` from it.

Nothing here changes the primary analysis. Adding MICS would be a new primary version.
